// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "CSPrimitiveIndex.h"
#include "CSPrimitives.h"
#include "ContinuousStructure.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

CSPrimitiveIndex::CSPrimitiveIndex(ContinuousStructure& structure, const Config& config)
	: m_Config(config), m_Primitives(structure.GetAllPrimitives(true, config.type))
{
	std::size_t buckets = 1;
	for (unsigned int axis=0; axis<3; ++axis)
	{
		const double low = config.bounds[2*axis];
		const double high = config.bounds[2*axis+1];
		const auto count = config.bins[axis];
		if (!std::isfinite(low) || !std::isfinite(high) || high <= low ||
			!std::isfinite(high-low) || count == 0)
		{
			throw std::invalid_argument("Spatial index needs finite increasing bounds and nonzero bins");
		}
		if (buckets > config.maxBuckets / count)
		{
			throw std::length_error("Spatial index bucket limit exceeded");
		}
		buckets *= count;
		auto& edges = m_Edges[axis];
		edges.reserve(static_cast<std::size_t>(count)+1);
		for (std::size_t i=0; i<=count; ++i)
		{
			const double edge = i == count ? high : low+(high-low)*(static_cast<double>(i)/count);
			if (!edges.empty() && edge <= edges.back())
			{
				throw std::invalid_argument("Spatial index bin edges are not distinguishable at this precision");
			}
			edges.push_back(edge);
		}
	}
	if (m_Primitives.size() > std::numeric_limits<Rank>::max() ||
		m_Primitives.size() > config.maxReferences)
	{
		throw std::length_error("Spatial index primitive limit exceeded");
	}
	m_Buckets.resize(buckets);
	m_Bounds.resize(m_Primitives.size());
	const double tolerance = structure.GetDrawingTolerance();
	if (!std::isfinite(tolerance) || tolerance < 0)
	{
		throw std::invalid_argument("Spatial index needs a finite nonnegative drawing tolerance");
	}
	auto append = [&](std::vector<Rank>& list, Rank rank) {
		if (m_References >= config.maxReferences)
		{
			throw std::length_error("Spatial index reference limit exceeded");
		}
		list.push_back(rank);
		++m_References;
	};
	// Iterating the priority-sorted source preserves its exact tie order in
	// every bucket. Unknown bounds are stored once, not once per bucket.
	for (std::size_t i=0; i<m_Primitives.size(); ++i)
	{
		auto* primitive = m_Primitives[i];
		const Rank rank = static_cast<Rank>(i);
		auto& bounds = m_Bounds[i];
		bounds.known = primitive->GetConservativeBoundBox(bounds.values.data());
		if (bounds.known)
		{
			for (int axis = 0; axis < 3; ++axis)
			{
				double& low = bounds.values[2 * axis];
				double& high = bounds.values[2 * axis + 1];
				low = std::nextafter(low - tolerance, -std::numeric_limits<double>::infinity());
				high = std::nextafter(high + tolerance, std::numeric_limits<double>::infinity());
				bounds.known = bounds.known && std::isfinite(low) && std::isfinite(high);
			}
		}
		const int containment = primitive->IsInsideBox(config.bounds.data());
		if (!Intersects(rank, config.bounds.data()))
		{
			continue;
		}
		if (containment == 0 && !bounds.known)
		{
			append(m_Global, rank);
			continue;
		}
		for (unsigned int z=0; z<config.bins[2]; ++z)
		{
			for (unsigned int y=0; y<config.bins[1]; ++y)
			{
				for (unsigned int x=0; x<config.bins[0]; ++x)
				{
					const double box[6] = {m_Edges[0][x], m_Edges[0][x+1],
						m_Edges[1][y], m_Edges[1][y+1], m_Edges[2][z], m_Edges[2][z+1]};
					if (Intersects(rank, box))
					{
						append(m_Buckets[Offset(x, y, z)], rank);
					}
				}
			}
		}
	}
	if (config.prepareCandidateViews)
	{
		m_CandidateBuckets.resize(m_Buckets.size());
		for (std::size_t i = 0; i < m_Buckets.size(); ++i)
		{
			const auto& local = m_Buckets[i];
			const std::size_t count = local.size() + m_Global.size();
			if (count > config.maxReferences - m_References)
			{
				throw std::length_error("Spatial index prepared-candidate limit exceeded");
			}
			std::vector<Rank> merged;
			merged.reserve(count);
			std::merge(local.begin(), local.end(), m_Global.begin(), m_Global.end(),
				std::back_inserter(merged));
			auto& candidates = m_CandidateBuckets[i];
			candidates.reserve(count);
			for (Rank rank : merged)
			{
				candidates.push_back(m_Primitives[rank]);
			}
			m_References += count;
		}
	}
}

std::size_t CSPrimitiveIndex::Offset(unsigned int x, unsigned int y, unsigned int z) const
{
	return (static_cast<std::size_t>(z)*m_Config.bins[1]+y)*m_Config.bins[0]+x;
}

bool CSPrimitiveIndex::Intersects(Rank rank, const double* bounds) const
{
	const auto& box = m_Bounds[rank];
	if (!box.known)
	{
		return m_Primitives[rank]->IsInsideBox(bounds) >= 0;
	}
	for (int axis = 0; axis < 3; ++axis)
	{
		const double low = std::min(bounds[2 * axis], bounds[2 * axis + 1]);
		const double high = std::max(bounds[2 * axis], bounds[2 * axis + 1]);
		if (high < box.values[2 * axis] || low > box.values[2 * axis + 1])
		{
			return false;
		}
	}
	return true;
}

void CSPrimitiveIndex::Filter(const double* bounds, const std::vector<Rank>& ranks,
	std::vector<CSPrimitives*>& result) const
{
	for (Rank rank : ranks)
	{
		if (Intersects(rank, bounds))
		{
			result.push_back(m_Primitives[rank]);
		}
	}
}

std::vector<CSPrimitives*> CSPrimitiveIndex::Query(const double* bounds) const
{
	std::vector<CSPrimitives*> result;
	Query(bounds, result);
	return result;
}

bool CSPrimitiveIndex::BucketRange(const double* bounds, std::array<unsigned int, 3>& first,
	std::array<unsigned int, 3>& last) const
{
	if (!bounds)
	{
		throw std::invalid_argument("Spatial index needs a query box");
	}
	first.fill(0);
	last.fill(0);
	bool outside = false;
	for (unsigned int axis=0; axis<3; ++axis)
	{
		if (!std::isfinite(bounds[2*axis]) || !std::isfinite(bounds[2*axis+1]))
		{
			throw std::invalid_argument("Spatial index query bounds must be finite");
		}
		const double low = std::min(bounds[2*axis], bounds[2*axis+1]);
		const double high = std::max(bounds[2*axis], bounds[2*axis+1]);
		if (low < m_Config.bounds[2*axis] || high > m_Config.bounds[2*axis+1])
		{
			outside = true;
			continue;
		}
		const auto& edges = m_Edges[axis];
		// Both closed buckets adjacent to an exact boundary must participate.
		const auto begin = std::lower_bound(edges.begin(), edges.end(), low)-edges.begin();
		const auto end = std::upper_bound(edges.begin(), edges.end(), high)-edges.begin();
		first[axis] = static_cast<unsigned int>(begin == 0 ? 0 : begin-1);
		last[axis] = static_cast<unsigned int>(std::min<std::size_t>(end-1, m_Config.bins[axis]-1));
	}
	return !outside;
}

const std::vector<CSPrimitives*>& CSPrimitiveIndex::GetCandidates(const double* bounds,
	std::vector<CSPrimitives*>& scratch) const
{
	std::array<unsigned int, 3> first{}, last{};
	if (BucketRange(bounds, first, last) && first == last && !m_CandidateBuckets.empty())
	{
		return m_CandidateBuckets[Offset(first[0], first[1], first[2])];
	}
	Query(bounds, scratch);
	return scratch;
}

void CSPrimitiveIndex::Query(const double* bounds, std::vector<CSPrimitives*>& result) const
{
	std::array<unsigned int, 3> first{}, last{};
	const bool outside = !BucketRange(bounds, first, last);
	result.clear();
	std::size_t touched = 1;
	for (unsigned int axis=0; axis<3; ++axis)
	{
		touched *= last[axis]-first[axis]+1;
	}
	if (outside || touched >= (m_Buckets.size()+1)/2)
	{
		// Unbounded and wide queries reuse the sorted snapshot, not CSXCAD's
		// collection/sort path. This also covers geometry outside the index.
		for (std::size_t rank = 0; rank < m_Primitives.size(); ++rank)
		{
			if (Intersects(static_cast<Rank>(rank), bounds))
			{
				result.push_back(m_Primitives[rank]);
			}
		}
		return;
	}
	if (touched == 1 && m_Global.empty())
	{
		Filter(bounds, m_Buckets[Offset(first[0], first[1], first[2])], result);
		return;
	}
	struct Head
	{
		const std::vector<Rank>* list;
		std::size_t position;
		Rank Value() const
		{
			return (*list)[position];
		}
	};
	std::vector<Head> heads;
	heads.reserve(touched+1);
	if (!m_Global.empty())
	{
		heads.push_back({&m_Global, 0});
	}
	for (unsigned int z=first[2]; z<=last[2]; ++z)
	{
		for (unsigned int y=first[1]; y<=last[1]; ++y)
		{
			for (unsigned int x=first[0]; x<=last[0]; ++x)
			{
				const auto& bucket = m_Buckets[Offset(x, y, z)];
				if (!bucket.empty())
				{
					heads.push_back({&bucket, 0});
				}
			}
		}
	}
	if (heads.size() == 1)
	{
		Filter(bounds, *heads.front().list, result);
		return;
	}
	// Merge already sorted rank streams. No query-time geometry sort or
	// mutable shared scratch; duplicate figures from adjacent buckets vanish.
	auto later = [](const Head& a, const Head& b) {return a.Value() > b.Value();};
	std::make_heap(heads.begin(), heads.end(), later);
	Rank previous = std::numeric_limits<Rank>::max();
	while (!heads.empty())
	{
		std::pop_heap(heads.begin(), heads.end(), later);
		Head head = heads.back();
		heads.pop_back();
		const Rank rank = head.Value();
		if (rank != previous && Intersects(rank, bounds))
		{
			result.push_back(m_Primitives[rank]);
		}
		previous = rank;
		if (++head.position < head.list->size())
		{
			heads.push_back(head);
			std::push_heap(heads.begin(), heads.end(), later);
		}
	}
}
