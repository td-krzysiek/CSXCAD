// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#pragma once

#include "CSXCAD_Global.h"
#include "CSProperties.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#define CSXCAD_HAS_PREPARED_PRIMITIVE_GRID 1

class CSPrimitives;
class ContinuousStructure;

//! Immutable spatial lookup for a fixed, updated geometry (non-owning).
/*!
Construct from ContinuousStructure after Update(). Keep the structure and
its primitives alive and unchanged until the index is destroyed. Rebuild after
any geometry, parameter, membership, coordinate-system or priority change.
Queries retain priority/tie order and conservatively cull supported polygons.
Exact point containment remains GetPropertyByCoordPriority()'s job.
 */
class CSXCAD_EXPORT CSPrimitiveIndex
{
public:
	class PreparedGrid;
	using QueryRanges = std::array<std::vector<std::array<double, 2>>, 3>;

	struct Config
	{
		std::array<double, 6> bounds{{0, 1, 0, 1, 0, 1}};
		std::array<unsigned int, 3> bins{{16, 16, 1}};
		CSProperties::PropertyType type = CSProperties::ANY;
		std::size_t maxBuckets = 1048576;
		std::size_t maxReferences = 8388608;
		bool prepareCandidateViews = false;
	};

	//! Collect and priority-sort the requested figures from a fixed, updated structure.
	CSPrimitiveIndex(ContinuousStructure& structure, const Config& config);

	CSPrimitiveIndex(CSPrimitiveIndex&&) noexcept = default;
	CSPrimitiveIndex& operator=(CSPrimitiveIndex&&) noexcept = default;
	CSPrimitiveIndex(const CSPrimitiveIndex&) = delete;
	CSPrimitiveIndex& operator=(const CSPrimitiveIndex&) = delete;

	//! Return priority-ordered candidates for an axis-aligned query box.
	std::vector<CSPrimitives*> Query(const double* bounds) const;
	//! Same query, reusing caller-owned output storage (one buffer per thread).
	void Query(const double* bounds, std::vector<CSPrimitives*>& result) const;
	//! Borrow a precomputed, possibly broader tile list; crossing queries use caller scratch.
	//! The returned reference lives until index destruction or the next write to scratch.
	const std::vector<CSPrimitives*>& GetCandidates(const double* bounds,
		std::vector<CSPrimitives*>& scratch) const;
	std::size_t GetBucketCount() const {return m_Buckets.size();}
	std::size_t GetReferenceCount() const {return m_References;}
	std::size_t GetPrimitiveCount() const {return m_Primitives.size();}
	std::size_t GetGlobalCount() const {return m_Global.size();}

private:
	using Rank = std::uint32_t;
	struct Bounds
	{
		std::array<double, 6> values{};
		bool known = false;
	};
	std::size_t Offset(unsigned int x, unsigned int y, unsigned int z) const;
	bool Intersects(Rank rank, const double* bounds) const;
	bool BucketRange(const double* bounds, std::array<unsigned int, 3>& first,
		std::array<unsigned int, 3>& last) const;
	void Filter(const double* bounds, const std::vector<Rank>& ranks,
		std::vector<CSPrimitives*>& result) const;

	Config m_Config;
	std::array<std::vector<double>, 3> m_Edges;
	std::vector<CSPrimitives*> m_Primitives;
	std::vector<Bounds> m_Bounds;
	std::vector<std::vector<Rank>> m_Buckets;
	std::vector<Rank> m_Global;
	std::vector<std::vector<CSPrimitives*>> m_CandidateBuckets;
	std::size_t m_References = 0;
};

//! Exact candidate lists prepared once for a Cartesian product of query ranges.
/*!
Identical priority-ordered lists share one immutable vector. Integer lookup
does no geometry query, merge, filtering or allocation. The source index need
not survive construction; its geometry must remain alive and unchanged.
 */
class CSXCAD_EXPORT CSPrimitiveIndex::PreparedGrid
{
public:
	PreparedGrid(const CSPrimitiveIndex& index, const QueryRanges& ranges,
		std::size_t maxQueries = 8388608, std::size_t maxReferences = 8388608);
	const std::vector<CSPrimitives*>& GetCandidates(std::size_t x,
		std::size_t y, std::size_t z) const;
	std::size_t GetQueryCount() const {return m_QueryLists.size();}
	std::size_t GetListCount() const {return m_Lists.size();}
	std::size_t GetReferenceCount() const {return m_References;}
	std::size_t GetMappingBytes() const {return m_QueryLists.size() * sizeof(std::uint32_t);}

private:
	std::array<std::size_t, 3> m_Counts{};
	std::vector<std::uint32_t> m_QueryLists;
	std::vector<std::vector<CSPrimitives*>> m_Lists;
	std::size_t m_References = 0;
};
