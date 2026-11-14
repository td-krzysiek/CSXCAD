// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#pragma once

#include "CSXCAD_Global.h"
#include "CSProperties.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

class CSPrimitives;
class ContinuousStructure;

//! Immutable spatial lookup for a fixed, updated geometry (non-owning).
/*!
Construct from ContinuousStructure after Update(). Keep the structure and
its primitives alive and unchanged until the index is destroyed. Rebuild after
any geometry, parameter, membership, coordinate-system or priority change.
Queries retain the original priority/tie order and conservative IsInsideBox()
semantics; exact point containment remains GetPropertyByCoordPriority()'s job.
 */
class CSXCAD_EXPORT CSPrimitiveIndex
{
public:
	struct Config
	{
		std::array<double, 6> bounds{{0, 1, 0, 1, 0, 1}};
		std::array<unsigned int, 3> bins{{16, 16, 1}};
		CSProperties::PropertyType type = CSProperties::ANY;
		std::size_t maxBuckets = 1048576;
		std::size_t maxReferences = 8388608;
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
	std::size_t GetBucketCount() const {return m_Buckets.size();}
	std::size_t GetReferenceCount() const {return m_References;}
	std::size_t GetPrimitiveCount() const {return m_Primitives.size();}

private:
	using Rank = std::uint32_t;
	std::size_t Offset(unsigned int x, unsigned int y, unsigned int z) const;
	void Filter(const double* bounds, const std::vector<Rank>& ranks,
		std::vector<CSPrimitives*>& result) const;

	Config m_Config;
	std::array<std::vector<double>, 3> m_Edges;
	std::vector<CSPrimitives*> m_Primitives;
	std::vector<std::vector<Rank>> m_Buckets;
	std::vector<Rank> m_Global;
	std::size_t m_References = 0;
};
