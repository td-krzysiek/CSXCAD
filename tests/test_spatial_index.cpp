// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "ContinuousStructure.h"
#include "CSPrimitiveIndex.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"

#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>

namespace
{
void Require(bool condition, const char* message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

template<class Error, class Action>
void Throws(Action action)
{
	try
	{
		action();
	}
	catch (const Error&)
	{
		return;
	}
	throw std::runtime_error("expected configuration/query error was not raised");
}

CSPrimBox* Box(ContinuousStructure& csx, CSProperties* property, int priority,
	const std::array<double, 6>& bounds)
{
	auto* box = new CSPrimBox(csx.GetParameterSet(), property);
	box->SetPriority(priority);
	for (int i=0; i<6; ++i)
	{
		box->SetCoord(i, bounds[i]);
	}
	return box;
}

class UnknownBoundsBox : public CSPrimBox
{
public:
	using CSPrimBox::CSPrimBox;
	int IsInsideBox(const double*) override
	{
		return 0;
	}
};

void CheckGeometry(CoordinateSystem coordinates)
{
	ContinuousStructure csx;
	csx.SetCoordInputType(coordinates);
	auto* material = new CSPropMaterial(csx.GetParameterSet());
	auto* metal = new CSPropMetal(csx.GetParameterSet());
	csx.AddProperty(material);
	csx.AddProperty(metal);
	// Equal priorities deliberately test the existing, type-specific sort order.
	for (int i=0; i<80; ++i)
	{
		const double x = (i%8)*.5;
		const double y = ((i/8)%5)*.5;
		const double z = (i/40)*.5;
		Box(csx, i%3 ? static_cast<CSProperties*>(metal) : material,
			i%7, {{x, x+.5, y, y+.5, z, z+.5}});
	}
	// High-priority dielectric void, boundary contact, outside and transformed objects.
	Box(csx, material, 100, {{0, 1, 0, 1, 0, 1}});
	Box(csx, metal, 5, {{4, 4, 2, 3, 1, 2}});
	Box(csx, metal, 9, {{20, 21, 20, 21, 20, 21}});
	auto* shifted = Box(csx, metal, 3, {{0, .5, 0, .5, 0, .5}});
	const double translation[3] = {1, 1, 1};
	shifted->GetTransform()->Translate(translation);
	auto* unknown = new UnknownBoundsBox(csx.GetParameterSet(), material);
	unknown->SetPriority(-10);
	auto* cylinder = new CSPrimCylinder(csx.GetParameterSet(), metal);
	const double ends[6] = {2, 2, 1, 1, 0, 3};
	for (int i=0; i<6; ++i)
	{
		cylinder->SetCoord(i, ends[i]);
	}
	cylinder->SetRadius(.2);
	Require(csx.Update().empty(), "fixture must update before indexing");
	std::vector<std::array<double, 6>> queries = {
		{{0, 0, 0, 0, 0, 0}}, {{4, 4, 3, 3, 3, 3}},
		{{.5, .5, .5, .5, .5, .5}}, {{1, 0, 1, 0, 1, 0}},
		{{-1, 0, -1, 0, -1, 0}}, {{0, 4, 0, 3, 0, 3}},
		{{20, 21, 20, 21, 20, 21}}, {{1, 1.1, 1, 1.1, 1, 1.1}}};
	std::mt19937 random(42);
	std::uniform_real_distribution<double> position(-.5, 4.5);
	for (int i=0; i<1000; ++i)
	{
		std::array<double, 6> query;
		for (int axis=0; axis<3; ++axis)
		{
			query[2*axis] = position(random);
			query[2*axis+1] = query[2*axis] + (i%11)*.07;
		}
		queries.push_back(query);
	}
	const CSProperties::PropertyType masks[] = {CSProperties::ANY,
		CSProperties::MATERIAL, CSProperties::METAL,
		static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL)};
	for (const auto bins : {std::array<unsigned int, 3>{{8, 6, 1}},
		std::array<unsigned int, 3>{{8, 6, 6}}, std::array<unsigned int, 3>{{1, 1, 1}}})
	{
		for (auto mask : masks)
		{
			CSPrimitiveIndex::Config config;
			config.bounds = {{0, 4, 0, 3, 0, 3}};
			config.bins = bins;
			config.type = mask;
			CSPrimitiveIndex index(csx, config);
			Require(index.GetBucketCount() == bins[0]*bins[1]*bins[2], "wrong bucket count");
			std::vector<std::vector<CSPrimitives*>> expected;
			for (const auto& query : queries)
			{
				expected.push_back(csx.GetPrimitivesByBoundBox(query.data(), true, mask));
				Require(index.Query(query.data()) == expected.back(),
					"indexed candidates/order differ from unchanged legacy API");
			}
			for (auto* primitive : csx.GetAllPrimitives())
			{
				Require(!primitive->GetPrimitiveUsed(), "bounds queries must not mark Used");
			}
			std::atomic<bool> same{true};
			std::vector<std::thread> workers;
			for (int thread=0; thread<4; ++thread)
			{
				workers.emplace_back([&]() {
					std::vector<CSPrimitives*> result;
					for (std::size_t i=0; i<queries.size(); ++i)
					{
						index.Query(queries[i].data(), result);
						if (result != expected[i])
						{
							same = false;
						}
					}
				});
			}
			for (auto& worker : workers)
			{
				worker.join();
			}
			Require(same, "parallel queries differ");
			// Exact point lookup, including a higher-priority material over metal.
			const double point[3] = {.25, .25, .25};
			const auto candidates = index.Query(queries[5].data());
			const auto legacy = csx.GetPrimitivesByBoundBox(queries[5].data(), true, mask);
			Require(csx.GetPropertyByCoordPriority(point, candidates, false) ==
				csx.GetPropertyByCoordPriority(point, legacy, false), "point winner differs");
		}
	}
	// Explicit lifetime: discard the index before membership/priority mutations.
	csx.DeletePrimitive(shifted);
	Box(csx, material, 200, {{0, 1, 0, 1, 0, 1}});
	Require(csx.Update().empty(), "changed fixture must update");
	const CSPrimitiveIndex rebuilt(csx, CSPrimitiveIndex::Config{});
	Require(rebuilt.Query(queries[0].data()) ==
		csx.GetPrimitivesByBoundBox(queries[0].data(), true), "rebuild is stale");
}

void CheckErrors()
{
	ContinuousStructure csx;
	auto* material = new CSPropMaterial(csx.GetParameterSet());
	csx.AddProperty(material);
	Box(csx, material, 1, {{0, 1, 0, 1, 0, 1}});
	Require(csx.Update().empty(), "error fixture must update");
	CSPrimitiveIndex::Config config;
	config.bins[0] = 0;
	Throws<std::invalid_argument>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	config.bounds[0] = std::numeric_limits<double>::quiet_NaN();
	Throws<std::invalid_argument>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	config.maxBuckets = 1;
	Throws<std::length_error>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	config.maxReferences = 1;
	Throws<std::length_error>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	config.bounds[0] = 1e16;
	config.bounds[1] = std::nextafter(1e16, std::numeric_limits<double>::infinity());
	Throws<std::invalid_argument>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	config.bins = {{std::numeric_limits<unsigned int>::max(), 2, 2}};
	Throws<std::length_error>([&]() {CSPrimitiveIndex invalid(csx, config);});
	config = {};
	CSPrimitiveIndex index(csx, config);
	Throws<std::invalid_argument>([&]() {index.Query(nullptr);});
	double invalid[6] = {0, 1, 0, 1, 0, std::numeric_limits<double>::infinity()};
	Throws<std::invalid_argument>([&]() {index.Query(invalid);});
	ContinuousStructure empty;
	const CSPrimitiveIndex emptyIndex(empty, config);
	Require(emptyIndex.Query(config.bounds.data()).empty(), "empty index must stay empty");
}
}

int main()
{
	try
	{
		CheckGeometry(CARTESIAN);
		CheckGeometry(CYLINDRICAL);
		CheckErrors();
	}
	catch (const std::exception& error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	std::cout << "Spatial index tests passed\n";
}
