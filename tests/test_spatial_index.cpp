// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "ContinuousStructure.h"
#include "CSPrimitiveIndex.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPrimLinPoly.h"
#include "CSPrimRotPoly.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
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
	config.prepareCandidateViews = true;
	config.bins = {{1, 1, 1}};
	config.maxReferences = 1;
	Throws<std::length_error>([&]() {CSPrimitiveIndex invalid(csx, config);});
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

void PolygonCoords(CSPrimPolygon& polygon, double x, double y)
{
	// Concave contour: the missing upper-right corner must stay outside.
	for (const auto& point : {std::array<double, 2>{{0, 0}}, {{1, 0}}, {{1, .4}},
		{{.4, .4}}, {{.4, 1}}, {{0, 1}}})
	{
		polygon.AddCoord(x + point[0]);
		polygon.AddCoord(y + point[1]);
	}
}

void CheckPolygonBounds()
{
	ContinuousStructure csx;
	csx.SetDrawingTolerance(.01);
	auto* metal = new CSPropMetal(csx.GetParameterSet());
	auto* dielectric = new CSPropMaterial(csx.GetParameterSet());
	csx.AddProperty(metal);
	csx.AddProperty(dielectric);
	Box(csx, dielectric, -1, {{-1, 17, -1, 17, -1, 2}});
	std::vector<CSPrimLinPoly*> polygons;
	for (int i = 0; i < 64; ++i)
	{
		auto* polygon = new CSPrimLinPoly(csx.GetParameterSet(), metal);
		polygon->SetNormDir(2);
		polygon->SetElevation(.25);
		polygon->SetLength(.1);
		polygon->SetPriority(5);
		PolygonCoords(*polygon, 2 * (i % 8), 2 * (i / 8));
		polygons.push_back(polygon);
	}
	// Higher-priority dielectric and tied-priority metal exercise winner ordering.
	Box(csx, dielectric, 10, {{0, .2, 0, .2, .2, .4}});
	Box(csx, metal, 5, {{0, .3, 0, .3, .2, .4}});
	auto* shifted = new CSPrimLinPoly(csx.GetParameterSet(), metal);
	PolygonCoords(*shifted, 0, 0);
	shifted->SetNormDir(2);
	shifted->SetLength(.2);
	const double move[3] = {5, 5, .8};
	shifted->GetTransform()->Translate(move);
	Require(csx.Update().empty(), "polygon fixture must update");
	double bounds[6];
	Require(!polygons[0]->GetBoundBox(bounds), "legacy polygon bounds contract changed");
	Require(polygons[0]->IsInsideBox(bounds) == 0, "legacy unknown-box result changed");
	Require(polygons[0]->GetConservativeBoundBox(bounds), "polygon lacks conservative bounds");
	Require(!shifted->GetConservativeBoundBox(bounds), "transformed polygon must remain conservative");
	CSPrimitiveIndex::Config config;
	config.bounds = {{-1, 17, -1, 17, -1, 2}};
	config.prepareCandidateViews = true;
	const auto all = csx.GetAllPrimitives(true);
	for (const auto bins : {std::array<unsigned int, 3>{{8, 8, 1}}, {{8, 8, 4}}})
	{
		config.bins = bins;
		CSPrimitiveIndex index(csx, config);
		Require(index.GetGlobalCount() == 1, "untransformed polygons escaped the partition");
		const double query[6] = {.1, .2, .1, .2, .26, .30};
		Require(index.Query(query).size() < 8, "polygon candidates were not reduced");
		std::vector<CSPrimitives*> scratch;
		const auto& prepared = index.GetCandidates(query, scratch);
		Require(&prepared != &scratch, "single-tile candidate list was not borrowed");
		Require(&prepared == &index.GetCandidates(query, scratch), "prepared view address changed");
		const double crossing[6] = {1.2, 1.3, .1, .2, .26, .30};
		const auto& merged = index.GetCandidates(crossing, scratch);
		Require(&merged == &scratch, "crossing query did not use caller scratch");
		Require(merged == index.Query(crossing), "crossing candidate list changed");
		std::mt19937 random(37);
		std::uniform_real_distribution<double> position(-.5, 15.5);
		std::vector<std::array<double, 3>> points;
		for (int i = 0; i < 2000; ++i)
		{
			points.push_back({{position(random), position(random), .245 + (i % 12) * .01}});
		}
		for (const auto point : {std::array<double, 3>{{-.005, .1, .3}}, {{.1, -.005, .3}},
			{{.1, .1, .245}}, {{.1, .1, .355}}, {{.405, .405, .3}}, {{.8, .8, .3}},
			{{5.1, 5.1, .9}}})
		{
			points.push_back(point);
		}
		for (const auto& point : points)
		{
			// A point query specifically checks the drawing-tolerance expansion.
			const double box[6] = {point[0], point[0], point[1], point[1], point[2], point[2]};
			const auto found = index.Query(box);
			for (auto* polygon : polygons)
			{
				if (polygon->IsInside(point.data(), csx.GetDrawingTolerance()))
				{
					Require(std::find(found.begin(), found.end(), polygon) != found.end(),
						"index dropped a true polygon hit, including tolerance");
				}
			}
			// Pad the query for legacy boxes, whose old box contract is unchanged.
			const double padded[6] = {point[0] - .02, point[0] + .02, point[1] - .02,
				point[1] + .02, point[2] - .02, point[2] + .02};
			Require(csx.GetPropertyByCoordPriority(point.data(), index.Query(padded), false) ==
				csx.GetPropertyByCoordPriority(point.data(), all, false), "polygon material winner differs");
			Require(csx.GetPropertyByCoordPriority(point.data(), index.GetCandidates(padded, scratch), false) ==
				csx.GetPropertyByCoordPriority(point.data(), all, false), "prepared-view material winner differs");
		}
		std::atomic<bool> parallelSame{true};
		std::vector<std::thread> workers;
		for (int thread = 0; thread < 4; ++thread)
		{
			workers.emplace_back([&]() {
				std::vector<CSPrimitives*> local;
				for (const auto& point : points)
				{
					const double padded[6] = {point[0] - .02, point[0] + .02, point[1] - .02,
						point[1] + .02, point[2] - .02, point[2] + .02};
					if (csx.GetPropertyByCoordPriority(point.data(), index.GetCandidates(padded, local), false) !=
						csx.GetPropertyByCoordPriority(point.data(), all, false))
					{
						parallelSame = false;
					}
				}
			});
		}
		for (auto& worker : workers)
		{
			worker.join();
		}
		Require(parallelSame, "parallel borrowed/scratch candidate queries differ");
	}
	// Normal-axis orientation and signed extrusion are part of the bounds contract.
	for (int axis = 0; axis < 3; ++axis)
	{
		for (double length : {-.5, .5})
		{
			polygons[0]->SetNormDir(axis);
			polygons[0]->SetLength(length);
			Require(csx.Update().empty(), "signed polygon fixture must update");
			Require(polygons[0]->GetConservativeBoundBox(bounds), "signed extrusion bounds missing");
			Require(bounds[2 * axis] == std::min(.25, .25 + length) &&
				bounds[2 * axis + 1] == std::max(.25, .25 + length), "wrong extrusion interval");
		}
	}
	CSPrimRotPoly rotated(csx.GetParameterSet(), nullptr);
	Require(!rotated.GetConservativeBoundBox(bounds), "RotPoly must not use planar bounds");
	polygons[0]->SetCoordInputType(CYLINDRICAL, false);
	Require(!polygons[0]->GetConservativeBoundBox(bounds), "non-Cartesian polygon bounds accepted");
}

void CheckPreparedGrid()
{
	ContinuousStructure csx;
	auto* material = new CSPropMaterial(csx.GetParameterSet());
	csx.AddProperty(material);
	Box(csx, material, 2, {{0, 1, 0, 1, 0, 1}});
	Box(csx, material, 4, {{.5, 2, .5, 2, .5, 2}});
	Require(csx.Update().empty(), "prepared fixture must update");
	CSPrimitiveIndex::Config config;
	config.bounds = {{-1, 3, -1, 3, -1, 3}};
	config.bins = {{4, 4, 4}};
	CSPrimitiveIndex index(csx, config);
	CSPrimitiveIndex::QueryRanges ranges;
	for (auto& axis : ranges)
	{
		axis = {{{-.25, .25}}, {{.75, 1.25}}, {{2, 1}}, {{-1, 3}}, {{-.25, .25}}};
	}
	std::unique_ptr<CSPrimitiveIndex::PreparedGrid> grid;
	{
		CSPrimitiveIndex temporary(csx, config);
		grid = std::make_unique<CSPrimitiveIndex::PreparedGrid>(temporary, ranges);
	}
	Require(grid->GetQueryCount() == 125 && grid->GetMappingBytes() == 500,
		"prepared grid dimensions or mapping bytes differ");
	Require(grid->GetListCount() < grid->GetQueryCount(), "prepared lists were not shared");
	for (std::size_t z = 0; z < 5; ++z)
	{
		for (std::size_t y = 0; y < 5; ++y)
		{
			for (std::size_t x = 0; x < 5; ++x)
			{
				const double bounds[6] = {ranges[0][x][0], ranges[0][x][1],
					ranges[1][y][0], ranges[1][y][1], ranges[2][z][0], ranges[2][z][1]};
				Require(grid->GetCandidates(x, y, z) == index.Query(bounds),
					"prepared ordered candidates differ from exact query");
			}
		}
	}
	Require(&grid->GetCandidates(0, 0, 0) == &grid->GetCandidates(4, 4, 4),
		"identical query lists do not share a vector");
	std::atomic<bool> same{true};
	std::vector<std::thread> readers;
	for (int i = 0; i < 4; ++i)
	{
		readers.emplace_back([&]() {
			for (int j = 0; j < 1000; ++j)
			{
				if (&grid->GetCandidates(0, 0, 0) != &grid->GetCandidates(4, 4, 4))
				{
					same = false;
				}
			}
		});
	}
	for (auto& reader : readers)
	{
		reader.join();
	}
	Require(same, "parallel prepared readers differ");
	Throws<std::out_of_range>([&]() {grid->GetCandidates(5, 0, 0);});
	Throws<std::length_error>([&]() {CSPrimitiveIndex::PreparedGrid tooMany(index, ranges, 124);});
	Throws<std::length_error>([&]() {CSPrimitiveIndex::PreparedGrid tooMany(index, ranges, 125, 0);});
	ranges[0][0][0] = std::numeric_limits<double>::infinity();
	Throws<std::invalid_argument>([&]() {CSPrimitiveIndex::PreparedGrid invalid(index, ranges);});
	ranges[0].clear();
	Throws<std::length_error>([&]() {CSPrimitiveIndex::PreparedGrid empty(index, ranges);});
}
}

int main()
{
	try
	{
		CheckGeometry(CARTESIAN);
		CheckGeometry(CYLINDRICAL);
		CheckErrors();
		CheckPolygonBounds();
		CheckPreparedGrid();
	}
	catch (const std::exception& error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	std::cout << "Spatial index tests passed\n";
}
