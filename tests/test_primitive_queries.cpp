// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"

#include <iostream>
#include <atomic>
#include <stdexcept>
#include <thread>
#include <vector>

namespace
{
void Require(bool condition, const char* message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

CSPrimBox* AddBox(ContinuousStructure& csx, CSProperties* property,
	int priority, double offset=0)
{
	auto* box = new CSPrimBox(csx.GetParameterSet(), property);
	box->SetPriority(priority);
	for (int axis=0; axis<3; ++axis)
	{
		box->SetCoord(2*axis, offset);
		box->SetCoord(2*axis+1, offset+1);
	}
	return box;
}

// The broad phase must retain shapes whose bounding test is inconclusive.
class UnknownBoundsBox : public CSPrimBox
{
public:
	using CSPrimBox::CSPrimBox;
	int IsInsideBox(const double*) override
	{
		return 0;
	}
};

void CheckQueries()
{
	ContinuousStructure csx;
	auto* material = new CSPropMaterial(csx.GetParameterSet());
	auto* metal = new CSPropMetal(csx.GetParameterSet());
	csx.AddProperty(material);
	csx.AddProperty(metal);
	auto* lower = AddBox(csx, metal, 2);
	auto* upper = AddBox(csx, material, 7);
	auto* tied = AddBox(csx, metal, 7);
	AddBox(csx, metal, 1, 10);
	auto* unknown = new UnknownBoundsBox(csx.GetParameterSet(), material);
	unknown->SetPriority(-10);
	for (int axis=0; axis<3; ++axis)
	{
		unknown->SetCoord(2*axis, 30.0);
		unknown->SetCoord(2*axis+1, 31.0);
	}
	const double translation[3] = {3, 0, 0};
	auto* shifted = AddBox(csx, metal, 4);
	shifted->GetTransform()->Translate(translation);
	auto* cylinder = new CSPrimCylinder(csx.GetParameterSet(), metal);
	cylinder->SetPriority(5);
	cylinder->SetCoord(0, -2.0);
	cylinder->SetCoord(1, -2.0);
	cylinder->SetCoord(2, 0.0);
	cylinder->SetCoord(3, 0.0);
	cylinder->SetCoord(4, 0.0);
	cylinder->SetCoord(5, 1.0);
	cylinder->SetRadius(.5);
	Require(csx.Update().empty(), "fixture geometry must update");

	const CSProperties::PropertyType masks[] = {CSProperties::ANY,
		CSProperties::MATERIAL, CSProperties::METAL,
		static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL)};
	const double bounds[][6] = {{-.1, .5, -.1, .5, -.1, .5},
		{1, 1, 0, 1, 0, 1}, {3.1, 3.2, .1, .2, .1, .2},
		{-2.1, -1.9, -.1, .1, .2, .8}, {20, 21, 20, 21, 20, 21}};
	for (auto mask : masks)
	{
		for (bool sorted : {false, true})
		{
			const auto all = csx.GetAllPrimitives(sorted, mask);
			const auto original = all;
			for (const auto& bound : bounds)
			{
				Require(csx.GetPrimitivesByBoundBox(bound, all) ==
					csx.GetPrimitivesByBoundBox(bound, sorted, mask),
					"reused query must preserve legacy candidates and tie order");
			}
			Require(all == original, "queries must not modify the borrowed list");
		}
	}
	const auto all = csx.GetAllPrimitives(true);
	const auto local = csx.GetPrimitivesByBoundBox(bounds[0], all);
	Require(local.size() == 5 && local.back() == unknown,
		"unknown bounds must remain candidates; distant boxes must be excluded");
	Require(csx.GetPrimitivesByBoundBox(bounds[2], all) ==
		std::vector<CSPrimitives*>({shifted, unknown}),
		"transformed box must appear at its transformed location");
	Require(csx.GetPrimitivesByBoundBox(bounds[3], all) ==
		std::vector<CSPrimitives*>({cylinder, shifted, unknown}),
		"cylinder bounds and inconclusive transformed bounds must be respected");
	Require(csx.GetPrimitivesByBoundBox(bounds[4], all) ==
		std::vector<CSPrimitives*>({shifted, unknown}),
		"transformed and unknown bounds must remain in a distant query");
	Require(csx.GetPrimitivesByBoundBox(bounds[0], csx.GetAllPrimitives(true, CSProperties::METAL)) ==
		std::vector<CSPrimitives*>({tied, shifted, lower}),
		"a property-filtered snapshot must retain only the requested property type");
	const std::vector<CSPrimitives*> empty;
	Require(csx.GetPrimitivesByBoundBox(bounds[0], empty).empty(), "empty list must stay empty");

	// Also lock the upstream const-reference point-query API and its observable
	// priority, used-flag and no-hit semantics. No copy may be reintroduced.
	using Lookup = CSProperties* (ContinuousStructure::*)(const double*,
		const std::vector<CSPrimitives*>&, bool, CSPrimitives**);
	Lookup lookup = static_cast<Lookup>(&ContinuousStructure::GetPropertyByCoordPriority);
	const double point[3] = {.5, .5, .5};
	CSPrimitives* found = nullptr;
	auto* winner = (csx.*lookup)(point, local, true, &found);
	Require(found == local.front() && winner == found->GetProperty(),
		"highest-priority candidate must win with unchanged tie order");
	Require(found->GetPrimitiveUsed(), "winner must be marked used");
	const double transformedPoint[3] = {3.2, .2, .2};
	(csx.*lookup)(transformedPoint, csx.GetPrimitivesByBoundBox(bounds[2], all), false, &found);
	Require(found == shifted, "point containment must evaluate the transformation");
	Require(!lower->GetPrimitiveUsed(), "lower-priority primitive must remain unused");
	upper->SetPrimitiveUsed(false);
	tied->SetPrimitiveUsed(false);
	(csx.*lookup)(point, local, false, &found);
	Require(!upper->GetPrimitiveUsed() && !tied->GetPrimitiveUsed(),
		"read-only query must not mark primitives used");
	const double outside[3] = {100, 100, 100};
	Require((csx.*lookup)(outside, local, false, nullptr) == nullptr,
		"point query outside every candidate must return no property");
	csx.SetDrawingTolerance(1e-4);
	const double nearBoundary[3] = {1+5e-5, .5, .5};
	Require((csx.*lookup)(nearBoundary, local, false, nullptr) != nullptr,
		"point lookup must preserve drawing tolerance");

	// Independent output vectors can share one immutable priority list. The
	// existing atomic used flags are the only intentional shared writes here.
	std::atomic<bool> equivalent{true};
	std::vector<std::thread> workers;
	for (int thread=0; thread<4; ++thread)
	{
		workers.emplace_back([&]() {
			for (int query=0; query<100; ++query)
			{
				const auto candidates = csx.GetPrimitivesByBoundBox(bounds[0], all);
				CSPrimitives* selected = nullptr;
				(csx.*lookup)(point, candidates, true, &selected);
				if (candidates != local || selected != local.front())
				{
					equivalent.store(false);
				}
			}
		});
	}
	for (auto& worker : workers)
	{
		worker.join();
	}
	Require(equivalent.load(), "parallel queries must preserve candidates and winners");

	// Membership and priority changes require a fresh caller snapshot; there is
	// no hidden cache that could retain deleted pointers or old priorities.
	csx.DeletePrimitive(shifted);
	auto* newBox = AddBox(csx, material, 20);
	lower->SetPriority(30);
	Require(csx.Update().empty(), "changed geometry must update");
	const auto fresh = csx.GetAllPrimitives(true);
	Require(fresh.front() == lower, "refreshed list must reflect changed priorities");
	const auto refreshed = csx.GetPrimitivesByBoundBox(bounds[0], fresh);
	Require(refreshed == csx.GetPrimitivesByBoundBox(bounds[0], true),
		"refreshed membership must match a fresh legacy query");
	Require(newBox->GetProperty() == material, "new primitive ownership must be unchanged");
}
}

int main()
{
	try
	{
		CheckQueries();
	}
	catch (const std::exception& error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	std::cout << "Primitive query tests passed\n";
	return 0;
}
