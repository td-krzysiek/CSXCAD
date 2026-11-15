// SPDX-License-Identifier: LGPL-3.0-or-later
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include <cstdlib>
#include <iostream>
#include <new>
#include <vector>

namespace
{
bool tracking = false;
std::size_t allocations = 0;
}

void* operator new(std::size_t bytes)
{
	if (tracking) {
		++allocations;
	}
	void* pointer = std::malloc(bytes ? bytes : 1);
	if (!pointer) {
		throw std::bad_alloc();
	}
	return pointer;
}
void operator delete(void* pointer) noexcept
{
	std::free(pointer);
}
void operator delete(void* pointer, std::size_t) noexcept
{
	std::free(pointer);
}

int main()
{
	ContinuousStructure geometry;
	auto* low = new CSPropMaterial(geometry.GetParameterSet());
	auto* high = new CSPropMaterial(geometry.GetParameterSet());
	geometry.AddProperty(low);
	geometry.AddProperty(high);
	auto makeBox = [&](CSProperties* property, int priority) {
		auto* box = new CSPrimBox(geometry.GetParameterSet(), property);
		box->SetPriority(priority);
		for (int i = 0; i < 6; ++i) {
			box->SetCoord(i, i % 2 ? 1.0 : -1.0);
		}
		return box;
	};
	auto* lowBox = makeBox(low, 1);
	auto* highBox = makeBox(high, 2);
	if (!geometry.Update().empty()) {
		return 2;
	}
	const std::vector<CSPrimitives*> primitives = {highBox, lowBox};
	const std::vector<CSPrimitives*> empty;
	const double inside[3] = {0, 0, 0}, outside[3] = {3, 3, 3};
	CSPrimitives* found = nullptr;
	// Warm native shape-coordinate caches before measuring list ownership.
	geometry.GetPropertyByCoordPriority(inside, primitives, false);
	allocations = 0;
	tracking = true;
	auto* first = geometry.GetPropertyByCoordPriority(inside, primitives, false, &found);
	auto* miss = geometry.GetPropertyByCoordPriority(outside, primitives, false);
	auto* none = geometry.GetPropertyByCoordPriority(inside, empty, false);
	tracking = false;
	if (allocations != 0 || first != high || found != highBox || miss || none ||
		highBox->GetPrimitiveUsed() || lowBox->GetPrimitiveUsed() || primitives[0] != highBox) {
		std::cerr << "Borrowed lookup failed: allocations=" << allocations << '\n';
		return 1;
	}
	geometry.GetPropertyByCoordPriority(inside, primitives, true);
	if (!highBox->GetPrimitiveUsed() || lowBox->GetPrimitiveUsed()) {
		return 1;
	}
	return 0;
}
