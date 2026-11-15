// SPDX-License-Identifier: LGPL-3.0-or-later
#include "ContinuousStructure.h"
#include "CSPrimLinPoly.h"
#include "CSPrimRotPoly.h"
#include "CSPropMetal.h"
#include "CSTransform.h"
#include <algorithm>
#include <array>
#include <iostream>

int main()
{
	ContinuousStructure geometry;
	auto* metal = new CSPropMetal(geometry.GetParameterSet());
	geometry.AddProperty(metal);
	auto* polygon = new CSPrimLinPoly(geometry.GetParameterSet(), metal);
	for (double coordinate : {0., 0., 1., 0., .2, 1.}) {
		polygon->AddCoord(coordinate);
	}
	double bounds[6];
	if (polygon->GetConservativeBoundBox(bounds)) {
		return 1;
	}
	for (int axis = 0; axis < 3; ++axis) {
		for (double length : {-.5, .5}) {
			polygon->SetNormDir(axis);
			polygon->SetElevation(.25);
			polygon->SetLength(length);
			if (!geometry.Update().empty() || !polygon->GetConservativeBoundBox(bounds)) {
				return 2;
			}
			if (bounds[2*axis] != std::min(.25, .25+length) ||
				bounds[2*axis+1] != std::max(.25, .25+length)) {
				return 3;
			}
			if (polygon->GetBoundBox(bounds) || polygon->IsInsideBox(bounds) != 0) {
				return 4; // Legacy contracts must not change.
			}
			polygon->GetConservativeBoundBox(bounds);
			const double tolerance = .01;
			for (int x = -2; x <= 22; ++x) {
				for (int y = -2; y <= 22; ++y) {
					for (int z = -12; z <= 22; ++z) {
						const double point[3] = {x*.05, y*.05, z*.05};
						if (!polygon->IsInside(point, tolerance)) {
							continue;
						}
						for (int coordinate = 0; coordinate < 3; ++coordinate) {
							if (point[coordinate] < bounds[2*coordinate]-tolerance ||
								point[coordinate] > bounds[2*coordinate+1]+tolerance) {
								return 5;
							}
						}
					}
				}
			}
		}
	}
	polygon->SetLength(.1);
	if (polygon->GetConservativeBoundBox(bounds)) {
		return 6; // Invalidated caches must never be used.
	}
	geometry.Update();
	const double translation[3] = {1, 2, 3};
	polygon->GetTransform()->Translate(translation);
	geometry.Update();
	if (polygon->GetConservativeBoundBox(bounds)) {
		return 7;
	}
	CSPrimRotPoly rotated(geometry.GetParameterSet(), nullptr);
	if (rotated.GetConservativeBoundBox(bounds)) {
		return 8;
	}
	polygon->SetCoordInputType(CYLINDRICAL, false);
	if (polygon->GetConservativeBoundBox(bounds) || polygon->GetConservativeBoundBox(nullptr)) {
		return 9;
	}
	std::cout << "Conservative polygon bounds: PASS\n";
	return 0;
}
