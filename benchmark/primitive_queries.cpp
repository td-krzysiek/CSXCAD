// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include "CSPropMetal.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
	try
	{
		const unsigned int count = argc > 1 ? std::stoul(argv[1]) : 1024;
		const unsigned int queries = argc > 2 ? std::stoul(argv[2]) : 20000;
		const unsigned int repeats = argc > 3 ? std::stoul(argv[3]) : 5;
		if (count == 0 || queries == 0 || repeats == 0)
		{
			throw std::invalid_argument("counts must be positive");
		}
		ContinuousStructure csx;
		auto* material = new CSPropMaterial(csx.GetParameterSet());
		auto* metal = new CSPropMetal(csx.GetParameterSet());
		csx.AddProperty(material);
		csx.AddProperty(metal);
		const unsigned int side = static_cast<unsigned int>(std::ceil(std::sqrt(count)));
		for (unsigned int i=0; i<count; ++i)
		{
			auto* box = new CSPrimBox(csx.GetParameterSet(), i % 3 == 0 ?
				static_cast<CSProperties*>(material) : static_cast<CSProperties*>(metal));
			const double x = i % side;
			const double y = i / side;
			box->SetPriority(i % 13);
			box->SetCoord(0, x);
			box->SetCoord(1, x+.8);
			box->SetCoord(2, y);
			box->SetCoord(3, y+.8);
			box->SetCoord(4, 0.0);
			box->SetCoord(5, .1);
		}
		if (!csx.Update().empty())
		{
			throw std::runtime_error("fixture update failed");
		}
		const auto sorted = csx.GetAllPrimitives(true);
		std::uint64_t reference = 0;
		std::cout << "repeat,mode,primitives,queries,seconds,checksum\n";
		for (unsigned int repeat=0; repeat<repeats; ++repeat)
		{
			// Reverse order on alternate trials to reduce systematic warm-up bias.
			for (unsigned int order=0; order<2; ++order)
			{
				const bool reuse = (order + repeat) % 2 != 0;
				std::uint64_t checksum = 0;
				const auto start = std::chrono::steady_clock::now();
				for (unsigned int i=0; i<queries; ++i)
				{
					const unsigned int index = (i * 101ULL) % count;
					const double x = index % side + .4;
					const double y = index / side + .4;
					const double bounds[6] = {x-.1, x+.1, y-.1, y+.1, -.1, .2};
					const auto candidates = reuse ? csx.GetPrimitivesByBoundBox(bounds, sorted) :
						csx.GetPrimitivesByBoundBox(bounds, true);
					const double point[3] = {x, y, .05};
					CSPrimitives* found = nullptr;
					csx.GetPropertyByCoordPriority(point, candidates, false, &found);
					if (!found)
					{
						throw std::runtime_error("missing fixture primitive");
					}
					checksum += found->GetID();
				}
				const double seconds = std::chrono::duration<double>(
					std::chrono::steady_clock::now() - start).count();
				if (repeat == 0 && order == 0)
				{
					reference = checksum;
				}
				if (checksum != reference)
				{
					throw std::runtime_error("lookup checksum changed");
				}
				std::cout << repeat << ',' << (reuse ? "reused" : "legacy") << ',' << count << ','
					<< queries << ',' << std::setprecision(9) << seconds << ',' << checksum << '\n';
			}
		}
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}
