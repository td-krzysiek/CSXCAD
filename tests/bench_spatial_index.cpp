// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Krzysztof Rosinski

#include "ContinuousStructure.h"
#include "CSPrimitiveIndex.h"
#include "CSPrimBox.h"
#include "CSPropMaterial.h"

#include <chrono>
#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char** argv)
{
	const bool indexed = argc > 1 && std::string(argv[1]) == "indexed";
	ContinuousStructure csx;
	auto* material = new CSPropMaterial(csx.GetParameterSet());
	csx.AddProperty(material);
	for (int i=0; i<1024; ++i)
	{
		auto* box = new CSPrimBox(csx.GetParameterSet(), material);
		const double bounds[6] = {double(i%32), double(i%32)+.8,
			double(i/32), double(i/32)+.8, 0, 1};
		for (int n=0; n<6; ++n)
		{
			box->SetCoord(n, bounds[n]);
		}
		box->SetPriority(i%13);
	}
	if (!csx.Update().empty())
	{
		throw std::runtime_error("fixture update failed");
	}
	using Clock = std::chrono::steady_clock;
	const auto start = Clock::now();
	CSPrimitiveIndex::Config config;
	config.bounds = {{0, 32, 0, 32, 0, 1}};
	config.bins = {{16, 16, 1}};
	auto index = indexed ? std::make_unique<CSPrimitiveIndex>(csx, config) : nullptr;
	const auto built = Clock::now();
	std::size_t checksum = 0;
	for (int i=0; i<200000; ++i)
	{
		const double x = (i%997)/997.0*31.7;
		const double y = ((i*17)%991)/991.0*31.7;
		const double box[6] = {x, x+.25, y, y+.25, 0, 1};
		const auto result = index ? index->Query(box) : csx.GetPrimitivesByBoundBox(box, true);
		checksum += result.size();
		for (auto* primitive : result)
		{
			checksum += primitive->GetPriority()+1;
		}
	}
	const auto end = Clock::now();
	std::cout << (indexed ? "indexed" : "legacy") << ','
		<< std::chrono::duration<double>(built-start).count() << ','
		<< std::chrono::duration<double>(end-built).count() << ',' << checksum << ','
		<< (index ? index->GetReferenceCount() : 0) << '\n';
}
