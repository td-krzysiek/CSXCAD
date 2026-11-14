# Repeated primitive queries

Build with the normal library dependencies and `-DCSXCAD_BUILD_BENCHMARKS=ON`,
then run `benchmark/primitive_queries [primitive-count] [query-count] [repeats]`
from the build directory. Defaults are 1,024 primitives, 20,000 queries and 5 repeats.
The fixture is synthetic; no solver or private geometry is needed.

The CSV compares collecting/priority-sorting the geometry for each bounding-box
query with reusing one list obtained before the queries. Trials alternate order
and verify the same winning-primitive checksum. Report medians and the build,
compiler and machine. These timings measure a geometry-query workload, not total
openEMS startup or FDTD throughput.

```cpp
// Refresh after adding/deleting primitives or changing property types/priorities.
const auto primitives = csx.GetAllPrimitives(true, CSProperties::MATERIAL);
for (const auto& bounds : queryBounds)
{
    const auto candidates = csx.GetPrimitivesByBoundBox(bounds.data(), primitives);
    // GetPropertyByCoordPriority(point, candidates) preserves their priority order.
}
```

The list owns no primitives. Keep the structure/primitives alive, update geometry
before querying changed coordinates, and discard the list before deleting its
objects. The new overload does not install a cache in CSXCAD. Existing callers
continue to use the original overload until they explicitly reuse a list.
