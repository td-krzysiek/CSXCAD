# Spatial primitive index (C++ setup API)

Include `CSPrimitiveIndex.h` and construct it from `ContinuousStructure` and its
configuration. No factory method or friendship is required.

For many bounding-box queries against a fixed geometry, build one index rather
than repeatedly collecting and sorting the structure's primitives:

```cpp
// Populate csx, then check the error string returned by csx.Update().
CSPrimitiveIndex::Config config;
config.bounds = {{xmin, xmax, ymin, ymax, zmin, zmax}};
config.bins = {{16, 16, 1}};
config.type = CSProperties::MATERIAL;
CSPrimitiveIndex index(csx, config);

std::vector<CSPrimitives*> candidates; // one output buffer per query thread
index.Query(queryBounds, candidates);
// Pass candidates to the existing exact material-at-point lookup.
```

The caller chooses logical XYZ tile counts in the geometry's input coordinate
system; these tiles do not change the solver mesh. Full-Z column queries can
use `Nz=1`, while spatial queries can use bins in all three directions. The
index is a flat array internally, addressed by logical X/Y/Z tile coordinates.
Each bucket stores ranks into a single priority-sorted snapshot. Queries merge
already ordered lists and remove duplicate ranks. They do not collect or
re-sort the entire geometry.
Broad queries and boxes outside the indexed bounds scan the sorted snapshot.
Unknown bounds and transformed primitives remain candidates even at a distance.
The exact point-containment/material-priority operation is unchanged.

Untransformed Cartesian `Polygon` and `LinPoly` use the new C++
`GetConservativeBoundBox()` contract. It supplies a conservative enclosure after
`Update()` without changing legacy `GetBoundBox()`/`IsInsideBox()` behaviour.
The index expands it by the structure's drawing tolerance, with outward rounding.
Signed extrusion and any polygon normal axis are supported. Transforms,
rotational polygons and non-Cartesian polygon input remain unknown. For supported
polygons, indexed lists can therefore be shorter than legacy candidate lists;
exact material winners must agree, rather than the lists themselves.

For a zero-copy, possibly broader candidate list, opt in at construction:

```cpp
config.prepareCandidateViews = true;
CSPrimitiveIndex index(csx, config);
std::vector<CSPrimitives*> scratch; // one per query thread
const auto& candidates = index.GetCandidates(queryBounds, scratch);
```

For a query contained in one tile, this returns its premerged priority-ordered
list by reference, without a heap or output-vector fill. Crossing/outside queries
use the caller's reusable scratch and ordinary filtering. A borrowed reference
lives until index destruction; a scratch reference until the caller changes the
buffer. Prepared tile lists include the global candidates and may be broader
than `Query()` results, so downstream exact tests can cost more. Measure the full
lookup-plus-containment workload before choosing this API. openEMS currently
retains filtered `Query()` results for that reason.

Ownership stays with `ContinuousStructure`. Keep it and its primitives alive
and unchanged, including drawing tolerance, while the index exists. Call
`Update()` first, and discard/rebuild
the index after any parameter, primitive membership, geometry, coordinate-system,
property type or priority change. There is no automatic invalidation and no
hidden persistent cache. Concurrent read-only queries use separate output
vectors; callers must also respect primitive-specific thread-safety contracts.

Build separate indices for separate property masks to preserve the existing
mask-specific equal-priority tie order. Index construction and bounds queries
do not mark primitives as used. The move-only index cannot be copied by mistake.

`maxBuckets` (default 1,048,576) bounds the product of XYZ tile counts;
`maxReferences` (default 8,388,608) bounds stored bucket/global ranks and the
snapshot primitive count. Invalid/nonfinite bounds, zero bins, indistinguishable
floating-point bin edges, and exhausted limits throw. These are count limits,
not exact RAM limits: vector capacity, pointers and query scratch also consume
memory. Optional prepared-view pointers also count against `maxReferences`.
Cached conservative bounds add a small per-primitive record. Construction
currently costs O(primitives × buckets); it is intended
for repeated queries, not rebuilding at each cell. No new dependencies or XML
fields are introduced. This internal C++ query API has no scripting geometry
semantics to serialize.

```bash
cmake -S . -B build -DCSXCAD_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
# Public synthetic microbenchmark: 1024 boxes, 200000 localized queries.
build/tests/bench_spatial_index legacy
build/tests/bench_spatial_index indexed
```

The benchmark prints mode, build seconds, query seconds, checksum and stored
reference count. Alternate the modes on one machine and report construction
separately. Query speed-up is not a claim about total openEMS setup or FDTD time.

The added virtual bounds method changes the C++ ABI: rebuild CSXCAD and its
dependent native code together. No XML or scripting geometry setting changes.
