#pragma once

#include <cstdint>
#include <vector>
#include <wiScene.h>

// ============================================================================
// SPATIAL HASHING IMPLEMENTATION AND KNOWN ISSUES
// ============================================================================
//
// This class implements a spatial hashing data structure for efficient
// nearest-neighbor queries in 3D space.  The main idea is to divide space into
// a uniform grid of cells, hash each cell to a bucket in a hash table, and
// store object indices in those buckets.
//
// ALGORITHM OVERVIEW:
// -------------------
// 1. The 3D space is conceptually divided into a uniform grid with cell size
//    equal to 'spacing'.
// 2. Each object is assigned to exactly ONE cell based on its position.
// 3. Cell coordinates (xi, yi, zi) are hashed to an integer in range
//    [0, tableSize-1] using the HashCoords() function.
// 4. All objects that hash to the same value are stored contiguously in the
//    cellEntries array.
// 5. The cellStart array (computed via partial sums) stores the starting
//    index in cellEntries for each hash bucket.
//
// QUERY PROCESS:
// --------------
// When querying for nearby objects around a position with radius maxDist:
// 1. Compute the bounding box of cells that could contain nearby objects
// 2. For each cell (xi, yi, zi) in this bounding box:
//    a. Compute hash h = HashCoords(xi, yi, zi)
//    b. Iterate from cellStart[h] to cellStart[h+1] in cellEntries
//    c.  Collect all object indices found in this range
//
// ISSUE:  HASH COLLISIONS CAUSE DUPLICATE IDs IN QUERY RESULTS
// ---------------------------------------------------------------------
// The current implementation suffers from a fundamental problem:  when multiple
// cells hash to the same bucket (hash collision), and those cells fall within
// the same query bounding box, object IDs will be returned MULTIPLE TIMES.
//
// WHY THIS HAPPENS:
// -----------------
// Consider this concrete scenario:
//
// Object ID=3 at position (0.5, 0.5, 0.5) → cell (0,0,0) → HashCoords(0,0,0) = 42
// Object ID=5 at position (1.5, 0.5, 0.5) → cell (1,0,0) → HashCoords(1,0,0) = 42
//
// Both cells (0,0,0) and (1,0,0) hash to the same value:  42 (COLLISION!)
//
// During Hash::Create():
// ----------------------
// First pass (counting objects per hash bucket):
//   cellStart[42]++ for ID=3  →  cellStart[42] = 1
//   cellStart[42]++ for ID=5  →  cellStart[42] = 2
//
// After partial sums:
//   cellStart[42] points to the START of the region for hash value 42
//   Let's say this is position X in cellEntries
//
// Second pass (filling cellEntries):
//   For ID=3:
//     cellStart[42]--          →  cellStart[42] now points to position X
//     cellEntries[X] = 3       →  ID=3 stored at position X
//
//   For ID=5:
//     cellStart[42]--          →  cellStart[42] now points to position X+1
//     cellEntries[X+1] = 5     →  ID=5 stored at position X+1
//
// Result: cellEntries contains both IDs in the SAME contiguous region:
//   cellEntries = [... , 3, 5, ...]
//                      ↑  ↑
//                      X  X+1
//
// During Hash::Query():
// ---------------------
// Now suppose we query around a position where the bounding box includes
// BOTH cells (0,0,0) and (1,0,0). The query will iterate over all cells
// in the bounding box:
//
//   for xi from 0 to 1:      // Covers both cells (0,0,0) and (1,0,0)
//     for yi from 0 to 0:
//       for zi from 0 to 0:
//
//         Iteration 1 - Cell (0,0,0):
//           h = HashCoords(0,0,0)  →  h = 42
//           start = cellStart[42]   →  start = X
//           end = cellStart[43]     →  end = X+2
//           Loop from X to X+2:
//             queryIds[querySize++] = cellEntries[X]    →  adds ID=3
//             queryIds[querySize++] = cellEntries[X+1]  →  adds ID=5
//
//         Iteration 2 - Cell (1,0,0):
//           h = HashCoords(1,0,0)  →  h = 42  (SAME HASH!)
//           start = cellStart[42]   →  start = X  (SAME START!)
//           end = cellStart[43]     →  end = X+2  (SAME END!)
//           Loop from X to X+2:
//             queryIds[querySize++] = cellEntries[X]    →  adds ID=3 AGAIN!
//             queryIds[querySize++] = cellEntries[X+1]  →  adds ID=5 AGAIN!
//
// RESULT: queryIds = [3, 5, 3, 5]  →  DUPLICATES!
//
// ROOT CAUSE ANALYSIS:
// --------------------
// 1. Objects in DIFFERENT cells can hash to the SAME bucket (hash collision)
// 2. These objects are stored together in the SAME region of cellEntries
// 3. When a query bounding box includes multiple cells with the same hash
// 4. The query iterates over the SAME region of cellEntries MULTIPLE TIMES
// 5. Therefore, the same object IDs are added to queryIds multiple times
//
// IMPACT:
// -------
// In collision detection (as in Balls:: Simulate()), duplicate IDs cause:
// - The same collision to be processed multiple times
// - Redundant computations and velocity updates
// - Potential numerical instability from repeated corrections
// - Performance degradation
//
// POSSIBLE SOLUTIONS:
// -------------------
// 1. CHECK FOR DUPLICATES IN QUERY (Simple but adds O(n²) or O(n) overhead):
//    - Linear search through queryIds before adding (slow)
//    - Use a temporary flag array to mark already-seen IDs (faster, O(1) check)
//
// 2. REDUCE HASH COLLISIONS:
//    - Increase tableSize (e.g., from 2*maxNumObjects to 4* or 8*)
//    - Use a better hash function (e.g., Morton codes / Z-order curve)
//    - Both reduce but don't eliminate collisions
//
// 3. USE PERFECT HASHING OR COLLISION-FREE STRUCTURE:
//    - Explicit 3D grid with std::unordered_map<CellCoord, vector<int>>
//    - Guaranteed no duplicates but higher memory overhead
//    - More complex implementation
//
// 4. POST-PROCESS QUERY RESULTS:
//    - Sort queryIds and remove duplicates (std::sort + std::unique)
//    - Adds O(n log n) overhead but guarantees unique results
//
// CURRENT STATUS:
// ---------------
// This implementation follows the approach from the "Ten Minute Physics"
// tutorial series and demonstrates the spatial hashing concept.  However,
// it does NOT handle hash collisions, leading to duplicate IDs in query
// results. For production use, one of the solutions above should be
// implemented based on performance requirements.
//
// PERFORMANCE TRADE-OFFS:
// -----------------------
// - Leaving duplicates:  Fast query, but caller must handle duplicates
// - Flag array check: O(1) duplicate detection, requires extra memory
// - Larger table: Reduces collisions but increases memory usage
// - Perfect hashing: No duplicates guaranteed, but more complex
//
// The choice depends on:
// - Number of objects (affects collision probability)
// - Query frequency vs. creation frequency
// - Available memory
// - Acceptable computational overhead
//
// ============================================================================

namespace SpatialHashing
{
class Hash
{
  public:
    Hash(float spacing, int32_t maxNumObjects);

    // Hash a 3D position
    int32_t HashPos(const XMFLOAT3 &pos) const;

    // Create the spatial hashing structures from positions
    void Create(const std::vector<XMFLOAT3> &positions);

    // Query nearby objects within maxDist of a position
    void Query(const std::vector<XMFLOAT3> &positions, int32_t objectIndex,
               float maxDist);

    // Get query results
    const std::vector<int32_t> &GetQueryIds() const { return queryIds; }
    int32_t GetQuerySize() const { return querySize; }

  private:
    // Hash integer coordinates
    int32_t HashCoords(int32_t xi, int32_t yi, int32_t zi) const;

    // Convert float coordinate to integer grid coordinate
    int32_t IntCoord(float coord) const;

    // Interleave the bits of x, y, z (up to 10 bits per coordinate for 30 bits total)
    // static uint32_t Part1By2(uint32_t n) {
    //     n &= 0x3ff; // 10 bits
    //     n = (n | (n << 16)) & 0x30000ff;
    //     n = (n | (n << 8))  & 0x300f00f;
    //     n = (n | (n << 4))  & 0x30c30c3;
    //     n = (n | (n << 2))  & 0x9249249;
    //     return n;
    // }
    //
    // static uint32_t Morton3D(uint32_t x, uint32_t y, uint32_t z) {
    //     return (Part1By2(z) << 2) | (Part1By2(y) << 1) | Part1By2(x);
    // }

    float spacing;						 //  cell size
    int32_t tableSize;					 // number of hash table entries
    std::vector<int32_t> cellStart;		 // partial sums for cell entries
    std::vector<int32_t> cellEntries;	 // object indices per cell (are stored contiguously)
    std::vector<int32_t> queryIds;		 // nearby object indices from last query
    int32_t querySize;;					 // number of nearby objects from last query
};
} // namespace SpatialHashing
