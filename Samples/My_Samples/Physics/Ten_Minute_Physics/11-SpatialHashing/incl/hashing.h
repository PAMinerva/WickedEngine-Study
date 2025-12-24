#pragma once

#include <cstdint>
#include <vector>
#include <wiScene.h>

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
