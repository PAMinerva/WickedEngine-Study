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

    // Query nearby objects within maxDist of a position (object) specified by objectIndex, and store results in queryIds
    void Query(const std::vector<XMFLOAT3> &positions, int32_t objectIndex, float maxDist);

    // Query all nearby objects for every position\object and store results in adjacency list format (for self-collision)
    void QueryAll(const std::vector<XMFLOAT3> &positions, float maxDist);

    // Get query results (for Query)
    const std::vector<int32_t> &GetQueryIds() const { return queryIds; }
    int32_t GetQuerySize() const { return querySize; }

    // Get adjacency results (for QueryAll)
    int32_t GetFirstAdjId(int32_t index) const { return firstAdjId[index]; }
    int32_t GetAdjId(int32_t index) const { return adjIds[index]; }

  private:
    // Hash integer coordinates
    int32_t HashCoords(int32_t xi, int32_t yi, int32_t zi) const;

    // Convert float coordinate to integer grid coordinate
    int32_t IntCoord(float coord) const;

    float spacing;                       // cell size
    int32_t tableSize;                   // number of hash table entries
    int32_t maxNumObjects;               // maximum number of objects
    std::vector<int32_t> cellStart;      // partial sums for cell entries
    std::vector<int32_t> cellEntries;    // object indices per cell
    std::vector<int32_t> queryIds;       // nearby object indices from last query
    int32_t querySize;                   // number of nearby objects from last query

    // Adjacency data for QueryAll (self-collision)
    std::vector<int32_t> firstAdjId;     // start index in adjIds for each object
    std::vector<int32_t> adjIds;         // adjacent object indices
};
} // namespace SpatialHashing
