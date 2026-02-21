#include "hashing.h"
#include <cmath>
#include <algorithm>

namespace SpatialHashing
{
    Hash::Hash(float spacing, int32_t maxNumObjects)
        : spacing(spacing)
        , tableSize(2 * maxNumObjects)
        , maxNumObjects(maxNumObjects)
        , cellStart(tableSize + 1, 0)
        , cellEntries(maxNumObjects, 0)
        , queryIds(maxNumObjects, 0)
        , querySize(0)
        , firstAdjId(maxNumObjects + 1, 0)
        , adjIds(10 * maxNumObjects, 0)  // Initial size, will grow if needed
    {
    }

    int32_t Hash::HashCoords(int32_t xi, int32_t yi, int32_t zi) const
    {
        // Fantasy function from the original code
        int32_t h = (xi * 92837111) ^ (yi * 689287499) ^ (zi * 283923481);
        return std::abs(h) % tableSize;
    }

	// int32_t Hash::HashCoords(int32_t xi, int32_t yi, int32_t zi) const
	// {
	// 	uint32_t morton = Morton3D(static_cast<uint32_t>(xi), static_cast<uint32_t>(yi), static_cast<uint32_t>(zi));
	// 	return static_cast<int32_t>(morton % tableSize);
	// }

    int32_t Hash::IntCoord(float coord) const
    {
        return static_cast<int32_t>(std::floor(coord / spacing));
    }

    int32_t Hash::HashPos(const XMFLOAT3& pos) const
    {
        return HashCoords(
            IntCoord(pos.x),
            IntCoord(pos.y),
            IntCoord(pos.z)
        );
    }

    void Hash::Create(const std::vector<XMFLOAT3>& positions)
    {
        int32_t numObjects = std::min(
            static_cast<int32_t>(positions.size()),
            static_cast<int32_t>(cellEntries.size())
        );

        // Clear previous data
        std::fill(cellStart.begin(), cellStart.end(), 0);
        std::fill(cellEntries.begin(), cellEntries.end(), 0);

		// Compute the number of objects in each cell of the conceptual 3D grid.
		// This populates cellStart with counts
        for (int32_t i = 0; i < numObjects; i++)
        {
			// Each cell is defined by integer coordinates (xi, yi, zi), obtained by dividing
			// the object's position by the grid cell spacing and flooring the result.
			// The hash function maps these integer coordinates to an index in cellStart,
			// which counts how many objects fall into each grid cell.
            int32_t h = HashPos(positions[i]);
            cellStart[h]++;
        }

        // Compute the partial sums
        int32_t start = 0;
        for (int32_t i = 0; i < tableSize; i++)
        {
            start += cellStart[i];
            cellStart[i] = start;
        }
        cellStart[tableSize] = start; // guard; the value should be equal to numObjects

        // Fill in object ids
        for (int32_t i = 0; i < numObjects; i++)
        {
			// Retrieve the cellStart index for this object
			// Decrement the count in cellStart to get the index within cellEntries
			// and store the object index in cellEntries.
			// See the lesson notes for details:
			// https://matthias-research.github.io/pages/tenMinutePhysics/11-hashing.pdf
            int32_t h = HashPos(positions[i]);
            cellStart[h]--;
            cellEntries[cellStart[h]] = i;
        }
    }

	void Hash::Query(const std::vector<XMFLOAT3>& positions, int32_t objectIndex, float maxDist)
	{
		const XMFLOAT3& pos = positions[objectIndex];

		int32_t x0 = IntCoord(pos.x - maxDist);
		int32_t y0 = IntCoord(pos.y - maxDist);
		int32_t z0 = IntCoord(pos.z - maxDist);

		int32_t x1 = IntCoord(pos.x + maxDist);
		int32_t y1 = IntCoord(pos.y + maxDist);
		int32_t z1 = IntCoord(pos.z + maxDist);

		querySize = 0;

		// flag array to avoid duplicates
		std::vector<bool> alreadyAdded(queryIds.size(), false);

		for (int32_t xi = x0; xi <= x1; xi++)
		{
			for (int32_t yi = y0; yi <= y1; yi++)
			{
				for (int32_t zi = z0; zi <= z1; zi++)
				{
					int32_t h = HashCoords(xi, yi, zi);
					int32_t start = cellStart[h];
					int32_t end = cellStart[h + 1];

					for (int32_t i = start; i < end; i++)
					{
						int32_t candidateId = cellEntries[i];

						// Skip if already added
						if (alreadyAdded[candidateId])
							continue;

						// Add ID to the query result and mark as added
						if (querySize < static_cast<int32_t>(queryIds.size()))
						{
							queryIds[querySize++] = candidateId;
							alreadyAdded[candidateId] = true;
						}
					}
				}
			}
		}
	}

    void Hash::QueryAll(const std::vector<XMFLOAT3>& positions, float maxDist)
    {
        int32_t num = 0;
        float maxDist2 = maxDist * maxDist;

		// For each object...
        for (int32_t i = 0; i < maxNumObjects; i++)
        {
			// Store num as the starting index in adjIds for the adjacency list corresponding to this object (indexed by id0).
            int32_t id0 = i;
            firstAdjId[id0] = num;

            // Query nearby objects for this particle
            Query(positions, id0, maxDist);

            for (int32_t j = 0; j < querySize; j++)
            {
                int32_t id1 = queryIds[j];

                // Only consider pairs where id1 < id0 to avoid duplicates
                if (id1 >= id0)
                    continue;

                // Check actual distance
                float dx = positions[id0].x - positions[id1].x;
                float dy = positions[id0].y - positions[id1].y;
                float dz = positions[id0].z - positions[id1].z;
                float dist2 = dx * dx + dy * dy + dz * dz;

				// Only consider objects within maxDist.
                if (dist2 > maxDist2)
                    continue;

                // Grow adjIds array if needed
                if (num >= static_cast<int32_t>(adjIds.size()))
                {
                    adjIds.resize(2 * num);
                }

				// Store the id of the adjacent object in adjIds and increment the count.
				// All adjacent IDs for object id0 will be stored contiguously starting from firstAdjId[id0].
				// Example: if firstAdjId[3] = 10 and the next firstAdjId[4] = 15,
				// then adjIds[10] to adjIds[14] are the adjacent IDs for object 3.
				// So, to access the adiacent IDs for object 3, we can retrieve the first index in adjIds
				// from firstAdjId[3] and the last one from firstAdjId[4] (exclusive) and then iterate
				// through adjIds in that range to get all adjacent IDs.
                adjIds[num++] = id1;
            }
        }

		// Store the total number of adiacent objects at the end of firstAdjId as a guard
		// (not strictly necessary, but can be useful for debugging)
        firstAdjId[maxNumObjects] = num;
    }
}
