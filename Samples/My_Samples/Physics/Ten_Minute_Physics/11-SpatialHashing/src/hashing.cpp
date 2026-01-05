#include "hashing.h"
#include <wiBacklog.h>

namespace SpatialHashing
{
    Hash::Hash(float spacing, int32_t maxNumObjects)
        : spacing(spacing)
        , tableSize(2 * maxNumObjects)
        , cellStart(tableSize + 1, 0)
        , cellEntries(maxNumObjects, 0)
        , queryIds(maxNumObjects, 0)
        , querySize(0)
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
		// Get the position of the object to query around
        const XMFLOAT3& pos = positions[objectIndex];

		// Compute the integer coordinates of the bounding box surrounding
		// the query position with the given maxDist radius,
		// which in this case is 2 * radius of the balls.
        int32_t x0 = IntCoord(pos.x - maxDist);
        int32_t y0 = IntCoord(pos.y - maxDist);
        int32_t z0 = IntCoord(pos.z - maxDist);

        int32_t x1 = IntCoord(pos.x + maxDist);
        int32_t y1 = IntCoord(pos.y + maxDist);
        int32_t z1 = IntCoord(pos.z + maxDist);

        querySize = 0;

		// Iterate over all cells in the bounding box
		// and collect the object indices from those cells.
		// At the end, queryIds will contain the list of nearby object indices,
		// and querySize will indicate how many were found.
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
                        queryIds[querySize] = cellEntries[i];
                        querySize++;
                    }
                }
            }
        }
    }
}
