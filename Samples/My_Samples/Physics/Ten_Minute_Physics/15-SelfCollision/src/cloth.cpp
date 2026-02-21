#include <Utility/DirectXMath/DirectXMath.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <wiScene.h>
#include "cloth.h"
#include "simulation_utils.h"

ClothMesh::ClothMesh(const SimulationParams& params)
    : numX(params.numX)
    , numY(params.numY)
    , thickness(params.thickness)
    , handleCollisions(true)
{
    float spacing = params.spacing;
    float jitter = 0.001f * spacing;

    // Random number generator for jitter
    std::random_device rd;
    std::mt19937 gen(rd());
    // std::uniform_real_distribution<float> dis(-jitter, jitter);
	std::uniform_real_distribution<float> unit(0.0f, 1.0f);

    // ========== PARTICLES ==========
    numParticles = numX * numY;

    pos.resize(numParticles);
    prevPos.resize(numParticles);
    restPos.resize(numParticles);
    vel.resize(numParticles);
    invMass.resize(numParticles);

    // Initialize particle positions in a grid
    // Grid is in XY plane, centered at origin
	// Initialize positions column by column starting from the bottom-left corner and going upwards,
	// then moving to the next column to the right.
	// id = i * numY + j follows the same convention: 0 refers to (0,0), 1 to (0,1), ..., numY-1 to (0,numY-1), numY to (1,0), etc.
    for (int i = 0; i < numX; i++)
    {
        for (int j = 0; j < numY; j++)
        {
            int id = i * numY + j;

            float x = -numX * spacing * 0.5f + i * spacing;
            float y = 0.2f + j * spacing;  // Start slightly above ground
            float z = 0.0f;

            // Add small jitter to break symmetry
            // pos[id] = XMFLOAT3(x + dis(gen), y + dis(gen), z + dis(gen));
            float ox = -2.0f * jitter * jitter * unit(gen);
            float oy = -2.0f * jitter * jitter * unit(gen);
            float oz = -2.0f * jitter * jitter * unit(gen);
            pos[id] = XMFLOAT3(x + ox, y + oy, z + oz);

            invMass[id] = 1.0f;

            // Optional: attach top corners (uncomment to enable)
            // if (j == numY - 1 && (i == 0 || i == numX - 1))
            //     invMass[id] = 0.0f;
        }
    }

    // Copy initial positions to rest and prev
    restPos = pos;
    prevPos = pos;

    // Initialize velocities to zero
    for (int i = 0; i < numParticles; i++)
    {
        vel[i] = XMFLOAT3(0.0f, 0.0f, 0.0f);
    }

    // ========== SPATIAL HASHING ==========
    hash = std::make_unique<SpatialHashing::Hash>(spacing, numParticles);

    // ========== CONSTRAINTS ==========
    // 6 constraint types:
    // 0: Stretch vertical   (0,0)-(0,1)
    // 1: Stretch horizontal (0,0)-(1,0)
    // 2: Shear diagonal 1   (0,0)-(1,1)
    // 3: Shear diagonal 2   (0,1)-(1,0)
    // 4: Bending vertical   (0,0)-(0,2)
    // 5: Bending horizontal (0,0)-(2,0)

    const int numConstraintTypes = 6;
    int offsets[numConstraintTypes * 4] = {
        0, 0, 0, 1,  // Stretch vertical
        0, 0, 1, 0,  // Stretch horizontal
        0, 0, 1, 1,  // Shear diagonal 1
        0, 1, 1, 0,  // Shear diagonal 2
        0, 0, 0, 2,  // Bending vertical
        0, 0, 2, 0   // Bending horizontal
    };

    float complianceValues[numConstraintTypes] = {
        params.stretchingCompliance,  // Stretch vertical
        params.stretchingCompliance,  // Stretch horizontal
        params.shearCompliance,       // Shear diagonal 1
        params.shearCompliance,       // Shear diagonal 2
        params.bendingCompliance,     // Bending vertical
        params.bendingCompliance      // Bending horizontal
    };

    // Temporary vectors to collect constraints
    std::vector<XMINT2> tempIds;
    std::vector<float> tempCompliances;

	// This loop generates all constraint pairs for the cloth grid.
	// For each constraint type (stretch, shear, bending), it connects two particles
	// according to the offsets defined above. Each constraint is stored as a pair of indices.
	//
	// Example for numX = 4, numY = 4 (particle indices: id = i * numY + j):
	//
	// Grid layout (indices):
	// (i,j):
	// (0,0)=0   (1,0)=4   (2,0)=8   (3,0)=12
	// (0,1)=1   (1,1)=5   (2,1)=9   (3,1)=13
	// (0,2)=2   (1,2)=6   (2,2)=10  (3,2)=14
	// (0,3)=3   (1,3)=7   (2,3)=11  (3,3)=15
	//
	// Example: Stretch vertical: (offset: 0,0,0,1)
	//   Connects (i,j) with (i,j+1) if j+1 < numY
	//   Pairs: (0,0)-(0,1): 0-1, (0,1)-(0,2): 1-2, (0,2)-(0,3): 2-3, ...
	//
	// Example: Stretch horizontal (offset: 0,0,1,0)
	//   Connects (i,j) with (i+1,j) if i+1 < numX
	//   Pairs: (0,0)-(1,0): 0-4, (1,0)-(2,0): 4-8, (2,0)-(3,0): 8-12, ...
	//
	// Example: Shear diagonal 1 (offset: 0,0,1,1)
	//   Connects (i,j) with (i+1,j+1) if i+1 < numX && j+1 < numY
	//   Pairs: (0,0)-(1,1): 0-5, (1,0)-(2,1): 4-9, (2,0)-(3,1): 8-13, ...
	//
	// Example: Shear diagonal 2 (offset: 0,1,1,0)
	//   Connects (i,j+1) with (i+1,j) if i+1 < numX && j+1 < numY
	//   Pairs: (0,1)-(1,0): 1-4, (1,1)-(2,0): 5-8, (2,1)-(3,0): 9-12, ...
	//
	// Example: Bending vertical (offset: 0,0,0,2)
	//   Connects (i,j) with (i,j+2) if j+2 < numY
	//   Pairs: (0,0)-(0,2): 0-2, (0,1)-(0,3): 1-3, ...
	//
	// Example: Bending horizontal (offset: 0,0,2,0)
	//   Connects (i,j) with (i+2,j) if i+2 < numX
	//   Pairs: (0,0)-(2,0): 0-8, (1,0)-(3,0): 4-12, (0,1)-(2,1): 1-9, ...
    for (int constType = 0; constType < numConstraintTypes; constType++)
    {
        int p = constType * 4;

        for (int i = 0; i < numX; i++)
        {
            for (int j = 0; j < numY; j++)
            {
                int i0 = i + offsets[p + 0];
                int j0 = j + offsets[p + 1];
                int i1 = i + offsets[p + 2];
                int j1 = j + offsets[p + 3];

                // Check bounds
                if (i0 < numX && j0 < numY && i1 < numX && j1 < numY)
                {
                    int id0 = i0 * numY + j0;
                    int id1 = i1 * numY + j1;

                    tempIds.push_back(XMINT2(id0, id1));
                    tempCompliances.push_back(complianceValues[constType]);
                }
            }
        }
    }

    numConstraints = static_cast<int>(tempIds.size());
    constraintIds = std::move(tempIds);
    compliances = std::move(tempCompliances);

	// REMOVE THIS!
	// Count bending constraints (last 2 types: vertical and horizontal bending)
	// These are constraints of type 4 and 5
	numBendingConstraints = 0;
	for (int i = 0; i < numX; i++)
	{
		for (int j = 0; j < numY; j++)
		{
			// Bending vertical: (0,0)-(0,2)
			if (i < numX && (j + 2) < numY)
				numBendingConstraints++;

			// Bending horizontal: (0,0)-(2,0)
			if ((i + 2) < numX && j < numY)
				numBendingConstraints++;
		}
	}

	// Store compliance values
	stretchingCompliance = params.stretchingCompliance;
	shearCompliance = params.shearCompliance;
	bendingCompliance = params.bendingCompliance;

    // Compute rest lengths
    restLens.resize(numConstraints);
    for (int i = 0; i < numConstraints; i++)
    {
        int id0 = constraintIds[i].x;
        int id1 = constraintIds[i].y;

        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);
        XMVECTOR diff = XMVectorSubtract(p0, p1);

        restLens[i] = XMVectorGetX(XMVector3Length(diff));
    }

    // ========== TRIANGLE INDICES (for rendering) ==========
	// It generates two triangles for each quad in the grid. Each quad is formed by four particles:
	// (i,j), (i+1,j), (i,j+1), (i+1,j+1) with indices id, id+numY, id+1, id+1+numY respectively.
	// The two triangles counter-clockwise are:
	// 1) (id+1, id, id+1+numY) -> (i,j+1), (i,j), (i+1,j+1)
	// 2) (id+1+numY, id, id+numY) -> (i+1,j+1), (i,j), (i+1,j)
    for (int i = 0; i < numX - 1; i++)
    {
        for (int j = 0; j < numY - 1; j++)
        {
            int id = i * numY + j;

            // First triangle: (id+1, id, id+1+numY)
            triIds.push_back(id + 1);
            triIds.push_back(id);
            triIds.push_back(id + 1 + numY);

            // Second triangle: (id+1+numY, id, id+numY)
            triIds.push_back(id + 1 + numY);
            triIds.push_back(id);
            triIds.push_back(id + numY);
        }
    }
    numTris = static_cast<int>(triIds.size()) / 3;

	// REMOVE THIS!
    // ========== EDGE INDICES (for wireframe) ==========
    for (int i = 0; i < numX; i++)
    {
        for (int j = 0; j < numY; j++)
        {
            int id = i * numY + j;

            // Horizontal edge
            if (i < numX - 1)
            {
                edgeIds.push_back(id);
                edgeIds.push_back(id + numY);
            }

            // Vertical edge
            if (j < numY - 1)
            {
                edgeIds.push_back(id);
                edgeIds.push_back(id + 1);
            }
        }
    }

	// REMOVE THIS!
    // ========== TEMPORARY VECTORS FOR CALCULATIONS ==========
    vecs.resize(4);  // For intermediate calculations in solvers
}

void ClothMesh::Simulate(float frameDt, int numSubSteps, const XMFLOAT3& gravity)
{
	// dt is the timestep for each substep.
	// maxVelocity is computed based on the thickness and dt to avoid tunneling issues in collisions:
	// Tunneling occurs when particles move so fast that they pass through each other without detecting a collision.
	// By clamping the velocity, we can avoid collisions being missed due to large time steps or thin cloth.
	// It ensures that particles don't move more than a fraction of the thickness in one substep.
    float dt = frameDt / static_cast<float>(numSubSteps);
    float maxVelocity = 0.2f * thickness / dt;

    if (handleCollisions)
    {
	    // Build spatial hash once per frame (before substeps)
        hash->Create(pos);
		// maxTravelDist is the maximum distance a particle can move in one frame (not substep) based on the max velocity and frame time.
		// Remember that s = v * t.
        float maxTravelDist = maxVelocity * frameDt;
		// Query neighbors for each particle within maxTravelDist.
        hash->QueryAll(pos, maxTravelDist);
    }

    XMVECTOR gravityVec = XMLoadFloat3(&gravity);

	// For each substep...
    for (int step = 0; step < numSubSteps; step++)
    {
		//
        // ========== INTEGRATE ==========
		//
		// For each particle...
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] > 0.0f)
            {
				// Load velocity and position into XMVECTOR for calculations
                XMVECTOR velVec = XMLoadFloat3(&vel[i]);
                XMVECTOR posVec = XMLoadFloat3(&pos[i]);

                // Apply gravity:
				// vel += gravity * dt
                velVec = XMVectorAdd(velVec, XMVectorScale(gravityVec, dt));

                // Clamp velocity to avoid tunneling issues in collisions
                float v = XMVectorGetX(XMVector3Length(velVec));
                float maxV = 0.2f * thickness / dt;
                if (v > maxV)
                {
                    velVec = XMVectorScale(velVec, maxV / v);
                }

                XMStoreFloat3(&vel[i], velVec);

                // Save previous position
                prevPos[i] = pos[i];

                // Update position: pos += vel * dt
                posVec = XMVectorAdd(posVec, XMVectorScale(velVec, dt));
                XMStoreFloat3(&pos[i], posVec);
            }
        }

        // ========== SOLVE ==========
        SolveGroundCollisions();
        SolveConstraints(dt);

        if (handleCollisions)
            SolveCollisions(dt);

        // ========== UPDATE VELOCITIES ==========
        float invDt = 1.0f / dt;
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] > 0.0f)
            {
                XMVECTOR posVec = XMLoadFloat3(&pos[i]);
                XMVECTOR prevPosVec = XMLoadFloat3(&prevPos[i]);

                // vel = (pos - prevPos) / dt
                XMVECTOR velVec = XMVectorScale(XMVectorSubtract(posVec, prevPosVec), invDt);
                XMStoreFloat3(&vel[i], velVec);
            }
        }
    }
}

void ClothMesh::SolveGroundCollisions()
{
	// The ground level is set to half the thickness of the cloth.
	// This allows the cloth to rest on the ground without penetrating it, while still allowing for some visual thickness.
	// Also sets a damping factor to -1 to stop particles that collide with the ground.
    float groundLevel = 0.5f * thickness;
    float damping = 1.0f;

	// For each particle...
    for (int i = 0; i < numParticles; i++)
    {
		// If the particle is fixed (invMass == 0), skip it since it cannot move or collide.
        if (invMass[i] == 0.0f)
            continue;

		// Check if the particle is below the ground level.
		// If that's the case, we need to resolve the collision by moving it back to
		// the ground level and stopping it from moving further in every direction.
        float y = pos[i].y;
        if (y < groundLevel)
        {
            // Compute the difference between the current position and
			// the previous position to determine how much the particle has
			// moved since the last frame.
            XMVECTOR posVec = XMLoadFloat3(&pos[i]);
            XMVECTOR prevPosVec = XMLoadFloat3(&prevPos[i]);
            XMVECTOR diff = XMVectorSubtract(posVec, prevPosVec);

            // Stop the particle from moving further in every direction after
			// collision with the ground by scaling the difference vector by
			// -1, which effectively negates the movement that occurred in the
			// integreation step of the Simulate function.
            posVec = XMVectorAdd(posVec, XMVectorScale(diff, -damping));

            // Clamp position to ground level
            XMStoreFloat3(&pos[i], posVec);
            pos[i].y = groundLevel;
        }
    }
}

void ClothMesh::SolveConstraints(float dt)
{
    float dtSq = dt * dt;

	// For each constraint...
    for (int i = 0; i < numConstraints; i++)
    {
		// Get the indices of the two particles connected by this constraint.
		// Regardless of the constraint type (stretch, shear, bending), we can always treat it as
		// a distance constraint that tries to keep the two particles at a certain rest length.
        int id0 = constraintIds[i].x;
        int id1 = constraintIds[i].y;

		// Get the inverse masses of the two particles and compute their sum.
        float w0 = invMass[id0];
        float w1 = invMass[id1];
        float w = w0 + w1;

		// If both are zero, it means they are fixed and we can skip this constraint.
        if (w == 0.0f)
            continue;

		// Get the current positions of the two particles.
        XMVECTOR p0 = XMLoadFloat3(&pos[id0]);
        XMVECTOR p1 = XMLoadFloat3(&pos[id1]);

		// Compute the difference vector between the two particles and its length.
        XMVECTOR diff = XMVectorSubtract(p0, p1);
        float len = XMVectorGetX(XMVector3Length(diff));

		// If the length is zero, it means the particles are at the same position, so we can skip this constraint to avoid division by zero.
        if (len == 0.0f)
            continue;

        // Normalize the difference vector to get the direction of the constraint force.
        XMVECTOR dir = XMVectorScale(diff, 1.0f / len);

		// Retrieve the rest length of the constraint and compute the constraint violation (C) as the
		// difference between the current length and the rest length.
		// Compute the scaling factor (s) for the constraint correction:
		// See 10-SoftBodies sample for further explanation.
        float restLen = restLens[i];
        float C = len - restLen;
        float alpha = compliances[i] / dtSq;
        float s = -C / (w + alpha);

        // Apply corrections
        p0 = XMVectorAdd(p0, XMVectorScale(dir, s * w0));
        p1 = XMVectorSubtract(p1, XMVectorScale(dir, s * w1));

		// Store the corrected positions back to the pos array.
        XMStoreFloat3(&pos[id0], p0);
        XMStoreFloat3(&pos[id1], p1);
    }
}

void ClothMesh::SolveCollisions(float dt)
{
    float thickness2 = thickness * thickness;

	// For each particle...
    for (int i = 0; i < numParticles; i++)
    {
		// If the particle is fixed (invMass == 0), skip it since it cannot move or collide.
        if (invMass[i] == 0.0f)
            continue;

		// Get the range of adjacent particle indices from the spatial hash for the current particle.
        int first = hash->GetFirstAdjId(i);
        int last = hash->GetFirstAdjId(i + 1);

		// Get the current position of the particle.
        XMVECTOR p0 = XMLoadFloat3(&pos[i]);

		// For each index in the adjacency list for this particle...
        for (int j = first; j < last; j++)
        {
			// Get the index of the adjacent particle.
            int id1 = hash->GetAdjId(j);

			// If the adjacent particle is fixed (invMass == 0), skip it since it cannot move or collide.
            if (invMass[id1] == 0.0f)
                continue;

			// Get the current position of the adjacent particle.
            XMVECTOR p1 = XMLoadFloat3(&pos[id1]);

            // Compute difference vector and squared distance between the two particles.
            XMVECTOR diff = XMVectorSubtract(p0, p1);
            float dist2 = XMVectorGetX(XMVector3LengthSq(diff));

			// If the squared distance is greater than the squared thickness,
			// it means the particles are not colliding, so we can skip to the next adjacent particle.
            if (dist2 > thickness2 || dist2 == 0.0f)
                continue;

			// Compute the actual distance by taking the square root of the squared distance.
            float dist = std::sqrt(dist2);

			// Normalize the difference vector to get the direction of the collision response.
            XMVECTOR dir = XMVectorScale(diff, 1.0f / dist);

			// Set the rest distance to the thickness of the cloth, which is the minimum allowed distance between particles.
			// Then proceed to resolve the collision by treating it as a distance constraint.
			// Note that no compliance is used for collision constraints since we want to resolve collisions rigidly.
            float restDist = thickness;
            float C = dist - restDist;

            float w0 = invMass[i];
            float w1 = invMass[id1];
            float w = w0 + w1;

            float s = -C / w;

            // Apply corrections
            p0 = XMVectorAdd(p0, XMVectorScale(dir, s * w0));
            p1 = XMVectorSubtract(p1, XMVectorScale(dir, s * w1));

            // Store p1 immediately since we might not see it again
            XMStoreFloat3(&pos[id1], p1);
        }

        // Store p0 after processing all its neighbors
        XMStoreFloat3(&pos[i], p0);
    }
}

void ClothMesh::StartGrab(const wi::scene::PickResult& pick)
{
    // Find the closest vertex to pick.position
    XMVECTOR pickPos = XMVectorSet(pick.position.x, pick.position.y, pick.position.z, 0.0f);
    float minDistSq = std::numeric_limits<float>::max();
    grabId = -1;

    for (int i = 0; i < numParticles; i++)
    {
        XMVECTOR p = XMLoadFloat3(&pos[i]);
        float d2 = XMVectorGetX(XMVector3LengthSq(XMVectorSubtract(p, pickPos)));
        if (d2 < minDistSq) {
            minDistSq = d2;
            grabId = i;
        }
    }

    // If a valid vertex is found, fix it at pick.position
    if (grabId >= 0)
    {
        grabInvMass = invMass[grabId];
        invMass[grabId] = 0.0f;
        XMStoreFloat3(&pos[grabId], pickPos);
    }
}

void ClothMesh::MoveGrabbed(const std::vector<float>& pos, const std::vector<float>& vel)
{
    if (grabId >= 0)
    {
        // Copy pos to this->pos[grabId]
        XMVECTOR p = XMVectorSet(pos[0], pos[1], pos[2], 0.0f);
        XMStoreFloat3(&this->pos[grabId], p);
    }
}

void ClothMesh::EndGrab(const std::vector<float>& pos, std::vector<float>& vel)
{
    if (grabId >= 0)
    {
        // Restore inverse mass
        invMass[grabId] = grabInvMass;

        // Copy vel to this->vel[grabId]
        XMVECTOR v = XMVectorSet(vel[0], vel[1], vel[2], 0.0f);
        XMStoreFloat3(&this->vel[grabId], v);
    }
    grabId = -1;
}
