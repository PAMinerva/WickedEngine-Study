#include "balls.h"
#include "simulation_utils.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <wiScene.h>

Balls::Balls(float radius_, const std::vector<XMFLOAT3> &initialPositions,
             const std::vector<XMFLOAT3> &initialVelocities,
             const std::array<XMFLOAT3, 2> &worldBounds_)
    : radius(radius_), numBalls(static_cast<int>(initialPositions.size())),
      pos(initialPositions), prevPos(initialPositions), vel(initialVelocities),
      worldBounds(worldBounds_)
{
    hash = std::make_unique<SpatialHashing::Hash>(radius * 2.0f, // cell size
												  numBalls		 // max number of objects
												  );

    CreateMaterials();
    CreateBalls();
}

void Balls::CreateBalls()
{
    using namespace wi::ecs;
    using namespace wi::scene;

    // Create a temporary mesh+object
    wi::ecs::Entity tempEntity =
        GetScene().Entity_CreateSphere("SharedSphereMesh", radius, 8, 8);

    // Extract the meshID
    wi::scene::ObjectComponent *tempObj =
        GetScene().objects.GetComponent(tempEntity);
    wi::ecs::Entity meshID = tempObj->meshID;

    // Store the shared mesh entity
    mesh_entity = meshID;

    // Remove object, transform and layer components (we only need the mesh)
    GetScene().objects.Remove(tempEntity);
    GetScene().transforms.Remove(tempEntity);
    GetScene().layers.Remove(tempEntity);

    // Set the white material to the mesh
    wi::scene::MeshComponent *mesh = GetScene().meshes.GetComponent(meshID);
    if (mesh && !mesh->subsets.empty())
    {
        mesh->subsets[0].materialID = white_mat_entity;
    }

    // Create instances of the sphere mesh
    for (size_t i = 0; i < pos.size(); ++i)
    {
        wi::ecs::Entity sphereEntity =
            GetScene().Entity_CreateObject("Sphere_" + std::to_string(i));

        wi::scene::ObjectComponent *obj =
            GetScene().objects.GetComponent(sphereEntity);
        obj->meshID = meshID; // Use shared mesh for all sphere instances

        // Set initial color to red
        obj->color = XMFLOAT4{0.9f, 0.0f, 0.0f, 1.0f};

        // Use initial position (see init_positions_and_velocities in
        // sample.cpp)
        wi::scene::TransformComponent *transform =
            GetScene().transforms.GetComponent(sphereEntity);
        transform->translation_local = pos[i];
        transform->SetDirty();

        // Store object (instance) entity
        // It will be used later to update position and material of sphere
        // instances during simulation see Balls::UpdateSpherePosition and
        // Balls::SetSphereMaterial
        sphereEntities.push_back(sphereEntity);
    }
}

void Balls::CreateMaterials()
{
    using namespace wi::ecs;
    using namespace wi::scene;
    // Create a material that returns white color when enlighted
    if (white_mat_entity == INVALID_ENTITY)
    {
        white_mat_entity = CreateEntity();
        GetScene().materials.Create(white_mat_entity);
        MaterialComponent &material =
            *GetScene().materials.GetComponent(white_mat_entity);
        material.shaderType = MaterialComponent::SHADERTYPE_PBR;
        material.baseColor = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
    }
}

void Balls::SetSphereMaterial(size_t sphereIndex, bool inCollision)
{
    wi::scene::ObjectComponent *obj =
        wi::scene::GetScene().objects.GetComponent(sphereEntities[sphereIndex]);

    if (obj)
    {
        if (inCollision)
        {
            obj->color = XMFLOAT4(0.9f, 0.5f, 0.0f, 1.0f);
        }
        else
        {
            obj->color = XMFLOAT4(0.9f, 0.0f, 0.0f, 1.0f);
        }
    }
}

void Balls::UpdateSpherePosition(int32_t i)
{
    wi::scene::TransformComponent *transform =
        wi::scene::GetScene().transforms.GetComponent(sphereEntities[i]);

    if (transform)
    {
        transform->translation_local = pos[i];
        transform->SetDirty();
    }
}

void Balls::Simulate(float dt, const XMFLOAT3 &gravity)
{
    float minDist = 2.0f * radius;
    XMVECTOR gravityVec = XMLoadFloat3(&gravity);

    // Integrate
    for (int32_t i = 0; i < numBalls; i++)
    {
        XMVECTOR velVec = XMLoadFloat3(&vel[i]);
        XMVECTOR posVec = XMLoadFloat3(&pos[i]);

        // vel += gravity * dt
        velVec = XMVectorAdd(velVec, XMVectorScale(gravityVec, dt));

        // prevPos = pos
        prevPos[i] = pos[i];

        // pos += vel * dt
        posVec = XMVectorAdd(posVec, XMVectorScale(velVec, dt));

        // Store back to XMFLOAT3
        // Now velocity and position are updated for all balls
        XMStoreFloat3(&vel[i], velVec);
        XMStoreFloat3(&pos[i], posVec);
    }

    // Create and initialize spatial hashing structures
	// (see Hash::Create and the lesson notes for details:
	// https://matthias-research.github.io/pages/tenMinutePhysics/11-hashing.pdf)
    hash->Create(pos);

    // Handle collisions
    for (int32_t i = 0; i < numBalls; i++)
    {
        bool inCollision = false;
        XMVECTOR posVec = XMLoadFloat3(&pos[i]);
        XMVECTOR velVec = XMLoadFloat3(&vel[i]);

        // World bounds collision
        if (XMVectorGetX(posVec) < worldBounds[0].x + radius)
        {
            posVec = XMVectorSetX(posVec, worldBounds[0].x + radius);
            velVec = XMVectorSetX(velVec, -XMVectorGetX(velVec));
            inCollision = true;
        }
        else if (XMVectorGetX(posVec) > worldBounds[1].x - radius)
        {
            posVec = XMVectorSetX(posVec, worldBounds[1].x - radius);
            velVec = XMVectorSetX(velVec, -XMVectorGetX(velVec));
            inCollision = true;
        }

        if (XMVectorGetY(posVec) < worldBounds[0].y + radius)
        {
            posVec = XMVectorSetY(posVec, worldBounds[0].y + radius);
            velVec = XMVectorSetY(velVec, -XMVectorGetY(velVec));
            inCollision = true;
        }
        else if (XMVectorGetY(posVec) > worldBounds[1].y - radius)
        {
            posVec = XMVectorSetY(posVec, worldBounds[1].y - radius);
            velVec = XMVectorSetY(velVec, -XMVectorGetY(velVec));
            inCollision = true;
        }

        if (XMVectorGetZ(posVec) < worldBounds[0].z + radius)
        {
            posVec = XMVectorSetZ(posVec, worldBounds[0].z + radius);
            velVec = XMVectorSetZ(velVec, -XMVectorGetZ(velVec));
            inCollision = true;
        }
        else if (XMVectorGetZ(posVec) > worldBounds[1].z - radius)
        {
            posVec = XMVectorSetZ(posVec, worldBounds[1].z - radius);
            velVec = XMVectorSetZ(velVec, -XMVectorGetZ(velVec));
            inCollision = true;
        }

        // Store back after world collision
        // We need to store them before inter-ball collision detection and
        // response because we need the updated position and velocity
        XMStoreFloat3(&pos[i], posVec);
        XMStoreFloat3(&vel[i], velVec);

        // Query nearby balls
        hash->Query(pos, i, 2.0f * radius);

        // Inter-ball collision detection and response between
        // ball i and all nearby balls returned by the spatial hash query.
        for (int32_t nr = 0; nr < hash->GetQuerySize(); nr++)
        {
            int32_t j = hash->GetQueryIds()[nr];

            // Reload position and velocity of ball i (might have been modified
            // by world bounds collision)
            posVec = XMLoadFloat3(&pos[i]);
            velVec = XMLoadFloat3(&vel[i]);

            // Load position and velocity of nearby ball j
            XMVECTOR otherPosVec = XMLoadFloat3(&pos[j]);
            XMVECTOR velOtherVec = XMLoadFloat3(&vel[j]);

            // Calculate squared lenght of pos[i] - pos[j]
            XMVECTOR normalVec = XMVectorSubtract(posVec, otherPosVec);
            float d2 = XMVectorGetX(XMVector3LengthSq(normalVec));

            // Are balls overlapping?
            if (d2 > 0.0f && d2 < minDist * minDist)
            {
                // Compute distance between ball centers
                float d = std::sqrt(d2);

                // Normalize
                normalVec = XMVector3Normalize(normalVec);

                // Separate balls to avoid interpenetration:
				// This correction moves each ball away from the other along the collision normal.
				// (minDist - d) computes the overlap distance between the two balls.
				// Multiplying by 0.5f splits the correction equally between the two balls,
				// so each ball is moved half the overlap distance in opposite directions along 
				// the normal vector. This ensures both balls are separated fairly and prevents
				// them from remaining interpenetrated.
                float corr = (minDist - d) * 0.5f;
                posVec = XMVectorAdd(posVec, XMVectorScale(normalVec, corr));
                otherPosVec = XMVectorSubtract(otherPosVec, XMVectorScale(normalVec, corr));

				// Project each velocity vector onto the collision normal to obtain the scalar
				// component of the velocity in the direction of the collision. This represents
				// how fast each ball is moving towards or away from the other along the normal.
                float vi = XMVectorGetX(XMVector3Dot(velVec, normalVec));
                float vj = XMVectorGetX(XMVector3Dot(velOtherVec, normalVec));

				// Swap the normal components of the velocities between the two balls. This simulates
				// an elastic collision by exchanging velocity components along the collision normal, while
				// leaving the tangential components unchanged.
				//
				// CONCRETE EXAMPLE:
				//
				//   normalVec     = (1, 0, 0)
				//   velVec        = ( 3, 2, 0)
				//   velOtherVec   = (-1, 5, 0)
				//
				// Projection onto the collision normal:
				//   vi = dot(( 3, 2, 0), (1, 0, 0)) =  3
				//   vj = dot((-1, 5, 0), (1, 0, 0)) = -1
				//
				// Elastic collision (equal masses):
				//   the first sphere must end up with normal component vj = -1
				//   the second sphere must end up with normal component vi = 3
				//
				// Correction applied to the first sphere:
				//   vj - vi = -1 - 3 = -4
				//   velVec_new = (3,2,0) + (-4)*(1,0,0) = (-1,2,0)
				//
				// Correction applied to the second sphere:
				//   vi - vj = 3 - (-1) = 4
				//   velOther_new = (-1,5,0) + 4*(1,0,0) = (3,5,0)
				//
				// Result:
				//   - the normal (X) velocity components are exchanged
				//   - the tangential (Y,Z) components remain unchanged
                velVec = XMVectorAdd(velVec, XMVectorScale(normalVec, vj - vi));
                velOtherVec = XMVectorAdd(velOtherVec, XMVectorScale(normalVec, vi - vj));

                // Store back updated positions and velocities
                XMStoreFloat3(&pos[i], posVec);
                XMStoreFloat3(&pos[j], otherPosVec);
                XMStoreFloat3(&vel[i], velVec);
                XMStoreFloat3(&vel[j], velOtherVec);

                inCollision = true;
            }
        }

        // Update material based on collision
        if (showCollisions && inCollision)
            SetSphereMaterial(i, true);
        else
            SetSphereMaterial(i, false);

        // Update transform
        UpdateSpherePosition(i);
    }
}
