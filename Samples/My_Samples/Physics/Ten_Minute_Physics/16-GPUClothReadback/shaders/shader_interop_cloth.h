#pragma once

// This file is shared between C++ and HLSL.
// Keep types compatible with both languages.

#ifdef __cplusplus
#include <cstdint>
#endif

// Thread group size for compute shaders
#define CLOTH_THREAD_GROUP_SIZE     256

// UAV buffer slots
#define CLOTH_SLOT_POS              0
#define CLOTH_SLOT_PREV_POS         1
#define CLOTH_SLOT_VEL              2
#define CLOTH_SLOT_INV_MASS         3
#define CLOTH_SLOT_CONST_IDS        4
#define CLOTH_SLOT_REST_LENGTHS     5
#define CLOTH_SLOT_CORRECTIONS      6
#define CLOTH_SLOT_NORMALS          7
#define CLOTH_SLOT_TRI_IDS          8
#define CLOTH_SLOT_NORM_ACCUM       9
#define CLOTH_SLOT_TRI_DIST         10

// Constant buffer slots
#define CLOTH_CB_SLOT               0
#define CLOTH_CB_RAYCAST_SLOT       1

// SimParams struct
struct ClothSimConstants
{
    float dt;
    float gravX;
    float gravY;
    float gravZ;

    float sphereCX;
    float sphereCY;
    float sphereCZ;
    float sphereR;

    float jacobiScale;
    int numParticles;
    int firstConstraint;
    int numConstraintsInPass;

    int dragParticleNr;       // -1 if no drag active
    float dragPosX;
    float dragPosY;
    float dragPosZ;
};

// RaycastConstants
struct RaycastConstants
{
    float origX;
    float origY;
    float origZ;
    float _pad0;

    float dirX;
    float dirY;
    float dirZ;
    float _pad1;
};
