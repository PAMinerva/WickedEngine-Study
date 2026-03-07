#pragma once

// This file is shared between C++ and HLSL.
// Keep types compatible with both languages.

#ifdef __cplusplus
#include <cstdint>
#endif

// Thread group size for compute shaders (matches HTML: 64)
#define CLOTH_THREAD_GROUP_SIZE     256

// UAV buffer slots
#define CLOTH_SLOT_POS              0
#define CLOTH_SLOT_PREV_POS         1
#define CLOTH_SLOT_VEL              2
#define CLOTH_SLOT_INV_MASS         3
#define CLOTH_SLOT_CONST_IDS        4
#define CLOTH_SLOT_REST_LENGTHS     5
#define CLOTH_SLOT_CORRECTIONS      6   // uint[3*N], stores float bits via CAS atomics
#define CLOTH_SLOT_NORMALS          7
#define CLOTH_SLOT_TRI_IDS          8
#define CLOTH_SLOT_NORM_ACCUM       9   // uint[3*N], stores float bits via CAS atomics
#define CLOTH_SLOT_TRI_DIST         10

// Constant buffer slots
#define CLOTH_CB_SLOT               0
#define CLOTH_CB_RAYCAST_SLOT       1

// SimParams struct — matches HTML layout exactly (16 floats = 64 bytes)
//
// HTML layout:
//   [0]  dt            : f32
//   [1]  gravX         : f32
//   [2]  gravY         : f32
//   [3]  gravZ         : f32
//   [4]  sphereCX      : f32
//   [5]  sphereCY      : f32
//   [6]  sphereCZ      : f32
//   [7]  sphereR       : f32
//   [8]  jacobiScale   : f32
//   [9]  numParticles  : u32
//   [10] firstConstraint    : u32
//   [11] numConstraintsInPass : u32
//   [12] dragParticleNr: i32   (-1 = no drag)
//   [13] dragPosX      : f32
//   [14] dragPosY      : f32
//   [15] dragPosZ      : f32
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

// RaycastConstants — matches HTML RayParams layout (8 floats = 32 bytes)
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
