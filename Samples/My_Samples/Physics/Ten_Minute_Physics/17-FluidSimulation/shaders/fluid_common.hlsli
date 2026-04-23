// ============================================================================
// fluid_common.hlsli
// Shared declarations for all fluid simulation compute shaders.
// ============================================================================

#include "shader_interop_fluid.h"

// ---------------------------------------------------------------------------
// Constant buffer
// ---------------------------------------------------------------------------
cbuffer FluidCB : register(b0)
{
    FluidConstants cb;
};

// ---------------------------------------------------------------------------
// UAV bindings (3D textures) — simulation fields
// ---------------------------------------------------------------------------
RWTexture3D<float>  texU          : register(u0);  // velocity X (staggered: gridX+1, gridY, gridZ)
RWTexture3D<float>  texV          : register(u1);  // velocity Y (staggered: gridX, gridY+1, gridZ)
RWTexture3D<float>  texW          : register(u2);  // velocity Z (staggered: gridX, gridY, gridZ+1)
RWTexture3D<float>  texPressure   : register(u3);  // pressure (cell-centered)
RWTexture3D<float>  texDivergence : register(u4);  // divergence (cell-centered)
RWTexture3D<float>  texSmoke      : register(u5);  // smoke density (cell-centered)
RWTexture3D<uint>   texSolid      : register(u6);  // 0=solid, 1=fluid

// Temp buffers for advection ping-pong
RWTexture3D<float>  texU_temp     : register(u7);
RWTexture3D<float>  texV_temp     : register(u8);
RWTexture3D<float>  texW_temp     : register(u9);
RWTexture3D<float>  texSmoke_temp : register(u10);

// Temperature field (cell-centered)
RWTexture3D<float>  texTemp       : register(u11);
RWTexture3D<float>  texTemp_temp  : register(u12);

// Curl / vorticity vector (cell-centered, xyz = curl components)
RWTexture3D<float4> texCurl       : register(u13);

// ---------------------------------------------------------------------------
// Helper: check if cell (i,j,k) is within fluid domain and is not solid
// ---------------------------------------------------------------------------
bool isFluid(int i, int j, int k)
{
    if (i < 0 || i >= cb.gridX || j < 0 || j >= cb.gridY || k < 0 || k >= cb.gridZ)
        return false;
    return texSolid[int3(i, j, k)] != 0;
}

bool isInBounds(int i, int j, int k)
{
    return (i >= 0 && i < cb.gridX && j >= 0 && j < cb.gridY && k >= 0 && k < cb.gridZ);
}

// ---------------------------------------------------------------------------
// Helper: number of non-solid neighbors (for pressure solver)
// Each neighbor contributes 1 if it's fluid, 0 if solid/boundary.
// ---------------------------------------------------------------------------
float countFluidNeighbors(int i, int j, int k)
{
    float s = 0.0;
    s += (i > 0            && texSolid[int3(i-1, j, k)] != 0) ? 1.0 : 0.0;  // left
    s += (i < cb.gridX - 1 && texSolid[int3(i+1, j, k)] != 0) ? 1.0 : 0.0;  // right
    s += (j > 0            && texSolid[int3(i, j-1, k)] != 0) ? 1.0 : 0.0;  // bottom
    s += (j < cb.gridY - 1 && texSolid[int3(i, j+1, k)] != 0) ? 1.0 : 0.0;  // top
    s += (k > 0            && texSolid[int3(i, j, k-1)] != 0) ? 1.0 : 0.0;  // back
    s += (k < cb.gridZ - 1 && texSolid[int3(i, j, k+1)] != 0) ? 1.0 : 0.0;  // front
    return s;
}

// ============================================================================
// HELPER: sampleVelocity
// Manual trilinear interpolation of one staggered velocity component at an
// arbitrary continuous position in grid-index space.
//
// ── What this helper does ────────────────────────────────────────────────────
// This function returns a single velocity component (u, v, or w) evaluated
// at a continuous position `pos`.
// The requested component is selected by `component`:
//   0 = u
//   1 = v
//   2 = w
//
// Important: this helper does NOT reconstruct the full velocity vector.
// It samples only one scalar staggered field: texU, texV, or texW.
//
// ── Why interpolation is needed ──────────────────────────────────────────────
// During Semi-Lagrangian advection, the backtraced departure point almost
// never lands exactly on a stored grid sample.  It usually lies somewhere
// between neighboring samples, so we must interpolate.
//
// Because the simulation uses a MAC grid, each component lives on a
// different staggered layout:
//   u lives on X faces: (i,     j+0.5, k+0.5)
//   v lives on Y faces: (i+0.5, j,     k+0.5)
//   w lives on Z faces: (i+0.5, j+0.5, k)
//
// Therefore, the first thing this helper does is choose the dimensions of
// the correct staggered texture for the requested component.
//
// ── What happens step by step ────────────────────────────────────────────────
// 1. Select the correct source texture and its dimensions.
// 2. Clamp `pos` so that interpolation never reads outside valid bounds.
// 3. Compute:
//      i0 = floor(pos)    -> lower integer corner
//      f  = pos - i0      -> fractional offset inside that local cell
//      i1 = i0 + 1        -> upper corner (clamped to texture bounds)
// 4. Read the 8 scalar samples at the corners surrounding `pos`.
// 5. Blend them with standard trilinear interpolation.
//
// The 8 values are conventionally named:
//   c000, c100, c010, c110, c001, c101, c011, c111
// where each bit tells whether we use the lower (0) or upper (1) corner
// along x, y, and z.
//
// ── Numerical example ────────────────────────────────────────────────────────
// Suppose we call:
//
//   sampleVelocity(float3(10.2, 5.7, 3.1), 1)
//
// This asks for the v component at a continuous point in grid-index space.
//
// Then:
//   i0 = (10, 5, 3)
//   f  = (0.2, 0.7, 0.1)
//   i1 = (11, 6, 4)
//
// So the function reads these 8 neighboring v samples:
//
//   texV[10,5,3], texV[11,5,3],
//   texV[10,6,3], texV[11,6,3],
//   texV[10,5,4], texV[11,5,4],
//   texV[10,6,4], texV[11,6,4]
//
// and blends them according to the fractional offsets:
//   20% along x, 70% along y, 10% along z.
//
// So this helper is simply:
//   "sample one staggered scalar velocity field continuously by trilinear
//    interpolation."
//
// ── Relation to advectVelCS ────────────────────────────────────────────────
// advectVelCS decides which position and which component are needed.
// sampleVelocity only performs the actual continuous lookup of that one
// component at that one position.
// ============================================================================
float sampleVelocity(float3 pos, int component)
{
    // Determine the texture dimensions for this component
    int3 dims;
    if (component == 0)
        dims = int3(cb.gridX + 1, cb.gridY, cb.gridZ);
    else if (component == 1)
        dims = int3(cb.gridX, cb.gridY + 1, cb.gridZ);
    else
        dims = int3(cb.gridX, cb.gridY, cb.gridZ + 1);

    // Clamp pos to valid range for interpolation (so we don't read outside the texture)
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    // Compute lower and upper corners of the cell containing pos
    int3 i0 = int3(floor(pos));
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));

	// Compute fractional offset of pos within the cell
    float3 f = pos - float3(i0);

    // Clamp i0 and i1
    i0 = max(i0, int3(0, 0, 0));

	// Read the 8 corner samples for the requested component
    float c000, c100, c010, c110, c001, c101, c011, c111;

    if (component == 0)
    {
        c000 = texU[int3(i0.x, i0.y, i0.z)];
        c100 = texU[int3(i1.x, i0.y, i0.z)];
        c010 = texU[int3(i0.x, i1.y, i0.z)];
        c110 = texU[int3(i1.x, i1.y, i0.z)];
        c001 = texU[int3(i0.x, i0.y, i1.z)];
        c101 = texU[int3(i1.x, i0.y, i1.z)];
        c011 = texU[int3(i0.x, i1.y, i1.z)];
        c111 = texU[int3(i1.x, i1.y, i1.z)];
    }
    else if (component == 1)
    {
        c000 = texV[int3(i0.x, i0.y, i0.z)];
        c100 = texV[int3(i1.x, i0.y, i0.z)];
        c010 = texV[int3(i0.x, i1.y, i0.z)];
        c110 = texV[int3(i1.x, i1.y, i0.z)];
        c001 = texV[int3(i0.x, i0.y, i1.z)];
        c101 = texV[int3(i1.x, i0.y, i1.z)];
        c011 = texV[int3(i0.x, i1.y, i1.z)];
        c111 = texV[int3(i1.x, i1.y, i1.z)];
    }
    else
    {
        c000 = texW[int3(i0.x, i0.y, i0.z)];
        c100 = texW[int3(i1.x, i0.y, i0.z)];
        c010 = texW[int3(i0.x, i1.y, i0.z)];
        c110 = texW[int3(i1.x, i1.y, i0.z)];
        c001 = texW[int3(i0.x, i0.y, i1.z)];
        c101 = texW[int3(i1.x, i0.y, i1.z)];
        c011 = texW[int3(i0.x, i1.y, i1.z)];
        c111 = texW[int3(i1.x, i1.y, i1.z)];
    }

    // Trilinear blend
    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ============================================================================
// Helper: sample smoke at an arbitrary position (cell-centered scalar field)
//
// Smoke is conceptually stored at cell centers:
//   cell (i,j,k) has its smoke value at physical position
//   (i+0.5, j+0.5, k+0.5)
// in grid-index space.
//
// However, this helper uses the simpler "texture index" convention:
//   texSmoke[i,j,k] is sampled as if it lived at coordinate (i,j,k).
//
// That means the caller must do the half-cell conversion when needed:
//   center-space position  (bx, by, bz)
// becomes
//   sampleSmoke position   (bx-0.5, by-0.5, bz-0.5)
//
// Example in 1D:
//   cell 3 has center at x = 3.5
//   if backtracing gives bx = 2.8, then the point lies between
//   cell-center 2.5 (cell 2) and cell-center 3.5 (cell 3).
//
// To sample correctly with this helper we pass:
//   x = 2.8 - 0.5 = 2.3
//
// Then:
//   i0 = floor(2.3) = 2
//   i1 = 3
//   f  = 0.3
//
// so the result is:
//   0.7 * smoke[2] + 0.3 * smoke[3]
//
// which is exactly the interpolation we want at physical position 2.8.
// Without the -0.5 shift, the helper would interpolate at the wrong place,
// introducing a systematic half-cell offset.
//
// As with sampleVelocity, the interpolation itself is standard trilinear:
// 1. clamp pos to valid bounds,
// 2. compute i0 = floor(pos), i1 = i0 + 1, and fractional offset f,
// 3. read the 8 neighboring smoke samples,
// 4. blend them according to f.
//
// So this helper is simply:
//   "sample the smoke scalar field continuously by trilinear interpolation,
//    using texture-index coordinates."
// ============================================================================
float sampleSmoke(float3 pos)
{
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    int3 i0 = int3(floor(pos));
    float3 f = pos - float3(i0);
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));
    i0 = max(i0, int3(0, 0, 0));

    float c000 = texSmoke[int3(i0.x, i0.y, i0.z)];
    float c100 = texSmoke[int3(i1.x, i0.y, i0.z)];
    float c010 = texSmoke[int3(i0.x, i1.y, i0.z)];
    float c110 = texSmoke[int3(i1.x, i1.y, i0.z)];
    float c001 = texSmoke[int3(i0.x, i0.y, i1.z)];
    float c101 = texSmoke[int3(i1.x, i0.y, i1.z)];
    float c011 = texSmoke[int3(i0.x, i1.y, i1.z)];
    float c111 = texSmoke[int3(i1.x, i1.y, i1.z)];

    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ---------------------------------------------------------------------------
// Helper: sample temperature at an arbitrary position (cell-centered field)
//
// This uses the same coordinate convention as sampleSmoke():
// texTemp[i,j,k] is sampled as if it lived at coordinate (i,j,k), so callers
// that work in cell-center coordinates (i+0.5,j+0.5,k+0.5) must shift by
// -0.5 before calling this helper.
// ---------------------------------------------------------------------------
float sampleTemperature(float3 pos)
{
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    int3 i0 = int3(floor(pos));
    float3 f = pos - float3(i0);
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));
    i0 = max(i0, int3(0, 0, 0));

    float c000 = texTemp[int3(i0.x, i0.y, i0.z)];
    float c100 = texTemp[int3(i1.x, i0.y, i0.z)];
    float c010 = texTemp[int3(i0.x, i1.y, i0.z)];
    float c110 = texTemp[int3(i1.x, i1.y, i0.z)];
    float c001 = texTemp[int3(i0.x, i0.y, i1.z)];
    float c101 = texTemp[int3(i1.x, i0.y, i1.z)];
    float c011 = texTemp[int3(i0.x, i1.y, i1.z)];
    float c111 = texTemp[int3(i1.x, i1.y, i1.z)];

    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ---------------------------------------------------------------------------
// Helper: average velocity at cell center (i,j,k)
// For staggered grid, the velocity at cell center is the average of the
// two face values in each direction.
// ---------------------------------------------------------------------------
float3 avgVelAtCenter(int i, int j, int k)
{
    float u = 0.5 * (texU[int3(i, j, k)] + texU[int3(i + 1, j, k)]); // average of left and right faces
    float v = 0.5 * (texV[int3(i, j, k)] + texV[int3(i, j + 1, k)]); // average of bottom and top faces
    float w = 0.5 * (texW[int3(i, j, k)] + texW[int3(i, j, k + 1)]); // average of back and front faces

	// Return the average velocity
    return float3(u, v, w);
}

// ---------------------------------------------------------------------------
// Helper: sample velocity from _temp textures (for MacCormack correction).
// Same logic as sampleVelocity but reads from ping-pong temp buffers.
// ---------------------------------------------------------------------------
float sampleVelocityTemp(float3 pos, int component)
{
    int3 dims;
    if (component == 0)
        dims = int3(cb.gridX + 1, cb.gridY, cb.gridZ);
    else if (component == 1)
        dims = int3(cb.gridX, cb.gridY + 1, cb.gridZ);
    else
        dims = int3(cb.gridX, cb.gridY, cb.gridZ + 1);

    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    int3 i0 = int3(floor(pos));
    float3 f = pos - float3(i0);
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));
    i0 = max(i0, int3(0, 0, 0));

    float c000, c100, c010, c110, c001, c101, c011, c111;

    if (component == 0)
    {
        c000 = texU_temp[int3(i0.x, i0.y, i0.z)];
        c100 = texU_temp[int3(i1.x, i0.y, i0.z)];
        c010 = texU_temp[int3(i0.x, i1.y, i0.z)];
        c110 = texU_temp[int3(i1.x, i1.y, i0.z)];
        c001 = texU_temp[int3(i0.x, i0.y, i1.z)];
        c101 = texU_temp[int3(i1.x, i0.y, i1.z)];
        c011 = texU_temp[int3(i0.x, i1.y, i1.z)];
        c111 = texU_temp[int3(i1.x, i1.y, i1.z)];
    }
    else if (component == 1)
    {
        c000 = texV_temp[int3(i0.x, i0.y, i0.z)];
        c100 = texV_temp[int3(i1.x, i0.y, i0.z)];
        c010 = texV_temp[int3(i0.x, i1.y, i0.z)];
        c110 = texV_temp[int3(i1.x, i1.y, i0.z)];
        c001 = texV_temp[int3(i0.x, i0.y, i1.z)];
        c101 = texV_temp[int3(i1.x, i0.y, i1.z)];
        c011 = texV_temp[int3(i0.x, i1.y, i1.z)];
        c111 = texV_temp[int3(i1.x, i1.y, i1.z)];
    }
    else
    {
        c000 = texW_temp[int3(i0.x, i0.y, i0.z)];
        c100 = texW_temp[int3(i1.x, i0.y, i0.z)];
        c010 = texW_temp[int3(i0.x, i1.y, i0.z)];
        c110 = texW_temp[int3(i1.x, i1.y, i0.z)];
        c001 = texW_temp[int3(i0.x, i0.y, i1.z)];
        c101 = texW_temp[int3(i1.x, i0.y, i1.z)];
        c011 = texW_temp[int3(i0.x, i1.y, i1.z)];
        c111 = texW_temp[int3(i1.x, i1.y, i1.z)];
    }

    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ---------------------------------------------------------------------------
// Helper: sample smoke from _temp texture (for MacCormack correction)
// ---------------------------------------------------------------------------
float sampleSmokeTemp(float3 pos)
{
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    int3 i0 = int3(floor(pos));
    float3 f = pos - float3(i0);
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));
    i0 = max(i0, int3(0, 0, 0));

    float c000 = texSmoke_temp[int3(i0.x, i0.y, i0.z)];
    float c100 = texSmoke_temp[int3(i1.x, i0.y, i0.z)];
    float c010 = texSmoke_temp[int3(i0.x, i1.y, i0.z)];
    float c110 = texSmoke_temp[int3(i1.x, i1.y, i0.z)];
    float c001 = texSmoke_temp[int3(i0.x, i0.y, i1.z)];
    float c101 = texSmoke_temp[int3(i1.x, i0.y, i1.z)];
    float c011 = texSmoke_temp[int3(i0.x, i1.y, i1.z)];
    float c111 = texSmoke_temp[int3(i1.x, i1.y, i1.z)];

    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ---------------------------------------------------------------------------
// Helper: sample temperature from _temp texture (for MacCormack correction)
// ---------------------------------------------------------------------------
float sampleTemperatureTemp(float3 pos)
{
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));

    int3 i0 = int3(floor(pos));
    float3 f = pos - float3(i0);
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));
    i0 = max(i0, int3(0, 0, 0));

    float c000 = texTemp_temp[int3(i0.x, i0.y, i0.z)];
    float c100 = texTemp_temp[int3(i1.x, i0.y, i0.z)];
    float c010 = texTemp_temp[int3(i0.x, i1.y, i0.z)];
    float c110 = texTemp_temp[int3(i1.x, i1.y, i0.z)];
    float c001 = texTemp_temp[int3(i0.x, i0.y, i1.z)];
    float c101 = texTemp_temp[int3(i1.x, i0.y, i1.z)];
    float c011 = texTemp_temp[int3(i0.x, i1.y, i1.z)];
    float c111 = texTemp_temp[int3(i1.x, i1.y, i1.z)];

    float c00 = lerp(c000, c100, f.x);
    float c10 = lerp(c010, c110, f.x);
    float c01 = lerp(c001, c101, f.x);
    float c11 = lerp(c011, c111, f.x);
    float c0 = lerp(c00, c10, f.y);
    float c1 = lerp(c01, c11, f.y);
    return lerp(c0, c1, f.z);
}

// ---------------------------------------------------------------------------
// Helper: get min/max of 8 corners of a cell-centered field at given position.
// Used for MacCormack clamping to prevent oscillations.
// ---------------------------------------------------------------------------
void getMinMaxSmoke(float3 pos, out float mn, out float mx)
{
    // Read from _temp (forward SL results) — read-only during MacCormack, no race.
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));
    int3 i0 = max(int3(floor(pos)), int3(0, 0, 0));
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));

    float v0 = texSmoke_temp[int3(i0.x, i0.y, i0.z)];
    float v1 = texSmoke_temp[int3(i1.x, i0.y, i0.z)];
    float v2 = texSmoke_temp[int3(i0.x, i1.y, i0.z)];
    float v3 = texSmoke_temp[int3(i1.x, i1.y, i0.z)];
    float v4 = texSmoke_temp[int3(i0.x, i0.y, i1.z)];
    float v5 = texSmoke_temp[int3(i1.x, i0.y, i1.z)];
    float v6 = texSmoke_temp[int3(i0.x, i1.y, i1.z)];
    float v7 = texSmoke_temp[int3(i1.x, i1.y, i1.z)];

    mn = min(min(min(v0, v1), min(v2, v3)), min(min(v4, v5), min(v6, v7)));
    mx = max(max(max(v0, v1), max(v2, v3)), max(max(v4, v5), max(v6, v7)));
}

void getMinMaxTemp(float3 pos, out float mn, out float mx)
{
    // Read from _temp (forward SL results) — read-only during MacCormack, no race.
    int3 dims = int3(cb.gridX, cb.gridY, cb.gridZ);
    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));
    int3 i0 = max(int3(floor(pos)), int3(0, 0, 0));
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));

    float v0 = texTemp_temp[int3(i0.x, i0.y, i0.z)];
    float v1 = texTemp_temp[int3(i1.x, i0.y, i0.z)];
    float v2 = texTemp_temp[int3(i0.x, i1.y, i0.z)];
    float v3 = texTemp_temp[int3(i1.x, i1.y, i0.z)];
    float v4 = texTemp_temp[int3(i0.x, i0.y, i1.z)];
    float v5 = texTemp_temp[int3(i1.x, i0.y, i1.z)];
    float v6 = texTemp_temp[int3(i0.x, i1.y, i1.z)];
    float v7 = texTemp_temp[int3(i1.x, i1.y, i1.z)];

    mn = min(min(min(v0, v1), min(v2, v3)), min(min(v4, v5), min(v6, v7)));
    mx = max(max(max(v0, v1), max(v2, v3)), max(max(v4, v5), max(v6, v7)));
}

void getMinMaxVel(float3 pos, int component, out float mn, out float mx)
{
    // Read from _temp (forward SL results) — read-only during MacCormack, no race.
    int3 dims;
    if (component == 0)      dims = int3(cb.gridX + 1, cb.gridY, cb.gridZ);
    else if (component == 1) dims = int3(cb.gridX, cb.gridY + 1, cb.gridZ);
    else                     dims = int3(cb.gridX, cb.gridY, cb.gridZ + 1);

    pos = clamp(pos, float3(0.0, 0.0, 0.0), float3(dims) - float3(1.0, 1.0, 1.0));
    int3 i0 = max(int3(floor(pos)), int3(0, 0, 0));
    int3 i1 = min(i0 + int3(1, 1, 1), dims - int3(1, 1, 1));

    float v0, v1, v2, v3, v4, v5, v6, v7;
    if (component == 0) {
        v0 = texU_temp[int3(i0.x,i0.y,i0.z)]; v1 = texU_temp[int3(i1.x,i0.y,i0.z)];
        v2 = texU_temp[int3(i0.x,i1.y,i0.z)]; v3 = texU_temp[int3(i1.x,i1.y,i0.z)];
        v4 = texU_temp[int3(i0.x,i0.y,i1.z)]; v5 = texU_temp[int3(i1.x,i0.y,i1.z)];
        v6 = texU_temp[int3(i0.x,i1.y,i1.z)]; v7 = texU_temp[int3(i1.x,i1.y,i1.z)];
    } else if (component == 1) {
        v0 = texV_temp[int3(i0.x,i0.y,i0.z)]; v1 = texV_temp[int3(i1.x,i0.y,i0.z)];
        v2 = texV_temp[int3(i0.x,i1.y,i0.z)]; v3 = texV_temp[int3(i1.x,i1.y,i0.z)];
        v4 = texV_temp[int3(i0.x,i0.y,i1.z)]; v5 = texV_temp[int3(i1.x,i0.y,i1.z)];
        v6 = texV_temp[int3(i0.x,i1.y,i1.z)]; v7 = texV_temp[int3(i1.x,i1.y,i1.z)];
    } else {
        v0 = texW_temp[int3(i0.x,i0.y,i0.z)]; v1 = texW_temp[int3(i1.x,i0.y,i0.z)];
        v2 = texW_temp[int3(i0.x,i1.y,i0.z)]; v3 = texW_temp[int3(i1.x,i1.y,i0.z)];
        v4 = texW_temp[int3(i0.x,i0.y,i1.z)]; v5 = texW_temp[int3(i1.x,i0.y,i1.z)];
        v6 = texW_temp[int3(i0.x,i1.y,i1.z)]; v7 = texW_temp[int3(i1.x,i1.y,i1.z)];
    }

    mn = min(min(min(v0, v1), min(v2, v3)), min(min(v4, v5), min(v6, v7)));
    mx = max(max(max(v0, v1), max(v2, v3)), max(max(v4, v5), max(v6, v7)));
}

