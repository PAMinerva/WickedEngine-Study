// cloth_common.hlsli
// Shared buffer and constant declarations for cloth simulation compute shaders.

#include "shader_interop_cloth.h"

// Constant buffer (matches ClothSimConstants, 64 bytes)
cbuffer ClothCB : register(b0)
{
    ClothSimConstants cb;
};

// ── Particle buffers (UAVs) ──
RWStructuredBuffer<float4> pos       : register(u0);   // xyz position, w unused
RWStructuredBuffer<float4> prevPos   : register(u1);
RWStructuredBuffer<float4> vel       : register(u2);
RWStructuredBuffer<float>  invMass   : register(u3);

// ── Constraint buffers ──
RWStructuredBuffer<uint>  constIds    : register(u4);   // flat pairs: [2*cNr] = id0, [2*cNr+1] = id1
RWStructuredBuffer<float> restLengths : register(u5);

// ── Jacobi corrections (CAS-loop float atomics, 3 floats per particle stored as uint bits) ──
RWStructuredBuffer<uint> corrections : register(u6);    // [3*pNr+0]=x, [3*pNr+1]=y, [3*pNr+2]=z

// ── Normals buffer ──
RWStructuredBuffer<float4> normals : register(u7);      // xyz normal, w unused

// ── Triangle buffers ──
RWStructuredBuffer<uint> triIds : register(u8);         // flat: [3*triNr+0/1/2]

// ── Normal accumulation (CAS-loop float atomics, 3 floats per particle stored as uint bits) ──
RWStructuredBuffer<uint> normAccum : register(u9);      // [3*pNr+0]=x, [3*pNr+1]=y, [3*pNr+2]=z

// ── Raycast output ──
RWStructuredBuffer<float> triDist : register(u10);

// ── CAS-loop float atomic add ──
// Atomically adds a float 'val' to buffer[idx] using InterlockedCompareExchange.
// The buffer stores float values as uint bits (asuint/asfloat reinterpretation).
#define ATOMIC_ADD_FLOAT(buf, idx, val)                          \
{                                                                \
    uint _cas_exp = (buf)[(idx)];                                \
    [loop] while (true)                                          \
    {                                                            \
        uint _cas_des = asuint(asfloat(_cas_exp) + (val));       \
        uint _cas_orig;                                          \
        InterlockedCompareExchange((buf)[(idx)],                 \
            _cas_exp, _cas_des, _cas_orig);                      \
        if (_cas_orig == _cas_exp) break;                        \
        _cas_exp = _cas_orig;                                    \
    }                                                            \
}
