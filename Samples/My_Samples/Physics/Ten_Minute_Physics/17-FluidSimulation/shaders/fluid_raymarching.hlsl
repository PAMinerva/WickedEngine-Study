// ============================================================================
// fluid_raymarching.hlsl
//
// Screen-sized compute shader for volumetric rendering of the fluid simulation.
//
// This shader does not draw triangle geometry. Instead, it takes the 3D data
// produced by the simulation (smoke density, pressure, velocity), looks at that
// data from the current camera point of view, and converts it into a 2D RGBA
// image that can later be composited over the already rendered scene.
//
// The basic idea is:
//   1. each thread is responsible for one output pixel
//   2. for that pixel, reconstruct the corresponding camera ray in world space
//   3. find where that ray enters and exits the fluid domain box
//   4. advance through the volume in small steps
//   5. at each step, sample the 3D textures and convert the sample into
//      color + opacity
//   6. accumulate those contributions into the final pixel color
//
// Front-to-back alpha compositing:
// This shader accumulates samples starting from the part of the volume nearest
// to the camera and then moving farther away. Each sample contributes some
// color and some opacity. If the accumulated opacity becomes high enough, the
// ray can stop early because farther samples would no longer be visible.
// In other words, if the first few samples contain dense smoke, the pixel
// quickly becomes more opaque, so smoke deeper in the volume contributes less
// and less. This is the volumetric equivalent of saying "nearer smoke hides
// farther smoke".
//
// Per-pixel jitter for anti-banding:
// "Banding" means visible stripes or layers in places where the volume should
// look smooth. Instead of seeing a gradual change of density, the eye notices
// repeated steps, as if the smoke were built from stacked slices.
// Why it happens without jitter:
// if every pixel started at exactly the same ray distance and used exactly the
// same step spacing, neighboring pixels would tend to sample the same density
// layers at almost the same positions. That makes many adjacent pixels produce
// very similar accumulated values, then suddenly switch together when the next
// sample layer is crossed. The result is usually a visible pattern of bands
// aligned with the sampling pattern.
// Concrete example:
// suppose the step size is 0.10 and there is a thin dense region around
// distances [1.05, 1.09] along the ray.
// Without jitter:
//   pixel A samples at 1.00, 1.10, 1.20, ...
//   pixel B samples at 1.00, 1.10, 1.20, ...
// Both pixels completely miss that dense region, so both become too dark.
// A nearby region might instead be hit by many neighboring pixels at once,
// making them all brighten together. That collective miss/hit behavior is what
// creates visible bands rather than a smooth transition.
// With jitter:
// each pixel keeps the same ray direction and the same step spacing, but the
// first sample is shifted by a small per-pixel offset inside the first step.
// For example:
//   pixel A samples at 1.03, 1.13, 1.23, ...
//   pixel B samples at 1.08, 1.18, 1.28, ...
// Now one pixel might partially hit the dense region while the other misses it,
// and across many pixels those small differences break up the repeated pattern.
// The image becomes noisier at a tiny scale, but much less visibly layered.
// Note: this tradeoff can also show as a temporary blocky/grainy edge look in smoke
// boundaries when some neighboring pixels hit density and others miss it.
//
// Adaptive step count proportional to ray length:
// A short ray crossing only a small corner of the domain does not need as many
// samples as a long ray crossing almost the whole box. This shader estimates
// the ray length inside the domain and chooses the number of samples from that.
// This keeps the result detailed enough without paying the maximum cost for
// every pixel.
//
// Configurable smoke color:
// In smoke mode, density mainly controls opacity, while the visible tint comes
// from the smoke color passed by C++. That lets the same simulation be shown
// with different visual styles without changing the simulation data itself.
//
// Why ray marching here:
// The simulation produces 3D fields (smoke density, pressure, velocity), not a
// mesh surface that can be drawn directly with the normal raster pipeline. To
// render smoke from the camera point of view, we reconstruct one viewing ray
// per screen pixel, march through the domain volume, sample the 3D textures,
// and integrate their contribution into a final 2D color/alpha result.
//
// The output of this shader is therefore not "the scene" itself, but a 2D
// overlay image of the fluid that will later be composited over the already
// rendered 3D scene in SampleRenderPath::Compose().
//
// Rendering approach inspired by "3D Fluid Simulation and Rendering" (Iacchini)
// and Wicked Engine's volumetric cloud rendering pattern.
// ============================================================================

#include "shader_interop_fluid.h"

// ---------------------------------------------------------------------------
// Constant buffer
// ---------------------------------------------------------------------------
cbuffer RenderCB : register(b0)
{
    FluidRenderConstants rc;
};

// ---------------------------------------------------------------------------
// SRV bindings (read-only 3D textures from simulation)
// ---------------------------------------------------------------------------
Texture3D<float>  smokeTex     : register(t0);
Texture3D<float>  pressureTex  : register(t1);
Texture3D<float>  velUTex      : register(t2);
Texture3D<float>  velVTex      : register(t3);
Texture3D<float>  velWTex      : register(t4);

// Sampler for trilinear interpolation (Wicked Engine global sampler slot)
SamplerState sampler_linear_clamp : register(s100);

// UAV: output render target
RWTexture2D<float4> outputTex : register(u0);

// ---------------------------------------------------------------------------
// Ray / AABB intersection using the slab method
// ---------------------------------------------------------------------------
//
// AABB = Axis-Aligned Bounding Box:
// a box whose faces are parallel to the X, Y, and Z axes.
//
// Goal
// -----
// Given a ray
//
//   P(t) = rayOrigin + t * rayDir
//
// and an AABB defined by [boxMin, boxMax], we want to find for which values
// of t the ray is inside the box.
//
// If such values exist, they form an interval:
//
//   [tNear, tFar]
//
// where:
// - tNear is where the ray enters the box
// - tFar  is where the ray exits  the box
//
// So the whole idea is:
// each axis tells us when the ray is inside that axis-aligned slab,
// and the box hit is where all 3 slab conditions are true together.
//
// ---------------------------------------------------------------------------
// Why it is called the "slab" method
// ---------------------------------------------------------------------------
//
// A slab is the infinite region between two parallel planes.
//
// In 2D, the analogous idea is the strip between two parallel lines.
// This is easier to picture before going back to 3D.
//
// Vertical strip:
//
//        |       |
//        |       |
//        |       |
//        |       |
//
// Horizontal strip:
//
//        -----------
//
//        -----------
//
// If we call the first one the X slab and the second one the Y slab,
// their intersection is the finite rectangle where both conditions hold:
//
//        |       |
//      --|-------|--
//        |       |
//      --|-------|--
//        |       |
//
// In words:
// - inside the vertical strip   => x is between 2 limits
// - inside the horizontal strip => y is between 2 limits
// - inside both at once         => point is inside the rectangle
//
// In 3D the idea is exactly the same:
//
// - X slab: between x = boxMin.x and x = boxMax.x
// - Y slab: between y = boxMin.y and y = boxMax.y
// - Z slab: between z = boxMin.z and z = boxMax.z
//
// Their intersection is exactly the box.
//
// Being inside the box means satisfying all 3 coordinate constraints
// simultaneously:
//
//   boxMin.x <= P.x <= boxMax.x
//   boxMin.y <= P.y <= boxMax.y
//   boxMin.z <= P.z <= boxMax.z
//
// ---------------------------------------------------------------------------
// 1) Solve the crossing times on one axis
// ---------------------------------------------------------------------------
//
// Start with one axis, for example X.
//
// Along the ray, the X coordinate is
//
//   P.x(t) = rayOrigin.x + t * rayDir.x
//
// The ray crosses the two X planes when
//
//   rayOrigin.x + t * rayDir.x = boxMin.x
//   rayOrigin.x + t * rayDir.x = boxMax.x
//
// Solving for t gives
//
//   t0.x = (boxMin.x - rayOrigin.x) / rayDir.x
//   t1.x = (boxMax.x - rayOrigin.x) / rayDir.x
//
// These are simply the two times at which the ray reaches the two X planes.
//
// The same reasoning applies to Y and Z.
//
// In vector form we compute all 3 axis pairs at once:
//
//   invDir = 1 / rayDir
//   t0 = (boxMin - rayOrigin) * invDir
//   t1 = (boxMax - rayOrigin) * invDir
//
// Using invDir replaces division with multiplication and keeps the code short.
//
// ---------------------------------------------------------------------------
// 2) For each axis, determine entry and exit
// ---------------------------------------------------------------------------
//
// On a single axis, the ray crosses the 2 planes that bound that slab.
// One crossing is where the ray enters the slab, the other is where it exits.
//
// The important detail is that the order depends on the ray direction.
//
// Think only about the X axis.
//
// Case A: rayDir.x > 0
// Increasing t moves toward larger x, so the ray reaches boxMin.x first
// and boxMax.x second:
//
//   x -------------------------------------------------------------------->
//
//                 boxMin.x                    boxMax.x
//                    |---------------------------|
//                         X slab interior
//
//   ray motion as t increases:  ---------------------------->
//
//   first crossing  = enter slab  = t0.x
//   second crossing = exit slab   = t1.x
//
// Case B: rayDir.x < 0
// Increasing t moves toward smaller x, so the ray reaches boxMax.x first
// and boxMin.x second:
//
//   x -------------------------------------------------------------------->
//
//                 boxMin.x                    boxMax.x
//                    |---------------------------|
//                         X slab interior
//
//   ray motion as t increases:  <----------------------------
//
//   first crossing  = enter slab  = t1.x
//   second crossing = exit slab   = t0.x
//
// Therefore we cannot assume that:
// - t0 is always entry
// - t1 is always exit
//
// We fix that by reordering component-wise:
//
//   tmin = min(t0, t1)
//   tmax = max(t0, t1)
//
// Here, "Component-wise" means:
// - compare t0.x with t1.x
// - compare t0.y with t1.y
// - compare t0.z with t1.z
//
// In other words, min/max are applied separately to X, Y, and Z:
//
//   tmin.x = min(t0.x, t1.x)
//   tmin.y = min(t0.y, t1.y)
//   tmin.z = min(t0.z, t1.z)
//
//   tmax.x = max(t0.x, t1.x)
//   tmax.y = max(t0.y, t1.y)
//   tmax.z = max(t0.z, t1.z)
//
// We do this because, on a given axis, t0 and t1 are just the two plane
// crossing times. Depending on the ray direction, either one may happen first.
//
// After reordering:
// - tmin is the per-axis entry time
// - tmax is the per-axis exit time
//
// For example, on the X axis:
// - tmin.x is the entry time for the X slab
// - tmax.x is the exit  time for the X slab
//
// and likewise for Y and Z.
//
// ---------------------------------------------------------------------------
// 3) Intersect the 3 slab intervals
// ---------------------------------------------------------------------------
//
// Each slab tells us for which t values the ray is inside that slab:
//
//   X slab: [tmin.x ------------------------------- tmax.x]
//   Y slab:        [tmin.y ------------------------------- tmax.y]
//   Z slab:                 [tmin.z ------- tmax.z]
//
// Put differently, each slab says:
//
//   "the ray is inside me only for this interval of t"
//
// If we draw all 3 intervals on the same t axis, the box hit is just the
// part where all 3 intervals overlap:
//
//   t ------------------------------------------------------------------->
//
//   X:      [===============================]
//           tmin.x                     tmax.x
//
//   Y:             [===============================]
//                   tmin.y                     tmax.y
//
//   Z:                      [=======]
//                           tmin.z tmax.z
//
//   box:                    [=======]
//                           tNear   tFar
//
// Notice what happened:
// - the overlap starts at the LAST entry among the 3 slabs
// - the overlap ends   at the FIRST exit  among the 3 slabs
//
// Why?
// - Before that last entry, the ray is still missing at least one slab,
//   so it cannot yet be inside the box.
// - After that first exit, the ray has already left at least one slab,
//   so it is no longer inside the box.
//
// Therefore:
//
//   tNear = max(tmin.x, tmin.y, tmin.z)
//   tFar  = min(tmax.x, tmax.y, tmax.z)
//
// Good mnemonic:
//
//   latest entry wins
//   earliest exit wins
//
// ---------------------------------------------------------------------------
// 4) Decide whether the overlap exists
// ---------------------------------------------------------------------------
//
// There are 2 possible cases.
//
// Case A: the intervals overlap
//
//   t --------------------------------------------------------------->
//
//   X:      [===========================]
//   Y:            [=============================]
//   Z:                 [===========]
//
//   box:               [===========]
//                     tNear      tFar
//
// Here:
//
//   tNear <= tFar
//
// so there is a common overlap.
// That means there exists a portion of the ray that is inside all 3 slabs
// at the same time, so the ray intersects the box.
//
// Case B: the intervals do not overlap
//
//   t --------------------------------------------------------------->
//
//   X:      [===========]
//                       ^
//   Y:                  |  [================]
//                       |  ^
//   Z:            [================]
//                       |  |
//                       |  |
//                       |  |
//                earliest  latest
//                  exit    entry
//                 = tFar  = tNear
//
//   box:           --- no overlap ---
//
// Since the earliest exit is before the latest entry,
//
//   tNear > tFar
//
// so the 3 slab intervals never overlap all at once.
// There is no t for which the ray is inside all 3 slabs together,
// so the ray misses the box.
//
// In this particular function we do not explicitly test hit/miss here.
// We simply return the candidate interval
//
//   [max(tNear, 0), tFar]
//
// and leave it to the caller to decide whether that interval is usable.
//
// Typical caller-side tests are:
//
// - returned.x <= returned.y   -> valid hit interval
// - returned.x <  returned.y   -> valid interval with positive length
//
// The second form is often used in ray marching, where a zero-length interval
// (tangent touch) is not useful.
//
// ---------------------------------------------------------------------------
// 5) Why we clamp tNear to 0
// ---------------------------------------------------------------------------
//
// For ray marching we usually care only about the forward half-ray:
//
//   t >= 0
//
// If the ray origin is already inside the box, the box entry happens
// behind the ray origin, so tNear is negative.
//
// In that case, the useful starting point is not "behind the camera",
// but the ray origin itself. So we clamp:
//
//   max(tNear, 0)
//
// and keep tFar as the forward exit distance.
//
// ---------------------------------------------------------------------------
// 6) Small numerical example
// ---------------------------------------------------------------------------
//
// Let
//
//   boxMin = (1, 2, 3)
//   boxMax = (5, 6, 7)
//
// and
//
//   rayOrigin = (0, 4, 4)
//   rayDir    = (1, 0.5, -1)
//
// Then
//
//   t0 = (boxMin - rayOrigin) / rayDir = ( 1, -4,  1)
//   t1 = (boxMax - rayOrigin) / rayDir = ( 5,  4, -3)
//
// Reordering per axis:
//
//   tmin = (1, -4, -3)
//   tmax = (5,  4,  1)
//
// Therefore:
//
//   tNear = max(1, -4, -3) = 1
//   tFar  = min(5,  4,  1) = 1
//
// The overlap collapses to a single t value.
// So the ray only touches the box at one point.
//
// This function would return:
//
//   float2(max(1, 0), 1) = float2(1, 1)
//
// A caller that accepts tNear <= tFar would treat this as a hit.
// A caller that requires a strictly positive interval (tNear < tFar)
// would reject it as a tangent touch.
//
// ---------------------------------------------------------------------------
// Summary
// ---------------------------------------------------------------------------
//
// 1. Compute the two ray/plane crossing times on each axis
// 2. Reorder them into per-axis entry/exit times
// 3. Take the latest entry  => tNear = max(all tmin)
// 4. Take the earliest exit => tFar  = min(all tmax)
// 5. Clamp the start to the forward ray: max(tNear, 0)
// 6. Let the caller decide whether the returned interval is valid
// ---------------------------------------------------------------------------
float2 rayBoxIntersect(float3 rayOrigin, float3 rayDir, float3 boxMin, float3 boxMax)
{
    float3 invDir = 1.0 / rayDir;
    float3 t0 = (boxMin - rayOrigin) * invDir;
    float3 t1 = (boxMax - rayOrigin) * invDir;
    float3 tmin = min(t0, t1);
    float3 tmax = max(t0, t1);
    float tNear = max(max(tmin.x, tmin.y), tmin.z);
    float tFar  = min(min(tmax.x, tmax.y), tmax.z);
    return float2(max(tNear, 0.0), tFar);
}

// ---------------------------------------------------------------------------
// Pseudo-random hash for per-pixel jitter (prevents banding artifacts)
// ---------------------------------------------------------------------------
float hash21(float2 p)
{
    p = frac(p * float2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return frac(p.x * p.y);
}

// ---------------------------------------------------------------------------
// Scientific color map
//
// This helper converts a scalar value in [0, 1] into a false-color RGB value.
// It is useful for diagnostic rendering of quantities such as pressure or
// speed, where we want the viewer to distinguish low, medium, and high values
// at a glance.
//
// The mapping is piecewise-linear and goes through these color bands:
//
//   0.00  -> blue
//   0.25  -> cyan
//   0.50  -> green
//   0.75  -> yellow
//   1.00  -> red
//
// In each quarter of the [0, 1] range, one channel stays fixed while another
// channel ramps up or down linearly:
//
// - [0.00, 0.25): blue stays at 1, green ramps up      -> blue to cyan
// - [0.25, 0.50): green stays at 1, blue ramps down    -> cyan to green
// - [0.50, 0.75): green stays at 1, red ramps up       -> green to yellow
// - [0.75, 1.00]: red stays at 1, green ramps down     -> yellow to red
//
// Values outside [0, 1] are clamped first, so the function always returns a
// valid color on that scale.
// ---------------------------------------------------------------------------
float3 getSciColor(float val)
{
    val = clamp(val, 0.0, 1.0);
    float r = 0.0, g = 0.0, b = 0.0;

    if (val < 0.25)
    {
        b = 1.0;
        g = val / 0.25;
    }
    else if (val < 0.5)
    {
        g = 1.0;
        b = 1.0 - (val - 0.25) / 0.25;
    }
    else if (val < 0.75)
    {
        g = 1.0;
        r = (val - 0.5) / 0.25;
    }
    else
    {
        r = 1.0;
        g = 1.0 - (val - 0.75) / 0.25;
    }

    return float3(r, g, b);
}

// ---------------------------------------------------------------------------
// Front-to-back alpha compositing helper
//
// This function merges one newly sampled volumetric contribution into the
// running result for the current pixel.
//
// Meaning of the inputs:
// - sampleColor = RGB color contributed by the current sample
// - sampleAlpha = opacity contributed by the current sample
// - weight      = optional extra multiplier used to scale that sample
//
// Meaning of the inout values:
// - color = RGB accumulated so far for this pixel
// - alpha = total opacity accumulated so far for this pixel
//
// At the start of ray marching we have
//
//   color = (0, 0, 0)
//   alpha = 0
//
// which means:
// - no light/color has been accumulated yet
// - the pixel is still fully transparent
//
// As we march from the front of the volume toward the back, each new sample is
// composited over what has already been accumulated.
//
// First, the sample opacity is optionally scaled by weight and clamped:
//
//   sampleAlpha = min(sampleAlpha * weight, 1)
//
// Then we compute
//
//   t = sampleAlpha * (1 - alpha)
//
// The factor (1 - alpha) is the fraction of the pixel that is still
// transparent. So:
// - if alpha is small, much of the pixel is still uncovered and the new sample
//   can contribute strongly
// - if alpha is already close to 1, most of the pixel is already opaque and the
//   new sample contributes very little
//
// The new sample is then merged by
//
//   color += sampleColor * t
//   alpha += t
//
// So this helper implements the idea:
// "nearer volumetric samples are composited first, and farther ones only fill
// in whatever transparency is still left."
// ---------------------------------------------------------------------------
void accumSample(float3 sampleColor, float sampleAlpha, float weight,
                 inout float3 color, inout float alpha)
{
    sampleAlpha = min(sampleAlpha * weight, 1.0);
    float t = sampleAlpha * (1.0 - alpha);
    color += sampleColor * t;
    alpha += t;
}

// ---------------------------------------------------------------------------
//                            Ray-marching
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Reconstructing A World-Space Ray From A Pixel
// ---------------------------------------------------------------------------
//
// Each compute thread corresponds to one pixel of the output texture, so
// dtid.xy are the integer pixel coordinates in screen space:
//
//   pixel = (dtid.x, dtid.y)
//
// Screen-space pixel coordinates are discrete integer indices, but for ray
// reconstruction we usually want the ray through the pixel center, not its
// top-left corner. That is why we add 0.5:
//
//   pixelCenter = (dtid.xy + 0.5)
//
// We then convert that pixel position to UV coordinates:
//
//   uv = pixelCenter / float2(screenWidth, screenHeight)
//
// UV means normalized 2D texture/screen coordinates:
// - U is the horizontal coordinate
// - V is the vertical coordinate
// - both are typically in the range [0, 1]
//
// So:
// - uv = (0, 0) is the top-left corner of the screen
// - uv = (1, 1) is the bottom-right corner
// - uv = (0.5, 0.5) is the screen center
//
// Next we remap UVs to NDC (Normalized Device Coordinates), whose X and Y
// range is [-1, 1]:
//
//   ndc.x = uv.x * 2 - 1
//   ndc.y = 1 - uv.y * 2 = (1 - uv.y) * 2 - 1
//
// The Y flip is needed because screen/texture coordinates usually grow
// downward, while NDC Y grows upward.
//
// At this point we have a homogeneous clip/NDC-style position:
//
//   ndc = float4(ndc.x, ndc.y, z, 1)
//
// Here z = 0.5 is just a point somewhere between near and far planes.
// We do not need an exact surface depth here: any point along the camera
// frustum for this pixel is enough, because the ray direction is obtained
// from the camera position toward that point.
//
// To move from clip/NDC-style coordinates back into world space, we multiply
// by the inverse view-projection matrix:
//
//   worldPt = mul(ndc, invVP)
//
// This gives a homogeneous world-space position, so we must divide by W:
//
//   worldPt.xyz = worldPt.xyz / worldPt.w
//
// After perspective divide, worldPt.xyz is a 3D point in world space lying
// on the camera ray corresponding to this pixel.
//
// The ray origin is simply the camera world position:
//
//   rayOrigin = cameraPosition
//
// The ray direction is the normalized vector from the camera to that
// unprojected world-space point:
//
//   rayDir = normalize(worldPt.xyz - rayOrigin)
//
// In short:
//   pixel coords -> UV [0,1] -> NDC [-1,1] -> unproject with invVP
//   -> perspective divide -> world-space point -> normalized camera-to-point vector
//
// Web reference for the same derivation in native article form:
// Scratchapixel, "Generating Camera Rays":
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays.html
// It walks through the same pixel-center -> normalized coordinates -> screen/NDC -> camera/world-space ray construction path used here.
// Naming note:
// Scratchapixel uses "NDC" for the intermediate [0, 1] normalized pixel
// coordinates. Here uv is that intermediate quantity, while ndc is the
// already remapped [-1, 1] value used for unprojection.
//
// ---------------------------------------------------------------------------
// Intersecting The Camera Ray With The Fluid Domain
// ---------------------------------------------------------------------------
//
// Once the world-space ray for the current pixel has been reconstructed, the
// next question is: which part of that ray actually passes through the fluid?
//
// The simulation lives inside an axis-aligned box (the fluid domain), defined
// by boxMin and boxMax. So before marching through the volume, we intersect the
// camera ray with that box.
//
// Conceptually, rayBoxIntersect() returns the interval
//
//   [tNear, tFar]
//
// along the ray
//
//   P(t) = rayOrigin + t * rayDir
//
// for which the ray is inside the domain box.
//
// This interval has a very important meaning:
// - tNear is the first useful point where the forward ray is inside the box
// - tFar  is the point where the ray leaves the box
//
// So if the interval is valid, the segment
//
//   P(t),  t in [tNear, tFar]
//
// is exactly the portion of the viewing ray that crosses the simulated volume.
//
// If the interval is empty or collapses to a point, there is no useful segment
// to march through. In that case the current pixel receives no fluid
// contribution, and the shader writes transparent black.
//
// In other words:
// - ray misses box        -> no volumetric contribution
// - ray touches box only  -> still no useful marching segment here
// - ray crosses box       -> march only inside that bounded segment
//
// This is why the ray-box test is done before choosing sample count or taking
// any volume samples: it tells us whether there is anything to render for this
// pixel at all, and if so, where the marching interval begins and ends.
//
// ---------------------------------------------------------------------------
// Choosing How Densely To March Along The Ray
// ---------------------------------------------------------------------------
//
// Once we know the ray segment inside the fluid domain, the next step is to
// decide how many samples to take along that segment.
//
// A ray that crosses only a tiny corner of the domain does not need as many
// samples as a ray that travels through a large portion of the box. So instead
// of always using a fixed number of steps for every pixel, this shader adapts
// the sample count to the ray length inside the domain.
//
// First we measure the length of the valid marching interval:
//
//   rayLength = tFar - tNear
//
// This is the distance, in world space, that the ray travels inside the fluid.
//
// Next we estimate a representative grid cell size in world units. The domain
// has a world-space size
//
//   domainSize = boxMax - boxMin
//
// and the simulation grid has resolution
//
//   (gridX, gridY, gridZ)
//
// The shader takes the largest domain dimension and divides it by the largest
// grid resolution:
//
//   cellWorldSize = max(domainSize) / max(gridResolution)
//
// This gives an approximate voxel size in world space.
//
// Once we have that estimate, we can ask:
// "roughly how many cells does this ray cross?"
//
//   rayLength / cellWorldSize
//
// This shader then multiplies that value by 2, so it aims for about two
// samples per voxel crossed by the ray:
//
//   fSamples = rayLength / cellWorldSize * 2
//
// The result is clamped to the range [1, rc.maxSteps], so:
// - very short rays still take at least one sample
// - very long rays do not exceed the configured maximum cost
//
// Then we compute:
//
//   nSamples = floor(fSamples)
//   stepSize = rayLength / fSamples
//   stepVec  = rayDir * stepSize
//
// Their roles are:
// - nSamples = number of regularly spaced samples taken by the main marching loop
// - stepSize = world-space distance between consecutive samples
// - stepVec  = 3D offset added each iteration to move along the ray
//
// In other words:
// - short in-box segment -> fewer samples
// - long in-box segment  -> more samples
//
// This keeps the marching work roughly proportional to how much volume the
// pixel ray actually traverses, instead of spending the same cost on every
// pixel regardless of path length.
//
// ---------------------------------------------------------------------------
// Handling Fractional Sample Counts
// ---------------------------------------------------------------------------
//
// The adaptive sample count is first computed as a real-valued quantity:
//
//   fSamples
//
// but the main marching loop can only execute an integer number of full steps:
//
//   nSamples = floor(fSamples)
//
// So if fSamples is not an integer, the loop covers only the integer part and
// leaves a fractional remainder:
//
//   frac(fSamples)
//
// This matters because nearby pixels often get slightly different fSamples.
// Their rays are similar, but not identical, so they can have slightly
// different in-box ray lengths and therefore slightly different adaptive
// sample counts.
//
// For example, one pixel might get
//
//   fSamples = 10.99
//
// while a neighboring pixel gets
//
//   fSamples = 11.01
//
// Without any correction, the first pixel would take 10 full samples and the
// second would take 11, even though their rays are almost the same. That
// sudden one-sample jump can create visible discontinuities.
//
// To reduce that effect, after the main loop the shader may take one extra
// tail sample and scale its opacity contribution by
//
//   weight = frac(fSamples)
//
// This acts like a partial final step:
// - weight = 0.0 -> no extra contribution
// - weight = 0.5 -> half of one last sample
// - weight = 0.9 -> almost a full extra sample
//
// So the total contribution changes more smoothly as fSamples varies from
// pixel to pixel, reducing banding or popping caused by discrete step-count
// changes.
//
// In this implementation, that fractional tail sample is only applied in smoke
// mode, where opacity continuity is the most visually important.
//
// ---------------------------------------------------------------------------
// Choosing The Initial Sample Position With Per-Pixel Jitter
// ---------------------------------------------------------------------------
//
// Once the step spacing has been chosen, we still need to decide where the
// first sample should be placed along the valid marching interval.
//
// A naive approach would be to always start exactly at
//
//   P(tNear) = rayOrigin + tNear * rayDir
//
// and then advance by one fixed step each iteration.
//
// That works, but it tends to produce visible banding artifacts. The reason is
// that neighboring pixels often have similar rays, so if every pixel starts
// sampling at exactly the start of its marching interval and then advances with
// its own regular step spacing, many nearby pixels end up sampling similar
// density layers in a similarly aligned pattern.
//
// Instead, this shader applies a small per-pixel pseudo-random offset, called
// jitter, to the starting position. The offset is generated from the integer
// pixel coordinates:
//
//   jitter = hash21(dtid.xy)
//
// so each pixel gets a stable pseudo-random value in [0, 1).
//
// That jitter is then converted into an offset inside the first marching step:
//
//   startT = tNear + jitter * stepSize
//
// and the initial sample position becomes
//
//   pos = rayOrigin + startT * rayDir
//
// Different pixels generally have different rays and different step sizes.
// Jitter does not change those per-pixel values; it only offsets where each
// pixel begins sampling along its own ray.
//
// The result is that the regular sampling pattern becomes less visually aligned
// across neighboring pixels. Instead of large coherent bands, the error is
// broken up into much finer noise, which is usually much less objectionable.
//
// So the purpose of jitter here is not to change the physical meaning of the
// ray march, but to improve visual quality by reducing structured sampling
// artifacts.
//
// ---------------------------------------------------------------------------
// Marching Through The Volume And Building The Output Pixel
// ---------------------------------------------------------------------------
//
// Once the starting point and the marching step have been chosen, the shader
// advances through the valid ray segment one sample at a time.
//
// At each iteration of the main loop:
// 1. the current world-space position is examined
// 2. special regions such as obstacle or emitter spheres may contribute
// 3. the position is converted from world space to texture-space UVW
// 4. the relevant 3D field is sampled
// 5. that sample is converted into color and opacity
// 6. the result is composited into the running pixel result
//
// The accumulation is done front-to-back, meaning from the part of the volume
// nearest to the camera toward the part farther away.
//
// This ordering matters because nearer volumetric samples should reduce the
// visibility of farther ones. The exact compositing rule used here is described
// in the dedicated comment above accumSample().
//
// So:
// - if little opacity has been accumulated yet, a new sample can contribute a lot
// - if the ray is already nearly opaque, a new sample contributes very little
//
// This is the volumetric analogue of the idea that nearer smoke hides smoke
// behind it.
//
// In this shader, the running result is stored in
//
//   color
//   alpha
//
// These two values represent the RGBA result being built for the current pixel
// of the 2D output texture. Once the ray march is finished, that pixel value is
// written to outputTex and later composited over the already rendered 3D scene.
//
// After each iteration, the sample position is advanced by
//
//   pos += stepVec
//
// so the shader walks forward through the volume until all regularly spaced
// samples have been processed, or until the accumulated opacity is high enough
// that further work is no longer useful.
//
[numthreads(FLUID_THREADS_2D, FLUID_THREADS_2D, 1)]
void main(uint3 dtid : SV_DispatchThreadID)
{
    if ((int)dtid.x >= rc.screenWidth || (int)dtid.y >= rc.screenHeight)
        return;

    // Reconstruct the world-space ray for this pixel.
    // See the detailed note above for the full screen -> UV -> NDC -> world path.
    float3 rayOrigin = float3(rc.cameraPosX, rc.cameraPosY, rc.cameraPosZ);
    float2 uv = (float2(dtid.xy) + 0.5) / float2(rc.screenWidth, rc.screenHeight);
    float4 ndc = float4(uv.x * 2.0 - 1.0, (1.0 - uv.y) * 2.0 - 1.0, 0.5, 1.0);
    float4 worldPt = mul(ndc, rc.invVP);
    worldPt.xyz /= worldPt.w;
    float3 rayDir = normalize(worldPt.xyz - rayOrigin);

    // Compute the forward ray segment that lies inside the simulation box.
    // rayBoxIntersect() returns [tNear, tFar], that is, the entry/exit interval
    // after clamping the start to the forward ray. If that interval is empty or
    // collapses to a point (tNear >= tFar), there is no useful segment to march
    // through, so this pixel gets no fluid contribution.
	// See the detailed note above raBoxIntersect() for the full derivation of
	// how ray/box intersection works and why it returns the interval it does.
    float3 boxMin = float3(rc.domainMinX, rc.domainMinY, rc.domainMinZ);
    float3 boxMax = float3(rc.domainMaxX, rc.domainMaxY, rc.domainMaxZ);

    float2 tHit = rayBoxIntersect(rayOrigin, rayDir, boxMin, boxMax);

    if (tHit.x >= tHit.y)
    {
        outputTex[dtid.xy] = float4(0, 0, 0, 0);
        return;
    }

    // Adaptive step count: take roughly two samples per voxel crossed by the
    // ray, capped by rc.maxSteps. This keeps short rays cheap and long rays
    // detailed enough for smooth volumetric accumulation.
    float3 domainSize = boxMax - boxMin;
    float rayLength = tHit.y - tHit.x;
    float maxGridDim = (float)max(rc.gridX, max(rc.gridY, rc.gridZ));
    float cellWorldSize = max(domainSize.x, max(domainSize.y, domainSize.z)) / maxGridDim;
    float fSamples = rayLength / cellWorldSize * 2.0;
    fSamples = clamp(fSamples, 1.0, rc.maxSteps);

    int nSamples = (int)floor(fSamples);
    float stepSize = rayLength / fSamples;
    float3 stepVec = rayDir * stepSize;

    // Initial sample position along the valid in-box ray segment.
    // jitter shifts that starting point within the first step so neighboring
    // pixels do not all sample the same regularly spaced layers.
    float jitter = hash21(float2(dtid.xy));
    float3 pos = rayOrigin + rayDir * (tHit.x + jitter * stepSize);

    // Obstacle sphere
    float3 obsCenter = float3(rc.obsCenterX, rc.obsCenterY, rc.obsCenterZ);
    float obsR2 = rc.obsRadius * rc.obsRadius;

    // Smoke emitter sphere
    float3 srcCenter = float3(rc.smokeSrcRenderX, rc.smokeSrcRenderY, rc.smokeSrcRenderZ);
    float srcR2 = rc.smokeSrcRenderRadius * rc.smokeSrcRenderRadius;

    // Smoke color from constants
    float3 smokeColor = float3(rc.smokeColorR, rc.smokeColorG, rc.smokeColorB);

    // Front-to-back accumulation. Each step samples the relevant field from the
    // 3D textures, converts that sample into color/opacity, and composites it
    // into the running RGBA result for this pixel.
    float3 color = float3(0, 0, 0);
    float alpha = 0.0;

    for (int s = 0; s < nSamples; s++)
    {
        // Obstacle sphere test:
        // for the current sample position, compute the vector from the sphere
        // center to the sample. If its squared length is smaller than the
        // squared radius, this sample lies inside the obstacle sphere.
        float3 toObs = pos - obsCenter;
        if (dot(toObs, toObs) < obsR2)
        {
            // Use the center-to-surface direction as an approximate sphere
            // normal and apply a simple directional Lambert term so the sphere
            // appears shaded instead of flat.
            // L is the light direction: a fixed world-space direction pointing
            // from the surface toward an imagined distant light source.
            // We hard-code it as (0.3, 1.0, -0.5) just to get a pleasing
            // top/front-side lighting direction for this debug-style obstacle
            // rendering. It is normalized because the Lambert dot product
            // assumes unit-length directions, so dot(N, L) measures only the
            // cosine of the angle between normal and light, not their lengths.
            float3 N = normalize(toObs);
            float3 L = normalize(float3(0.3, 1.0, -0.5));
            float NdL = max(dot(N, L), 0.0);
            float3 sphereCol = float3(0.35, 0.35, 0.4) * (0.3 + 0.7 * NdL);

            // Render the obstacle as fully opaque at the first hit sample.
            // After compositing it, stop marching because anything farther
            // along the ray is hidden behind this solid sphere.
            accumSample(sphereCol, 1.0, 1.0, color, alpha);
            break;
        }

        // Emitter sphere test:
        // this is a small visual glow around the smoke source position, used
        // to make the emitter region visible in the final render.
        float3 toSrc = pos - srcCenter;
        float srcDist2 = dot(toSrc, toSrc);
        if (srcDist2 < srcR2)
        {
            // Convert squared distance to a normalized [0, 1] radial factor:
            // - fade = 0 at the center of the sphere
            // - fade = 1 at the sphere boundary
            //
            // The emitted color is strongest near the center and becomes
            // slightly dimmer toward the edge.
            float fade = srcDist2 / srcR2;
            float3 emitCol = float3(1.0, 0.5, 0.1) * (1.0 - 0.5 * fade);

            // Add only a small opacity contribution so the emitter appears as
            // a soft warm glow rather than an opaque solid object.
            accumSample(emitCol, 0.02, 1.0, color, alpha);
        }

        // Convert the current world-space sample position into texture-space
        // UVW coordinates for the simulation 3D textures.
        //
        // This works because the fluid domain in world space is exactly the
        // axis-aligned box [boxMin, boxMax], while the 3D textures are sampled
        // in normalized coordinates [0, 1]^3 over that same box.
        //
        // So:
        // - pos = boxMin maps to uvw = (0, 0, 0)
        // - pos = boxMax maps to uvw = (1, 1, 1)
        // - points in between map linearly in each axis
        //
        // In other words, this is just the box-to-unit-cube normalization
        //
        //   uvw = (pos - boxMin) / (boxMax - boxMin)
        //
        // The clamp keeps sampling slightly away from the exact texture border
        // to avoid edge/boundary artifacts during filtered 3D texture lookups.
        float3 uvw = (pos - boxMin) / domainSize;
        uvw = clamp(uvw, float3(0.001, 0.001, 0.001), float3(0.999, 0.999, 0.999));

        if (rc.renderMode == 0)
        {
            // SMOKE mode:
            // sample the smoke density field at the current UVW position.
            // This density says how much smoke is present locally in the volume.
            //
            // The sampled density is not used here to change hue; instead it is
            // converted into opacity. The visible tint comes from smokeColor,
            // while density controls how strongly this step contributes.
            //
            // rc.smokeAbsorption controls how aggressively smoke density turns
            // into opacity. stepSize matters too: a longer marching step means
            // this single sample stands in for a thicker segment of volume
            // along the ray, so it should contribute more opacity than a very
            // short step through the same local density.
            //
            // The final 0.1 has no special geometric meaning by itself.
            // It simply scales the opacity down after density, absorption, and
            // step length have been combined, so the smoke does not become too
            // opaque too quickly for the chosen rendering setup.
            float density = smokeTex.SampleLevel(sampler_linear_clamp, uvw, 0);
            float sAlpha = density * rc.smokeAbsorption * stepSize * 0.1;
            accumSample(smokeColor, sAlpha, 1.0, color, alpha);
        }
        else if (rc.renderMode == 1)
        {
            // PRESSURE mode: diagnostic visualization of the scalar field.
            // This is not physically rendered smoke, just a color-mapped view
            // of pressure values in the volume.
            float pressure = pressureTex.SampleLevel(sampler_linear_clamp, uvw, 0);
            float pNorm = clamp(pressure * 0.001 + 0.5, 0.0, 1.0);
            accumSample(getSciColor(pNorm), 0.15, 1.0, color, alpha);
        }
        else
        {
            // VELOCITY mode: diagnostic visualization of speed magnitude.
            float u = velUTex.SampleLevel(sampler_linear_clamp, uvw, 0);
            float v = velVTex.SampleLevel(sampler_linear_clamp, uvw, 0);
            float w = velWTex.SampleLevel(sampler_linear_clamp, uvw, 0);
            float speed = length(float3(u, v, w));
            float sNorm = clamp(speed * 0.2, 0.0, 1.0);
            accumSample(getSciColor(sNorm), 0.15, 1.0, color, alpha);
        }

        // Early exit: once the accumulated opacity is almost 1, the pixel is
        // already nearly fully covered, so farther samples would contribute
        // negligibly and are not worth evaluating.
        if (alpha > 0.99)
            break;

        // Advance from the current sample position to the next regularly
        // spaced sample along the ray.
        pos += stepVec;
    }

    // Optional fractional tail sample for non-integer adaptive sample counts.
    if (alpha < 0.99)
    {
        float weight = frac(fSamples);
        if (weight > 0.001)
        {
            float3 uvw = (pos - boxMin) / domainSize;
            uvw = clamp(uvw, float3(0.001, 0.001, 0.001), float3(0.999, 0.999, 0.999));

            if (rc.renderMode == 0)
            {
                float density = smokeTex.SampleLevel(sampler_linear_clamp, uvw, 0);
                float sAlpha = density * rc.smokeAbsorption * stepSize * 0.1;
                accumSample(smokeColor, sAlpha, weight, color, alpha);
            }
        }
    }

    // Write the final RGBA result for this pixel into the 2D output texture.
    // This texture is later composited over the already rendered scene.
    outputTex[dtid.xy] = float4(color, alpha);
}
