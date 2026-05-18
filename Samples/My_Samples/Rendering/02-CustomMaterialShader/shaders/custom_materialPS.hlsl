#define OBJECTSHADER_LAYOUT_COMMON
#include "../../../../../WickedEngine/shaders/objectHF.hlsli"

[earlydepthstencil]
float4 main(PixelInput input) : SV_TARGET
{
	ShaderMaterial material = GetMaterial();
	float4 uvsets = input.GetUVSets();

	write_mipmap_feedback(push.materialIndex, ddx_coarse(uvsets), ddy_coarse(uvsets));

	float stripeCount = max(1.0, asfloat(material.userdata.x));
	float warpAmount = asfloat(material.userdata.y);
	float scrollSpeed = asfloat(material.userdata.z);

	float time = GetFrame().time * scrollSpeed;
	float2 uv = uvsets.xy;
	float3 pos = input.GetPos3D();
	float3 normal = normalize(input.nor);
	float3 view = normalize(input.GetViewVector());

	float bridgeWave = sin((uv.x + sin(pos.y * 2.8 + time) * warpAmount) * stripeCount * 6.28318530718 + time * 2.0);
	float band = smoothstep(-0.05, 0.05, bridgeWave);

	float plankMask = smoothstep(0.46, 0.50, abs(frac(uv.y * 6.0) - 0.5));
	float3 white = float3(0.95, 0.95, 0.90);
	float3 red = float3(0.95, 0.03, 0.02);
	float3 green = float3(0.0, 0.80, 0.20);

	float nDotL = saturate(dot(normal, normalize(float3(-0.45, 0.85, -0.25))));
	float rim = pow(1.0 - saturate(dot(normal, -view)), 2.5);

	float3 color = lerp(white, red, band);
	color = lerp(color, green, plankMask * 0.85);
	color *= 0.30 + nDotL * 0.95;
	color += rim * float3(0.35, 0.55, 1.0);

	return saturateMediump(float4(color, 1.0));
}
