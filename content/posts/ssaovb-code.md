---
title: "SSAO using Visibility Bitmasks"
date: 2023-05-11T09:48:54-04:00
draft: false
_build:
  render: true
  list: true
author: therrieno 
---

<!---
```csharp
// dirSign: 1 -> sampling right side of the slice
//         -1 -> sampling left side of the slice
// pos: position of current pixel
// angleVtoN: angle from V to N
// return: GI (RGB), AO (A)
float4 SampleDir(float dirSign, float3 pos, float3 viewDir, float angleVtoN, float thickness)
{
    uint bitmaskSlice = 0;
    float4 GIAO = 0;

    for (uint i = 0; i < sampleCount; i++)
    {
        float3 frontPos = GetSamplePos(pos, i);
        float3 frontDir = frontPos - pos;

        // Find backface dir by moving of thickness along -viewDir
        float3 backDir = frontDir - viewDir * thickness;

        // Project sample onto the unit circle and compute the angle relative to viewDir
        float2 angles = float2(dot(normalize(frontDir), viewDir), dot(normalize(backDir), viewDir));
        angles = FastAcos(angles);

        // Shift slice from viewDir to normal, clamp in [0..1]
        angles = saturate(((dirSign * -angles) - angleVtoN + HALF_PI) / PI);

        // Sampling direction inverts min/max angles
        angles = dirSign > 0 ? angles.yx : angles.xy;
        uint bitmask = OccludeSectors(angles.x, angles.y);

        GIAO.rgb += ComputeGI(bitmask, bitmaskSlice, frontPos);
        bitmaskSlice |= bitmask;
    }

    GIAO.rgb += ComputeAmbient(bitmaskSlice, normal, thickness);
    GIAO.a = ComputeAO(bitmaskSlice);

    return GIAO;
}
```

```csharp
uint OccludeSectors(float minAngle, float maxAngle){
  uint startAngle = minAngle * NUM_SECTOR;
  float angle = round((maxAngle-minAngle) * NUM_SECTOR);
  return uint(exp2(angle)-1) << startAngle;
}
```

```csharp
float ComputeAO(uint bf){
  return 1.0 - float(countbits(bf)) / NUM_SECTOR;
}
```

```csharp
float3 ComputeAmbient(uint bf, float3 T, float3 N){
  uint subRegionSize = NUM_SECTOR / NUM_SAMPLE;
  uint mask = exp2(subRegionSize) - 1;
  float3 light = 0;
  for (float i = 0; i < NUM_SAMPLE; i++){
    float hits = subRegionSize - countbits(bf & mask);
    // Generate a ray centered in the subregion
    float a = lerp(0, PI, (i + 0.5) / NUM_SAMPLE);
    float3 dir = normalize(T * cos(a) + N * sin(a));
    light += SampleAmbient(dir) * hits / NUM_SECTOR;
    bf >>= size;
  }
  return light;
```

```csharp
float3 ComputeGI(uint bf, uint bfSlice, float3 samplePos, 
  float3 pixelToSample, float3 N){
  // Compute the number of newly occluded sectors only
  float hits = countbits(bf & (~bfSlice));
  if(hits > 0){
    float3 l = SampleLightTexture(samplePos.xy);
    float3 L = normalize(pixelToSample);
    float NDotL = saturate(dot(N, L));
    // Continue if light is facing surface normal
    if(NDotL > 0.0){
      float3 Ln = GetNormal(samplePos.xy);
      float LnDotL = saturate(dot(Ln, -L));
      return hits / NUM_SECTOR * l * NDotL * LnDotL;
    }
  }
  return 0;
}
```
-->

We recently released the paper [Screen Space Indirect Lighting with Visibility Bitmask](https://arxiv.org/abs/2301.11376), a rendering technique based on GTAO that replaces the two horizon angles by a bit field representing the binary state (occluded / un-occluded) of N sectors uniformly distributed around the hemisphere slice. It allows light to pass behind surfaces of constant thickness while keeping the efficiency of horizon-based methods.

We wanted to share some HLSL code to make it easier to implement in a 3D applications and clarify some details from the paper. Here we take the Unity GTAO code from HDRP as a starting point, to show that our method can be integrated in an existing implementation with only a few modifications.

The GTAO algorithm starts by taking a fixed number (_dirCount_) of random directions per pixel (_dir_). Then a 2D circular slice is defined with the following parametrization:

- _V_: The view direction (direction from current pixel to camera)
- _normalVS_: The view-space normal at the current pixel
- _sliceN_: Unit vector that is perpendicular to the slice plane
- _projN_: The normal projected onto the slice plane (the normal is almost never aligned with the slice plane)
- _T_: The slice tangent (perpendicular to _V_ and _sliceN_)
- _N_: The angle in radians from _V_ to _projN_

```hlsl
float integral = 0;
float3 V = normalize(-positionVS);

for (int i = 0; i < dirCount; ++i)
{
    float2 dir = GetDirection(dispatchThreadId.xy, i);
    float3 normalVS = GetNormalVS(normalBufferData);
    float3 sliceN = normalize(cross(float3(dir.xy, 0.0f), V.xyz));
    float3 projN = normalVS - sliceN * dot(normalVS, sliceN);
    float projNLen = length(projN);
    float cosN = dot(projN / projNLen, V);

    float3 T = cross(V, sliceN);
    float N = -sign(dot(projN, T)) * GTAOFastAcos(cosN);

    // Find horizons
    float2 maxHorizons;
    maxHorizons.x = HorizonLoop(positionVS, V, rayStart, dir, offset, step, 0);
    maxHorizons.y = HorizonLoop(positionVS, V, rayStart, negDir, offset, step, 0);

    // Now we find the actual horizon angles
    maxHorizons.x = -GTAOFastAcos(maxHorizons.x);
    maxHorizons.y = GTAOFastAcos(maxHorizons.y);
    maxHorizons.x = N + max(maxHorizons.x - N, -HALF_PI);
    maxHorizons.y = N + min(maxHorizons.y - N, HALF_PI);
    integral += AnyIsNaN(maxHorizons) ? 1 : IntegrateArcCosWeighted(maxHorizons.x, maxHorizons.y, N, cosN);
}

integral /= dirCount;
```

For a given direction _i_, the _HorizonLoop()_ function finds the maximum elevation for the right part of the hemisphere (_maxHorizon.x_) and for the left part (_maxHorizon.y_).

With those maximum elevations, the algorithm then clamps the horizons to the hemisphere centered on the projected normal, and then computes the integral of the un-occluded regions on each side of the view vector using the _IntegrateArcCosWeighted()_ function.

To use the visibility bitmask approach, we slightly modify the _HorizonLoop()_ function to pass an _inout uint globalOccludedBitfield_ that will be marked with occluded sectors.

We can then compute the integral by counting the number of occluded sectors and dividing by the total amount of sectors. Note that we do not take the cosine weight into account in this case. 

```hlsl
// Find horizons
float2 maxHorizons;
uint globalOccludedBitfield = 0;
maxHorizons.x = HorizonLoop(positionVS, V, rayStart, dir, offset, step, 0, globalOccludedBitfield, 1, N);
maxHorizons.y = HorizonLoop(positionVS, V, rayStart, negDir, offset, step, 0, globalOccludedBitfield, -1, N);

#ifdef VISIBILITY_BITMASK
    integral += 1.0 - float(countbits(globalOccludedBitfield)) / float(SECTOR_COUNT);
#else
    // Now we find the actual horizon angles
    ...
#endif
```

Most of the work happens in the HorizonLoop() function. The GTAO implementation take a fixed number of samples (__AOStepCount_) along the current direction (_rayDir_).

For each sample found, the view space position is reconstructed (_samplePosVS_), and the direction from the current pixel to the current sample is computed (_deltaPos_). This is used to compute the elevation (_currHorizon_). The _UpdateHorizon()_ function is used to apply a falloff over the distance and keep the highest found elevation. 

```hlsl
float HorizonLoop(float3 positionVS, float3 V, float2 rayStart, float2 rayDir, 
    float rayOffset, float rayStep, int mipModifier)
{
    float maxHorizon = -1.0f;  // cos(pi)
    float t = rayOffset * rayStep + rayStep;

    const uint startWithLowerRes = min(max(0, _AOStepCount / 2 - 2), 3);
    for (uint i = 0; i < _AOStepCount; i++)
    {
        float2 samplePos = max(2, min(rayStart + t * rayDir, _AOBufferSize.xy - 2));

        // Find horizons at these steps:
        float sampleDepth = GetDepthSample(samplePos, i > startWithLowerRes);
        float3 samplePosVS = GetPositionVS(samplePos.xy, sampleDepth);

        float3 deltaPos = samplePosVS - positionVS;

        float deltaLenSq = dot(deltaPos, deltaPos);
        float currHorizon = dot(deltaPos, V) * rsqrt(deltaLenSq);
        maxHorizon = UpdateHorizon(maxHorizon, currHorizon, deltaLenSq);

        t += rayStep;
    }

    return maxHorizon;
}
```

With visibility bitmasks things are a bit different; We need to compute not only the front-face sample (_deltaPos_) but a back-face one too (_deltaPosBackface_) that is determined by the constant thickness value used (__Thickness_). 

The horizon angles for front and back are computed, shifted from viewDir to normal  and clamped in [0, 1] where 0 is the left side of the hemisphere, and 1 the right side. Then the _UpdatePartitions()_ function is used to set the bits that lay between the font and back angles in _globalOccludedBitfield_.

```hlsl
...

float3 deltaPos = samplePosVS - positionVS;

#ifdef VISIBILITY_BITMASK
    float2 frontBackHorizon;
    float3 deltaPosBackface = deltaPos - V * _Thickness;

    // Project sample onto the unit circle and compute the angle relative to V
    frontBackHorizon = float2(dot(normalize(deltaPos), V), dot(normalize(deltaPosBackface), V));
    frontBackHorizon = GTAOFastAcos(frontBackHorizon);

    // Shift sample from V to normal, clamp in [0..1]
    frontBackHorizon = saturate(((samplingDirection * -frontBackHorizon) - N + HALF_PI) / PI);

    // Sampling direction inverts min/max angles
    frontBackHorizon = samplingDirection >= 0 ? frontBackHorizon.yx : frontBackHorizon.xy;

    globalOccludedBitfield = UpdatePartitions(frontBackHorizon.x, frontBackHorizon.y, globalOccludedBitfield);
#else
    float deltaLenSq = dot(deltaPos, deltaPos);
    float currHorizon = dot(deltaPos, V) * rsqrt(deltaLenSq);
    maxHorizon = UpdateHorizon(maxHorizon, currHorizon, deltaLenSq);
#endif
```

With visibility bitmasks the _UpdateHorizon()_ function from GTAO is not needed anymore, because we don't need to apply any falloff! The constant thickness and the bitmask is enough to reproduce the attenuation over the distance in a plausible manner. It's similar to ray tracing AO that doesn't need any falloff heuristic either.

The _UpdateParitions()_ function takes as input minHorizon and maxHorizon that represent the position in the hemisphere with a value between 0 and 1. 

```hlsl
#define SECTOR_COUNT 32

uint UpdatePartitions(float minHorizon, float maxHorizon, uint globalOccludedBitfield)
{
    uint startHorizonInt = minHorizon * SECTOR_COUNT;
    uint angleHorizonInt = ceil((maxHorizon-minHorizon) * SECTOR_COUNT);
    uint angleHorizonBitfield = angleHorizonInt > 0 ? (0xFFFFFFFF >> (SECTOR_COUNT-angleHorizonInt)) : 0;
    uint currentOccludedBitfield = angleHorizonBitfield << startHorizonInt;
    return globalOccludedBitfield | currentOccludedBitfield;
}
```

Sectors get activated from minHorizon to maxHorizon depending on the chosen rounding function:

- ceil: Sample needs to at least touch a sector to activate it
- round: Sample needs to cover at least half a sector to activate it
- floor: Sample needs to cover the entire sector to activate it 

That's it! That pretty much covers the essential parts related to ambient occlusion using visibility bimasks. We left out the GI and the ambient sampling parts, maybe for a later post!


