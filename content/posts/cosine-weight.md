---
title: "Cosine Weight"
date: 2022-08-25T09:48:54-04:00
draft: true
---

When doing ray tracing it's important to sample the hemisphere around the shading point with the correct distribution to avoid bias. 

**Spherical sampling (hemisphere)**

The naive solution to generate points on a sphere is to input 2 random numbers between 0 and 1 directly into the sphere equation.

``` csharp
float theta = 2 * Mathf.PI * rnd2;
float phi = Mathf.PI * 0.5f * rnd1;

p = new Vector3(Mathf.Sin(phi) * Mathf.Cos(theta), Mathf.Sin(phi) * Mathf.Sin(theta), Mathf.Cos(phi));
```

![image.png](/cosine-weight/.attachments/image-617e6190-9f28-442c-8664-72c08b3cc3fb.png)

This method is incorrect because the area covered by samples is not constant. Samples clump to the poles, and have a much smaller area there compared to samples near the hemisphere.

