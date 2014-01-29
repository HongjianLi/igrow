#pragma once
#ifndef IGROW_MAT3_HPP
#define IGROW_MAT3_HPP

#include <cassert>
#include "vec3.hpp"

/// Constructs a rotation matrix from a normalized axis and the cosine value of an angle.
array<double, 9> vec3_to_mat3(const array<double, 3>& a, const double c)
{
	if (zero(a))
	{
		assert(eq(c, 1) || eq(c, -1));
		return
		{
			1, 0, 0,
			0, 1, 0,
			0, 0, 1,
		};
	}
	else
	{
		assert(normalized(a));
		assert(c >= -1);
		assert(c <=  1);
		const double t = 1 - c;
		const double ta0a0 = t * a[0] * a[0];
		const double ta1a1 = t * a[1] * a[1];
		const double ta2a2 = t * a[2] * a[2];
		const double ta0a1 = t * a[0] * a[1];
		const double ta0a2 = t * a[0] * a[2];
		const double ta1a2 = t * a[1] * a[2];
		const double s = sqrt(1 - c * c); // s = sin(acos(c))
		const double sa0 = s * a[0];
		const double sa1 = s * a[1];
		const double sa2 = s * a[2];
		return
		{
			ta0a0 + c, ta0a1 - sa2, ta0a2 + sa1,
			ta0a1 + sa2, ta1a1 + c, ta1a2 - sa0,
			ta0a2 - sa1, ta1a2 + sa0, ta2a2 + c,
		};
	}
}

/// Transforms a vector by current 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
	};
}

#endif
