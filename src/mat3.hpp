#pragma once
#ifndef IGROW_MAT3_HPP
#define IGROW_MAT3_HPP

#include "vec3.hpp"

/// Represents a row-major 3x3 matrix for vector transformation.
class mat3 : public std::array<fl, 9>
{
public:
	/// Constructs an empty 3x3 matrix.
	mat3() {}

	/// Constructs a matrix with specified values.
	mat3(const fl d00, const fl d01, const fl d02,
		 const fl d10, const fl d11, const fl d12,
		 const fl d20, const fl d21, const fl d22)
	{
		(*this)[0] = d00; (*this)[1] = d01; (*this)[2] = d02;
		(*this)[3] = d10; (*this)[4] = d11; (*this)[5] = d12;
		(*this)[6] = d20; (*this)[7] = d21; (*this)[8] = d22;
	}

	/// Constructs a rotation matrix from a normalized axis and the cosine value of an angle.
	mat3(const vec3& a, const fl c)
	{
		if (a.zero())
		{
			assert(eq(c, 1) || eq(c, -1));
			(*this)[0] = 1;
			(*this)[1] = 0;
			(*this)[2] = 0;
			(*this)[3] = 0;
			(*this)[4] = 1;
			(*this)[5] = 0;
			(*this)[6] = 0;
			(*this)[7] = 0;
			(*this)[8] = 1;
		}
		else
		{
			assert(a.normalized());
			assert(c >= -1);
			assert(c <=  1);
			const fl t = 1 - c;
			const fl ta0a0 = t * a[0] * a[0];
			const fl ta1a1 = t * a[1] * a[1];
			const fl ta2a2 = t * a[2] * a[2];
			const fl ta0a1 = t * a[0] * a[1];
			const fl ta0a2 = t * a[0] * a[2];
			const fl ta1a2 = t * a[1] * a[2];
			const fl s = sqrt(1 - c * c); // s = sin(acos(c))
			const fl sa0 = s * a[0];
			const fl sa1 = s * a[1];
			const fl sa2 = s * a[2];
			(*this)[0] = ta0a0 + c;
			(*this)[1] = ta0a1 - sa2;
			(*this)[2] = ta0a2 + sa1;
			(*this)[3] = ta0a1 + sa2;
			(*this)[4] = ta1a1 + c;
			(*this)[5] = ta1a2 - sa0;
			(*this)[6] = ta0a2 - sa1;
			(*this)[7] = ta1a2 + sa0;
			(*this)[8] = ta2a2 + c;
		}
	}

	/// Returns the value at index (i, j) where j is the lowest dimension.
	fl operator()(const size_t i, const size_t j) const
	{
		assert(i < 3);
		assert(j < 3);
		return (*this)[3 * i + j];
	}

	/// Transforms a vector by current 3x3 matrix.
	vec3 operator*(const vec3& v) const
	{
		return vec3
		(
			(*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2],
			(*this)[3] * v[0] + (*this)[4] * v[1] + (*this)[5] * v[2],
			(*this)[6] * v[0] + (*this)[7] * v[1] + (*this)[8] * v[2]
		);
	}
};

#endif
