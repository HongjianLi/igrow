#pragma once
#ifndef IGROW_VEC3_HPP
#define IGROW_VEC3_HPP

#include <cmath>
#include <array>
using namespace std;

const double epsilon = 0.00001; ///< Tolerance for equality comparison of two floating point values.

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const double a, const double b)
{
	return fabs(a - b) < epsilon;
}

/// Returns the square of a generic value.
template<typename T>
inline T sqr(const T x)
{
	return x * x;
}

/// Returns true is the vector is (0, 0, 0).
inline bool zero(const array<double, 3>& v)
{
	return (eq(v[0], 0) && eq(v[1], 0) && eq(v[2], 0));
}

/// Returns the square norm.
inline double norm_sqr(const array<double, 3>& v)
{
	return sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
}

/// Returns the norm.
inline double norm(const array<double, 3>& v)
{
	return sqrt(norm_sqr(v));
}

/// Returns true if the norm equals 1.
inline bool normalized(const array<double, 3>& v)
{
	return eq(norm_sqr(v), 1);
}

/// Normalize the vector.
inline array<double, 3> normalize(const array<double, 3>& v)
{
	const double nrm = norm(v);
	if (eq(nrm, 0)) return v;
	const double f = 1 / nrm;
	return
	{
		f * v[0],
		f * v[1],
		f * v[2],
	};
}

/// Returns the dot product of the current vector and the given vector.
inline double operator*(const array<double, 3>& t, const array<double, 3>& v)
{
	return t[0] * v[0] + t[1] * v[1] + t[2] * v[2];
}

/// Returns the result of pairwise addition of the current vector and the given vector.
inline array<double, 3> operator+(const array<double, 3>& t, const array<double, 3>& v)
{
	return
	{
		t[0] + v[0],
		t[1] + v[1],
		t[2] + v[2],
	};
}

/// Returns the result of pairwise subtraction of the current vector and the given vector.
inline array<double, 3> operator-(const array<double, 3>& t, const array<double, 3>& v)
{
	return
	{
		t[0] - v[0],
		t[1] - v[1],
		t[2] - v[2],
	};
}

/// Pairwise multiply a constant to the current vector.
inline array<double, 3> operator*(const double s, const array<double, 3>& v)
{
	return
	{
		s * v[0],
		s * v[1],
		s * v[2],
	};
}

/// Returns the cross product of two vectors.
inline array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b)
{
	return
	{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};
}

/// Returns the square distance between two vectors.
inline double distance_sqr(const array<double, 3>& a, const array<double, 3>& b)
{
	return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

#endif
