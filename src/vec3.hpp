#pragma once
#ifndef IGROW_VEC3_HPP
#define IGROW_VEC3_HPP

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

/// Represents a vector of 3 floating point elements.
class vec3 : public array<double, 3>
{
public:
	/// Constructs a vector with specified values.
	vec3(const double d0, const double d1, const double d2)
	{
		(*this)[0] = d0;
		(*this)[1] = d1;
		(*this)[2] = d2;
	}

	/// Returns true is the vector is (0, 0, 0).
	bool zero() const
	{
		return (eq((*this)[0], 0) && eq((*this)[1], 0) && eq((*this)[2], 0));
	}

	/// Returns the square norm.
	double norm_sqr() const
	{
		return sqr((*this)[0]) + sqr((*this)[1]) + sqr((*this)[2]);
	}

	/// Returns the norm.
	double norm() const
	{
		return sqrt(norm_sqr());
	}

	/// Returns true if the norm equals 1.
	bool normalized() const
	{
		return eq(norm_sqr(), 1);
	}

	/// Normalize the vector.
	vec3 normalize() const
	{
		if (zero()) return *this;
		const double factor = 1 / norm();
		return vec3(factor * (*this)[0], factor * (*this)[1], factor * (*this)[2]);
	}

	/// Returns the dot product of the current vector and the given vector.
	double operator*(const vec3& v) const
	{
		return (*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2];
	}

	/// Returns the result of pairwise multiplication of the current vector and the given vector.
	vec3 operator*(const array<size_t, 3>& v) const
	{
		return vec3((*this)[0] * v[0], (*this)[1] * v[1], (*this)[2] * v[2]);
	}

	/// Returns the result of pairwise addition of the current vector and the given vector.
	vec3 operator+(const vec3& v) const
	{
		return vec3((*this)[0] + v[0], (*this)[1] + v[1], (*this)[2] + v[2]);
	}

	/// Returns the result of pairwise subtraction of the current vector and the given vector.
	vec3 operator-(const vec3& v) const
	{
		return vec3((*this)[0] - v[0], (*this)[1] - v[1], (*this)[2] - v[2]);
	}
};

/// Pairwise multiply a constant to the current vector.
inline vec3 operator*(const double s, const vec3& v)
{
	return vec3(s * v[0], s * v[1], s * v[2]);
}

/// Returns the normalized vector of a vector.
inline vec3 normalize(const vec3& v)
{
	return (1 / v.norm()) * v;
}

/// Returns the cross product of two vectors.
inline vec3 cross_product(const vec3& a, const vec3& b)
{
	return vec3(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Returns the square distance between two vectors.
inline double distance_sqr(const vec3& a, const vec3& b)
{
	return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

#endif
