#include <cassert>
#include <cmath>
#include "array.hpp"

const double epsilon = 1e-5; //!< Tolerance for equality comparison of two floating point values.

//! Returns the square of a generic value.
template<typename T>
T sqr(const T x)
{
	return x * x;
}

//! Returns true if the absolute difference between two floating point values is within the constant tolerance.
bool zero(const double a)
{
	return fabs(a) < epsilon;
}

//! Returns true is the vector is (0, 0, 0).
bool zero(const array<double, 3>& a)
{
	return zero(a[0]) && zero(a[1]) && zero(a[2]);
}

//! Returns the square norm.
double norm_sqr(const array<double, 3>& a)
{
	return sqr(a[0]) + sqr(a[1]) + sqr(a[2]);
}

//! Returns the norm.
double norm(const array<double, 3>& a)
{
	return sqrt(norm_sqr(a));
}

//! Returns true if the norm equals 1.
bool normalized(const array<double, 3>& a)
{
	return zero(norm_sqr(a) - 1);
}

//! Normalize the vector.
array<double, 3> normalize(const array<double, 3>& a)
{
	const double nrm = norm(a);
	if (zero(nrm)) return a;
	return (1 / nrm) * a;
}

//! Returns the dot product of the current vector and the given vector.
double operator*(const array<double, 3>& a, const array<double, 3>& b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

//! Returns the result of pairwise addition of the current vector and the given vector.
array<double, 3> operator+(const array<double, 3>& a, const array<double, 3>& b)
{
	return
	{
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2],
	};
}

//! Returns the result of pairwise subtraction of the current vector and the given vector.
array<double, 3> operator-(const array<double, 3>& a, const array<double, 3>& b)
{
	return
	{
		a[0] - b[0],
		a[1] - b[1],
		a[2] - b[2],
	};
}

//! Pairwise multiply a constant to the current vector.
array<double, 3> operator*(const double s, const array<double, 3>& a)
{
	return
	{
		s * a[0],
		s * a[1],
		s * a[2],
	};
}

//! Returns the cross product of two vectors.
array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b)
{
	return
	{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};
}

//! Returns the square distance between two vectors.
double distance_sqr(const array<double, 3>& a, const array<double, 3>& b)
{
	return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

//! Constructs a rotation matrix from a normalized axis and the cosine value of an angle.
array<double, 9> vec3_to_mat3(const array<double, 3>& a, const double c)
{
	if (zero(a))
	{
		assert(zero(c - 1) || zero(c + 1));
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

//! Transforms a vector by current 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
	};
}
