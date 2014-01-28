#pragma once
#ifndef IGROW_COMMON_HPP
#define IGROW_COMMON_HPP

/// igrow uses double precision floating point computation by default.
/// This could possible be demoted to single precision for better performance.
typedef double fl;

const fl epsilon = static_cast<fl>(0.00001); ///< Tolerance for equality comparison of two floating point values.

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
inline bool eq(const fl a, const fl b)
{
	return fabs(a - b) < epsilon;
}

/// Returns the square of a generic value.
template<typename T>
inline T sqr(const T x)
{
	return x * x;
}

#endif
