#pragma once
#ifndef IGROW_ARRAY_HPP
#define IGROW_ARRAY_HPP

#include <array>
using namespace std;

/// Returns true if the absolute difference between two floating point values is within the constant tolerance.
bool eq(const double a, const double b);

/// Returns true is the vector is (0, 0, 0).
bool zero(const array<double, 3>& v);

/// Returns the square norm.
double norm_sqr(const array<double, 3>& v);

/// Returns the norm.
double norm(const array<double, 3>& v);

/// Returns true if the norm equals 1.
bool normalized(const array<double, 3>& v);

/// Normalize the vector.
array<double, 3> normalize(const array<double, 3>& v);

/// Returns the dot product of the current vector and the given vector.
double operator*(const array<double, 3>& t, const array<double, 3>& v);

/// Returns the result of pairwise addition of the current vector and the given vector.
array<double, 3> operator+(const array<double, 3>& t, const array<double, 3>& v);

/// Returns the result of pairwise subtraction of the current vector and the given vector.
array<double, 3> operator-(const array<double, 3>& t, const array<double, 3>& v);

/// Pairwise multiply a constant to the current vector.
array<double, 3> operator*(const double s, const array<double, 3>& v);

/// Returns the cross product of two vectors.
array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b);

/// Returns the square distance between two vectors.
double distance_sqr(const array<double, 3>& a, const array<double, 3>& b);

/// Constructs a rotation matrix from a normalized axis and the cosine value of an angle.
array<double, 9> vec3_to_mat3(const array<double, 3>& a, const double c);

/// Transforms a vector by current 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v);

#endif
