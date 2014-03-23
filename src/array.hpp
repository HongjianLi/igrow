#pragma once
#ifndef IGROW_ARRAY_HPP
#define IGROW_ARRAY_HPP

#include <array>
using namespace std;

//! Returns true if the absolute difference between a scalar and zero is within the constant tolerance.
bool zero(const double a);

//! Returns true is the vector is approximately (0, 0, 0).
bool zero(const array<double, 3>& a);

//! Returns the square norm of a vector.
double norm_sqr(const array<double, 3>& a);

//! Returns the norm of a vector.
double norm(const array<double, 3>& a);

//! Returns true if the norm of a vector is approximately 1.
bool normalized(const array<double, 3>& a);

//! Normalize the vector.
array<double, 3> normalize(const array<double, 3>& a);

//! Returns the dot product of the two vectors.
double operator*(const array<double, 3>& a, const array<double, 3>& b);

//! Returns the result of pairwise addition of two vectors.
array<double, 3> operator+(const array<double, 3>& a, const array<double, 3>& b);

//! Returns the result of pairwise subtraction of two vectors.
array<double, 3> operator-(const array<double, 3>& a, const array<double, 3>& b);

//! Pairwise multiply a constant to a vector.
array<double, 3> operator*(const double s, const array<double, 3>& a);

//! Returns the cross product of two vectors.
array<double, 3> cross_product(const array<double, 3>& a, const array<double, 3>& b);

//! Returns the square distance between two vectors.
double distance_sqr(const array<double, 3>& a, const array<double, 3>& b);

//! Constructs a rotation matrix from a normalized axis and the cosine value of an angle.
array<double, 9> vec3_to_mat3(const array<double, 3>& a, const double c);

//! Transforms a vector by a 3x3 matrix.
array<double, 3> operator*(const array<double, 9>& m, const array<double, 3>& v);

#endif
