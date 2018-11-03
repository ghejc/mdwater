#ifndef MATRIX3_H_
#define MATRIX3_H_

#include <cassert>
#include "OpenMM.h"

using OpenMM::Vec3; // so we can just say "Vec3" below

/**
 * This class represents a 3x3 matrix.  It is used for rotation matrices
 */

class Matrix3 {
public:
    /**
     * Create a Matrix3 whose elements are all 0.
     */
    Matrix3() {
    }
    /**
     * Create a Matrix3 with specified x, y, and z components as diagonal elements.
     */
    Matrix3(double x, double y, double z) {
        data[0][0] = x;
        data[1][1] = y;
        data[2][2] = z;
    }

    /**
     * Create a Matrix3 from a single vector.
     */
    Matrix3(const Vec3 &x) {
        data[0][0] = x[0];
        data[1][1] = x[1];
        data[2][2] = x[2];
    }

    /**
     * Create a Matrix3 from another matrix.
     */
    Matrix3(const Matrix3 &x) {
        data[0] = x[0];
        data[1] = x[1];
        data[2] = x[2];
    }

    /**
     * Create a Matrix3 from an outer product of two vectors.
     */
    Matrix3(const Vec3 &x, const Vec3 &y) {
        data[0] = y * x[0];
        data[1] = y * x[1];
        data[2] = y * x[2];
    }

    /**
     * Create a Matrix3 with specified x, y and z vectors.
     */
    Matrix3(const Vec3 &x, const Vec3 &y, const Vec3 &z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    Vec3 operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    Vec3& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    bool operator==(const Matrix3& rhs) const {
        return (data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
    }

    bool operator!=(const Matrix3& rhs) const {
        return (data[0] != rhs[0] || data[1] != rhs[1] || data[2] != rhs[2]);
    }

    // Arithmetic operators

    // unary plus
    Matrix3 operator+() const {
        return Matrix3(*this);
    }

    // plus
    Matrix3 operator+(const Matrix3& rhs) const {
        const Matrix3& lhs = *this;
        return Matrix3(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }

    Matrix3& operator+=(const Matrix3& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    // unary minus
    Matrix3 operator-() const {
        const Matrix3& lhs = *this;
        return Matrix3(-lhs[0], -lhs[1], -lhs[2]);
    }

    // minus
    Matrix3 operator-(const Matrix3& rhs) const {
        const Matrix3& lhs = *this;
        return Matrix3(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
    }

    Matrix3& operator-=(const Matrix3& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

    // transpose the matrix
    Matrix3 transpose() const {
        const Matrix3& lhs = *this;
        return Matrix3(Vec3(lhs[0][0], lhs[1][0], lhs[2][0]),
                       Vec3(lhs[0][1], lhs[1][1], lhs[0][1]),
                       Vec3(lhs[0][2], lhs[1][2], lhs[2][2]));
    }

    // matrix multiplication
    Matrix3 operator*(const Matrix3 &rhs) const {
        const Matrix3& lhs = *this;
        const Matrix3 m = rhs.transpose();
        return Matrix3(lhs.dot(m[0]), lhs.dot(m[1]), lhs.dot(m[2]));
    }

    // scalar product
    Matrix3 operator*(double rhs) const {
        const Matrix3& lhs = *this;
        return Matrix3(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
    }

    Matrix3& operator*=(double rhs) {
        data[0] *= rhs;
        data[1] *= rhs;
        data[2] *= rhs;
        return *this;
    }

    // scalar division
    Matrix3 operator/(double rhs) const {
        const Matrix3& lhs = *this;
        double scale = 1.0/rhs;
        return Matrix3(lhs[0]*scale, lhs[1]*scale, lhs[2]*scale);
    }

    Matrix3& operator/=(double rhs) {
        double scale = 1.0/rhs;
        data[0] *= scale;
        data[1] *= scale;
        data[2] *= scale;
        return *this;
    }

    // dot product
    Vec3 dot(const Vec3& rhs) const {
        const Matrix3& lhs = *this;
        return Vec3(lhs[0].dot(rhs), lhs[1].dot(rhs), lhs[2].dot(rhs));
    }

    // cross product
    Matrix3 cross(const Vec3& rhs) const {
        const Matrix3& lhs = *this;
        return Matrix3(lhs[0].cross(rhs), lhs[1].cross(rhs), lhs[2].cross(rhs));
    }

protected:
    Vec3 data[3]; // row vectors
};


#endif /*MATRIX3_H_*/
