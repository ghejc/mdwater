#ifndef ROTATIONMATRIX_H_
#define ROTATIONMATRIX_H_

#include <cassert>
#include <cmath>
#include "Matrix3.h"

using OpenMM::Vec3; // so we can just say "Vec3" below

/**
 * This class represents a rotation matrix.
 */

class RotationMatrix: public Matrix3 {
public:
    RotationMatrix() : Matrix3(1,1,1) {
    }

    RotationMatrix(double angle, const Vec3 & rotationAxis) {
        double sina = sin(angle);
        double cosa = cos(angle);
        double len = sqrt(rotationAxis.dot(rotationAxis));
        Vec3 rotationAxisUnitVector;
        if (len > 0.0) {
            rotationAxisUnitVector = rotationAxis / len;
        } else {
            rotationAxisUnitVector = Vec3(0,0,1);
        }
        // rotation matrix around unit vector
        Matrix3 R(cosa, cosa, cosa);
        Matrix3 O(rotationAxisUnitVector, rotationAxisUnitVector);
        R += O * (1.0 - cosa);
        rotationAxisUnitVector *= sina;
        Matrix3 P(Vec3(0.0, - rotationAxisUnitVector[2], rotationAxisUnitVector[1]),
                  Vec3(rotationAxisUnitVector[2], 0.0, - rotationAxisUnitVector[0]),
                  Vec3(- rotationAxisUnitVector[1], rotationAxisUnitVector[0], 0.0)
                  );
        R += P;
        data[0] = R[0];
        data[1] = R[1];
        data[2] = R[2];
    }
};


#endif /*ROTATIONMATRIX_H_*/
