#ifndef FORCE_FIELD_H_
#define FORCE_FIELD_H_

#include "OpenMM.h"

using OpenMM::Vec3;

class ForceField {
public:
    virtual void eval(int particleIndex, double t, const Vec3 &x, const Vec3 &v, Vec3 &f) = 0;
};

#endif // FORCE_FIELD_H_
