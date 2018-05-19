#ifndef FORCE_FIELD_H_
#define FORCE_FIELD_H_

#include "OpenMM.h"

using OpenMM::Vec3;

class ForceField {
public:
    virtual void eval(double t, const std::vector<Vec3> &x, const std::vector<Vec3> &v, std::vector<Vec3> &f) = 0;
};

#endif // FORCE_FIELD_H_
