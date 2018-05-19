#ifndef FORCE_FIELD_INTEGRATOR_H_
#define FORCE_FIELD_INTEGRATOR_H_

#include "OpenMM.h"
#include "ForceField.h"
#include <vector>

class ForceFieldIntegrator: public OpenMM::CustomIntegrator {
public:
    ForceFieldIntegrator(double stepSize);
    ~ForceFieldIntegrator() {
        for(std::vector<ForceField *>::iterator it = forceFields.begin(); it != forceFields.end(); it++) {
            if (*it != nullptr) {
                delete *it;
                *it = nullptr;
            }
        }
    };
    int addForceField(ForceField *ff) {
        forceFields.push_back(ff);
        return forceFields.size() - 1;
    };
    void step(int steps);
private:
    std::vector<Vec3> f;
    std::vector<ForceField *> forceFields;
};

#endif // FORCE_FIELD_INTEGRATOR_H_
