#ifndef FORCE_FIELD_INTEGRATOR_H_
#define FORCE_FIELD_INTEGRATOR_H_

#include "OpenMM.h"
#include "ForceField.h"
#include <vector>

class ForceFieldIntegrator: public OpenMM::CustomIntegrator {
public:
    ForceFieldIntegrator(double stepSize) : CustomIntegrator(stepSize) {}
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

    void distributeVirtualForces(std::vector<Vec3> &f) {
        /* Check all forces applied to virtual particles and distribute it to the connected real particles.
           Only ThreeParticleAverageSite and TwoParticleAverageSite are supported.
        */
        for(int i = 0; i < f.size(); i++) {
            if (owner->getSystem().isVirtualSite(i)) {
                const OpenMM::VirtualSite &v = owner->getSystem().getVirtualSite(i);
                if (dynamic_cast<const OpenMM::ThreeParticleAverageSite*>(&v) != NULL) {
                    const OpenMM::ThreeParticleAverageSite &site = dynamic_cast<const OpenMM::ThreeParticleAverageSite &>(v);
                    for(int j = 0; j < site.getNumParticles(); j++) {
                        int particleIndex = site.getParticle(j);
                        f[particleIndex] += f[i] * site.getWeight(j);
                    }
                    f[i] = f_zero;
                } else if (dynamic_cast<const OpenMM::TwoParticleAverageSite*>(&v) != NULL) {
                    const OpenMM::TwoParticleAverageSite &site = dynamic_cast<const OpenMM::TwoParticleAverageSite &>(v);
                    for(int j = 0; j < site.getNumParticles(); j++) {
                        int particleIndex = site.getParticle(j);
                        f[particleIndex] += f[i] * site.getWeight(j);
                    }
                    f[i] = f_zero;
                }
            }
        }
    }
    virtual void step(int steps) {
            CustomIntegrator::step(steps);
    };
protected:
    Vec3 f_zero;
    std::vector<ForceField *> forceFields;
};

#endif // FORCE_FIELD_INTEGRATOR_H_
