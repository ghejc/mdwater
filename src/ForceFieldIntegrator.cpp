#include "ForceFieldIntegrator.h"

LorenzForceFieldIntegrator::LorenzForceFieldIntegrator(double stepSize) : ForceFieldIntegrator(stepSize) {

    addForceField(new ElectricForceField());
    addForceField(new MagneticForceField());

    addPerDofVariable("ff", 0);
    addPerDofVariable("xold", 0);

    // TODO: Velocity Verlet not ideal for Lorentz forces.
    addUpdateContextState();
    addComputePerDof("v", "v + 0.5*dt*(f+ff)/m");
    addComputePerDof("x", "x + dt*v");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    // TODO: ff should be recalculated here with new positions and velocities
    addComputePerDof("v", "v + 0.5*dt*(f+ff)/m + (x-xold)/dt");
    addConstrainVelocities();
}

void LorenzForceFieldIntegrator::step(int steps) {
    std::vector<Vec3> f;
    if (owner != nullptr) {
        int numParticles = owner->getSystem().getNumParticles();
        if (f.size() != numParticles) {
            f.resize(numParticles);
        }
    }

    for (int n = 0; n < steps; n++) {
        /* get the current positions and velocities */
        const OpenMM::State state = owner->getState(OpenMM::State::Positions|OpenMM::State::Velocities);
        const std::vector<Vec3> &x = state.getPositions();
        const std::vector<Vec3> &v = state.getVelocities();

        /* reset forcefields */
        std::fill(f.begin(), f.end(), f_zero);

        /* evaluate all force fields */
        for(std::vector<ForceField *>::iterator it = forceFields.begin(); it != forceFields.end(); it++) {
            if (*it != nullptr) {
                (*it)->eval(particleIndex, state.getTime(),x,v,f);
            }
            particleIndex++;
        }

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
        setPerDofVariableByName("fe", f);
        CustomIntegrator::step(1);
    }
}


