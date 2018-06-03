#include "ForceFieldIntegrator.h"

ForceFieldIntegrator::ForceFieldIntegrator(double stepSize) : CustomIntegrator(stepSize) {

    addPerDofVariable("ff", 0);

    // TODO: Velocity Verlet not ideal for Lorentz forces. Replace by an integrator,
    //       where x and v are at equal times
    addUpdateContextState();
    addPerDofVariable("x1", 0);
    addComputePerDof("v", "v + 0.5*dt*(f+ff)/m");
    addComputePerDof("x", "x + dt*v");
    addComputePerDof("x1", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + 0.5*dt*(f+ff)/m + (x-x1)/dt");
    addConstrainVelocities();
}

void ForceFieldIntegrator::step(int steps) {
    if (owner != nullptr) {
        int numParticles = owner->getSystem().getNumParticles();
        if (f.size() != numParticles) {
            f.resize(numParticles);
        }
    }

    for (int i = 0; i < steps; i++) {
        /* get the current positions and velocities */
        const OpenMM::State state = owner->getState(OpenMM::State::Positions|OpenMM::State::Velocities);
        const std::vector<Vec3> &x = state.getPositions();
        const std::vector<Vec3> &v = state.getVelocities();

        /* reset forcefields */
        std::fill(f.begin(), f.end(), f_zero);

        /* evaluate all force fields */
        for(std::vector<ForceField *>::iterator it = forceFields.begin(); it != forceFields.end(); it++) {
            if (*it != nullptr) {
                (*it)->eval(state.getTime(),x,v,f);
            }
        }

        /* Check all forces applied to virtual particles and distribute it to the connected real particles.
           Only ThreeParticleAverageSite are supported.
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
                }
            }
        }
        setPerDofVariableByName("ff", f);
        CustomIntegrator::step(1);
    }
}


