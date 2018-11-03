#include "ForceFieldIntegrator.h"

/*
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
    // TODO: addConstrainVelocities removes non-zero velocity of center of gravity
    // addConstrainVelocities();
}
*/
ForceFieldIntegrator::ForceFieldIntegrator(double stepSize) : CustomIntegrator(stepSize) {

    addPerDofVariable("ff", 0);
    addPerDofVariable("x1", 0); /* x(t) */
    addPerDofVariable("x2", 0); /* x(t-dt/2) */
    addPerDofVariable("v1", 0); /* v(t) */
    addPerDofVariable("v2", 0); /* v(t-dt/2) */
    addPerDofVariable("xold", 0);
    addGlobalVariable("steps", 0);

    /* Euler integrator for the first 3 steps */
    beginIfBlock("steps <= 2");

    addUpdateContextState();
    addComputePerDof("v", "v + dt*(f+ff)/m");
    addComputePerDof("x", "x + dt*v");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + (x-xold)/dt");
    addComputePerDof("x2", "x1");
    addComputePerDof("x1", "x");
    addComputePerDof("v2", "v1");
    addComputePerDof("v1", "v");

    endBlock();

    /* Two-Step Leap-frog Velocity Verlet integrator for the remaining steps */
    beginIfBlock("steps > 2");

    /* calculate x(t+dt/2),v(t+dt/2) */
    addUpdateContextState();
    addComputePerDof("x", "2*x1 - x2 + 0.25*dt*dt*(f+ff)/m");
    addComputePerDof("v", "v2 + dt*(f+ff)/m");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + (x-xold)/dt");
    addComputePerDof("x2", "x1");
    addComputePerDof("x1", "x");
    addComputePerDof("v2", "v1");
    addComputePerDof("v1", "v");

    /* calculate x(t+dt),v(t+dt) */
    addUpdateContextState();
    addComputePerDof("x", "2*x1 - x2 + 0.25*dt*dt*(f+ff)/m");
    addComputePerDof("v", "v2 + dt*(f+ff)/m");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + (x-xold)/dt");
    addComputePerDof("x2", "x1");
    addComputePerDof("x1", "x");
    addComputePerDof("v2", "v1");
    addComputePerDof("v1", "v");

    endBlock();
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
        setGlobalVariableByName("steps", i);
        CustomIntegrator::step(1);
    }
}


