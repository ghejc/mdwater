#include "LorenzForceFieldIntegrator.h"

LorenzForceFieldIntegrator::LorenzForceFieldIntegrator(double stepSize, std::vector<double> &charges) : ForceFieldIntegrator(stepSize) {

    addForceField(new ElectricForceField(charges));
    addForceField(new MagneticForceField(charges));

    addPerDofVariable("fe", 0);
    addPerDofVariable("fm", 0);
    addPerDofVariable("xold", 0);

    // TODO: Velocity Verlet not ideal for Lorentz forces.
    addUpdateContextState();
    addComputePerDof("v", "v + 0.5*dt*(f+fe+fm)/m");
    addComputePerDof("x", "x + dt*v");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + 0.5*dt*(f+fe+fm)/m + (x-xold)/dt");
    addConstrainVelocities();
}

void LorenzForceFieldIntegrator::step(int steps) {
    if (owner != nullptr) {
        int numParticles = owner->getSystem().getNumParticles();
        if (fe.size() != numParticles) {
            fe.resize(numParticles);
        }
        if (fm.size() != numParticles) {
            fm.resize(numParticles);
        }
    }

    for (int n = 0; n < steps; n++) {
        /* get the current positions and velocities */
        const OpenMM::State state = owner->getState(OpenMM::State::Positions|OpenMM::State::Velocities);
        const std::vector<Vec3> &x = state.getPositions();
        const std::vector<Vec3> &v = state.getVelocities();

        /* evaluate all force fields */
        std::vector<ForceField *>::iterator it = forceFields.begin();
        if (*it != nullptr) {
            for (int i; i < fe.size(); i++) {
                (*it)->eval(i, state.getTime(),x[i],v[i],fe[i]);
            }
        }
        it++;
        if (*it != nullptr) {
            for (int i; i < fm.size(); i++) {
                (*it)->eval(i, state.getTime(),x[i],v[i],fm[i]);
            }
        }

        distributeVirtualForces(fe);
        distributeVirtualForces(fm);

        setPerDofVariableByName("fe", fe);
        setPerDofVariableByName("fm", fm);
        CustomIntegrator::step(1);
    }
}


