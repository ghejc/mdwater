#include "LorentzForceIntegrator.h"

LorentzForceIntegrator::LorentzForceIntegrator(double stepSize, ElectricField *E, MagneticField *B, const std::vector<double> &charges) :
    CustomIntegrator(stepSize), E(E), B(B), q(charges) {

    addPerDofVariable("fe", 0);
    addPerDofVariable("fm", 0);
    addPerDofVariable("xold", 0);
    addPerDofVariable("t", 0);
    addPerDofVariable("s", 0);

    // Boris integrator.
    addUpdateContextState();
    addComputePerDof("v", "v + 0.5*dt*(f+fe)/m");
    addComputePerDof("t", "0.5*dt*fm/m");
    addComputePerDof("s", "2*t/(1+dot(t,t))");
    addComputePerDof("v", "v + cross(v,s) + cross(cross(v,t),s)");
    addComputePerDof("v", "v + 0.5*dt*fe/m");
    addComputePerDof("x", "x + dt*v");
    addComputePerDof("xold", "x");
    addConstrainPositions();
    addComputePerDof("v", "v + 0.5*dt*f/m + (x-xold)/dt");
    addConstrainVelocities();
}

void LorentzForceIntegrator::step(int steps) {
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

        /* evaluate all fields */
        for (int i = 0; i < fe.size(); i++) {
            fe[i] = (*E)(state.getTime(),x[i]) * q[i];
        }
        for (int i = 0; i < fm.size(); i++) {
            // fm[i] = v[i].cross((*B)(state.getTime(),x[i])) * q[i];
            fm[i] = (*B)(state.getTime(),x[i]) * q[i];
        }

        distributeVirtualForces(fe);
        distributeVirtualForces(fm);

        setPerDofVariableByName("fe", fe);
        setPerDofVariableByName("fm", fm);
        CustomIntegrator::step(1);
    }
}


