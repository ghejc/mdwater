#include "LorentzForceIntegrator.h"

LorentzForceIntegrator::LorentzForceIntegrator(double stepSize, ElectricField *E, MagneticField *B, const std::vector<double> &charges) :
    CustomIntegrator(stepSize), E(E), B(B), q(charges) {

    addPerDofVariable("fe", 0);
    addPerDofVariable("fm", 0);
    addPerDofVariable("xold", 0);
    addPerDofVariable("s1", 0);
    addPerDofVariable("s2", 0);

    // Boris integrator.
    addUpdateContextState();
    addComputePerDof("v", "v + 0.5*dt*(f+fe)/m");
    addComputePerDof("s1", "0.5*dt*fm/m");
    addComputePerDof("s2", "2*s1/(1+dot(s1,s1))");
    addComputePerDof("v", "v + cross(v,s2) + cross(cross(v,s1),s2)");
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


