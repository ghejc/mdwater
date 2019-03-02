#ifndef LORENZ_FORCE_FIELD_INTEGRATOR_H_
#define LORENZ_FORCE_FIELD_INTEGRATOR_H_

#include "OpenMM.h"
#include "ForceFieldIntegrator.h"
#include <vector>

static const double ElectricFieldStrength = 9.648533290731905e-08;
static const double MagneticFieldStrength = ElectricFieldStrength * 1000;

class ElectricForceField : public ForceField {

public:

    ElectricForceField(std::vector<double> &charges):q(charges) {
        // electric field amplitude
        E0[0] = 1000;
        E0[1] = 0;
        E0[2] = 0;
        // convert to MD units
        E0 *= ElectricFieldStrength;
        // frequency
        freq = 2 * M_PI * 0.001; // in THz
        // phase difference between E and B
        phas = M_PI / 2.0;
        // time offset
        t0 = 0;

    };

    void eval(int particleIndex, double t, const Vec3 &x, const Vec3 &v, Vec3 &f) {
        f = E0 * q[particleIndex] * cos(freq * (t - t0));
    };

private:

    std::vector<double> q;
    Vec3 E0; // in V/m
    double freq; // in THz
    double phas; // in radians
    double t0; // in ps
};

class MagneticForceField : public ForceField {

public:

    MagneticForceField(std::vector<double> &charges):q(charges) {
        // magnetic field amplitude
        B0[0] = 0;
        B0[1] = 0;
        B0[2] = 1;
        // convert to MD units
        B0 *= MagneticFieldStrength;
        // frequency
        freq = 2 * M_PI * 0.001; // in THz
        // phase difference between E and B
        phas = M_PI / 2.0;
        // time offset
        t0 = 0;
    };

    void eval(int particleIndex, double t, const Vec3 &x, const Vec3 &v, Vec3 &f) {
        f = v.cross(B0) * q[particleIndex] * cos(freq * (t - t0) + phas);
    };

private:

    std::vector<double> q;
    Vec3 B0;  // in T
    double freq; // in THz
    double phas; // in radians
    double t0; // in ps
};

class LorenzForceFieldIntegrator: public ForceFieldIntegrator {
public:
    LorenzForceFieldIntegrator(double stepSize, std::vector<double> &charges);

    virtual void step(int steps);
private:
    std::vector<Vec3> fe;
    std::vector<Vec3> fm;
};

#endif // LORENZ_FORCE_FIELD_INTEGRATOR_H_
