#ifndef LORENZ_FORCE_INTEGRATOR_H_
#define LORENZ_FORCE_INTEGRATOR_H_

#include "OpenMM.h"
#include <vector>

using OpenMM::Vec3; // so we can just say "Vec3" below

static const double ElectricFieldStrength = 9.648533290731905e-08;
static const double MagneticFieldStrength = ElectricFieldStrength * 1000;

class ElectricField {

public:
    ElectricField() {
        // electric field amplitude
        E0[0] = 1000;
        E0[1] = 0;
        E0[2] = 0;
        // convert to MD units
        E0 *= ElectricFieldStrength;
        // frequency
        freq = 2 * M_PI * 0.001; // in THz
        // time offset
        t0 = 0;
    }

    virtual Vec3 operator()(double t, const Vec3 &x) {
        return E0 * cos(freq * (t - t0));
    }

private:
    Vec3 E0; // in V/m
    double freq; // in THz
    double t0; // in ps
};

class MagneticField {

public:
    MagneticField() {
        // magnetic field amplitude
        B0[0] = 0;
        B0[1] = 0;
        B0[2] = 1;
        // convert to MD units
        B0 *= MagneticFieldStrength;
        // frequency
        freq = 2 * M_PI * 0.001; // in THz
        // time offset
        t0 = 0;
    }

    virtual Vec3 operator()(double t, const Vec3 &x) {
        return B0 * cos(freq * (t - t0));
    }

private:
    Vec3 B0;  // in T
    double freq; // in THz
    double t0; // in ps
};

class LorentzForceIntegrator: public OpenMM::CustomIntegrator {
public:
    LorentzForceIntegrator(double stepSize, ElectricField E, MagneticField B, const std::vector<double> &charges);

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

    virtual void step(int steps);
private:
    std::vector<double> q;
    Vec3 f_zero;
    ElectricField E;
    MagneticField B;
    std::vector<Vec3> fe;
    std::vector<Vec3> fm;
};

#endif // LORENZ_FORCE_INTEGRATOR_H_
