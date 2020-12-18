#ifndef LORENZ_FORCE_INTEGRATOR_H_
#define LORENZ_FORCE_INTEGRATOR_H_

#include "OpenMM.h"
#include <vector>
#include <memory>

using OpenMM::Vec3; // so we can just say "Vec3" below

static const double ElectricFieldStrength = 9.648533290731905e-08;
static const double MagneticFieldStrength = ElectricFieldStrength * 1000;

//base class for oscillating fields
class OscillatingField {

public:
    OscillatingField(const Vec3 &F = Vec3(),
                     double f = 0 /* THz */,
                     double t = 0 /* ps */) : F0(F), freq(f), t0(t) {}

    virtual ~OscillatingField() {}

    Vec3 operator()(double t, const Vec3 &x) {
        return F0 * cos(2 * M_PI * freq * (t - t0));
    }

protected:
    Vec3 F0; // in MD units
    double freq; // in THz
    double t0; // in ps
};

// base class for an electric field
class ElectricField : public OscillatingField {

public:
    ElectricField(const Vec3 &E = Vec3() /* V/m */,
                  double f = 0 /* THz */,
                  double t = 0 /* ps */) : OscillatingField(E * ElectricFieldStrength, f, t) {}
};

// base class for a magnetic field
class MagneticField: public OscillatingField {

public:
    MagneticField(const Vec3 &B = Vec3() /* T */,
                  double f = 0 /* THz */,
                  double t = 0 /* ps */) : OscillatingField(B * MagneticFieldStrength, f, t) {}
};

class LorentzForceIntegrator: public OpenMM::CustomIntegrator {
public:
    LorentzForceIntegrator(double stepSize, ElectricField *E, MagneticField *B, const std::vector<double> &charges, bool applyThermostat = false);

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
                    f[i] = Vec3();
                } else if (dynamic_cast<const OpenMM::TwoParticleAverageSite*>(&v) != NULL) {
                    const OpenMM::TwoParticleAverageSite &site = dynamic_cast<const OpenMM::TwoParticleAverageSite &>(v);
                    for(int j = 0; j < site.getNumParticles(); j++) {
                        int particleIndex = site.getParticle(j);
                        f[particleIndex] += f[i] * site.getWeight(j);
                    }
                    f[i] = Vec3();
                }
            }
        }
    }

    virtual void step(int steps);
    virtual bool kineticEnergyRequiresForce() const;
private:
    std::vector<double> q;
    std::unique_ptr<ElectricField> E;
    std::unique_ptr<MagneticField> B;
    std::vector<Vec3> fe;
    std::vector<Vec3> fm;
    bool applyThermostat;
};

#endif // LORENZ_FORCE_INTEGRATOR_H_
