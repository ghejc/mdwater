#ifndef WATER_SIMULATOR_H_
#define WATER_SIMULATOR_H_

#include "OpenMM.h"
#include "ForceField.h"
#include "ForceFieldIntegrator.h"
#include <cmath>

#define PARTICLE_M_IS_VIRTUAL

class WaterSimulator {
public:

    static const double AtomicMassUnitsPerKg;
    static const double MeterPerNm;

    static const bool  UseConstraints = true;   // Should we constrain O-H bonds?
    static const double CutoffDistanceInAng;    // Angstroms
    static const double FrictionInPerPs; // collisions per picosecond

    static const double Coulomb14Scale;
    static const double LennardJones14Scale;

    // Values are from tip4pew.xml.
    // oxygen
    static const double O_mass;
    static const double O_charge;
    static const double O_sigma;
    static const double O_epsilon;

    // hydrogen
    static const double H_mass;
    static const double H_charge;
    static const double H_sigma;
    static const double H_epsilon;

    // Parameters for the O-H bonds.
    static const double OH_nominalLengthInAng;
    static const double OH_stiffnessInKcalPerAng2; // that is, e=k(x-x0)^2

    // Parameters for the O-M bonds.
    static const double OM_nominalLengthInAng;

    // Parameters for the H-O-H angle.
    static const double HOH_nominalAngleInDeg;
    static const double HOH_stiffnessInKcalPerRad2; // that is e=k(a-a0)^2

    // negative charge center
    static const double M_mass;
    static const double M_charge;
    static const double M_sigma;
    static const double M_epsilon;
#if defined(PARTICLE_M_IS_VIRTUAL)
    static const double O_weight;
    static const double H_weight;
#endif

    WaterSimulator(unsigned int numOfMolecules,
                      double temperature = 300,
                      double stepSizeInFs = 1);
    ~WaterSimulator() {
        delete context;
        delete system;
        delete integrator;
    };

    void initSystemState(double startTemperatureInK, double densityInKgPerM3);

    void getSystemState(double& timeInPs, std::vector<double>& atomPositionsInAng);

    void run(double SimulationTimeInPs);

private:

    void setRandomPositions(const double boxLengthInNm);

    void addWaterMolecule();

    void getCenterOfMassCoordinates(std::vector<OpenMM::Vec3> &coords);

    OpenMM::Context *context;
    OpenMM::System *system;
    ForceFieldIntegrator *integrator;

    int NonbondedForce_Index;
    int HarmonicBondForce_Index;
    int HarmonicAngleForce_Index;
    int AndersenThermostat_Index;
};

class ElectroMagneticForceField : public ForceField {

public:

    ElectroMagneticForceField(std::vector<double> &charges):q(charges) {
        const double ElectricFieldStrength = 9.648533290731905e-08;
        const double MagneticFieldStrength = ElectricFieldStrength * 1000;
        // electric field amplitude
        E0[0] = 1000;
        E0[1] = 0;
        E0[2] = 0;
        // magnetic field amplitude
        B0[0] = 0;
        B0[1] = 0;
        B0[2] = 1;
        // convert to MD units
        E0 *= ElectricFieldStrength;
        B0 *= MagneticFieldStrength;
        // frequency
        freq = 2 * M_PI * 0.001; // in THz
        // phase difference between E and B
        phas = M_PI / 2.0;
        // time offset
        t0 = 0;

    };

    virtual void eval(double t, const std::vector<Vec3> &x, const std::vector<Vec3> &v, std::vector<Vec3> &f) {
        for (size_t i = 0; i < v.size(); i++) {
             // apply Lorentz force on particle i
            applyForce(f, i, (E(t,x[i]) + v[i].cross(B(t,x[i]))) * q[i]);
        }
    };

    Vec3 E(double t, const Vec3 &x) {
        return E0 * cos(freq * (t - t0));
    };

    Vec3 B(double t, const Vec3 &x) {
        return B0 * cos(freq * (t - t0) + phas);
    };

private:

    std::vector<double> q;
    Vec3 E0; // in V/m
    Vec3 B0;  // in T
    double freq; // in THz
    double phas; // in radians
    double t0; // in ps
};

#endif // WATER_SIMULATOR_H_
