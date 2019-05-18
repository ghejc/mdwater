#ifndef WATER_SIMULATOR_H_
#define WATER_SIMULATOR_H_

#include "OpenMM.h"
#include "LorentzForceIntegrator.h"
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
                      double densityInKgPerM3 = 1000,
                      double stepSizeInFs = 1);
    ~WaterSimulator() {
        delete context;
        delete system;
        delete integrator;
    };

    void initSystemState(double temperatureInK);

    void getSystemState(double& timeInPs, std::vector<double>& atomPositionsInAng);

    void run(double SimulationTimeInPs);

private:

    double setRandomPositions(const double boxLengthInNm);

    void addWaterMolecule();

    void getCenterOfMassCoordinates(std::vector<OpenMM::Vec3> &coords);

    OpenMM::Context *context;
    OpenMM::System *system;
    LorentzForceIntegrator *integrator;

    int NonbondedForce_Index;
    int HarmonicBondForce_Index;
    int HarmonicAngleForce_Index;
    int AndersenThermostat_Index;
};

#endif // WATER_SIMULATOR_H_
