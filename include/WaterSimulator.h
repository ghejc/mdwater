#ifndef WATER_SIMULATOR_H
#define WATER_SIMULATOR_H

#include "OpenMM.h"

class WaterSimulator {
public:

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

    // Parameters for the H-O-H angle.
    static const double HOH_nominalAngleInDeg;
    static const double HOH_stiffnessInKcalPerRad2; // that is e=k(a-a0)^2

    // negative charge center
    static const double M_mass;
    static const double M_charge;
    static const double M_sigma;
    static const double M_epsilon;

    // center of mass
    static const double X_mass;
    static const double X_charge;
    static const double X_sigma;
    static const double X_epsilon;

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

    void addWaterMolecule();

    OpenMM::Context *context;
    OpenMM::System *system;
    OpenMM::Integrator *integrator;

    int NonbondedForce_Index;
    int HarmonicBondForce_Index;
    int HarmonicAngleForce_Index;
    int AndersenThermostat_Index;
};

#endif // WATER_SIMULATOR_H
