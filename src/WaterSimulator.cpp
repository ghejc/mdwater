#include "WaterSimulator.h"
#include <string>
#include <vector>
#include <cstdlib>

using OpenMM::Vec3; // so we can just say "Vec3" below

static void
myWritePDBFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng);

//                     SIMULATION PARAMETERS
const double Temperature         = 300;    // Kelvins
const double StepSizeInFs        = 2;      // integration step size (fs)
const double ReportIntervalInFs  = 100;    // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 10;

const double WaterSimulator::FrictionInPerPs     = 91.;    // collisions per picosecond
const double WaterSimulator::CutoffDistanceInAng = 10.;    // Angstroms

const double WaterSimulator::Coulomb14Scale = 0.833333;
const double WaterSimulator::LennardJones14Scale = 0.5;

// Values are from tip4pew.xml.
// oxygen
const double WaterSimulator::O_mass             = 15.99943;
const double WaterSimulator::O_charge           = 0.0;
const double WaterSimulator::O_sigma            = 1;
const double WaterSimulator::O_epsilon          = 0;

// hydrogen
const double WaterSimulator::H_mass             = 1.007947;
const double WaterSimulator::H_charge           = 0.52422;
const double WaterSimulator::H_sigma            = 1;
const double WaterSimulator::H_epsilon          = 0;

// Parameters for the O-H bonds.
const double WaterSimulator::OH_nominalLengthInAng      = 0.9572;
const double WaterSimulator::OH_stiffnessInKcalPerAng2  = 553.0; // that is, e=k(x-x0)^2

// Parameters for the H-O-H angle.
const double WaterSimulator::HOH_nominalAngleInDeg      = 104.52;
const double WaterSimulator::HOH_stiffnessInKcalPerRad2 = 100.; // that is e=k(a-a0)^2

// negative charge center
const double WaterSimulator::M_mass             = 0;
const double WaterSimulator::M_charge           = -2 * H_charge;
const double WaterSimulator::M_sigma            = 1;
const double WaterSimulator::M_epsilon          = 0;

// center of mass
const double WaterSimulator::X_mass             = 0;
const double WaterSimulator::X_charge           = 0;
const double WaterSimulator::X_sigma            = 0.316435;
const double WaterSimulator::X_epsilon          = 0.680946;

// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
// (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     in a manner which is opaque to the caller.
// (3) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (4) Create an Integrator and a Context associating the Integrator with
//     the System.
// (5) Select the OpenMM platform to be used.
// (6) Return the MyOpenMMData struct and the name of the Platform in use.
//
// Note that this function must understand the calling MD code's molecule and
// force field data structures so will need to be customized for each MD code.

WaterSimulator::WaterSimulator      ( unsigned int        numOfMolecules,
                    double              temperature,
                    double              stepSizeInFs)
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
       (OpenMM::Platform::getDefaultPluginsDirectory());

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    system      = new OpenMM::System();

    OpenMM::NonbondedForce *        nonbond     = new OpenMM::NonbondedForce();
    NonbondedForce_Index = system->addForce(nonbond);

    OpenMM::HarmonicBondForce *     bondStretch = new OpenMM::HarmonicBondForce();
    HarmonicBondForce_Index = system->addForce(bondStretch);

    OpenMM::HarmonicAngleForce *    bondBend    = new OpenMM::HarmonicAngleForce();
    HarmonicAngleForce_Index = system->addForce(bondBend);

    OpenMM::AndersenThermostat *    thermostat  = new OpenMM::AndersenThermostat(
            temperature,      // kelvins
            FrictionInPerPs); // collision frequency in 1/picoseconds
    AndersenThermostat_Index = system->addForce(thermostat);

    // Create periodic box
    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond->setCutoffDistance(CutoffDistanceInAng * OpenMM::NmPerAngstrom);

    for (unsigned int i = 0; i < numOfMolecules; ++i) {
        addWaterMolecule();
    }

    integrator = new OpenMM::VerletIntegrator(StepSizeInFs * OpenMM::PsPerFs);
    context    = new OpenMM::Context(*system, *integrator);
}

void WaterSimulator::initSystemState(double startTemperatureInK, double densityInKgPerM3)
{
    const double boxLengthInNm = 30 * OpenMM::NmPerAngstrom;
    system->setDefaultPeriodicBoxVectors(Vec3(boxLengthInNm,0,0),
                                         Vec3(0,boxLengthInNm,0),
                                         Vec3(0,0,boxLengthInNm));

}

void WaterSimulator::getSystemState(double& timeInPs,
                                    std::vector<double>& atomPositionsInAng)
{
    const OpenMM::State state = context->getState(OpenMM::State::Positions, true);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy OpenMM positions into output array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    atomPositionsInAng.resize(3*positionsInNm.size());
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atomPositionsInAng[3*i+j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
}

void WaterSimulator::run(double SimulationTimeInPs) {
    context->setVelocitiesToTemperature(Temperature);

    // Run the simulation:
    //  (1) Write the first line of the PDB file and the initial configuration.
    //  (2) Run silently entirely within OpenMM between reporting intervals.
    //  (3) Write a PDB frame when the time comes.
    printf("REMARK  Using OpenMM platform %s\n", context->getPlatform().getName().c_str());

    std::vector<double> atomPositionsInAng; // x,y,z,x,y,z, ...
    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    for (int frame=1; ; ++frame)
    {
        double time;
        getSystemState(time, atomPositionsInAng);
        myWritePDBFrame(frame, time, atomPositionsInAng);

        if (time >= SimulationTimeInPs)
            break;

        context->getIntegrator().step(NumSilentSteps);
    }
}

void WaterSimulator::addWaterMolecule() {

    // Add atom masses to system
    int  oIndex = system->addParticle(O_mass); // O
    int h1Index = system->addParticle(H_mass); // H1
    int h2Index = system->addParticle(H_mass); // H2
    int  mIndex = system->addParticle(M_mass); // M
    OpenMM::VirtualSite *site = new OpenMM::ThreeParticleAverageSite(oIndex, h1Index, h2Index, 0.786646558, 0.106676721, 0.106676721);
    system->setVirtualSite(mIndex, site);
    int  xIndex = system->addParticle(X_mass); // X
    const double M = O_mass + 2 * H_mass;
    site = new OpenMM::ThreeParticleAverageSite(oIndex, h1Index, h2Index, O_mass/M, H_mass/M, H_mass/M);
    system->setVirtualSite(xIndex, site);

    OpenMM::NonbondedForce &nonbond = (OpenMM::NonbondedForce &)system->getForce(NonbondedForce_Index);
    // Add atom charge, sigma, and stiffness to nonbonded force
    nonbond.addParticle( // Oxygen
        O_charge,
        O_sigma,
        O_epsilon
        );
    nonbond.addParticle( // Hydrogen1
        H_charge,
        H_sigma,
        H_epsilon);
    nonbond.addParticle( // Hydrogen2
        H_charge,
        H_sigma,
        H_epsilon);
    nonbond.addParticle( // Negative charge center
        M_charge,
        M_sigma,
        M_epsilon);
    nonbond.addParticle( // Center of mass
        X_charge,
        X_sigma,
        X_epsilon
        );

    // Constrain O-H bond lengths or use harmonic forces.
    if (UseConstraints)
    {
        system->addConstraint(oIndex, h1Index,
                             OH_nominalLengthInAng * OpenMM::NmPerAngstrom);
        system->addConstraint(oIndex, h2Index,
                             OH_nominalLengthInAng * OpenMM::NmPerAngstrom);
    }
    else {
        // Add stretch parameters for two covalent bonds
        // Note factor of 2 for stiffness below because Amber specifies the constant
        // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
        // it as used in the force term kx, with energy kx^2/2.
        OpenMM::HarmonicBondForce &bondStretch = (OpenMM::HarmonicBondForce &)system->getForce(HarmonicBondForce_Index);
        bondStretch.addBond(oIndex, h1Index,
                            OH_nominalLengthInAng     * OpenMM::NmPerAngstrom,
                            OH_stiffnessInKcalPerAng2 * 2 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
        bondStretch.addBond(oIndex, h2Index,
                            OH_nominalLengthInAng     * OpenMM::NmPerAngstrom,
                            OH_stiffnessInKcalPerAng2 * 2 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
    }

    std::vector< std::pair<int,int> >   bondPairs;
    // Store bonds for exclusion list
    bondPairs.push_back(std::make_pair(oIndex, h1Index));
    bondPairs.push_back(std::make_pair(oIndex, h2Index));

    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    // Add bond bend parameters for one angle.
    // See note under bond stretch above regarding the factor of 2 here.
    OpenMM::HarmonicAngleForce &bondBend = (OpenMM::HarmonicAngleForce &)system->getForce(HarmonicAngleForce_Index);
    bondBend.addAngle(h1Index, oIndex, h2Index,
                      HOH_nominalAngleInDeg      * OpenMM::RadiansPerDegree,
                      HOH_stiffnessInKcalPerRad2 * 2 * OpenMM::KJPerKcal);
}

//                               PDB FILE WRITER
// This is a PDB writer that only knows how to write out water molecules. It is
// just here for this example and has nothing to do with OpenMM!
static void
myWritePDBFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng)
{
    const char* atomNames[] = {" O  ", " H1 ", " H2 ", " M  ", " X  "}; // cycle through these
    int N = sizeof(atomNames)/sizeof(atomNames[0]);
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f picoseconds\n", timeInPs);
    for (int atom=0; atom < (int)atomPosInAng.size()/3; ++atom)
    {
        printf("HETATM%5d %4s HOH  %4d    ",        // start of pdb HETATM line
            atom+1, atomNames[atom%N], 1 + atom/N); // atom number, name, residue #
        printf("%8.3f%8.3f%8.3f",                   // middle of pdb HETATM line
            atomPosInAng[3*atom+0], atomPosInAng[3*atom+1], atomPosInAng[3*atom+2]);
        printf("  1.00  0.00            \n");       // end of pdb HETATM line
    }
    printf("ENDMDL\n"); // end of trajectory frame
}

// -----------------------------------------------------------------------------
//                           WATER SIMULATOR MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        // Set up OpenMM data structures; return handle and OpenMM Platform name.
        WaterSimulator simulator(10, Temperature, StepSizeInFs);// output
        simulator.initSystemState(Temperature, 997.0);
        simulator.run(SimulationTimeInPs);

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}
