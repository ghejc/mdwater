#include "WaterSimulator.h"
#include <string>
#include <vector>
#include <cstdlib>

using OpenMM::Vec3; // so we can just say "Vec3" below

static void
myWritePDBFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng);

const double WaterSimulator::AtomicMassUnitsPerKg = 1.66054E-27;
const double WaterSimulator::MeterPerNm = 1.0E9;

//                     SIMULATION PARAMETERS
const unsigned int NumberOfMolecules = 10;
const double Temperature         = 300;    // Kelvins
const double Density             = 997;    // kg/m^3
const double StepSizeInFs        = 2;      // integration step size (fs)
const double ReportIntervalInFs  = 100;    // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 10;

const double WaterSimulator::FrictionInPerPs     = 91.;    // collisions per picosecond
const double WaterSimulator::CutoffDistanceInAng = 10.;    // Angstroms

const double WaterSimulator::Coulomb14Scale = 0.833333;
const double WaterSimulator::LennardJones14Scale = 0.5;

// Values are from tip4pew.xml.
// Oxygen
const double WaterSimulator::O_mass             = 15.99943;
const double WaterSimulator::O_charge           = 0.0;
const double WaterSimulator::O_sigma            = 1;
const double WaterSimulator::O_epsilon          = 0;

// Hydrogen
const double WaterSimulator::H_mass             = 1.007947;
const double WaterSimulator::H_charge           = 0.52422;
const double WaterSimulator::H_sigma            = 1;
const double WaterSimulator::H_epsilon          = 0;

// Negative charge center
const double WaterSimulator::M_mass             = 0;
const double WaterSimulator::M_charge           = -2 * H_charge;
const double WaterSimulator::M_sigma            = 1;
const double WaterSimulator::M_epsilon          = 0;
const double WaterSimulator::O_weight           = 0.786646558;
const double WaterSimulator::H_weight           = 0.106676721;

// Center of mass
const double WaterSimulator::X_mass             = 0;
const double WaterSimulator::X_charge           = 0;
const double WaterSimulator::X_sigma            = 0.316435;
const double WaterSimulator::X_epsilon          = 0.680946;

// Parameters for the O-H bonds.
const double WaterSimulator::OH_nominalLengthInAng      = 0.9572;
const double WaterSimulator::OH_stiffnessInKcalPerAng2  = 553.0; // that is, e=k(x-x0)^2

// Parameters for the H-O-H angle.
const double WaterSimulator::HOH_nominalAngleInDeg      = 104.52;
const double WaterSimulator::HOH_stiffnessInKcalPerRad2 = 100.; // that is e=k(a-a0)^2

// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins
// (2) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (3) Create an Integrator and a Context associating the Integrator with
//     the System.
// (4) Select the OpenMM platform to be used.
// -----------------------------------------------------------------------------

WaterSimulator::WaterSimulator ( unsigned int        numOfMolecules,
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

    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond->setCutoffDistance(CutoffDistanceInAng * OpenMM::NmPerAngstrom);

    for (unsigned int i = 0; i < numOfMolecules; ++i) {
        addWaterMolecule();
    }

    integrator = new OpenMM::VerletIntegrator(StepSizeInFs * OpenMM::PsPerFs);
    context    = new OpenMM::Context(*system, *integrator);
}

void WaterSimulator::setRandomPositions(const double boxLengthInNm)
{
    std::vector<Vec3> positions;
    size_t N = context->getMolecules().size();
    std::vector<Vec3> coords;
    getCenterOfMassCoordinates(coords);

    if (N == 1) {
        context->setPositions(coords);
        return;
    }
    size_t i = 0;
    while( i < N ) {
        Vec3 x(rand()/(double)RAND_MAX * boxLengthInNm,rand()/(double)RAND_MAX * boxLengthInNm, rand()/(double)RAND_MAX * boxLengthInNm);
        i++;
        // rotate coords by a random angle around a random unit vector
        for (std::vector<Vec3>::iterator it = coords.begin(); it != coords.end(); it++) {
            positions.push_back(*it + x);
        }
    }
    context->setPositions(positions);

}

void WaterSimulator::initSystemState(double startTemperatureInK, double densityInKgPerM3)
{
    size_t N = context->getMolecules().size();
    double V = N * (O_mass + 2 * H_mass) * AtomicMassUnitsPerKg / densityInKgPerM3;
    const double boxLengthInNm = pow(V, 1./3) * MeterPerNm;

    // Create periodic box
    system->setDefaultPeriodicBoxVectors(Vec3(boxLengthInNm, 0, 0),
                                         Vec3(0, boxLengthInNm, 0),
                                         Vec3(0, 0, boxLengthInNm));
    if( system->usesPeriodicBoundaryConditions() ) {
        printf("REMARK  System uses periodic boundary conditions with a box of side length %5.3f nm\n", boxLengthInNm);
    } else {
        printf("REMARK  System does not use periodic boundary conditions\n");
    }
    printf("REMARK  The box contains %d water molecules\n", N);

    setRandomPositions(boxLengthInNm);
    context->setVelocitiesToTemperature(startTemperatureInK);
}

void WaterSimulator::getSystemState(double& timeInPs,
                                    std::vector<double>& atomPositionsInAng)
{
    const OpenMM::State state = context->getState(OpenMM::State::Positions, true);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy only non-virtual OpenMM positions into output array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    size_t N = context->getMolecules().size();
    atomPositionsInAng.resize(9 * N);
    int k = 0;
    for (int i = 0; i < (int)positionsInNm.size(); ++i) {
        if (!system->isVirtualSite(i)) {
            for (int j=0; j < 3; ++j)
                atomPositionsInAng[3*k+j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;
            k++;
        }
    }
}

void WaterSimulator::getCenterOfMassCoordinates(std::vector<Vec3> &coords)
{
    Vec3 r_O(0, 0, 0);
    double r_OH = OH_nominalLengthInAng * OpenMM::NmPerAngstrom;
    double angle_HO = HOH_nominalAngleInDeg * OpenMM::RadiansPerDegree / 2.0;
    Vec3 r_H1(-r_OH * sin(angle_HO), -r_OH * cos(angle_HO), 0);
    Vec3 r_H2(+r_OH * sin(angle_HO), -r_OH * cos(angle_HO), 0);
    double M = O_mass + 2 * H_mass;
    Vec3 r_X = r_O * O_mass/M + (r_H1 + r_H2) * H_mass/M;
    Vec3 r_M = r_O * O_weight + (r_H1 + r_H2) * H_weight;
    coords.push_back(r_O - r_X);
    coords.push_back(r_H1 - r_X);
    coords.push_back(r_H2 - r_X);
    coords.push_back(r_M - r_X);
    coords.push_back(r_X - r_X);
}

void WaterSimulator::run(double SimulationTimeInPs) {
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

    // Add virtual sites
    int  mIndex = system->addParticle(M_mass); // M
    OpenMM::VirtualSite *site = new OpenMM::ThreeParticleAverageSite(oIndex, h1Index, h2Index, O_weight, H_weight, H_weight);
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
    const char* atomNames[] = {" O  ", " H1 ", " H2 "}; // cycle through these
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
        WaterSimulator simulator(NumberOfMolecules, Temperature, StepSizeInFs);// output
        simulator.initSystemState(Temperature, Density);
        simulator.run(SimulationTimeInPs);

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}
