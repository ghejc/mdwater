#include "WaterSimulator.h"
#include "LorentzForceIntegrator.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include "RotationMatrix.h"

using OpenMM::Vec3;

static void
myWritePDBFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng);

static void
myWriteXYZFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng);

const double WaterSimulator::AtomicMassUnitsPerKg = 1.66054E-27;
const double WaterSimulator::MeterPerNm = 1.0E9;

//                     SIMULATION PARAMETERS
unsigned int NumberOfMolecules = 10;
double Temperature             = 300;   // Kelvins
double Density                 = 997;   // kg/m^3
double StepSizeInFs            = 0.2;   // integration step size (fs)
double ReportIntervalInFs      = 100;   // how often to generate an output frame (fs)
double SimulationTimeInPs      = 10;

double ElectricFieldAmplitude = 1000;       // electric field in V/m
double MagneticFieldAmplitude = 1;          // magnetic field in T
double ElectricFieldFrequencyInTHz = 0.001;        // frequency of the electric field in THz
double MagneticFieldFrequencyInTHz = 0.001;        // frequency of the magnetic field in THz
double ElectricFieldPhaseInPs = 0.0;        // phase offset of the electric field in ps
double MagneticFieldPhaseInPs = 0.0;        // phase offset of the magnetic field in ps

const double WaterSimulator::FrictionInPerPs     = 91.;    // collisions per picosecond
const double WaterSimulator::CutoffDistanceInAng = 10.;    // Angstroms

const double WaterSimulator::Coulomb14Scale = 0.833333;
const double WaterSimulator::LennardJones14Scale = 0.5;

// Values are from tip4pew.xml.

// Hydrogen
const double WaterSimulator::H_mass             = 1.007947;
const double WaterSimulator::H_charge           = 0.52422;
const double WaterSimulator::H_sigma            = 1;
const double WaterSimulator::H_epsilon          = 0;

// Negative charge center
#if defined(PARTICLE_M_IS_VIRTUAL)
const double WaterSimulator::M_mass             = 0;
const double WaterSimulator::H_weight           = 0.106676721;
const double WaterSimulator::O_weight           = 1 - 2 * H_weight;
#else
const double WaterSimulator::M_mass             = 0.5;
#endif
const double WaterSimulator::M_charge           = -2 * H_charge;
const double WaterSimulator::M_sigma            = 1;
const double WaterSimulator::M_epsilon          = 0;

// Oxygen
const double WaterSimulator::O_mass             = 15.99943 - M_mass;
const double WaterSimulator::O_charge           = 0.0;
const double WaterSimulator::O_sigma            = 0.316435;
const double WaterSimulator::O_epsilon          = 0.680946;

// Parameters for the O-M bonds.
const double WaterSimulator::OM_nominalLengthInAng      = 0.15;

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

    std::vector<double> charges;
    for (unsigned int i = 0; i < numOfMolecules; ++i) {
        addWaterMolecule();
        charges.push_back(O_charge);
        charges.push_back(H_charge);
        charges.push_back(H_charge);
        charges.push_back(M_charge);
    }

    // use a Lorentz force integrator with an electric field pointing in x direction
    // and a magnetic field pointing in y direction both with a given oscillator
    // frequency and phase offset
    integrator = new LorentzForceIntegrator(StepSizeInFs * OpenMM::PsPerFs,
                                            new ElectricField(Vec3(ElectricFieldAmplitude,0,0),
                                                              ElectricFieldFrequencyInTHz,
                                                              ElectricFieldPhaseInPs),
                                            new MagneticField(Vec3(0,MagneticFieldAmplitude,0),
                                                              MagneticFieldFrequencyInTHz,
                                                              MagneticFieldPhaseInPs),
                                            charges);
    context    = new OpenMM::Context(*system, *integrator);
}

static double random()
{
    return std::rand()/(double)RAND_MAX;
}

// puts all molecules at random positions and returns the minimum distance between
double WaterSimulator::setRandomPositions(const double boxLengthInNm)
{
    std::vector<Vec3> positions;
    size_t N = context->getMolecules().size();
    // relative coordinates of the atoms in the molecule
    std::vector<Vec3> coords;
    getCenterOfMassCoordinates(coords);
    double min_dist = boxLengthInNm;

    if (N == 1) {
        context->setPositions(coords);
        return min_dist;
    }
    // number of divisions of the box length
    size_t n = (size_t)std::max(1.0, ceil(pow(N, 1.0/3.0)));
    // cell size
    double cellLengthInNm = boxLengthInNm / (double)n;
    size_t i = 0;
    std::vector<size_t> cells;
    std::map<size_t,size_t> cell_map;
    while ( i < N ) {
        Vec3 randVec( floor(random()*n), floor(random()*n), floor(random()*n));
        size_t cell = (size_t)(randVec[0] + n * (randVec[1] + n * randVec[2]));
        if (std::find(cells.begin(), cells.end(), cell) != cells.end()) {
            continue;
        }
        cells.push_back(cell);
        randVec *= cellLengthInNm;
        double shift = 0.5 * cellLengthInNm;
        randVec[0] += shift;
        randVec[1] += shift;
        randVec[2] += shift;

        // rotate coordinates by a random angle around a random unit vector
        Vec3 randAxis(random(), random(), random());
        double randAngleInRad = random() * 2 * M_PI;
        RotationMatrix m(randAngleInRad, randAxis);
        for (std::vector<Vec3>::iterator it = coords.begin(); it != coords.end(); it++) {
           Vec3 x = m.dot(*it);
           positions.push_back(x + randVec);
        }
        cell_map.insert(std::make_pair(cell,i));
        i++;
    }
    // calculates the minimum distance between atoms
    for (size_t cell : cells) {
        std::map<size_t,size_t>::iterator it2;
        Vec3 pos;
        for( size_t  k = 0; k < coords.size(); k++) {
            if ((it2 = cell_map.find(cell)) != cell_map.end()) {
                pos = positions[it2->second + k];
            } else {
                continue;
            }

            std::vector<size_t> neighbor_cells;
            neighbor_cells.push_back(cell + 1);
            neighbor_cells.push_back(cell - 1);
            neighbor_cells.push_back(cell + n);
            neighbor_cells.push_back(cell - n);
            neighbor_cells.push_back(cell + n*n);
            neighbor_cells.push_back(cell - n*n);
            for (size_t neighbor : neighbor_cells) {
                if ((it2 = cell_map.find(neighbor)) != cell_map.end()) {
                    for (size_t j = 0; j < coords.size(); j++) {
                        Vec3 dist = pos - positions[it2->second + j];
                        double d = sqrt(dist.dot(dist));
                        if (d < min_dist) {
                            min_dist = d;
                        }
                    }
                } else {
                    continue;
                }
            }
        }
    }
    context->setPositions(positions);
    return min_dist;
}

void WaterSimulator::initSystemState(double startTemperatureInK, double densityInKgPerM3)
{
    size_t N = context->getMolecules().size();
    double V = N * (O_mass + 2 * H_mass + M_mass) * AtomicMassUnitsPerKg / densityInKgPerM3;
    const double boxLengthInNm = pow(V, 1./3) * MeterPerNm;

    // Create periodic box
    system->setDefaultPeriodicBoxVectors(Vec3(boxLengthInNm, 0, 0),
                                         Vec3(0, boxLengthInNm, 0),
                                         Vec3(0, 0, boxLengthInNm));
    /*
    context->setPeriodicBoxVectors(Vec3(boxLengthInNm, 0, 0),
                                   Vec3(0, boxLengthInNm, 0),
                                   Vec3(0, 0, boxLengthInNm));
    */
    /*
    if( system->usesPeriodicBoundaryConditions() ) {
        printf("System uses periodic boundary conditions with a box of side length %5.3f nm\n", boxLengthInNm);
    } else {
        printf("System does not use periodic boundary conditions\n");
    }
    printf("The box contains %d water molecules\n", N);
    */

    setRandomPositions(boxLengthInNm);
    context->setVelocitiesToTemperature(startTemperatureInK);
}

void WaterSimulator::getSystemState(double& timeInPs,
                                    std::vector<double>& atomPositionsInAng)
{
    const OpenMM::State state = context->getState(OpenMM::State::Positions, false);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy only non-virtual OpenMM positions into output array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    size_t N = context->getMolecules().size();
#if defined(PARTICLE_M_IS_VIRTUAL)
    const int numberOfAtomsPerMolecule = 3; // no virtual sites
#else
    const int numberOfAtomsPerMolecule = 4; // including virtual sites
#endif // PARTICLE_M_IS_VIRTUAL
    atomPositionsInAng.resize(3 * numberOfAtomsPerMolecule * N);
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
    double M = O_mass + 2 * H_mass + M_mass;
    double r_OM = OM_nominalLengthInAng * OpenMM::NmPerAngstrom;
#ifdef PARTICLE_M_IS_VIRTUAL
    Vec3 r_M = r_O * O_weight + (r_H1 + r_H2) * H_weight;
#else
    Vec3 r_M(0,-r_OM,0);
#endif // PARTICLE_M_IS_VIRTUAL
    Vec3 r_X = r_O * O_mass/M + (r_H1 + r_H2) * H_mass/M + r_M/M * M_mass;
    coords.push_back(r_O - r_X);
    coords.push_back(r_H1 - r_X);
    coords.push_back(r_H2 - r_X);
    coords.push_back(r_M - r_X);
}

void WaterSimulator::run(double SimulationTimeInPs) {
    // Run the simulation:
    //  (1) Write the first line of the PDB file and the initial configuration.
    //  (2) Run silently entirely within OpenMM between reporting intervals.
    //  (3) Write a PDB frame when the time comes.
    // printf("Using OpenMM platform %s\n", context->getPlatform().getName().c_str());

    std::vector<double> atomPositionsInAng; // x,y,z,x,y,z, ...
    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    for (int frame=1; ; ++frame)
    {
        double time;
        getSystemState(time, atomPositionsInAng);
        myWriteXYZFrame(frame, time, atomPositionsInAng);

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
    // Add virtual sites
#ifdef PARTICLE_M_IS_VIRTUAL
    OpenMM::VirtualSite *site = new OpenMM::ThreeParticleAverageSite(oIndex, h1Index, h2Index, O_weight, H_weight, H_weight);
    system->setVirtualSite(mIndex, site);
#endif // PARTICLE_M_IS_VIRTUAL

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

    // Constrain O-M bond length
#ifndef PARTICLE_M_IS_VIRTUAL
    system->addConstraint(oIndex, mIndex,
                          OM_nominalLengthInAng * OpenMM::NmPerAngstrom);
#endif // PARTICLE_M_IS_VIRTUAL

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
#ifdef PARTICLE_M_IS_VIRTUAL
    const char* atomNames[] = {" O  ", " H1 ", " H2 "}; // cycle through these
#else
    const char* atomNames[] = {" O  ", " H1 ", " H2 ", " M  "}; // cycle through these
#endif // PARTICLE_M_IS_VIRTUAL
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

//                               XYZ FILE WRITER
static void
myWriteXYZFrame(int frameNum, double timeInPs, const std::vector<double>& atomPosInAng)
{
#ifdef PARTICLE_M_IS_VIRTUAL
    const char* atomNames[] = {"O", "H", "H"}; // cycle through these
#else
    const char* atomNames[] = {"O", "H", "H", "M"}; // cycle through these
#endif // PARTICLE_M_IS_VIRTUAL
    int N = sizeof(atomNames)/sizeof(atomNames[0]);
    printf("%zd\n", atomPosInAng.size()/3 );
    printf("frame=%d time=%.3f ps\n", frameNum, timeInPs);
    for (int atom=0; atom < (int)atomPosInAng.size()/3; ++atom)
    {
        printf("%2s %8.3f%8.3f%8.3f\n",        // start of pdb HETATM line
               atomNames[atom%N],
               atomPosInAng[3*atom+0],
               atomPosInAng[3*atom+1],
               atomPosInAng[3*atom+2]);
    }
}

// -----------------------------------------------------------------------------
//                           WATER SIMULATOR MAIN PROGRAM
// -----------------------------------------------------------------------------
int main(int argc, char **argv) {
    if (argc > 1) {
        // parse command line arguments
        std::map<std::string,double> options;
        for(int i = 0; i < argc; i++) {
            std::string arg(argv[i]);
            std::size_t pos = arg.find("=");
            if (pos != std::string::npos) {
                double val = atof(arg.substr(pos+1, std::string::npos).c_str());
                options.insert(std::make_pair(arg.substr(0,pos),val));
            }
        }
        std::map<std::string,double>::iterator it;
        it = options.find("NumberOfMolecules");
        if (it != options.end()) {
            NumberOfMolecules = (int)(it->second);
        }
        it = options.find("Temperature");
        if (it != options.end()) {
            Temperature = it->second;
        }
        it = options.find("SimulationTimeInPs");
        if (it != options.end()) {
            SimulationTimeInPs = it->second;
        }
        it = options.find("StepSizeInFs");
        if (it != options.end()) {
            StepSizeInFs = it->second;
        }
        it = options.find("Density");
        if (it != options.end()) {
            Density = it->second;
        }
        it = options.find("ReportIntervalInFs");
        if (it != options.end()) {
            ReportIntervalInFs = it->second;
        }
        it = options.find("ElectricFieldAmplitude");
        if (it != options.end()) {
            ElectricFieldAmplitude = it->second;
        }
        it = options.find("MagneticFieldAmplitude");
        if (it != options.end()) {
            MagneticFieldAmplitude = it->second;
        }
        it = options.find("ElectricFieldFrequencyInTHz");
        if (it != options.end()) {
            ElectricFieldFrequencyInTHz = it->second;
        }
        it = options.find("MagneticFieldFrequencyInTHz");
        if (it != options.end()) {
            MagneticFieldFrequencyInTHz = it->second;
        }
        it = options.find("ElectricFieldPhaseInPs");
        if (it != options.end()) {
            ElectricFieldPhaseInPs = it->second;
        }
        it = options.find("MagneticFieldPhaseInPs");
        if (it != options.end()) {
            MagneticFieldPhaseInPs = it->second;
        }
    }

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
