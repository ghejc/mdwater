# mdwater simulation script
#
# mdwater simulation parameter:
#   NumberOfMolecules=10
#   Temperature=300          // Kelvins
#   Density=997              // kg/m^3
#   StepSizeInFs=0.2         // integration step size (fs)
#   ReportIntervalInFs=100   // how often to generate an output frame (fs)
#   SimulationTimeInPs=10
#   ElectricFieldAmplitude=1000         // electric field in V/m
#   MagneticFieldAmplitude=1            // magnetic field in T
#   ElectricFieldFrequencyInTHz=0.001   // frequency of the electric field in THz
#   MagneticFieldFrequencyInTHz=0.001   // frequency of the magnetic field in THz
#   ElectricFieldPhaseInPs=0.0          // phase offset of the electric field in ps
#   MagneticFieldPhaseInPs=0.0          // phase offset of the magnetic field in ps

import os
import subprocess

cwd = os.getcwd()
cmd = '{}\\bin\\Release\\mdwater.exe'.format(cwd)
with open('{}\\output.xyz'.format(cwd), "w") as outfile:
    subprocess.call([cmd,
    "NumberOfMolecules=10",
    "SimulationTimeInPs=10",
    "StepSizeInFs=0.2",
    "ReportIntervalInFs=2",
    "ElectricFieldAmplitude=0.0",
    "MagneticFieldAmplitude=0.0",
    "ElectricFieldFrequencyInTHz=0.0",
    "MagneticFieldFrequencyInTHz=0.0"], stdout=outfile)
