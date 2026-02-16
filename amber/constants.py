# **************************************************************************
# *
# * Authors: Aida Pinacho Pérez
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


AMBER_HOME = 'AMBER_HOME'

AMBER = 'amber'
V2020 = '2020.1'
V2025 = '2025.1'

AMBER_DEFAULT_VERSION = V2025

VMD_HOME = 'VMD_HOME'

AMBER_DIC = {'name': 'amber', 'version': AMBER_DEFAULT_VERSION, 'home': 'AMBER_HOME', 'pmemd_home': 'PMEMD_HOME'}

PROTEIN_RES = "ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL," \
              "HID,HIE,HIP,CYX,ASH,GLH,LYN,ARN,ACE,NME,NHE"
ENV_RES = "WAT,HOH,TIP3,SPC,SPCE,Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca2+,Zn2+"
RESTRAINS_DIC = {'Protein+Ligand': f':{PROTEIN_RES},LIG & !@H=', 'Protein only': f':{PROTEIN_RES} & !@H=',
                 'Ligand': ':LIG', 'Ligand+Backbone': ':LIG | @CA,C,N,O',
                 'Backbone': '@CA,C,N,O', 'CA': '@CA', 'Ligand+CA': ':LIG  | @CA', 'Everything except wat+ions': '!(:WAT,HOH,TIP3,Na+,Cl-,K+)'}

# Default workflows for different system types
# Protein workflow
PROTWORK = "{'MaxCycles': 500, 'SdCycles': 250, 'IntCutoff': 8.0,'Restraint': False, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
           "{'MDSteps': 10000, 'TimeStep': 0.002, 'Traj': 1000, 'InTemp': 0, 'FiTemp': 300, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Backbone', 'CustomIn': None, 'stepType': 'Heating'}\n" \
           "{'MDSteps': 10000, 'TimeStep': 0.002, 'TrajStep': 1000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': True, 'RestrAtoms': 'Backbone', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
           "{'MDSteps': 1000000, 'TimeStep': 0.002, 'TrajStep': 10000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': False, 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n"

# Protein + Ligand workflow
PROTLIGWORK = "{'MaxCycles': 500, 'SdCycles': 250, 'IntCutoff': 8.0,'Restraint': False, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
"{'MDSteps': 10000, 'TimeStep': 0.002, 'Traj': 1000, 'InTemp': 0, 'FiTemp': 300, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'CustomIn': None, 'stepType': 'Heating'}\n" \
"{'MDSteps': 10000, 'TimeStep': 0.002, 'TrajStep': 1000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 1000000, 'TimeStep': 0.002, 'TrajStep': 10000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': False, 'RestrAtoms': 'Protein', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n"

# Membrane Protein workflow
MEMPROTWORK = "{'MaxCycles': 2000, 'SdCycles': 2000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 500.0, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
"{'MaxCycles': 10000, 'SdCycles': 5000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
"{'MDSteps': 25000, 'TimeStep': 0.002, 'Traj': 1000, 'InTemp': 0, 'FiTemp': 303, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 5.0, 'CustomIn': None, 'stepType': 'Heating'}\n" \
"{'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'Protein', 'RestrForce': 2.5, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'CA', 'RestrForce': 1.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 50000000, 'TimeStep': 0.002, 'TrajStep': 50000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 125000000, 'TimeStep': 0.004, 'TrajStep': 125000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n"

# Membrane Protein + Ligand workflow
MEMPROTLIGWORK = "{'MaxCycles': 2000, 'SdCycles': 2000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 500.0, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
"{'MaxCycles': 10000, 'SdCycles': 5000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
"{'MDSteps': 25000, 'TimeStep': 0.002, 'Traj': 1000, 'InTemp': 0, 'FiTemp': 303, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 5.0, 'CustomIn': None, 'stepType': 'Heating'}\n" \
"{'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 2.5, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'Ligand+CA', 'RestrForce': 1.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 50000000, 'TimeStep': 0.002, 'TrajStep': 50000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n" \
"{'MDSteps': 125000000, 'TimeStep': 0.004, 'TrajStep': 125000, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Simulation'}\n"

