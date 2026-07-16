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
ENV_RES = "WAT,HOH,TIP3,SPC,SPCE,Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca2+,Zn2+,POT,CLA,POPC,OL,PA,PC"
WATER_RES = 'WAT,HOH,TIP3,TIP4,TIP5,SPC,SPCE,T4E,OPC'
ION_RES   = 'Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca2+,Zn2+,POT,CLA'

RESTRAINS_DIC = {'Protein+Ligand': f':{PROTEIN_RES},LIG & !@H=', 'Protein only': f':{PROTEIN_RES} & !@H=',
                 'Ligand': ':LIG', 'Ligand+Backbone': ':LIG | @CA,C,N,O',
                 'Backbone': '@CA,C,N,O', 'CA': '@CA', 'Ligand+CA': ':LIG  | @CA',
                 'Everything except wat+ions': '!(:WAT,HOH,TIP3,Na+,Cl-,K+)'}

# Default workflows for different system types
# Protein workflow
PROTWORK = "{'stepType': 'Minimization', 'MaxCycles': 50000, 'SdCycles': 35000, 'IntCutoff': 8.0,'Restraint': False, 'CustomIn': None}\n" \
           "{'stepType': 'Heating', 'MDSteps': 20000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 0, 'FiTemp': 300, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Backbone', 'RestrForce': 25.0, 'CustomIn': None}\n" \
           "{'stepType': 'Production', 'MDSteps': 10000, 'TimeStep': 0.002, 'TrajStep': 1000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': True, 'RestrAtoms': 'Backbone', 'RestrForce': 0.0, 'CustomIn': None}\n" \
           "{'stepType': 'Production', 'MDSteps': 1000000, 'TimeStep': 0.002, 'TrajStep': 10000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': False, 'RestrForce': 0.0, 'CustomIn': None}\n"

# Protein + ligand workflow
PROTLIGWORK = "{'stepType': 'Minimization', 'MaxCycles': 50000, 'SdCycles': 35000, 'IntCutoff': 8.0,'Restraint': False, 'CustomIn': None, 'stepType': 'Minimization'}\n" \
              "{'stepType': 'Heating', 'MDSteps': 20000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 0, 'FiTemp': 300, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 25.0, 'CustomIn': None, 'stepType': 'Heating'}\n" \
              "{'stepType': 'Production', 'MDSteps': 10000, 'TimeStep': 0.002, 'TrajStep': 1000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 0.0, 'CustomIn': None}\n" \
              "{'stepType': 'Production', 'MDSteps': 1000000, 'TimeStep': 0.002, 'TrajStep': 10000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': False, 'RestrAtoms': 'Protein', 'RestrForce': 0.0, 'CustomIn': None}\n"

# Membrane protein workflow
MEMPROTWORK = "{'stepType': 'Minimization', 'MaxCycles': 2000, 'SdCycles': 2000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None}\n" \
              "{'stepType': 'Minimization', 'MaxCycles': 10000, 'SdCycles': 5000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None}\n" \
              "{'stepType': 'Minimization', 'MaxCycles': 100000, 'SdCycles': 90000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Protein only', 'RestrForce': 50.0, 'CustomIn': None}\n" \
              "{'stepType': 'Heating', 'MDSteps': 100000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 0.0, 'FiTemp': 100.0, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0}\n" \
              "{'stepType': 'Heating', 'MDSteps': 1000000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 100.0, 'FiTemp': 300.0, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Protein only', 'RestrForce': 50.0}\n" \
              "{'stepType': 'Production', 'MDSteps': 100000, 'TimeStep': 0.002, 'TrajStep': 1000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Protein only', 'RestrForce': 50.0}\n" \
              "{'stepType': 'Production', 'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'CA', 'RestrForce': 1.0, 'CustomIn': None}\n" \
              "{'stepType': 'Production', 'MDSteps': 50000000, 'TimeStep': 0.002, 'TrajStep': 50000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None}\n" \
              "{'stepType': 'Production', 'MDSteps': 125000000, 'TimeStep': 0.002, 'TrajStep': 125000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None}\n"

# Membrane protein + ligand workflow
MEMPROTLIGWORK = "{'stepType': 'Minimization', 'MaxCycles': 2000, 'SdCycles': 2000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None}\n" \
                 "{'stepType': 'Minimization', 'MaxCycles': 10000, 'SdCycles': 5000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0, 'CustomIn': None}\n" \
                 "{'stepType': 'Minimization', 'MaxCycles': 100000, 'SdCycles': 90000, 'IntCutoff': 8.0, 'Restraint': True, 'RestrAtoms': 'Protein+Ligand', 'RestrForce': 50.0, 'CustomIn': None}\n" \
                 "{'stepType': 'Heating', 'MDSteps': 100000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 0.0, 'FiTemp': 100.0, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Everything except wat+ions', 'RestrForce': 50.0}\n" \
                 "{'stepType': 'Heating', 'MDSteps': 1000000, 'TimeStep': 0.001, 'Traj': 1000, 'SaveTrj': False, 'InTemp': 100.0, 'FiTemp': 300.0, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 50.0}\n" \
                 "{'stepType': 'Production', 'MDSteps': 100000, 'TimeStep': 0.002, 'TrajStep': 1000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'CustomIn': None, 'Restraint': True, 'RestrAtoms': 'Ligand+Backbone', 'RestrForce': 50.0}\n" \
                 "{'stepType': 'Production', 'MDSteps': 500000, 'TimeStep': 0.002, 'TrajStep': 5000, 'SaveTrj': False, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': True, 'RestrAtoms': 'CA', 'RestrForce': 1.0, 'CustomIn': None}\n" \
                 "{'stepType': 'Production', 'MDSteps': 50000000, 'TimeStep': 0.002, 'TrajStep': 50000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None}\n" \
                 "{'stepType': 'Production', 'MDSteps': 125000000, 'TimeStep': 0.002, 'TrajStep': 125000, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'semiisotropic', 'Restraint': False, 'RestrAtoms': '', 'RestrForce': 0.0, 'CustomIn': None}\n"
