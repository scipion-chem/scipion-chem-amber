# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Aida Pinacho Pérez
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
"""
This module will prepare the system for the simulation
"""
from os.path import abspath


from pwem.protocols import EMProtocol, ProtImportFiles
from pyworkflow.protocol import params
from pyworkflow.utils import Message

import amber
from pwchem.utils import *

import amber.objects as amberobj
from amber.objects import *

waterFFDic = {'tip3p': 'TIP3PBOX', 'tip4pew': 'TIP4PEWBOX', 'spce': 'SPCEBOX', 'opc': 'OPCBOX',
              'opc3': 'OPCBOX'}
chargeDic = {'RESP': 'resp', 'AM1-BCC': 'bcc', 'CM1': 'cm1', 'CM2': 'cm2', 'ESP': 'esp',
             'Mulliken': 'mul', 'Gasteiger': 'gas'}
statusDic = {'brief': 0, 'deafult': 1, 'verbose': 2}

ATOMSTRUCT, LIGAND = 0, 1

class AmberSystemPrep(EMProtocol):
    """
    This protocol will prepare a system for MD simulation. It will clean the input PDB for further analysis
    and generate topology and coordinate files necessary for MD simulation.

    """

    _label = 'system preparation'
    IMPORT_FROM_FILE = 1
    IMPORT_FROM_SCIPION = 1

    _ChargeModel = ['RESP', 'AM1-BCC', 'CM1', 'CM2', 'ESP', 'Mulliken', 'Gasteiger', ]
    _Status = ['brief', 'default', 'verbose']

    # -------------------------- DEFINE param functions ----------------------

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtImportFiles.__init__(self, **kwargs)


    def _defineParams(self, form):
        """
        Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('fromInput', params.EnumParam, choices=['AtomStructure', 'Ligand'], 
                      label='Input from: ', default=0,
                      help='Whether to use a simple atomic structure or a docked ligand together with the atom '
                           'structure it is docked to')

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure:", allowsNull=True,
                      important=True, pointerClass='AtomStruct', condition='fromInput==0',
                      help='Atom structure to convert to Amber system')
        form.addParam('inputLigands', params.PointerParam,
                      label="Import docked small molecules:", allowsNull=True, condition='fromInput==1',
                      important=True, pointerClass='SetOfSmallMolecules')
        form.addParam('inputLigandSelect', params.StringParam,
                      label="Select ligand structure:", condition='fromInput==1',
                      important=True, pointerClass='SetOfSmallMolecules')

        group = form.addGroup('target modification options')
        group.addParam('proteinResidues', params.BooleanParam, default=False,
                       label='Keep only protein residues: ')
        group.addParam('AmberCompatibleResidues', params.BooleanParam, default=False,
                       label='Keep only Amber compatible residues: ')
        group.addParam('phSimulation', params.BooleanParam, default=False,
                       label='Rename GLU, ASP, HIS for constant pH simulation: ')
        group.addParam('reduce', params.BooleanParam, default=False,
                       label='Run reduce first to add hydrogens: ')
        group.addParam('tleap', params.BooleanParam, default=False,
                       label='Use tleap to add missing atoms (EXPERIMENTAL): ')

        group = form.addGroup('ligand modifications options', condition='fromInput==1')
        group.addParam('proteinResidues', params.BooleanParam, default=False,
                       label='Keep only protein residues: ')
        group.addParam('AmberCompatibleResidues', params.BooleanParam, default=False,
                       label='Keep only Amber compatible residues: ')
        group.addParam('phSimulation', params.BooleanParam, default=False,
                       label='Rename GLU, ASP, HIS for constant pH simulation: ')
        group.addParam('reduce', params.BooleanParam, default=False,
                       label='Run reduce first to add hydrogens: ', help='The addition of hydrogen helps to find the '
                                                                         'hydrogen bond interactions and more favorable '
                                                                         'to us to find binding affinity of ligand '
                                                                         'against protein.')
        group.addParam('tleap', params.BooleanParam, default=False,
                       label='Use tleap to add missing atoms (EXPERIMENTAL): ')

        group = form.addGroup('Ligand parametrization', condition='fromInput==1')
        group.addParam('ChargeModel', params.EnumParam,
                       choices=self._ChargeModel,
                       label='Choose the charge model in order to calculate the atomic point charges: ')
        group.addParam('Status', params.EnumParam,
                       choices=self._Status,
                       label='Choose status information: ')

        form.addSection('MD prep')
        group = form.addGroup('Force field', help='Force field applied to the system. Force fields are sets of '
                                                  'potential functions and '
                                                  'parametrized interactions that can be used to study physical '
                                                  'systems. You should select as '
                                                  'many force fields as molecules in your system (i.e protein + ligand')

        group.addParam('ProteinForceField', params.BooleanParam,
                       label='Protein Force Field')
        group.addParam('ProteinForceFieldType', params.EnumParam,
                       label='Type',
                       choices=['ff14SB', 'ff19SB', 'ff14SBonlysc', 'ff15ipq', 'fb15', 'ff03.r1', 'ff03ua'],
                       condition='ProteinForceField')
        group.addParam('LigandForceField', params.BooleanParam, default= False, 
                       condition='fromInput==1', label='Ligand Force Field', 
                       help= 'if you have chosen to introduce a ligand, this force field is mandatory')
        group.addParam('DNAForceField', params.BooleanParam,
                       label='DNA Force Field')
        group.addParam('RNAForceField', params.BooleanParam,
                       label='RNA Force Field')
        group.addParam('RNAForceFieldType', params.EnumParam, condition='RNAForceField',
                       label='Type', choices=['OL3', 'LJbb', 'YIL', 'ROC', 'Shaw'])
        group.addParam('LipidForceField', params.BooleanParam,
                       label='Lipid Force Field')
        group.addParam('WaterForceField', params.EnumParam,
                       choices=['tip4pew', 'spce', 'spceb', 'opc', 'opc3', 'tip3p'],      
                       label='Water Force Field: ', help='Force field applied to the water')

        group = form.addGroup('Disulfide bridges')
        group.addParam('DisulfideBridges', params.BooleanParam,
                       label='Are there any S-S bridges?',
                       help='Residues involved must be renamed to CYX in the pdb file')
        group.addParam('DisulfideBridgesNumber', params.StringParam,
                       condition='DisulfideBridges',
                       label='Number of the residues involved in the disulfide bridge \n'
                             'with format 1º Residue - 2º Residue / 1º Residue - 2º Residue')

        group = form.addGroup('Solvate')
        group.addParam('SolvateStep', params.EnumParam,
                       choices=['Cubic', 'Octahedric'], defalult='Octahedric',
                       label='Solvation box',
                       help='Both solvation boxes will be isometric')
        line = group.addLine('Box size:')
        line.addParam('Distance', params.FloatParam,
                      default=20.0)

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        if self.fromInput == LIGAND:
            #self._insertFunctionStep('prepStep')
            self._insertFunctionStep('antechamberStep')
            self._insertFunctionStep('parmStep')
            self._insertFunctionStep('leapStep')
        self._insertFunctionStep('pdb4AmberStep')
        self._insertFunctionStep('forceFieldStep')
        self._insertFunctionStep('createOutputStep')

    def prepStep(self):
        myMolFile = self.getLigandFile()
        if not myMolFile.endswith('.pdb'):
            myMolFile = self.convertPDB(myMolFile)

        systemBasename = self.getSystemName()

        params = '{} > {}.LIG.pdb '.format(myMolFile, systemBasename)

        if self.proteinResidues:
            params += '-p '
        if self.AmberCompatibleResidues:
            params += '-a '
        if self.phSimulation:
            params += '--constantph '
        if self.reduce:
            params += '--reduce '
        if self.tleap:
            params += '--add-missing-atoms '

        amber.Plugin.runAmbertools(self, 'pdb4amber', params, cwd=self._getPath())

    def antechamberStep(self):
        myMolFile = self.getLigandFile()
        if not myMolFile.endswith('.pdb'):
            myMolFile = self.convertPDB(myMolFile)
        systemBasename = self.getSystemName()

        params = ' -i {} -fi pdb -o {}.LIG.mol2 -fo mol2 '.format(myMolFile, systemBasename)
        
        params += f"-c {chargeDic[self.getEnumText('ChargeModel')]} "
        params += f"-s {statusDic[self.getEnumText('Status')]} "

        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getPath())

    def parmStep(self):
        systemBasename = self.getSystemName()

        params = '-i {}.LIG.mol2 -o {}.LIG.frcmod -f mol2 '.format(*[systemBasename]*2)

        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self._getPath())

    def leapStep(self):
        systemBasename = self.getSystemName()

        params = 'source leaprc.gaff \n' \
                 'loadamberparams {}.LIG.frcmod \n' \
                 'LIG = loadmol2 {}.LIG.mol2 \n' \
                 'saveoff LIG {}.LIG.lib \n' \
                 'saveamberparm LIG {}.LIG.top {}.LIG.crd \n' \
                 'savepdb LIG {}.checkLIG.pdb \n' \
                 'quit'.format(*[systemBasename]*6)

        file = open(self._getExtraPath("leap_commandsLIG.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f extra/leap_commandsLIG.txt", cwd=self._getPath())
        os.rename(self._getPath('leap.log'), self._getPath('leapLIG.log'))

    def pdb4AmberStep(self):
        inputStructure = self.getProteinFile()
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)
        systemBasename = getBaseFileName(inputStructure)

        params = '{} > {}.amber.pdb -y '.format(inputStructure, systemBasename)

        if self.proteinResidues:
            params += '-p '
        if self.AmberCompatibleResidues:
            params += '-a '
        if self.phSimulation:
            params += '--constantph '
        if self.reduce:
            params += '--reduce '
        if self.tleap:
            params += '--add-missing-atoms '

        amber.Plugin.runAmbertools(self, 'pdb4amber', params, cwd=self._getPath())

    def forceFieldStep(self):
        systemBasename = self.getSystemName()

        params = '\n'

        if self.ProteinForceField:
            params += 'source leaprc.protein.{} \n'.format(self.getEnumText('ProteinForceFieldType'))
        if self.LigandForceField:
            params += 'source leaprc.gaff2 \n'
        if self.DNAForceField:
            params += 'source leaprc.DNA.OL15 \n'
        if self.RNAForceField:
            params += 'source leaprc.RNA.{} \n'.format(self.getEnumText('RNAForceFieldType'))
        if self.LipidForceField:
            params += 'source leaprc.lipid17 \n'
        if self.WaterForceField:
            params += 'source leaprc.water.{} \n'.format(self.getEnumText('WaterForceField'))

        params += 'APO = loadPdb {}.amber.pdb \n'.format(systemBasename)

        if self.fromInput == LIGAND:

            params += 'source leaprc.gaff \n' \
                      'loadamberparams {}.LIG.frcmod \n' \
                      'loadOff {}.LIG.lib \n' \
                      'LIG = loadmol2 {}.LIG.mol2 \n'.format(*[systemBasename]*3)

        if self.DisulfideBridges:

            for pair in self.DisulfideBridgesNumber.get().split('/'):
                first = pair.split('-')[0]
                second = pair.split('-')[1]
                params += 'bond APO.{}.SG APO.{}.SG \n'.format(first, second)

        if self.getEnumText('SolvateStep') == 'Cubic':
          Boxtype = 'SolvateBox'
        else:
          Boxtype = 'SolvateOct'

        if self.fromInput == LIGAND:
            complexName = 'COMPL'
            params += 'COMPL = combine { APO LIG } \n'
        else:
            complexName = 'APO'

        params += f"charge {complexName} \n {Boxtype} {complexName} " \
                  f"{waterFFDic[self.getEnumText('WaterForceField')]} {self.Distance.get()} iso \n"

        params += f'addIons {complexName} Cl- 0 \n addIons {complexName} Na+ 0 \n'
        params += f'saveAmberParm {complexName} {systemBasename}.top {systemBasename}.crd \n savepdb ' \
                  f'{complexName} {systemBasename}_check.pdb \nquit'

        file = open(self._getExtraPath("leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f extra/leap_commands.txt", cwd=self._getPath())

    def createOutputStep(self):
        systemBasename = self.getSystemName()

        topol_baseName = '{}.top'.format(systemBasename)
        crd_baseName = '{}.crd'.format(systemBasename)
        check_baseName = '{}_check.pdb'.format(systemBasename)


        topol_localPath = abspath(self._getPath(topol_baseName))
        crd_localPath = abspath(self._getPath(crd_baseName))
        check_localPath = abspath(self._getPath(check_baseName))


        amber_system = amberobj.AmberSystem(filename=crd_localPath, topoFile=topol_localPath,
                                           checkFile=check_localPath, ff=self.getEnumText('ProteinForceFieldType'),
                                           wff=self.getEnumText('WaterForceField'))

        self._defineOutputs(outputSystem=amber_system)
        if self.fromInput.get() == LIGAND:
            self._defineSourceRelation(self.inputLigands, amber_system)
        else:
            self._defineSourceRelation(self.inputStructure, amber_system)

    # --------------------------- INFO functions -----------------------------------

    def convertPDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
        runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        return oFile

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            summary.append(
                "This protocol has created a coordinate file, a topology file and a PDB file (for visualization)"
                "with the selected Main force fields and Water Force Field")

        else:
            summary.append("The protocol has not finished.")
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append("This protocol takes a clean pdb file and it uses the "
                           "AMBER software in order to transform the file into an amber format while applying to it "
                           'the force fields for the system and the water molecules.\n To do so, it calls the two main'
                           'preparation programs in AmberTools21: pdb4amber and LEaP. \n'
                           'Finally, the program LEap returns two files which will be necessary for the MD simulation'
                           '(.crd and .prmtop files) and a .pdb file to visualize the structure')

        return methods

    def getProteinFile(self):
        if self.fromInput == LIGAND:
            return abspath(self.inputLigands.get().getProteinFile())
        else:
            return abspath(self.inputStructure.get().getFileName())
        
    def getLigandFile(self):
      if self.fromInput == LIGAND:
          for mol in self.inputLigands.get():
            if mol.__str__() == self.inputLigandSelect.get():
              myMol = mol
              break
          
          molFile = mol.getPoseFile()
          if not molFile:
              molFile = mol.getFileName()
              print(f'Careful, using file {molFile} from input ligand, which might not be docked')
          return abspath(molFile)

    def getSystemName(self):
        return getBaseFileName(self.getProteinFile())