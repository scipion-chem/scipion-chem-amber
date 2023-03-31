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
from amber.constants import _cations, _anions

# Water FF from Ambertools21/dat/leap/lib/solvent.lib
waterFFDic = {'spce': 'SPCBOX', 'tip3p': 'TIP3PBOX', 'tip4p': 'TIP4PBOX', 'tip5p': 'TIP5PBOX',
              'opc': 'OPCBOX', 'opc3': 'OPC3BOX'}
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

    _ChargeModel = ['Gasteiger', 'RESP', 'AM1-BCC', 'CM1', 'CM2', 'ESP', 'Mulliken']

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

        group = form.addGroup('Target preparation')
        group.addParam('proteinResidues', params.BooleanParam, default=False,
                       label='Keep only protein residues: ', help='Keep only protein residues')
        group.addParam('AmberCompatibleResidues', params.BooleanParam, default=False,
                       label='Keep only Amber compatible residues: ', 
                       help='Keep only Amber compatible residues. https://ambermd.org/tutorials/basic/tutorial9/index.php')
        group.addParam('phSimulation', params.BooleanParam, default=False,
                       label='Rename GLU, ASP, HIS for constant pH simulation: ')
        group.addParam('reduce', params.BooleanParam, default=False,
                       label='Run reduce first to add hydrogens: ')
        group.addParam('tleap', params.BooleanParam, default=False,
                       label='Use tleap to add missing atoms (EXPERIMENTAL): ')
        group.addParam('DisulfideBridges', params.BooleanParam,
                       label='Are there any S-S bridges?', default=False,
                       help='Residues involved must be renamed to CYX in the pdb file')
        group.addParam('DisulfideBridgesNumber', params.StringParam, condition='DisulfideBridges',
                       label='Number of the residues involved in the disulfide bridge \n'
                             'with format 1º Residue - 2º Residue / 1º Residue - 2º Residue')

        group = form.addGroup('Ligand parametrization', condition='fromInput==1')
        group.addParam('addHydrogens', params.EnumParam, default=2, label='Add hydrogens: ',
                       choices=['Keep present', 'Use reduce (amber)', 'Use rdkit', 'Use obabel'])
        group.addParam('ChargeModel', params.EnumParam,
                       choices=self._ChargeModel, default=0,
                       label='Choose the charge model in order to calculate the atomic point charges: ')

        form.addSection('MD prep')
        group = form.addGroup('Force field', help='Force field applied to the system. Force fields are sets of '
                                                  'potential functions and '
                                                  'parametrized interactions that can be used to study physical '
                                                  'systems. You should select as '
                                                  'many force fields as molecules in your system (i.e protein + ligand')

        group.addParam('ProteinForceFieldType', params.EnumParam,
                       label='Protein Force Field: ', default=0, help='Force field applied to the Protein',
                       choices=['ff14SB', 'ff19SB', 'ff14SBonlysc', 'ff15ipq', 'fb15', 'ff03.r1', 'ff03ua'])
        group.addParam('WaterForceField', params.EnumParam, default=0,
                       choices=list(waterFFDic.keys()),
                       label='Water Force Field: ', help='Force field applied to the water')

        group.addParam('DNAForceField', params.BooleanParam, default=False,
                       label='DNA Force Field', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('RNAForceField', params.BooleanParam, default=False,
                       label='RNA Force Field', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('RNAForceFieldType', params.EnumParam, condition='RNAForceField', default=0,
                       label='Type', choices=['OL3', 'LJbb', 'YIL', 'ROC', 'Shaw'], expertLevel=params.LEVEL_ADVANCED)
        group.addParam('LipidForceField', params.BooleanParam, default=False,
                       label='Lipid Force Field', expertLevel=params.LEVEL_ADVANCED)

        group = form.addGroup('Solvate')
        group.addParam('SolvateStep', params.EnumParam, default=0, label='Solvation box',
                       choices=['Cubic', 'Octahedric'], help='Both solvation boxes will be isometric')
        line = group.addLine('Box size:')
        line.addParam('Distance', params.IntParam, default=20)

        group = form.addGroup('Ions')
        group.addParam('placeIons', params.EnumParam, default=1,
                       label='Add ions: ', choices=['None', 'Neutralize', 'Add number'],
                       help='Whether to add ions to the system.'
                            'https://manual.gromacs.org/documentation/2021.5/onlinehelp/gmx-genion.html')

        line = group.addLine('Cation type:', condition='placeIons!=0',
                             help='Type of the cations to add')
        line.addParam('cationType', params.EnumParam, condition='placeIons!=0',
                      label='Cation to add: ', choices=list(_cations.keys()), default=35,
                      help='Which anion to add in the system')
        line.addParam('cationNum', params.IntParam, condition='placeIons==2',
                      label='Number of cations to add: ',
                      help='Number of cations to add')

        line = group.addLine('Anion type:', condition='placeIons!=0',
                             help='Type of the anions to add')
        line.addParam('anionType', params.EnumParam, condition='placeIons!=0',
                      label='Anions to add: ', choices=list(_anions.keys()), default=1,
                      help='Which anion to add in the system')
        line.addParam('anionNum', params.IntParam, condition='placeIons==2',
                      label='Number of anions to add: ',
                      help='Number of anions to add')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        if self.fromInput == LIGAND:
            self._insertFunctionStep('ligandPrepStep')
            self._insertFunctionStep('leapStep')
        self._insertFunctionStep('pdb4AmberStep')
        self._insertFunctionStep('forceFieldStep')
        self._insertFunctionStep('createOutputStep')

    def ligandPrepStep(self):
        myMolFile = self.getLigandFile()
        systemBasename = self.getSystemName()
        outFile = '{}.LIG.pdb'.format(systemBasename)

        #  todo: problems to reduce and calculate charges, maybe using
        #   https://mmb.irbbarcelona.org/biobb/workflows/tutorials/amber_md_setup_lig
        # Convert into readable
        if myMolFile.endswith('.pdbqt'):
          myMolFile = self.convertPDB(myMolFile)

        if self.addHydrogens.get() == 0:
            outFile = myMolFile
        elif self.addHydrogens.get() == 1:
            args = '{} > {}'.format(myMolFile, outFile)
            amber.Plugin.runAmbertools(self, 'reduce', args, cwd=self._getPath())
        elif self.addHydrogens.get() == 2:
            outFile = self.rdkitAddHydrogens(myMolFile)
        else:
            args = ' -i{} {} -h -opdb -O {}'.format(os.path.splitext(myMolFile)[1][1:],
                                                     os.path.abspath(myMolFile), outFile)
            pwchemPlugin.runOPENBABEL(protocol=self, args=args, cwd=self._getPath(), popen=False)

        # Antechamber
        molFile = os.path.abspath(self._getPath('{}.LIG.mol2'.format(systemBasename)))
        params = ' -i {} -fi {} -o {} -fo mol2'.format(outFile, os.path.splitext(outFile)[1][1:], molFile)
        params += f" -c {chargeDic[self.getEnumText('ChargeModel')]} "
        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getPath())
        molFile = self.removeDuplicatedBondsMol2(molFile)

        # parmchk2: check if all parameters needed are available
        paramFile = '{}.LIG.frcmod'.format(systemBasename)
        params = '-i {} -o {} -f mol2 '.format(molFile, paramFile)
        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self._getPath())

        with open(self._getPath(paramFile)) as f:
            if 'ATTN: needs revision' in f.read():
                print('Antechamber was not able to parametrize the ligand, manual parametrization must be done')


    def leapStep(self):
        systemBasename = self.getSystemName()

        params = 'source leaprc.gaff \n' \
                 'loadamberparams {}.LIG.frcmod \n' \
                 'LIG = loadmol2 {}.LIG.mol2 \n' \
                 'saveoff LIG {}.LIG.lib \n' \
                 'saveamberparm LIG {}.LIG.prmtop {}.LIG.rst7 \n' \
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
        params += 'source leaprc.protein.{} \n'.format(self.getEnumText('ProteinForceFieldType'))

        if self.DNAForceField:
            params += 'source leaprc.DNA.OL15 \n'
        if self.RNAForceField:
            params += 'source leaprc.RNA.{} \n'.format(self.getEnumText('RNAForceFieldType'))
        if self.LipidForceField:
            params += 'source leaprc.lipid17 \n'

        params += 'source leaprc.water.{} \n'.format(self.getEnumText('WaterForceField'))

        params += 'APO = loadPdb {}.amber.pdb \n'.format(systemBasename)

        if self.fromInput == LIGAND:

            params += 'source leaprc.gaff2 \n' \
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

        params += f"charge {complexName} \n"

        if self.placeIons.get() != 0:
            if self.placeIons.get() == 1:
                numCat, numAni = 0, 0
            elif self.placeIons.get() == 2:
                numCat, numAni = self.cationNum.get(), self.anionNum.get()

            params += f'addIons {complexName} {_anions[self.getEnumText("anionType")]} {numCat} \n' \
                      f'addIons {complexName} {_cations[self.getEnumText("cationType")]} {numAni} \n'

        wFF = waterFFDic[self.getEnumText('WaterForceField')]
        params += f"{Boxtype} {complexName} {wFF} {self.Distance.get()} iso \n"
        params += f'saveAmberParm {complexName} {systemBasename}.prmtop {systemBasename}.rst7 \n savepdb ' \
                  f'{complexName} {systemBasename}_check.pdb \nquit'

        file = open(self._getExtraPath("leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f extra/leap_commands.txt", cwd=self._getPath())

    def createOutputStep(self):
        systemBasename = self.getSystemName()

        topol_baseName = '{}.prmtop'.format(systemBasename)
        rst_baseName = '{}.rst7'.format(systemBasename)
        check_baseName = '{}_check.pdb'.format(systemBasename)


        topol_localPath = abspath(self._getPath(topol_baseName))
        rst_localPath = abspath(self._getPath(rst_baseName))
        check_localPath = abspath(self._getPath(check_baseName))


        amber_system = amberobj.AmberSystem(filename=rst_localPath, topoFile=topol_localPath,
                                           checkFile=check_localPath, ff=self.getEnumText('ProteinForceFieldType'),
                                           wff=self.getEnumText('WaterForceField'))
        amber_system.setResNames(parse=True)

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

    def writeParamsFile(self, paramsFile, molFile):
        with open(paramsFile, 'w') as f:
            f.write('ligandFiles: {}\n'.format(molFile))

            f.write('outputDir: {}\n'.format(os.path.abspath(self._getExtraPath())))
            f.write('doHydrogens: True\n')
            f.write('doGasteiger: False\n')
            f.write('ffMethod: MMFF94\n')
            f.write('outFormat: pdb\n')
        return paramsFile

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
                           '(.rst7 and .prmtop files) and a .pdb file to visualize the structure')

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
          
          molFile = myMol.getPoseFile()
          if not molFile:
              molFile = myMol.getFileName()
              print(f'Careful, using file {molFile} from input ligand, which might not be docked')
          return abspath(molFile)

    def getSystemName(self):
        return getBaseFileName(self.getProteinFile())

    def removeDuplicatedBondsMol2(self, inFile, outFile=None):
        with open(inFile) as f:
            inStr = f.read()

        inConnect, outStr, bondedAtoms = False, '', []
        nBonds = 0
        for i, line in enumerate(inStr.split('\n')):
            includeLine = True
            if i == 2:
                infoLine = line
            if '@<TRIPOS>BOND' in line:
                inConnect = True
            elif inConnect and line.startswith('@'):
                inConnect = False

            elif inConnect:
                atoms = line.split()[1:3]
                if atoms in bondedAtoms:
                    includeLine = False
                else:
                    nBonds += 1
                    line = '{:>6}{:>6}{:>6}{:>2}   '.format(nBonds, atoms[0], atoms[1], line.split()[-1])
                    bondedAtoms.append(atoms)

            if includeLine:
                outStr += line + '\n'

        infoElems = infoLine.split()
        infoElems[1] = nBonds
        newInfoLine = '{:>5}{:>6}{:>6}{:>6}{:>6}'.format(*infoElems)

        if not outFile:
            outFile = inFile
        with open(outFile, 'w') as f:
            f.write(outStr.replace(infoLine, newInfoLine))

        return outFile


    def rdkitAddHydrogens(self, inFile):
      paramsFile = os.path.abspath(self._getTmpPath('addHydrogens.txt'))
      with open(paramsFile, 'w') as f:
          f.write('ligandFiles: {}\n'.format(inFile))

          f.write('outputDir: {}\n'.format(os.path.abspath(self._getExtraPath())))
          f.write('doHydrogens: True\n')
          f.write('doGasteiger: True\n')
          f.write('ffMethod: MMFF94\n')
          f.write('outFormat: pdb\n')

      pwchemPlugin.runScript(self, 'ligand_preparation_script.py', paramsFile, env='rdkit', cwd=self._getExtraPath())
      molBase = getBaseFileName(inFile)
      for file in os.listdir(self._getExtraPath()):
          if molBase in file:
              outFile = os.path.abspath(self._getExtraPath(file))
      return outFile