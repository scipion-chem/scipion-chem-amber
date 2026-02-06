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
from email.policy import default
from os.path import relpath, abspath

import os
import re
import numpy as np

from pwem.protocols import EMProtocol, ProtImportFiles
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC
from pwchem.utils import getBaseName, convertToSdf
from pwchem import Plugin as pwchemPlugin

import amber
from amber import Plugin as amberPlugin
from pwchem.utils import runOpenBabel

import amber.objects as amberobj
from amber.objects import *
from amber.constants import AMBER_DIC

scriptLigPrepName = 'rdkit_addHydrogens.py'

STRUCTURE, LIGAND = 0, 1
LIG_INPUT = f'inputFrom == {LIGAND}'

class AmberSystemPrep(EMProtocol):
    """
    This protocol will prepare a system for MD simulation. It will clean the input PDB for further analysis
    and generate topology and coordinate files necessary for MD simulation.
    """

    _label = 'system preparation'
    IMPORT_FROM_FILE = 1
    IMPORT_FROM_SCIPION = 1

    _ChargeModel = ['RESP', 'AM1-BCC', 'CM1', 'CM2', 'ESP', 'Mulliken', 'Gasteiger']
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
        iGroup = form.addGroup('Input')

        iGroup.addParam('inputFrom', params.EnumParam, default=STRUCTURE,
                        label='Input from: ', choices=['AtomStruct', 'SetOfSmallMolecules'],
                        help='Type of input you want to use')
        iGroup.addParam('inputStructure', params.PointerParam, pointerClass='AtomStruct',
                        label='Input structure to be prepared for MD:', condition='inputFrom==0', allowsNull=True,
                        help='Atomic structure to be prepared for MD by solvation, ions addition etc')
        iGroup.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input set of molecules:', condition=LIG_INPUT, allowsNull=True,
                        help='Input set of docked molecules. One of them will be prepared together with its target')
        iGroup.addParam('inputLigand', params.StringParam, condition=LIG_INPUT,
                        label='Ligand to prepare: ',
                        help='Specific ligand to prepare in the system')

        group = form.addGroup('Target modification options')
        group.addParam('targetProteinResidues', params.BooleanParam, default=False,
                       label='Keep only protein residues: ')
        group.addParam('targetAmberCompatibleResidues', params.BooleanParam, default=False,
                       label='Keep only Amber compatible residues: ')
        group.addParam('targetPhSimulation', params.BooleanParam, default=False,
                       label='Rename GLU, ASP, HIS for constant pH simulation: ')
        group.addParam('targetReduce', params.BooleanParam, default=True,
                       label='Run reduce first to add hydrogens: ')
        group.addParam('targetTleap', params.BooleanParam, default=False,
                       label='Use tleap to add missing atoms (EXPERIMENTAL): ')

        form.addParam('tMem', params.BooleanParam, default=False,
                      label='Model a transmembrane protein',
                      help='Embed the protein in a lipid bilayer')

        # Membrane params
        group = form.addGroup('Membrane options', condition='tMem')
        group.addParam('memLipids', params.StringParam, default='POPC',
                       label='Lipid composition:',
                       help='Lipids to embed the protein, use : for separating different lipids (e.g., "POPC:CHL1"). '
                            'Use "//" for different leaflets (e.g., "POPC//POPE"). '
                            'To see all available lipids run packmol-memgen --available_lipids')
        group.addParam('memRatio', params.StringParam, default='1',
                       label='Lipid ratio:',
                       help='Molar ratio matching the lipid string (e.g., "1:1"). '
                            'Set to 1 if single lipid.')

        #
        # group = form.addGroup('ligand modifications options', condition='ligand == True')
        # group.addParam('proteinResidues', params.BooleanParam, default=False,
        #                label='Keep only protein residues: ')
        # group.addParam('AmberCompatibleResidues', params.BooleanParam, default=False,
        #                label='Keep only Amber compatible residues: ')
        # group.addParam('phSimulation', params.BooleanParam, default=False,
        #                label='Rename GLU, ASP, HIS for constant pH simulation: ')
        # group.addParam('reduce', params.BooleanParam, default=False,
        #                label='Run reduce first to add hydrogens: ', help='The addition of hydrogen helps to find the '
        #                                                                  'hydrogen bond interactions and more favorable '
        #                                                                  'to us to find binding affinity of ligand '
        #                                                                  'against protein.')
        # group.addParam('tleap', params.BooleanParam, default=False,
        #                label='Use tleap to add missing atoms (EXPERIMENTAL): ')
        #
        # group = form.addGroup('Ligand parametrization', condition='ligand == True')
        # group.addParam('ChargeModel', params.EnumParam,
        #                choices=self._ChargeModel, allowsNull=True,
        #                label='Choose the charge model in order to calculate the atomic point charges: ')
        form.addParam('Status', params.EnumParam, allowsNull=True, default=1,
                       choices=self._Status,
                       label='Choose status information: ')


        form.addSection('MD prep')
        group = form.addGroup('Force field', help='Force field applied to the system. Force fields are sets of '
                                                  'potential functions and '
                                                  'parametrized interactions that can be used to study physical '
                                                  'systems. You should select as '
                                                  'many force fields as molecules in your system (i.e protein + ligand')

        group.addParam('ProteinForceField', params.BooleanParam, allowsNull=True, default=True,
                       label='Protein Force Field')
        group.addParam('ProteinFF', params.EnumParam,
                       label='Type',
                       choices=['ff14SB', 'ff19SB', 'ff14SBonlysc', 'ff15ipq', 'fb15', 'ff03.r1', 'ff03ua'],
                       default=0)
        group.addParam('protNetCharge', params.IntParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      label='Target net charge: ',
                      help="Enter the integer net charge of the molecule. \n")
        group.addParam('ligandCharge', params.EnumParam, default=2, choices=['AM1-BCC', 'Mulliken', 'Gasteiger'],
                       condition=LIG_INPUT, label="Small molecules charge method: ",
                       help='Small molecules charge method to use')
        group.addParam('ligandFF', params.EnumParam, default=1, choices=['gaff', 'gaff2', 'ESPALOMA'],
                      condition=LIG_INPUT, label="Small molecules force field: ",
                      help='Small molecules force field to use')

        group.addParam('ligNetCharge', params.IntParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                      label='Ligand net charge: ', condition=LIG_INPUT,
                      help="Enter the integer net charge of the molecule. \n"
                           "If antechamber reports an 'odd number of electrons', your charge is likely "
                           "mismatched with your structure's protonation state.")
        group.addParam('lipidFF', params.EnumParam, condition='tMem',
                       label='Type',choices=['lipid21','lipid17'], default=0)
        group.addParam('WaterForceField', params.EnumParam, default=1,
                       choices=['tip4pew', 'spce', 'spceb', 'opc', 'opc3', 'tip3p'],
                       allowsNull=True,
                       label='Water Force Field',
                       help='Force field applied to the water')
        group = form.addGroup('Disulfide bridges')
        group.addParam('DisulfideBridges', params.BooleanParam,
                       label='Are there any S-S bridges?', default=False,
                       help='Residues involved must be renamed to CYX in the pdb file')
        group.addParam('DisulfideBridgesNumber', params.StringParam,
                       condition='DisulfideBridges',
                       label='Number of the residues involved in the disulfide bridge \n'
                             'with format 1º Residue - 2º Residue / 1º Residue - 2º Residue')

        group = form.addGroup('Solvent box')
        group.addParam('SolvateStep', params.EnumParam, default=0, condition='tMem==False',
                       choices=['Cubic', 'Octahedric'], defalult='Cubic',
                       label='Solvation box',
                       help='Both solvation boxes will be isometric')
        group.addParam('minDist', params.FloatParam, condition='tMem==False', label='Padding distance:',
                      default=20.0, help='Minimum distance Å from the protein to the edge of the box')
        group.addParam('memDistXY', params.FloatParam, default=15.0, condition='tMem',
                       label='Min dist to XY boundary (Å):',
                       help='Minimum distance between the protein and the box boundaries in X/Y axes. ')
        group.addParam('memDistZ', params.FloatParam, default=17.5, condition='tMem',
                       label='Water layer width Z (Å):',
                       help='Thickness of the water layer above/below the membrane in Z axis.')

        ## AÑADIR CATION ANION Y MOLARIDAD
    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        recFile = self.getReceptorPDB()
        molFile = self.getSpecifiedMolFile() if self.inputFrom.get() == LIGAND else None

        print(recFile, molFile)

        if molFile:
            self._insertFunctionStep('AntechamberStep', molFile)
            self._insertFunctionStep('leapStep')
        self._insertFunctionStep('PDBAmberStep')
        if self.tMem.get():
            self._insertFunctionStep('membraneStep')
        self._insertFunctionStep('tleapStep')
        self._insertFunctionStep('createOutputStep')

    def PrepStep(self):
        for mol in self.inputLigands.get():
            if mol == self.inputLigandSelect:
                myMol = mol
                break
        myMolFile = os.path.abspath(myMol.getFileName())

        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

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

    def AntechamberStep(self, molFile):
        # inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        # systemBasename = os.path.basename(inputStructure.split(".")[0])
        molName = os.path.basename(molFile).split(".")[0]
        prepLigFile = f'{molName}_prep.mol2'

        nc = self.ligNetCharge.get()

        # params = ' -i {}.LIG.pdb -fi pdb -o {}.LIG.mol2 -fo mol2 '.format(*[systemBasename]*2)
        params = ' -i {} -fi sdf -o {} -fo mol2 -nc {} -rn LIG '.format(molFile, prepLigFile, nc)

        if self.getEnumText('ligandCharge') == 'RESP':
            params += '-c resp '
        if self.getEnumText('ligandCharge') == 'AM1-BCC':
            params += '-c bcc '
        if self.getEnumText('ligandCharge') == 'CM1':
            params += '-c cm1 '
        if self.getEnumText('ligandCharge') == 'CM2':
            params += '-c cm2 '
        if self.getEnumText('ligandCharge') == 'ESP':
            params += '-c esp '
        if self.getEnumText('ligandCharge') == 'Mulliken':
            params += '-c mul '
        if self.getEnumText('ligandCharge') == 'Gasteiger':
            params += '-c gas '

        if self.getEnumText('Status') == 'brief':
            params += '-s 0'
        if self.getEnumText('Status') == 'default':
            params += '-s 1'
        if self.getEnumText('Status') == 'verbose':
            params += '-s 2'

        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self.getLigandFileDir())

        params = '-i {} -o {}.frcmod -f mol2 '.format(prepLigFile, molName)

        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self.getLigandFileDir())


    def leapStep(self):
        # inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        # systemBasename = os.path.basename(inputStructure.split(".")[0])
        molFile = self.findFile(self.getLigandFileDir(),'.mol2')
        frcmodFile = self.findFile(self.getLigandFileDir(),'.frcmod')
        molName = os.path.basename(frcmodFile.split(".")[0])
        ligFF = self.getEnumText('ligandFF')

        params = 'source leaprc.{} \n' \
                 'LIG = loadmol2 {} \n' \
                 'loadamberparams {} \n' \
                 'saveoff LIG {}.lib \n' \
                 'saveamberparm LIG {}_LIG.prmtop {}_LIG.crd \n' \
                 'savepdb LIG {}.check.pdb \n' \
                 'quit'.format(ligFF, molFile, frcmodFile, *[molName]*4)

        file = open(os.path.join(self.getLigandFileDir(),"leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f leap_commands.txt", cwd=self.getLigandFileDir())

    def PDBAmberStep(self):
        # molFile = self.getSpecifiedMolFile() if self.inputFrom.get() == LIGAND else None
        inputStructure = self.getReceptorPDB()
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        params = '{} -o {}_amber.pdb --dry'.format(inputStructure, systemBasename)

        if self.targetProteinResidues:
            params += ' -p '
        if self.targetAmberCompatibleResidues.get():
            params += ' -a '
        if self.targetPhSimulation.get():
            params += ' --constantph '
        if self.targetReduce.get():
            params += ' --reduce '
        if self.targetTleap.get():
            params += ' --add-missing-atoms '

        amber.Plugin.runAmbertools(self, 'pdb4amber', params, cwd=self.getTargetFileDir())

    def tleapStep(self):
        inputStructure = self.findFile(self.getTargetFileDir(), '_amber.pdb')
        systemBasename = os.path.basename(inputStructure.split(".")[0])
        hasLigand = (self.inputFrom.get() == LIGAND)
        hasMembrane = self.tMem
        if self.inputFrom.get() == LIGAND:
            targetBasename = os.path.basename(self.findFile(self.getLigandFileDir(), '.sdf').split(".")[0])
        else:
            targetBasename = os.path.basename(self.getReceptorPDB().split(".")[0])

        cmdsTleap = []
        # load force fields
        cmdsTleap.append(f"source leaprc.protein.{self.getEnumText('ProteinFF')}")
        cmdsTleap.append(f"source leaprc.water.{self.getEnumText('WaterForceField')}")
        if hasMembrane:
            cmdsTleap.append(f"source leaprc.{self.getEnumText('lipidFF')}")
        if hasLigand:
            cmdsTleap.append(f"source leaprc.{self.getEnumText('ligandFF')}")

        # load components
        components = []
        cmdsTleap.append(f"PROT = loadPdb {inputStructure}")
        components.append('PROT')
        if self.DisulfideBridges:
            for pair in self.DisulfideBridgesNumber.get().split('/'):
                first = pair.split('-')[0]
                second = pair.split('-')[1]
                cmdsTleap.append(f"bond PROT.{first}.SG PROT.{second}.SG")

        # ligand
        if hasLigand:
            cmdsTleap.append(f"loadoff {self.findFile(self.getLigandFileDir(), '.lib')}")
            cmdsTleap.append(f"LIG = loadmol2 {self.findFile(self.getLigandFileDir(), '.mol2')}")
            components.append('LIG')
            cmdsTleap.append(f"loadamberparams {self.findFile(self.getLigandFileDir(), '.frcmod')}")

        # membrane
        if hasMembrane:
            memStructure = self.findFile(self.getTargetFileDir(), "_bilayer_aligned.pdb")
            cmdsTleap.append(f"MEMB = loadPdb {memStructure}")
            components.append('MEMB')

        if len(components) == 1:
            cmdsTleap.append('SYSTEM = PROT')
        else:
            cmdsTleap.append(f"SYSTEM = combine {{ {' '.join(components)} }}")

        # solvation non-membrane
        if not hasMembrane:
            boxtype = "SolvateBox" if self.getEnumText("SolvateStep") == "Cubic" else "SolvateOct"

            waterBoxes = {
                "tip3p": "TIP3PBOX",
                "tip4pew": "TIP4PEWBOX",
                "spece": "SPCEBOX",
                "opc": "OPCBOX",
                "opc3": "OPC3BOX"
            }
            waterModel = self.getEnumText("WaterForceField")
            wat = waterBoxes[waterModel]

            cmdsTleap.append("charge SYSTEM")
            cmdsTleap.append(f"{boxtype} SYSTEM {wat} {self.minDist.get()} iso")
            cmdsTleap.append("addIons SYSTEM Cl- 0")
            cmdsTleap.append("addIons SYSTEM Na+ 0")

        # solvation membrane
        if hasMembrane:
            x, y, z = self._memBox
            cmdsTleap.append(f"set SYSTEM box {{{x:.3f} {y:.3f} {z:.3f}}}")

        cmdsTleap.append(f"savepdb SYSTEM {targetBasename}.pdb")
        cmdsTleap.append(f"saveAmberParm SYSTEM {targetBasename}.prmtop {targetBasename}.crd")
        cmdsTleap.append(f"savepdb SYSTEM {targetBasename}_system.pdb")

        cmdsTleap.append("quit")

        # Write leap script + run
        leapFile = os.path.join(self.getTargetFileDir(), "leap_commands.txt")

        with open(leapFile, "w") as f:
            f.write("\n".join(cmdsTleap))

        amber.Plugin.runAmbertools(self, "tleap", "-f leap_commands.txt", cwd=self.getTargetFileDir())

    def membraneStep(self):
        inputStructure = self.findFile(self.getTargetFileDir(), '_amber.pdb')
        systemBasename = self.getSystemName()
        memOutput = os.path.join(f'{systemBasename}_bilayer.pdb')

        params = f'--lipids {self.memLipids.get()} --ratio {self.memRatio.get()} --dist {self.memDistXY.get()}' \
                 f' --dist_wat {self.memDistZ.get()} --pdb {inputStructure} --notprotonate --nottrim -o {memOutput}' \
                 f' --saltcon 0.15'

        amber.Plugin.runAmbertools(self, 'packmol-memgen ', params, cwd=self.getTargetFileDir())

        logFile = os.path.join(self.getTargetFileDir(), "packmol-memgen.log")
        x, y, z = self.extractBoxFromMemgenLog(logFile)
        print(f"Membrane box extracted: {x:.3f} {y:.3f} {z:.3f}")
        self._memBox = (x, y, z)

        # align the membrane to the input pdb to keep the coordinates
        memOutputAligned = (os.path.join(self.getTargetFileDir(), f'{systemBasename}_bilayer_aligned.pdb'))
        scriptParams = (
            f"-i {inputStructure} "
            f"-m {memOutput} "
            f"-o {memOutputAligned}"
        )
        amberPlugin.runScript(self, 'alignMembrane.py', args=scriptParams, env=AMBER_DIC,
                              cwd=self.getTargetFileDir())
        # amber.Plugin.runAmbertools(self, '{} -cq -d'.format(self.getPymolBin()), f'"{pymolCmd}"', cwd=self.getTargetFileDir())

    def createOutputStep(self):
        systemBasename = self.getSystemName()
        targetDir = self.getTargetFileDir()

        srcTop = self.findFile(targetDir, '.prmtop')
        srcCrd = self.findFile(targetDir, '.crd')
        srcSystemPdb = self.findFile(targetDir, '_system.pdb')

        destTop = abspath(self._getPath(f'{systemBasename}.parm7'))
        destCrd = abspath(self._getPath(f'{systemBasename}.rst7'))
        destSystemPdb = abspath(self._getPath(f'{systemBasename}_system.pdb'))

        shutil.copy(srcTop, destTop)
        shutil.copy(srcCrd, destCrd)
        shutil.copy(srcSystemPdb, destSystemPdb)

        createdSystem = amberobj.AmberSystem(filename=destSystemPdb, crdFile=destCrd, topoFile=destTop,
                                           ff=self.getEnumText('ProteinFF'),
                                           wff=self.getEnumText('WaterForceField'))

        if self.inputFrom.get() == LIGAND:
            molFile = self.findFile(self.getLigandFileDir(), '.mol2')
            createdSystem.setLigTopologyFile(molFile)

        self._defineOutputs(outputSystem=createdSystem)
        # self._defineSourceRelation(self.inputLigand, createdSystem)

    # --------------------------- INFO functions -----------------------------------
    def getReceptorPDB(self):
        recPDB = os.path.abspath(self._getExtraPath(f'{self.getSystemName()}.pdb'))
        if not os.path.exists(recPDB):
            recFile = self.getReceptorFilename()
            args = f'{recFile} --output {recPDB}'
            pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())
        return recPDB

    def getReceptorFilename(self):
      if self.inputFrom.get() == STRUCTURE:
          proteinFile = self.inputStructure.get().getFileName()
      elif self.inputFrom.get() == LIGAND:
          proteinFile = self.inputSetOfMols.get().getProteinFile()
      return os.path.abspath(proteinFile)

    def getSystemName(self):
      return getBaseName(self.getReceptorFilename())

    def getSpecifiedMolFile(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
          if mol.__str__() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            molFile = myMol.getPoseFile()
            sdfFile = convertToSdf(self, molFile)
            paramFile = self.writePrepParamsFile([sdfFile])
            pwchemPlugin.runScript(self, scriptLigPrepName, paramFile, env=RDKIT_DIC, cwd=self._getPath())
            return os.path.join(self.getLigandFileDir(), os.listdir(self.getLigandFileDir())[0])

    def writePrepParamsFile(self, molFiles):
        paramsFile = self.getLigParamFile()
        with open(paramsFile, 'w') as f:
            molFilesStr = ' '.join(molFiles)
            f.write(f"ligandFiles: {molFilesStr}\n")

            f.write(f'outputDir: {self.getLigandFileDir()}\n')
            f.write('doHydrogens: True\n')
            f.write('doGasteiger: False\n')
            f.write('sanitize: False\n')
        return paramsFile

    def getLigParamFile(self):
      return os.path.abspath(self._getExtraPath('addHydrogens.txt'))

    def getLigandFileDir(self):
      lDir = os.path.abspath(self._getExtraPath('ligand'))
      if not os.path.exists(lDir):
        os.mkdir(lDir)
      return lDir

    def getTargetFileDir(self):
      tDir = os.path.abspath(self._getExtraPath('target'))
      if not os.path.exists(tDir):
        os.mkdir(tDir)
      return tDir


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

    def findFile(self, directory, extension):
        if os.path.exists(directory):
            for f in os.listdir(directory):
                if f.endswith(extension):
                    return os.path.join(directory, f)
        return None

    def getSystemName(self):
      return getBaseName(self.getReceptorFilename())

    def extractBoxFromMemgenLog(self, logFile):
        """
        Reads packmol-memgen.log and extracts boxsize x_len, y_len, z_len.
        Returns tuple (x, y, z).
        """
        if not os.path.exists(logFile):
            raise FileNotFoundError(f"packmol-memgen log not found: {logFile}")

        with open(logFile, "r") as f:
            text = f.read()

        # Regex patterns
        x_match = re.search(r"x_len\s*=\s*([0-9.+-Ee]+)", text)
        y_match = re.search(r"y_len\s*=\s*([0-9.+-Ee]+)", text)
        z_match = re.search(r"z_len\s*=\s*([0-9.+-Ee]+)", text)

        if not (x_match and y_match and z_match):
            raise ValueError("Could not extract box dimensions from packmol-memgen.log")

        x_len = float(x_match.group(1))
        y_len = float(y_match.group(1))
        z_len = float(z_match.group(1))

        return x_len, y_len, z_len

    def getPymolBin(self):
        return pwchemPlugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')