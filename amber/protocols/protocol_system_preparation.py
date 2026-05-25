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
from os.path import abspath, relpath
import os, re, math, shutil

from pwem.protocols import EMProtocol, ProtImportFiles
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC
from pwchem.utils import getBaseName, convertToSdf, runOpenBabel
from pwchem import Plugin as pwchemPlugin

import amber
from amber import Plugin as amberPlugin

import amber.objects as amberobj
from amber.objects import *
from amber.constants import AMBER_DIC

from Bio import PDB

scriptLigPrepName = 'rdkit_addHydrogens.py'

STRUCTURE, LIGAND = 0, 1
LIG_INPUT = f'inputFrom == {LIGAND}'
GAPS_OPTIONS = ['No', 'Gaps termini', 'All termini']

class AmberSystemPrep(EMProtocol):
    """
    This protocol will prepare a system for MD simulation. It will clean the input PDB for further analysis
    and generate topology and coordinate files necessary for MD simulation.
    """

    _label = 'system preparation'
    IMPORT_FROM_FILE = 1
    IMPORT_FROM_SCIPION = 1

    _chargeModel = ['AM1-BCC', 'Mulliken', 'Gasteiger']
    _Status = ['brief', 'default', 'verbose']

    # -------------------------- DEFINE param functions ----------------------

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        # ProtImportFiles.__init__(self, **kwargs)

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
                       label='Add missing atoms: ')
        group.addParam('addCaps', params.EnumParam, choices=GAPS_OPTIONS, default=0,
                       label='Add ACE and NME caps: ',
                       help='Add acetyl (ACE) and N-methylamide (NME) capping groups to protein N-termini and C-termini respectively. '
                            'These caps neutralize terminal charges and are commonly used in MD simulations. '
                            '\n*None*: No caps added. '
                            '\n*Gaps termini*: Add caps only to missing loops (internal gaps in the structure), '
                            'preserving the real N- and C-termini uncapped. '
                            '\n*All termini*: Add caps to both gaps and real protein termini.')


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
        group.addParam('memPosition', params.EnumParam, choices=['Preoriented', 'MEMEMBED', 'PPM'],
                       label='Membrane orientation method', condition='tMem', default=1,
                       help='Select the method to orient the protein within the lipid bilayer:\n\n'
                            '"Preoriented": uses existing PDB coordinates (i.e. from OPM).\n'
                            '"MEMEMBED": performs a geometric search for the best embedding.\n'
                            '"PPM": uses an electrochemical model to calculate the optimal depth and tilt.')

        form.addParam('Status', params.EnumParam, allowsNull=True, default=1,
                      choices=self._Status,
                      label='Choose Status information: ')

        form.addSection('MD prep')
        group = form.addGroup('Force field', help='Force field applied to the system. Force fields are sets of '
                                                  'potential functions and '
                                                  'parametrized interactions that can be used to study physical '
                                                  'systems. You should select as '
                                                  'many force fields as molecules in your system (i.e protein + ligand')
        group.addParam('proteinFF', params.EnumParam,
                       label='Protein Force Field', help='Protein force filed to use',
                       choices=['ff14SB', 'ff19SB', 'ff14SBonlysc', 'ff15ipq', 'fb15', 'ff03.r1', 'ff03ua'],
                       default=0)

        group.addParam('ligandCharge', params.EnumParam, default=2, choices=self._chargeModel,
                       condition=LIG_INPUT, label="Small molecules charge method: ",
                       help='Small molecules charge method to use')
        group.addParam('ligandFF', params.EnumParam, default=1, choices=['gaff', 'gaff2', 'ESPALOMA'],
                       condition=LIG_INPUT, label="Small molecules force field: ",
                       help='Small molecules force field to use')

        group.addParam('lipidFF', params.EnumParam, condition='tMem',
                       label='Type', choices=['lipid21', 'lipid17'], default=0)
        group.addParam('waterForceField', params.EnumParam, default=0,
                       choices=['tip4pew', 'spce', 'spceb', 'opc', 'opc3', 'tip3p'],
                       allowsNull=True,
                       label='Water Force Field',
                       help='Force field applied to the water')
        group = form.addGroup('Disulfide bridges')
        group.addParam('disulfideBridges', params.BooleanParam,
                       label='Are there any S-S bridges?', default=False,
                       help='Residues involved must be renamed to CYX in the pdb file')
        group.addParam('disulfideBridgesNumber', params.StringParam,
                       condition='disulfideBridges',
                       label='Number of the residues involved in the disulfide bridge \n'
                             'with format 1º Residue - 2º Residue / 1º Residue - 2º Residue')

        group = form.addGroup('Solvent box')
        group.addParam('solvateStep', params.EnumParam, default=0, condition='tMem==False',
                       choices=['Cubic', 'Octahedric'], defalult='Cubic',
                       label='Solvation box',
                       help='Both solvation boxes will be isometric')
        group.addParam('minDist', params.FloatParam, condition='tMem==False', label='Padding distance:',
                       default=15.0, help='Minimum distance Å from the protein to the edge of the box')
        group.addParam('memDistXY', params.FloatParam, default=15.0, condition='tMem',
                       label='Min dist to XY boundary (Å):',
                       help='Minimum distance between the protein and the box boundaries in X/Y axes. ')
        group.addParam('memDistZ', params.FloatParam, default=17.5, condition='tMem',
                       label='Water layer width Z (Å):',
                       help='Thickness of the water layer above/below the membrane in Z axis.')

        group.addParam('addIons', params.BooleanParam, label='Add salt to the system: ',
                       help='Add a specific concentration of ions to the system.', default=True)
        line = group.addLine('Salt configuration: ', condition='addIons==True',
                             help='Ions to add to neutralize and to achive a desire salt concentration')
        line.addParam('cationType', params.EnumParam, label='Cation:',
                      default=1, help='Cation to add', choices=['Na+', 'K+'])
        line.addParam('anionType', params.EnumParam, label='Anion:',
                      default=0, help='Anion to add', choices=['Cl-'], defalut=0)
        line.addParam('ionConc', params.FloatParam, label='Concentration (M):',
                      default=0.15, help='Salt concentration of the system')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        molFile = self.getSpecifiedMolFile() if self.inputFrom.get() == LIGAND else None
        if molFile:
            self._insertFunctionStep(self.antechamberStep, molFile)
            self._insertFunctionStep(self.ligLeapStep)
        self._insertFunctionStep(self.prepPdb)
        if self.tMem.get():
            self._insertFunctionStep(self.membraneStep)
        self._insertFunctionStep(self.tleapStep)
        self._insertFunctionStep(self.createOutputStep)

    def antechamberStep(self, molFile):
        molName = os.path.basename(molFile).split(".")[0]
        prepLigFile = f'{molName}_prep.mol2'

        params = ' -i {} -fi sdf -o {} -fo mol2 -rn LIG '.format(molFile, prepLigFile)

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

    def ligLeapStep(self):
        molFile = self.findFile(self.getLigandFileDir(), '.mol2')
        frcmodFile = self.findFile(self.getLigandFileDir(), '.frcmod')
        molName = os.path.basename(frcmodFile.split(".")[0])
        ligFF = self.getEnumText('ligandFF')

        params = 'source leaprc.{} \n' \
                 'LIG = loadmol2 {} \n' \
                 'loadamberparams {} \n' \
                 'saveoff LIG {}.lib \n' \
                 'saveamberparm LIG {}_LIG.prmtop {}_LIG.crd \n' \
                 'savepdb LIG {}.check.pdb \n' \
                 'quit'.format(ligFF, molFile, frcmodFile, *[molName] * 4)

        file = open(os.path.join(self.getLigandFileDir(), "leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f leap_commands.txt", cwd=self.getLigandFileDir())

    def prepPdb(self):
        recPDB = self.getReceptorPDB()
        inputStructure = self.getInputReceptorFilename()
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)
        shutil.copy(inputStructure, recPDB)
        systemBasename = os.path.basename(recPDB.split(".")[0])
        addCapsMode = self.getEnumText('addCaps')
        if addCapsMode in GAPS_OPTIONS[1:]:
            mode = 'gaps' if addCapsMode == GAPS_OPTIONS[1]  else 'all'

            cappedPdb = os.path.abspath(os.path.join(self.getTargetFileDir(), f'{systemBasename}_capped.pdb'))
            pmlScript = self.addCapsPml(recPDB, cappedPdb, mode)

            self.runPymol(pmlScript, self.getTargetFileDir())
            # self.fixPdbTER(cappedPdb)
            recPDB = cappedPdb

        self.insertTERLines(recPDB)

        params = '{} -o {}_amber.pdb --dry'.format(recPDB, systemBasename)

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
        hasLigand = (self.inputFrom.get() == LIGAND)
        hasMembrane = self.tMem

        if hasMembrane:
            inputStructure = self.findFile(self.getTargetFileDir(), '_bilayer_aligned.pdb')
        else:
            inputStructure = self.findFile(self.getTargetFileDir(), '_amber.pdb')

        targetBasename = self.getTleapSystemName()
        if not hasMembrane:
            nCation, nAnion = self.calcIonConc(inputStructure, self.getEnumText('proteinFF'),
                                               self.getEnumText('waterForceField'), self.ionConc.get())
            print(f'{nCation} cations and {nAnion} anion will be added\n')

        cmdsTleap = []
        # load force fields
        cmdsTleap.append(f"source leaprc.protein.{self.getEnumText('proteinFF')}")
        cmdsTleap.append(f"source leaprc.water.{self.getEnumText('waterForceField')}")
        if hasMembrane:
            cmdsTleap.append(f"source leaprc.{self.getEnumText('lipidFF')}")
        if hasLigand:
            cmdsTleap.append(f"source leaprc.{self.getEnumText('ligandFF')}")

        # load components
        components = []
        cmdsTleap.append(f"PROT = loadPdb {inputStructure}")
        components.append('PROT')
        if self.disulfideBridges:
            for pair in self.disulfideBridgesNumber.get().split('/'):
                first = pair.split('-')[0]
                second = pair.split('-')[1]
                cmdsTleap.append(f"bond PROT.{first}.SG PROT.{second}.SG")

        # ligand
        if hasLigand and not hasMembrane:
            cmdsTleap.append(f"loadoff {self.findFile(self.getLigandFileDir(), '.lib')}")
            cmdsTleap.append(f"LIG = loadmol2 {self.findFile(self.getLigandFileDir(), '.mol2')}")
            components.append('LIG')
            cmdsTleap.append(f"loadamberparams {self.findFile(self.getLigandFileDir(), '.frcmod')}")

        # membrane
        if hasMembrane:
            memStructure = self.findFile(self.getTargetFileDir(),   "membrane.pdb")
            cmdsTleap.append(f"MEMB = loadPdb {memStructure}")
            components.append('MEMB')
            if hasLigand:
                cmdsTleap.append(f"loadoff {self.findFile(self.getLigandFileDir(), '.lib')}")
                cmdsTleap.append(f"LIG = loadmol2 {self.findFile(self.getLigandFileDir(), '_aligned.mol2')}")
                components.append('LIG')
                cmdsTleap.append(f"loadamberparams {self.findFile(self.getLigandFileDir(), '.frcmod')}")

        if len(components) == 1:
            cmdsTleap.append('SYSTEM = PROT')
        else:
            cmdsTleap.append(f"SYSTEM = combine {{ {' '.join(components)} }}")

        # solvation non-membrane
        if not hasMembrane:
            boxtype = "SolvateBox" if self.getEnumText("solvateStep") == "Cubic" else "SolvateOct"

            waterBoxes = {
                "tip3p": "TIP3PBOX",
                "tip4pew": "TIP4PEWBOX",
                "spce": "SPCBOX",
                "opc": "OPCBOX",
                "opc3": "OPC3BOX"
            }
            waterModel = self.getEnumText("waterForceField")
            wat = waterBoxes[waterModel]

            cmdsTleap.append("center SYSTEM")
            cmdsTleap.append("charge SYSTEM")
            cmdsTleap.append(f"{boxtype} SYSTEM {wat} {int(self.minDist.get())} iso")
            if self.addIons.get():
                cmdsTleap.append(
                    f"addIonsRand SYSTEM {self.getEnumText('cationType')} {nCation} {self.getEnumText('anionType')} {nAnion}")

        # solvation membrane
        if hasMembrane:
            x, y, z = self._memBox
            cmdsTleap.append(f"set SYSTEM box {{{x:.3f} {y:.3f} {z:.3f}}}")

        cmdsTleap.append(f"savepdb SYSTEM {targetBasename}.pdb")
        cmdsTleap.append(f"saveAmberParm SYSTEM {targetBasename}_woChains.parm7 {targetBasename}.crd")
        cmdsTleap.append(f"savepdb SYSTEM {targetBasename}_system.pdb")

        cmdsTleap.append("quit")
        leapFile = os.path.join(self.getTargetFileDir(), "leap_commands.txt")

        with open(leapFile, "w") as f:
            f.write("\n".join(cmdsTleap))

        amber.Plugin.runAmbertools(self, "tleap", "-f leap_commands.txt", cwd=self.getTargetFileDir())

        # Add chain information with parmed
        cmdsParmed = [
            f"parm {targetBasename}_woChains.parm7",
            f"addPDB {inputStructure}",
            f"outparm {targetBasename}.parm7",
            "quit"
        ]

        parmedFile = os.path.join(self.getTargetFileDir(), "parmed_commands.txt")
        with open(parmedFile, "w") as f:
            f.write("\n".join(cmdsParmed))

        amber.Plugin.runAmbertools(self, "parmed", "-i parmed_commands.txt", cwd=self.getTargetFileDir())

        # 2. Generate the final system PDB with Chain IDs using ambpdb
        # (ambpdb outputs to stdout, so we write it directly via Python's subprocess)

        finalPdbFile = os.path.join(self.getTargetFileDir(), f"{targetBasename}_system_chains.pdb")
        ambpdbCmd = f"-p {targetBasename}.parm7 -c {targetBasename}.crd -ext > {finalPdbFile}"
        amber.Plugin.runAmbertools(self, "ambpdb", ambpdbCmd, cwd=self.getTargetFileDir())

        print(f"System PDB with Chain IDs saved to: {finalPdbFile}")

    def membraneStep(self):
        inputStructure = self.findFile(self.getTargetFileDir(), '_amber.pdb')
        systemBasename = self.getSystemName()
        memOutput = os.path.join(f'{systemBasename}_bilayer.pdb')

        params = f'--lipids {self.memLipids.get()} --ratio {self.memRatio.get()} --dist {self.memDistXY.get()}' \
                 f' --dist_wat {self.memDistZ.get()} --pdb {inputStructure} --notprotonate --nottrim -o {memOutput} --salt' \
                 f' --salt_c {self.getEnumText("cationType")} --salt_a {self.getEnumText("anionType")} --saltcon {self.ionConc.get()}'

        if self.getEnumText('memPosition') == 'Preoriented':
            params += ' --preoriented'
        elif self.getEnumText('memPosition') == 'PPM':
            params += ' --ppm'

        params += ' 2>&1'

        amber.Plugin.runAmbertools(self, 'packmol-memgen ', params, cwd=self.getTargetFileDir())

        logFile = os.path.join(self.getTargetFileDir(), "packmol-memgen.log")
        x, y, z = self.extractBoxFromMemgenLog(logFile)
        print(f"Membrane box extracted: {x:.3f} {y:.3f} {z:.3f}")
        self._memBox = (x, y, z)

        # align the protein to the membrane system
        proteinOutputAligned = os.path.join(self.getTargetFileDir(), f'{systemBasename}_bilayer_aligned.pdb')

        scriptParams = (
            f"-i {inputStructure} "
            f"-m {memOutput} "
            f"-o {proteinOutputAligned}"
        )
        if self.inputLigand:
            molFile = self.findFile(self.getLigandFileDir(), '.mol2')
            if molFile and os.path.isfile(molFile):
                # Add ligand parameters to script
                ligandBasename = os.path.basename(molFile).split('.')[0]
                ligandOutputAligned = os.path.join(self.getLigandFileDir(),
                                                   f'{ligandBasename}_aligned.mol2')
                scriptParams += (f" -l {molFile} --ligand-out {ligandOutputAligned}")
                print(f"Ligand file found and will be aligned: {molFile}")

        amberPlugin.runScript(self, 'alignMembrane.py', args=scriptParams, env=AMBER_DIC,
                              cwd=self.getTargetFileDir())

    def createOutputStep(self):
        systemBasename = self.getTleapSystemName()
        targetDir = self.getTargetFileDir()

        srcTop = os.path.join(targetDir, f'{systemBasename}.parm7')
        srcCrd = os.path.join(targetDir, f'{systemBasename}.crd')
        srcSystemPdb = os.path.join(targetDir, f'{systemBasename}_system_chains.pdb')
        if not os.path.exists(srcSystemPdb):
            srcSystemPdb = os.path.join(targetDir, f'{systemBasename}_system.pdb')

        destTop = relpath(self._getPath(f'{systemBasename}.parm7'))
        destCrd = relpath(self._getPath(f'{systemBasename}.rst7'))
        destSystemPdb = relpath(self._getPath(f'{systemBasename}_system.pdb'))

        shutil.copy(srcTop, destTop)
        shutil.copy(srcCrd, destCrd)
        shutil.copy(srcSystemPdb, destSystemPdb)

        createdSystem = amberobj.AmberSystem(filename=destSystemPdb, crdFile=destCrd, topoFile=destTop,
                                             ff=self.getEnumText('proteinFF'),
                                             wff=self.getEnumText('waterForceField'))

        if self.inputFrom.get() == LIGAND:
            molFile = relpath(self.findFile(self.getLigandFileDir(), '.mol2'))
            createdSystem.setLigTopologyFile(molFile)
            createdSystem.setLigandID('LIG')

        self._defineOutputs(outputSystem=createdSystem)

    # --------------------------- INFO functions -----------------------------------
    def getReceptorPDB(self):
        recPDB = os.path.abspath(self._getExtraPath(f'{self.getSystemName()}.pdb'))
        return recPDB

    def getInputReceptorFilename(self):
        if self.inputFrom.get() == STRUCTURE:
            proteinFile = self.inputStructure.get().getFileName()
        elif self.inputFrom.get() == LIGAND:
            proteinFile = self.inputSetOfMols.get().getProteinFile()
        return os.path.abspath(proteinFile)

    def getSystemName(self):
        return getBaseName(self.getInputReceptorFilename())

    def getTleapSystemName(self):
        """Return the basename used by tleap/parmed for the final system files."""
        if self.inputFrom.get() == LIGAND:
            sdfFile = self.findFile(self.getLigandFileDir(), '.sdf')
            if sdfFile:
                return os.path.splitext(os.path.basename(sdfFile))[0]
        return os.path.basename(self.getReceptorPDB().split(".")[0])

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
            return os.path.join(self.getLigandFileDir(), os.path.basename(sdfFile))

    def getSpecifiedMol(self):
      myMol = None
      for mol in self.inputSetOfMols.get():
        if mol.__str__() == self.inputLigand.get():
          myMol = mol.clone()
          break
      if myMol == None:
        print('The input ligand is not found')
        return None
      else:
        return myMol

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

    def calcIonConc(self, proteinFile, proteinFF, waterFF, molar):
        '''calculate the number of ions to get a desired concentration based on
        the SPLIT method https://doi.org/10.1021/acs.jctc.9b00953'''
        ligandFile = None
        if self.inputFrom.get() == LIGAND:
            ligandFile = self.findFile(self.getLigandFileDir(), '.mol2')
        nWaters, totalQ = self.getNumberWaters(proteinFile, proteinFF, waterFF, ligandFile,
                                               self.getEnumText('ligandFF'))
        # 1. expected number of ions
        nIons = (nWaters * molar) / 56  # The constant 56 is the 'intuitive shortcut' for water molarity

        # 2. Calculate number of ions,  math.ceil to handle the "round up in case of odd Q" requirement
        nCation = int(math.ceil(nIons - (totalQ / 2)))
        nAnion = int(math.ceil(nIons + (totalQ / 2)))

        nCation = max(0, nCation)
        nAnion = max(0, nAnion)

        return nCation, nAnion

    def getNumberWaters(self, proteinFile, proteinFF, waterFF, ligFile=None, ligFF=None):
        '''runs short tleap to calculate nWaters of the box and total charge of the protein or protein+lig'''
        leap_input = os.path.join(self.getTargetFileDir(), 'get_info.in')
        leap_log = os.path.join(self.getTargetFileDir(), "get_info.log")

        commands = []
        commands.append(f"source leaprc.water.{waterFF}")
        commands.append(f"source leaprc.protein.{proteinFF}")

        if ligFF:
            commands.append(f"source leaprc.{ligFF}")

        # Load Solute
        commands.append(f"PROT = loadPdb {proteinFile}")

        if ligFile:
            commands.append(f"LIG = loadmol2 {ligFile}")
            commands.append("COMPLEX = combine {PROT LIG}")
            unit = "COMPLEX"
        else:
            unit = "PROT"

        commands.append(f"solvateBox {unit} TIP3PBOX {self.minDist.get()} iso")

        commands.append(f"charge {unit}")
        commands.append("quit")

        # Write and run tleap
        with open(leap_input, "w") as f:
            f.write("\n".join(commands))

        amber.Plugin.runAmbertools(self, 'tleap', f'-f get_info.in > {leap_log}', cwd=self.getTargetFileDir())

        # Parse the log for Nw and Q
        nWaters = 0
        totalQ = 0

        if os.path.exists(leap_log):
            with open(leap_log, "r") as f:
                log_content = f.read()

                # Find water count
                nw_match = re.search(r"Added\s+(\d+)\s+residues\.", log_content)
                if nw_match:
                    nWaters = int(nw_match.group(1))

                # Find charge
                q_match = re.search(r"Total\s+unperturbed\s+charge:\s+(-?\d+\.?\d*)", log_content)
                if q_match:
                    totalQ = float(q_match.group(1))

        return nWaters, totalQ

    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        if not hasattr(self, 'outputSystem'):
            return ['No output system produced yet.']

        outSystem = self.outputSystem
        summary = [
            f'System file    : {outSystem.getSystemFile()}',
            f'Protein FF     : {self.getEnumText("proteinFF")}',
            f'Water FF       : {self.getEnumText("waterForceField")}',
            f'Topology       : {outSystem.getTopologyFile()}',
            f'Coordinates    : {outSystem.getCrdFile()}',
        ]

        if self.inputFrom.get() == LIGAND:
            summary.extend([
                f'Ligand         : {self.inputLigand.get()}',
                f'Ligand FF      : {self.getEnumText("ligandFF")}',
            ])

        prepFlags = self._getPreparationFlagsSummary()
        if prepFlags:
            summary.append(f'Target prep    : {", ".join(prepFlags)}')

        if self.tMem.get():
            summary.append(f'Membrane       : {self.memLipids.get()}')
        else:
            summary.append(
                f'Solvation box  : {self.getEnumText("solvateStep")}'
            )

        return summary

    def _getPreparationFlagsSummary(self):
        flags = []

        options = [
            (self.targetProteinResidues.get(), 'protein-only residues'),
            (self.targetAmberCompatibleResidues.get(), 'Amber-compatible residues'),
            (self.targetReduce.get(), 'reduce (add H)'),
            (self.targetTleap.get(), 'tleap missing atoms'),
            (self.targetPhSimulation.get(), 'constant-pH renaming'),
        ]

        for enabled, label in options:
            if enabled:
                flags.append(label)

        return flags

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append(
                "This protocol prepares an Amber-compatible molecular dynamics system "
                "starting from a protein structure and, optionally, a ligand structure. "
                "The protocol uses AmberTools utilities such as pdb4amber, reduce, "
                "Antechamber, parmchk2 and tleap to clean and parameterize the system.\n"
                "The receptor structure can be processed to remove unwanted residues, "
                "add hydrogens, rebuild missing atoms and prepare the structure for "
                "constant-pH simulations. Optional ACE/NME capping groups can also be added.\n"
                "If a ligand is provided, the ligand is parameterized with the selected "
                "force field and charge method before being incorporated into the system.\n"
                "The final system is assembled with tleap using the selected protein, "
                "water and optional lipid force fields. The system can be solvated in "
                "a water box or embedded into a membrane environment, and ions can be "
                "added to neutralize the system and reach the desired salt concentration.\n"
                "Finally, Amber topology, coordinate and structure files required for "
                "molecular dynamics simulations are generated."
            )

        return methods

    # --------------------------- UTILS functions -----------------------------------

    def findFile(self, directory, extension):
        if os.path.exists(directory):
            for f in os.listdir(directory):
                if f.endswith(extension):
                    return os.path.join(directory, f)
        return None

    def extractBoxFromMemgenLog(self, logFile):
        """
        Reads packmol-memgen.log and extracts boxsize x_len, y_len, z_len.
        Returns tuple (x, y, z).
        """
        if not os.path.exists(logFile):
            raise FileNotFoundError(f"packmol-memgen log not found: {logFile}")

        with open(logFile, "r") as f:
            text = f.read()

        x_match = re.search(r"x_len\s*=\s*([0-9.+-Ee]+)", text)
        y_match = re.search(r"y_len\s*=\s*([0-9.+-Ee]+)", text)
        z_match = re.search(r"z_len\s*=\s*([0-9.+-Ee]+)", text)

        if not (x_match and y_match and z_match):
            raise ValueError("Could not extract box dimensions from packmol-memgen.log")

        x_len = float(x_match.group(1))
        y_len = float(y_match.group(1))
        z_len = float(z_match.group(1))

        return x_len, y_len, z_len

    def addCapsPml(self, inputPdb, outputPdb, mode='gaps'):
        data = self.identifyTermini(inputPdb)

        pmlLines = [
            "reinitialize",
            f"load {inputPdb}, protein",
            "remove hydro",
            "hide all",
            "show sticks, protein"
        ]

        for gap in data['gaps']:
            # Adds NME on the C-term and ACE on the N-term of gaps
            pmlLines.append(self.removeOXTCommand(gap['chain'], gap['c_term']))
            pmlLines.extend(self.addCapPmlCommand(gap['chain'], gap['c_term'], 'C', 'nme'))
            pmlLines.extend(self.addCapPmlCommand(gap['chain'], gap['n_term'], 'N', 'ace'))

        if mode == 'all':
            for term in data['protein_termini']:
                # Adds NME on the C-term and ACE on the N-term of chain termini
                pmlLines.append(self.removeOXTCommand(term['chain'], term['c_term']))
                pmlLines.extend(self.addCapPmlCommand(term['chain'], term['n_term'], 'N', 'ace'))
                pmlLines.extend(self.addCapPmlCommand(term['chain'], term['c_term'], 'C', 'nme'))

        pmlLines.extend([
            "remove hydro",
            "sort protein",
            f"save {outputPdb}, protein",
            "quit"
        ])

        scriptPath = os.path.join(self.getTargetFileDir(),"capping_script.pml")
        with open(scriptPath, "w") as f:
            f.write("\n".join(pmlLines))

        print(f"PML script for adding caps: {scriptPath} for mode: {mode}")
        return scriptPath

    def runPymol(self, pymolScript, workinDir):
        # run in the background
        self._log.info('Launching: ' + self._getPymol() + pymolScript)
        self.runJob(f'{self._getPymol()} -cq', pymolScript, cwd=workinDir)

    def _getPymol(self):
        return pwchemPlugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')

    def identifyTermini(self, inputPdb):
        """
        Parses a PDB file to identify chain N/C protein termini and
        internal gaps that require capping.
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", inputPdb)

        result = {
            'protein_termini': [],
            'gaps': []
        }

        for chain in structure.get_chains():
            residues = [r for r in chain if PDB.is_aa(r)]
            if not residues:
                continue

            # Extract global N and C termini from the first and last amino acids
            result['protein_termini'].append({
                'chain': chain.id,
                'n_term': residues[0].id[1],
                'c_term': residues[-1].id[1]
            })

            # Detect sequence breaks inline without building intermediate segment arrays
            for resCurr, resNext in zip(residues, residues[1:]):
                currId = resCurr.id[1]
                nextId = resNext.id[1]

                if nextId != currId + 1:
                    result['gaps'].append({
                        'chain': chain.id,
                        'c_term': currId,  # Needs NME
                        'n_term': nextId  # Needs ACE
                    })
        return result

    def addCapPmlCommand(self, chain, resi, atom, capType):
        return [
            f"select tmp_target, /protein//{chain}/{resi}/{atom}",
            "edit tmp_target",
            f"/editor.attach_amino_acid('pk1', '{capType}')"
        ]

    def removeOXTCommand(self, chain, resi):
        return f"remove /protein//{chain}/{resi}/OXT"

    def fixPdbTER(self, pdbPath):
        with open(pdbPath, 'r') as f:
            lines = f.readlines()

        cleanLines = [line for line in lines if not line.startswith("TER")]
        fixedLines = []
        numLines = len(cleanLines)

        for i, line in enumerate(cleanLines):
            fixedLines.append(line)

            # Guard Clause: Skip any line that isn't an NME ATOM record
            if not line.startswith("ATOM") or line[17:20].strip() != "NME":
                continue

            # If it's the very last line, it needs a TER
            if i + 1 >= numLines:
                fixedLines.append("TER\n")
                continue

            # Otherwise, check if the next atom belongs to a different residue
            nextLine = cleanLines[i + 1]
            if nextLine.startswith("ATOM") and nextLine[22:26].strip() != line[22:26].strip():
                fixedLines.append("TER\n")

        with open(pdbPath, 'w') as f:
            f.writelines(fixedLines)

    def getLigandName(self):
        return self.getSpecifiedMol().getMolName()

    def insertTERLines(self, inputPDB, outputPDB=None):
        """Insert TER between chains and gaps, fixing NME caps."""
        outputPDB = outputPDB or inputPDB
        lines = []

        prevChain = None
        prevResNum = 0
        prevResName = ""
        lastAtomSerial = 0

        with open(inputPDB, 'r') as f:
            for line in f:
                # Skip existing TER lines
                if line.startswith('TER'):
                    continue

                # Pass non-atom lines straight through and skip the rest of the loop
                if not line.startswith(('ATOM', 'HETATM')):
                    lines.append(line)
                    continue

                # --- Main ATOM/HETATM parsing ---
                atomSerial = int(line[6:11].strip())
                resName = line[17:20].strip()
                chainID = line[21:22].strip()
                resSeq = int(line[22:26].strip())

                # NME renaming: CH3 → C
                if resName == "NME" and "CH3" in line[12:16]:
                    line = line[:12] + " C  " + line[16:]

                # Detect chain break or gap to insert TER
                if prevChain and (chainID != prevChain or abs(resSeq - prevResNum) > 1):
                    terSerial = lastAtomSerial + 1
                    lines.append(f"TER   {terSerial:5d}      {prevResName:3s} {prevChain}{prevResNum:4d}\n")

                lines.append(line)

                # Update state variables for the next iteration
                prevChain = chainID
                prevResNum = resSeq
                prevResName = resName
                lastAtomSerial = atomSerial

        with open(outputPDB, 'w') as f:
            f.writelines(lines)

        return outputPDB