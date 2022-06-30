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
from os.path import relpath, abspath

import os

import numpy as np

from pwem.protocols import EMProtocol, ProtImportFiles
from pyworkflow.protocol import params
from pyworkflow.utils import Message

import amber
from amber import Plugin as amberPlugin
from pwchem.utils import runOpenBabel

import amber.objects as amberobj
from amber.objects import *


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

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure:", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure to convert to Amber system')
        form.addParam('ligand', params.BooleanParam,
                      label="Do you want to use a ligand?:", allowsNull=True,
                      help='Ligand structure to convert to Amber system, you can extract this structure by using'
                           ' the protocol "amber-Extract_ligand" or input a structure file')
        form.addParam('inputLigands', params.PointerParam,
                      label="Import extracted ligands:", allowsNull=True, condition='ligand == True',
                      important=True, pointerClass='SetOfSmallMolecules')
        form.addParam('inputLigandSelect', params.StringParam,
                      label="Select ligand structure:", allowsNull=True, condition='ligand == True',
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

        group = form.addGroup('ligand modifications options', condition='ligand == True')
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

        group = form.addGroup('Ligand parametrization', condition='ligand == True')
        group.addParam('ChargeModel', params.EnumParam,
                       choices=self._ChargeModel, allowsNull=True,
                       label='Choose the charge model in order to calculate the atomic point charges: ')
        group.addParam('Status', params.EnumParam, allowsNull=True,
                       choices=self._Status,
                       label='Choose status information: ')

        form.addSection('MD prep')
        group = form.addGroup('Force field', help='Force field applied to the system. Force fields are sets of '
                                                  'potential functions and '
                                                  'parametrized interactions that can be used to study physical '
                                                  'systems. You should select as '
                                                  'many force fields as molecules in your system (i.e protein + ligand')

        group.addParam('ProteinForceField', params.BooleanParam, allowsNull=True,
                       label='Protein Force Field')
        group.addParam('ProteinForceFieldType', params.EnumParam,
                       label='Type',
                       choices=['ff14SB', 'ff19SB', 'ff14SBonlysc', 'ff15ipq', 'fb15', 'ff03.r1', 'ff03ua'],
                       condition='ProteinForceField')
        group.addParam('LigandForceField', params.BooleanParam, allowsNull=True, default= False, condition='ligand == True',
                       label='Ligand Force Field', help= 'if you have chosen to introduce a ligand, this force field is mandatory')
        group.addParam('DNAForceField', params.BooleanParam, allowsNull=True,
                       label='DNA Force Field')
        group.addParam('RNAForceField', params.BooleanParam, allowsNull=True,
                       label='RNA Force Field')
        group.addParam('RNAForceFieldType', params.EnumParam, condition='RNAForceField',
                       label='Type', choices=['OL3', 'LJbb', 'YIL', 'ROC', 'Shaw'])
        group.addParam('LipidForceField', params.BooleanParam, allowsNull=True,
                       label='Lipid Force Field')
        group.addParam('WaterForceField', params.EnumParam,
                       choices=['tip4pew', 'spce', 'spceb', 'opc', 'opc3', 'tip3p'],
                       allowsNull=True,
                       label='Water Force Field: ',
                       help='Force field applied to the water')

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
        if self.ligand == True:
            self._insertFunctionStep('PrepStep')
            self._insertFunctionStep('AntechamberStep')
            self._insertFunctionStep('ParmStep')
            self._insertFunctionStep('LeapStep')
        self._insertFunctionStep('PDBAmberStep')
        self._insertFunctionStep('ForceFieldStep')
        self._insertFunctionStep('createOutputStep')

    def PrepStep(self):
        for mol in self.inputLigands.get():
            if mol.getUniqueName() == self.inputLigandSelect:
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

    def AntechamberStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        params = ' -i {}.LIG.pdb -fi pdb -o {}.LIG.mol2 -fo mol2 '.format(*[systemBasename]*2)

        if self.getEnumText('ChargeModel') == 'RESP':
            params += '-c resp '
        if self.getEnumText('ChargeModel') == 'AM1-BCC':
            params += '-c bcc '
        if self.getEnumText('ChargeModel') == 'CM1':
            params += '-c cm1 '
        if self.getEnumText('ChargeModel') == 'CM2':
            params += '-c cm2 '
        if self.getEnumText('ChargeModel') == 'ESP':
            params += '-c esp '
        if self.getEnumText('ChargeModel') == 'Mulliken':
            params += '-c mul '
        if self.getEnumText('ChargeModel') == 'Gasteiger':
            params += '-c gas '

        if self.getEnumText('Status') == 'brief':
            params += '-s 0'
        if self.getEnumText('Status') == 'default':
            params += '-s 1'
        if self.getEnumText('Status') == 'verbose':
            params += '-s 2'

        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getPath())

    def ParmStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        params = '-i {}.LIG.mol2 -o {}.LIG.frcmod -f mol2 '.format(*[systemBasename]*2)

        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self._getPath())

    def LeapStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        params = 'source leaprc.gaff \n' \
                 'loadamberparams {}.LIG.frcmod \n' \
                 'LIG = loadmol2 {}.LIG.mol2 \n' \
                 'saveoff LIG {}.LIG.lib \n' \
                 'saveamberparm LIG {}.LIG.top {}.LIG.crd \n' \
                 'savepdb LIG {}.checkLIG.pdb \n' \
                 'quit'.format(*[systemBasename]*6)

        file = open(self._getExtraPath("leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f extra/leap_commands.txt", cwd=self._getPath())


    def PDBAmberStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)
        systemBasename = os.path.basename(inputStructure.split(".")[0])

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

    def ForceFieldStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

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

        if self.ligand == True:

            params += 'loadamberparams {}.LIG.frcmod \n' \
                      'loadOff {}.LIG.lib \n' \
                      'LIG = loadmol2 {}.LIG.mol2 \n'.format(*[systemBasename]*3)

        if self.DisulfideBridges:

            for pair in self.DisulfideBridgesNumber.get().split('/'):
                first = pair.split('-')[0]
                second = pair.split('-')[1]
                params += 'bond APO.{}.SG APO.{}.SG \n'.format(first, second)

        if self.ligand == True:
            params += 'COMPL = combine { APO LIG } \n'

            if self.getEnumText('SolvateStep') == 'Cubic':
                Boxtype = 'SolvateBox'
            else:
                Boxtype = 'SolvateOct'

            if self.getEnumText('WaterForceField') == 'tip3p':
                params += 'charge COMPL \n {} COMPL TIP3PBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'tip4pew':
                params += 'charge COMPL \n {} COMPL TIP4PEWBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'spece':
                params += 'charge COMPL \n {} COMPL SPCEBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'opc':
                params += 'charge COMPL \n {} COMPL OPCBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'opc3':
                params += 'charge COMPL \n {} COMPL OPCBOX {} iso \n'.format(Boxtype, self.Distance.get())

            params += 'addIons COMPL Cl- 0 \n addIons COMPL Na+ 0 \n'
            params += 'saveAmberParm COMPL {}.top {}.crd \n savepdb COMPL {}_check.pdb \n' \
                      'quit'.format(systemBasename, systemBasename, systemBasename)

            file = open(self._getExtraPath("leap_commands.txt"), "w")
            file.write(params)
            file.close()

        else:
            if self.getEnumText('SolvateStep') == 'Cubic':
                Boxtype = 'SolvateBox'
            else:
                Boxtype = 'SolvateOct'

            if self.getEnumText('WaterForceField') == 'tip3p':
                params += 'charge APO \n {} APO TIP3PBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'tip4pew':
                params += 'charge APO \n {} APO TIP4PEWBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'spece':
                params += 'charge APO \n {} APO SPCEBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'opc':
                params += 'charge APO \n {} APO OPCBOX {} iso \n'.format(Boxtype, self.Distance.get())
            elif self.getEnumText('WaterForceField') == 'opc3':
                params += 'charge APO \n {} APO OPCBOX {} iso \n'.format(Boxtype, self.Distance.get())

            params += 'addIons APO Cl- 0 \n addIons APO Na+ 0 \n'
            params += 'saveAmberParm APO {}.top {}.crd \n savepdb APO {}_check.pdb \n' \
                      'quit'.format(systemBasename, systemBasename, systemBasename)

            file = open(self._getExtraPath("leap_commands.txt"), "w")
            file.write(params)
            file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f extra/leap_commands.txt", cwd=self._getPath())

    def createOutputStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        topol_baseName = '{}.top'.format(systemBasename)
        crd_baseName = '{}.crd'.format(systemBasename)
        check_baseName = '{}_check.pdb'.format(systemBasename)


        topol_localPath = abspath(self._getPath(topol_baseName))
        crd_localPath = abspath(self._getPath(crd_baseName))
        check_localPath = abspath(self._getPath(check_baseName))


        amber_files = amberobj.AmberSystem(filename=crd_localPath, topoFile=topol_localPath,
                                           checkFile=check_localPath, ff=self.getEnumText('ProteinForceFieldType'),
                                           wff=self.getEnumText('WaterForceField'))

        self._defineOutputs(outputSystem=amber_files)
        self._defineSourceRelation(self.inputStructure, amber_files)

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
