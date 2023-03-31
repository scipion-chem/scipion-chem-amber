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
This module will prepare the ligand for the simulation.
"""
from os.path import relpath, abspath

import os, shutil

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.utils import Message

import amber
from pwchem.utils import *
from pwchem import Plugin


import amber.objects as amberobj
from amber import Plugin as amberPlugin

anteChDic = {'RESP': 'resp', 'AM1-BCC': 'bcc', 'CM1': 'cm1', 'CM2': 'cm2', 'ESP': 'esp',
             'Mulliken': 'mul', 'Gasteiger': 'gas'}

statusDic = {'brief': 0, 'default': 1, 'verbose': 2}

class AmberLigandPrep(EMProtocol):
    """
With this protocol you will obtain coordinate and topology files from your ligand
using the pdb4amber and Antechamber programs from AMBERTOOLS
    """


    _label = 'Ligand preparation'

    _ChargeModel = ['RESP', 'AM1-BCC', 'CM1', 'CM2', 'ESP', 'Mulliken', 'Gasteiger',]
    _Status = ['brief', 'default', 'verbose']

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)


    # -------------------------- DEFINE param functions ----------------------


    def _defineParams(self, form):
        """
    Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Input set of molecules:', allowsNull=False,
                      help='Input set of docked molecules. One of them will be prepared together with its target')
        form.addParam('inputLigand', params.StringParam,
                      label='Ligand to prepare: ',
                      help='Specific ligand to prepare in the system')

        group = form.addGroup('pdb4amber options')
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

        group = form.addGroup('Ligand parametrization')
        group.addParam('ChargeModel', params.EnumParam,
                       choices=self._ChargeModel,
                       label='Choose the charge model in order to calculate the atomic point charges: ')
        group.addParam('Status', params.EnumParam,
                       choices=self._Status,
                       label='Choose status information: ')


    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('prepStep')
        self._insertFunctionStep('antechamberStep')
        self._insertFunctionStep('parmStep')
        self._insertFunctionStep('leapStep')
        self._insertFunctionStep('createOutputStep')


    def prepStep(self):
        ligandFile = self.getInputFile()
        inputStructure = self.getConvFile(getBaseFileName(ligandFile))
        if not ligandFile.endswith('.pdb'):
            inputStructure = self.convertPDB(ligandFile)
        else:
            shutil.copy(ligandFile, inputStructure)

        systemBasename = os.path.basename(inputStructure.split(".")[0])
        params = '{} > {}.pdb --no-conect '.format(inputStructure, systemBasename)

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

        amber.Plugin.runAmbertools(self, 'pdb4amber', params, cwd=self._getExtraPath())

    def antechamberStep(self):
        systemBasename = self.getInputBaseName()

        params = ' -i {}.pdb -fi pdb -o {}.mol2 -fo mol2 '.format(*[systemBasename]*2)
        params += '-c {} '.format(anteChDic[self.getEnumText('ChargeModel')])
        params += '-s {}'.format(statusDic[self.getEnumText('Status')])

        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getExtraPath())

    def parmStep(self):
        systemBasename = self.getInputBaseName()
        params = '-i {}.mol2 -o {}.frcmod -f mol2 '.format(*[systemBasename]*2)

        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self._getExtraPath())

    def leapStep(self):
        inputStructure = os.path.abspath(self.getConvFile(getBaseFileName(self.getInputFile())))
        systemBasename = os.path.basename(inputStructure.split(".")[0])
        params = 'source leaprc.gaff \n' \
                 'LIG = loadmol2 {}.mol2 \n' \
                 'loadamberparams {}.frcmod \n' \
                 'saveoff LIG {}.lib \n' \
                 'saveamberparm LIG {}.prmtop {}.rst7 \n' \
                 'savepdb LIG {}_check.pdb \n' \
                 'quit'.format(*[systemBasename]*6)

        file = open(self._getExtraPath("leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f leap_commands.txt", cwd=self._getExtraPath())


    def createOutputStep(self):
        systemBasename = self.getInputBaseName()

        topoFile = abspath(self._getPath('{}.prmtop'.format(systemBasename)))
        crdFile = abspath(self._getPath('{}.rst7'.format(systemBasename)))
        checkFile = abspath(self._getPath('{}_check.pdb'.format(systemBasename)))
        os.rename(abspath(self._getExtraPath('{}.prmtop'.format(systemBasename))), topoFile)
        os.rename(abspath(self._getExtraPath('{}.rst7'.format(systemBasename))), crdFile)
        os.rename(abspath(self._getExtraPath('{}_check.pdb'.format(systemBasename))), checkFile)

        libFile = abspath(self._getExtraPath('{}.lib'.format(systemBasename)))
        originFile = abspath(self._getExtraPath('{}.pdb'.format(systemBasename)))
        missingparamsFile = abspath(self._getExtraPath('{}.frcmod'.format(systemBasename)))

        amber_system = amberobj.AmberSystem(filename=crdFile, topoFile=topoFile,
                                            checkFile=checkFile)
                                            # libFile=lib_localPath, originLIGFile=origin_localPath,
                                            # missingFile=missingparams_localPath)

        self._defineOutputs(outputSystem=amber_system)
        self._defineSourceRelation(self.inputSetOfMols, amber_system)

    # --------------------------- INFO functions -----------------------------------

    def convertPDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = self.getConvFile(inName)

        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
        runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        return oFile

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            summary.append(
                "This protocol has created a coordinate file, a topology file and a PDB file (for visualization)"
                "from a ligand molecule")

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

    def _validate(self):
        vals = []
        if not self.inputLigand.get():
            vals.append('You must specify the ligand to prepare.')

    def getInputMolecule(self):
        molName = self.inputLigand.get()
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == molName:
                return mol

    def getInputFile(self):
        inputLigand = self.getInputMolecule()
        ligandFile = inputLigand.getPoseFile()
        if not ligandFile:
            ligandFile = inputLigand.getFileName()
        return ligandFile

    def getConvFile(self, inName):
        return os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

    def getInputBaseName(self):
        return getBaseFileName(self.getInputFile())















