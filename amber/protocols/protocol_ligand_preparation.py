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

import amber.objects as amberobj


class AmberLigandPrep(EMProtocol):
    """
With this protocol you will obtain coordinate and topology files from your ligand
using the pdb4amber and Antechamber programs from AMEBERTOOLS
    """


    _label = 'Ligand preparation'

    _ChargeModel = ['RESP', 'AM1-BCC', 'CM1', 'CM2', 'ESP', 'Mulliken', 'Gasteiger']
    _ChargeModelWorking = ['AM1-BCC', 'Mulliken', 'Gasteiger']
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
        form.addParam('netCharge', params.IntParam, default=0, expertLevel=params.LEVEL_ADVANCED,
                       label='Net Charge: ',
                       help="Enter the integer net charge of the molecule. \n"
                            "Common scenarios:\n"
                            "- 0: Neutral molecules (most drugs/ligands).\n"
                            "- -1: Deprotonated acids (e.g., carboxylates, phosphates).\n"
                            "- +1: Protonated bases (e.g., amines at physiological pH).\n\n"
                            "If antechamber reports an 'odd number of electrons', your charge is likely "
                            "mismatched with your structure's protonation state.")

        # group = form.addGroup('pdb4amber options')
        # group.addParam('proteinResidues', params.BooleanParam, default=False,
        #                label='Keep only protein residues: ')
        # group.addParam('AmberCompatibleResidues', params.BooleanParam, default=False,
        #                label='Keep only Amber compatible residues: ')
        # group.addParam('phSimulation', params.BooleanParam, default=False,
        #                label='Rename GLU, ASP, HIS for constant pH simulation: ')
        # group.addParam('reduce', params.BooleanParam, default=False,
        #                label='Run reduce first to add hydrogens: ')
        # group.addParam('tleap', params.BooleanParam, default=False,
        #                label='Use tleap to add missing atoms (EXPERIMENTAL): ')

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
        inputStructure = self.getConvFile(getBaseName(ligandFile))
        if not ligandFile.endswith('.pdb'):
            inputStructure = self.convertPDB(ligandFile)
        else:
            shutil.copy(ligandFile, inputStructure)

        systemBasename = os.path.basename(inputStructure.split(".")[0])
        params = '{} > {}_LIG.pdb '.format(inputStructure, systemBasename)

        amber.Plugin.runAmbertools(self, 'reduce', params, cwd=self._getExtraPath())

    def antechamberStep(self):
        systemBasename = self.getInputBaseName()
        nc = self.netCharge.get()

        params = ' -i {}_LIG.pdb -fi pdb -o {}_LIG.mol2 -fo mol2 -nc {} '.format(*[systemBasename]*2,nc)


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

        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getExtraPath())

    def parmStep(self):
        systemBasename = self.getInputBaseName()
        params = '-i {}_LIG.mol2 -o {}_LIG.frcmod -f mol2 '.format(*[systemBasename]*2)

        amber.Plugin.runAmbertools(self, 'parmchk2', params, cwd=self._getExtraPath())

    def leapStep(self):
        # inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        # systemBasename = os.path.basename(inputStructure.split(".")[0])
        systemBasename = self.getInputBaseName()
        params = 'source leaprc.gaff \n' \
                 'LIG = loadmol2 {}_LIG.mol2 \n' \
                 'loadamberparams {}_LIG.frcmod \n' \
                 'saveoff LIG {}_lig.lib \n' \
                 'saveamberparm LIG {}_LIG.top {}_LIG.crd \n' \
                 'savepdb LIG {}.check.pdb \n' \
                 'quit'.format(*[systemBasename]*6)

        file = open(self._getExtraPath("leap_commands.txt"), "w")
        file.write(params)
        file.close()

        amber.Plugin.runAmbertools(self, 'tleap ', "-f leap_commands.txt", cwd=self._getExtraPath())


    def createOutputStep(self):
        systemBasename = self.getInputBaseName()

        topol_baseName = '{}_LIG.top'.format(systemBasename)
        crd_baseName = '{}_LIG.crd'.format(systemBasename)
        lib_baseName = '{}_LIG.lib'.format(systemBasename)
        origin_baseName = '{}_LIG.pdb'.format(systemBasename)
        check_baseName = '{}_checkLIG.pdb'.format(systemBasename)
        missingparams_baseName = '{}_LIG.frcmod'.format(systemBasename)

        topol_localPath = abspath(self._getPath(topol_baseName))
        crd_localPath = abspath(self._getPath(crd_baseName))
        lib_localPath = abspath(self._getPath(lib_baseName))
        origin_localPath = abspath(self._getPath(origin_baseName))
        check_localPath = abspath(self._getPath(check_baseName))
        missingparams_localPath = abspath(self._getPath(missingparams_baseName))

        amber_system = amberobj.AmberSystem(crdLIGname=crd_localPath, topoLIGFile=topol_localPath,
                                           libFile=lib_localPath, originLIGFile=origin_localPath,
                                           checkLIGFile=check_localPath, missingFile=missingparams_localPath)

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
        return os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdb')))

    def getInputBaseName(self):
        return getBaseName(self.getInputFile())















