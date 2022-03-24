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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.utils import Message

import amber
from pwchem.utils import runOpenBabel

AMBERFF_FF19SB = 0
AMBERFF_0L15 = 1
AMBERFF_OL3 = 2
AMBERFF_GLYCAM_06J = 3
AMBERFF_LIPID17 = 4
AMBERFF_GAFF2 = 5

AMBER_MAINFF_NAME = dict()
AMBER_MAINFF_NAME[AMBERFF_FF19SB] = 'protein.ff19SB'
AMBER_MAINFF_NAME[AMBERFF_0L15] = 'DNA.OL15'
AMBER_MAINFF_NAME[AMBERFF_OL3] = 'RNA.OL3'
AMBER_MAINFF_NAME[AMBERFF_GLYCAM_06J] = 'GLYCAM_06J-1'
AMBER_MAINFF_NAME[AMBERFF_LIPID17] = 'lipid17'
AMBER_MAINFF_NAME[AMBERFF_GAFF2] = 'gaff2'

AMBERFF_LIST = [AMBER_MAINFF_NAME[AMBERFF_FF19SB], AMBER_MAINFF_NAME[AMBERFF_0L15],
                AMBER_MAINFF_NAME[AMBERFF_OL3], AMBER_MAINFF_NAME[AMBERFF_GLYCAM_06J],
                AMBER_MAINFF_NAME[AMBERFF_LIPID17], AMBER_MAINFF_NAME[AMBERFF_GAFF2]]

AMBER_TIP4 = 0
AMBER_TIP4PEW = 1
AMBER_TIP5P = 2
AMBER_SPCE = 3
AMBER_SPCEB = 4
AMBER_OPC = 5
AMBER_OPC3 = 6
AMBER_POL3 = 7
AMBER_TIP3PFB = 8
AMBER_TIP4PFB = 9

AMBER_WATERFF_NAME = dict()
AMBER_WATERFF_NAME[AMBER_TIP4] = 'tip4'
AMBER_WATERFF_NAME[AMBER_TIP4PEW] = 'tip4pew'
AMBER_WATERFF_NAME[AMBER_TIP5P] = 'tip5p'
AMBER_WATERFF_NAME[AMBER_SPCE] = 'spce'
AMBER_WATERFF_NAME[AMBER_SPCEB] = 'spceb'
AMBER_WATERFF_NAME[AMBER_OPC] = 'opc'
AMBER_WATERFF_NAME[AMBER_OPC3] = 'opc3'
AMBER_WATERFF_NAME[AMBER_POL3] = 'pol3'
AMBER_WATERFF_NAME[AMBER_TIP3PFB] = 'tip3pfb'
AMBER_WATERFF_NAME[AMBER_TIP4PFB] = 'tip4pfb'

AMBER_WATERS_LIST = [AMBER_WATERFF_NAME[AMBER_TIP4], AMBER_WATERFF_NAME[AMBER_TIP4PEW],
                     AMBER_WATERFF_NAME[AMBER_TIP5P], AMBER_WATERFF_NAME[AMBER_SPCE],
                     AMBER_WATERFF_NAME[AMBER_SPCEB], AMBER_WATERFF_NAME[AMBER_OPC],
                     AMBER_WATERFF_NAME[AMBER_OPC3], AMBER_WATERFF_NAME[AMBER_POL3],
                     AMBER_WATERFF_NAME[AMBER_TIP3PFB], AMBER_WATERFF_NAME[AMBER_TIP4PFB]]


class AmberSystemPrep(EMProtocol):
    """
    This protocol will prepare a system for MD simulation. It will clean the input PDB for further analysis
    and generate topology and coordinate files necessary for MD simulation.

    """

    _label = 'system preparation'
    IMPORT_FROM_FILE = 1
    IMPORT_FROM_SCIPION = 1

    # IMPORT_MDP_FILE = 0
    # IMPORT_MDP_SCIPION = 1

    # -------------------------- DEFINE param functions ----------------------

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self, form):
        """
        Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure:", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure to convert to Amber system')

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

        form.addSection('MD prep')
        group = form.addGroup('Force field')
        group.addParam('MainForceField', params.EnumParam,
                       choices=AMBERFF_LIST, default=AMBERFF_FF19SB,
                       label='Main Force Field',
                       help='Force field applied to the system. Force fields are sets of potential functions and '
                            'parametrized interactions that can be used to study physical systems')
        group.addParam('WaterForceField', params.EnumParam,
                       choices=AMBER_WATERS_LIST, default=AMBER_OPC,
                       label='Water Force Field: ',
                       help='Force field applied to the waters')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('PDBAmberStep')
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutputStep')

    def PDBAmberStep(self):

        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)

        systemBasename = os.path.basename(inputStructure.split(".")[0])
        params = '{} > {}.amber.pdb '.format(inputStructure, systemBasename)

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

    def preparationStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])
        Mainff = AMBER_MAINFF_NAME[self.MainForceField.get()]
        Waterff = AMBER_WATERFF_NAME[self.WaterForceField.get()]
        params = '\n source leaprc.{} \n' \
                 'source leaprc.water.{}\n' \
                 '{} = loadPdb "{}.amber.pdb"\n'\
                 'saveAmberParm {} {}.top {}.crd'.format(Mainff, Waterff, systemBasename, systemBasename, systemBasename, systemBasename, systemBasename)

        amber.Plugin.runAmbertools(self, 'tleap', params, cwd=self._getPath())

    # ------------------------POR AQUI VOY------------------------#

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions -----------------------------------

    def convertPDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
        runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        return oFile
