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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.utils import Message

import amber
from pwchem.utils import runOpenBabel

import amber.objects as amberobj


class AmberLigandPrep(EMProtocol):
    """
With this protocol you will obtain coordinate and topology files from your ligand
using the pdb4amber and Antechamber programs from AMEBERTOOLS
    """


    _label = 'Ligand preparation'
    IMPORT_FROM_FILE = 1
    IMPORT_FROM_SCIPION = 1

    _ChargeModel = ['']


    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    # -------------------------- DEFINE param functions ----------------------


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

        group = form.addGroup('Ligand parametrization')
        group.addParam('ChargeModel', params.EnumParam,
                       choices=self._ChargeModel, default='AM1-BCC',
                       label='Choose the charge model in order to calculate the atomic point charges: ')
        group.addParam('Verbosity', params.EnumParam,
                       choices=self._Verbosity, default='AM1-BCC',
                       label='Choose the charge model in order to calculate the atomic point charges: ')


    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('PrepStep')
        self._insertFunctionStep('AntechamberStep')
        self._insertFunctionStep('ParmStep')
        self._insertFunctionStep('createOutputStep')


    def PrepStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)

        systemBasename = os.path.basename(inputStructure.split(".")[0])
        params = '{} > {}.LIG.pdb --no-conect '.format(inputStructure, systemBasename)

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
        params = ' -i {}.LIG.pdb -fi pdb -o LIG.mol2 -c bbc -s 2 '


        amber.Plugin.runAmbertools(self, 'antechamber', params, cwd=self._getPath())
