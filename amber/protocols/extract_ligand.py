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
This module will extract the ligand from a complex pdb.
"""

import glob
from pyworkflow.utils.path import copyFile
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re
from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes
import csv
import pyworkflow.object as pwobj
from pwchem.utils import clean_PDB

scriptName = 'extract_ligand_script.py'

class ExtractStructures(EMProtocol):
    _label = 'Extract Ligand'


    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure:", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure from which you want to extract a ligand')
    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('ExtractStep')
        self._insertFunctionStep('CreateOutputStep')

    def extractStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        params = '{}.pdb'.format(systemBasename)
        Plugin.runRDKitScript(self, scriptName, params, cwd=self._getPath())

    def CreateOutputStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath())










