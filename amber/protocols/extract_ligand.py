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
from pwem.convert.atom_struct import toPdb
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

from os.path import relpath, abspath

import os

from pwem.protocols import EMProtocol, ProtImportFiles
from pyworkflow.protocol import params, protocol
from pyworkflow.utils import Message

import amber
from amber import Plugin as amberPlugin
from pwchem.utils import runOpenBabel

import amber.objects as amberobj
from amber.objects import *

import os
from Bio.PDB import PDBParser, PDBIO, Select

scriptName = 'extract_ligand_script.py'


def is_het(residue):
    res = residue.id[0]
    return res != " " and res != "W"


class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        """ Recognition of heteroatoms - Remove water molecules """
        return residue == self.residue and is_het(residue)


class ExtractStructures(EMProtocol):
    _label = 'Extract Ligand'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure:", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure from which you want to extract a ligand')

    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('CreateOutputStep')


    def CreateOutputStep(self):
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())
        if not inputStructure.endswith('.pdb'):
            inputStructure = self.convertPDB(inputStructure)
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        ligandFiles = self.extract_ligands(inputStructure)

        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
        for lFile in ligandFiles:
            oMol = SmallMolecule(smallMolFilename=lFile)
            outputSet.append(oMol)
        self._defineOutputs(outputSmallMolecules=outputSet)

    # --------------------------- INFO functions -----------------------------------

    def convertPDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdb')))
        oFile = toPdb(os.path.abspath(proteinFile), oFile)

        return oFile

    def extract_ligands(self, inputStructure):
        """ Extraction of the heteroatoms of .pdb files """
        SmallMolList = []
        i = 0
        pdb_code = inputStructure[:-4].split('/')[-1]
        print(pdb_code)
        print(inputStructure)
        pdb = PDBParser().get_structure(pdb_code, inputStructure)
        print(pdb)
        io = PDBIO()
        io.set_structure(pdb)
        for model in pdb:
            for chain in model:
                for residue in chain:
                    if not is_het(residue):
                        continue
                    print(f"saving {chain} {residue}")
                    molFile = self._getExtraPath(f"lig_{pdb_code}_{i}.pdb")
                    io.save(molFile, ResidueSelect(chain, residue))
                    i += 1

                    SmallMolList.append(molFile)

        return SmallMolList
