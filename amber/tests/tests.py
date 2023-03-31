# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from pwchem.protocols import ProtChemPrepareReceptor, ProtExtractLigands

from ..protocols import AmberSystemPrep


class TestAmberPrepareSystem(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareReceptor(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}')

        cls.launchProtocol(cls.protPrepareReceptor)

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            AmberSystemPrep,
            fromInput=0,
            inputStructure=cls.protPrepareReceptor.outputStructure)

        cls.launchProtocol(protPrepare)
        return protPrepare

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)

        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

class TestAmberPrepareSystemLigand(TestAmberPrepareSystem):
    @classmethod
    def _runExtractLigand(cls, inputProt):
        protExtLig = cls.newProtocol(
            ProtExtractLigands,
            cleanPDB=True, rchains=True, chain_name='{"model": 0, "chain": "C", "residues": 93}')

        protExtLig.inputStructure.set(inputProt)
        protExtLig.inputStructure.setExtended('outputPdb')

        cls.proj.launchProtocol(protExtLig, wait=True)
        cls.protExtLig = protExtLig
        return protExtLig

    @classmethod
    def _runPrepareSystem(cls, protExtract):
        protPrepare = cls.newProtocol(
            AmberSystemPrep,
            fromInput=1, inputLigands=protExtract.outputSmallMolecules,
            inputLigandSelect='SmallMolecule (g1_4erf_0R3-1_1 molecule)')

        cls.launchProtocol(protPrepare)
        return protPrepare

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)

        protPrepare = self._runPrepareSystem(protExtract)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

