# **************************************************************************
# *
# * Authors:     Joaquin Algorta Bove (joaquin.algorta@cnb.csic.es)
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

import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from pwchem.tests.tests_preparations import TestPrepareReceptor
from pwchem.tests.tests_docking import TestExtractLigand
from pwchem.protocols import ProtExtractLigands

from amber.protocols import *
from amber import Plugin as amberPlugin

BASETEST = """{'MaxCycles': 500, 'SdCycles': 250, 'IntCutoff': 8.0,'Restraint': False, 'CustomIn': None, 'stepType': 'Minimization'}
{'MDSteps': 500, 'TimeStep': 0.002, 'Traj': 100, 'SaveTrj': False, 'InTemp': 0, 'FiTemp': 300, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CustomIn': None, 'stepType': 'Heating'}
{'MDSteps': 500, 'TimeStep': 0.002, 'TrajStep': 100, 'SaveTrj': True, 'EnsemType': 'NPT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'Restraint': False, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'Restraint': False, 'RestrAtoms': 'Protein', 'RestrForce': 0.0, 'CustomIn': None, 'stepType': 'Production'}\n"""

LONGTEST = """{'MaxCycles': 10000, 'SdCycles': 5000, 'IntCutoff': 8.0, 'Restraint': False, 'RestrAtoms': 'Protein + Ligand', 'RestrForce': 50.0, 'CustomIn': None, 'stepType': 'Minimization'}
{'MDSteps': 5000, 'TimeStep': 0.002, 'Traj': 500, 'SaveTrj': False, 'InTemp': 0.0, 'FiTemp': 300.0, 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'CustomIn': None, 'Restraint': False, 'RestrAtoms': 'Backbone', 'RestrForce': 50.0, 'stepType': 'Heating'}
{'MDSteps': 5000, 'TimeStep': 0.002, 'TrajStep': 500, 'SaveTrj': True, 'EnsemType': 'NVT', 'Thermostat': 'Langevin', 'CollisFreq': 2.0, 'CoupConst': 2.0, 'FricConst': 2.0, 'Pressure': 1.0, 'Barostat': 'Monte Carlo', 'PressureScaling': 'isotropic', 'CustomIn': None, 'Restraint': False, 'RestrAtoms': 'Backbone', 'RestrForce': 50.0, 'stepType': 'Production'}"""

class TestAmberPrepareSystem(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
        cls.launchProtocol(cls.protImportPDB)  # synchronous, no wait=False

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            AmberSystemPrep,
            inputStructure=cls.protImportPDB.outputPdb,
            inputFrom=0, targetReduce=True, minDist=5)

        cls.launchProtocol(protPrepare)
        return protPrepare

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

class TestAmberPrepareSystemLig(TestPrepareReceptor, TestExtractLigand):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=0, pdbId='1uaz')
      cls.launchProtocol(protImportPDB)
      cls.protImportPDB = protImportPDB

    @classmethod
    def _runExtractLigand(cls, inputProt):
        protExtLig = cls.newProtocol(
            ProtExtractLigands,
            cleanPDB=True, rchains=True, chain_name='{"model": 0, "chain": "A", "residues": 236}')

        protExtLig.inputStructure.set(inputProt)
        protExtLig.inputStructure.setExtended('outputPdb')

        cls.launchProtocol(protExtLig)
        cls.protExtLig = protExtLig
        return protExtLig

    @classmethod
    def _runPrepareSystem(cls, protPrepare, inputFrom=STRUCTURE):
        protPrepareS = cls.newProtocol(
            AmberSystemPrep, inputFrom=inputFrom)

        if inputFrom == STRUCTURE:
            protPrepareS.inputStructure.set(protPrepare)
            protPrepareS.inputStructure.setExtended('outputStructure')
        else:
            protPrepareS.inputSetOfMols.set(protPrepare)
            protPrepareS.inputSetOfMols.setExtended('outputSmallMolecules')
            protPrepareS.inputLigand.set('SmallMolecule (g1_1uaz_RET_255-1_1 molecule)')

        cls.launchProtocol(protPrepareS)
        return protPrepareS

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

class TestAmberPrepareSystemMembrane(TestAmberPrepareSystemLig):

    @classmethod
    def _runPrepareSystem(cls, protPrepare, inputFrom=STRUCTURE):
        protPrepareM = cls.newProtocol(
            AmberSystemPrep, inputFrom=inputFrom, tMem=True, memDistXY=5, memDistZ=5)

        if inputFrom == STRUCTURE:
            protPrepareM.inputStructure.set(protPrepare)
            protPrepareM.inputStructure.setExtended('outputStructure')
        else:
            protPrepareM.inputSetOfMols.set(protPrepare)
            protPrepareM.inputSetOfMols.setExtended('outputSmallMolecules')
            protPrepareM.inputLigand.set('SmallMolecule (g1_1uaz_RET-1_1 molecule)')

        cls.launchProtocol(protPrepareM)
        return protPrepareM

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestAmberCpuSimulation(TestAmberPrepareSystem):

    def _runSimulation(self, protPrepare):
        protSim = self.newProtocol(
            AmberMDSimulation,
            amberSystem=protPrepare.outputSystem, workFlowSteps=BASETEST, useGpu=False)
        protSim.setObjLabel('amber - sander MD sim')

        self.launchProtocol(protSim)
        return protSim


    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

class TestAmberGpuSimulation(TestAmberPrepareSystem):

    def _runSimulation(self, protPrepare):
        protSim = self.newProtocol(
            AmberMDSimulation,
            amberSystem=protPrepare.outputSystem, workFlowSteps=BASETEST)
        protSim.setObjLabel('amber - sander MD sim')

        self.launchProtocol(protSim)
        return protSim


    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))

class TestAmberLigSimulation(TestAmberPrepareSystemLig):

    def _runSimulation(self, protPrepares):
        protSim = self.newProtocol(
            AmberMDSimulation,
            amberSystem=protPrepares.outputSystem, workFlowSteps=LONGTEST)
        protSim.setObjLabel('amber - pmemd MD sim')

        self.launchProtocol(protSim)
        return protSim

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))