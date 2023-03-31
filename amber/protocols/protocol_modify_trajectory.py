# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
This protocol modifies the trajectory of an amber system using cpptraj
"""
import os
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from amber import Plugin as amberPlugin
from amber.objects import AmberSystem


class AmberModifySystem(EMProtocol):
    """
    This protocol modifies a amber system trajectory and/or coordinates:
        - Cleans from waters and ions
        - Subsamples trajectory and applies filters
        - Fits trajectory to initial structure
    https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml#magicparlabel-1630
    """
    _label = 'system modification'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('amberSystem', params.PointerParam, label="Input Amber System: ",
                      pointerClass='AmberSystem',
                      help='Amber system to be modified')
        group = form.addGroup('Cleaning')
        group.addParam('cleaningW', params.BooleanParam, label="Clean water: ", default=False,
                       help='Remove waters from the system')
        group.addParam('cleaningI', params.BooleanParam, label="Clean ions: ", default=False,
                       help='Remove ions from the system')
        group = form.addGroup('Fitting')
        group.addParam('doFit', params.BooleanParam, label="Fit trajectory?: ", default=False,
                       help='Fit trajectory to initial structure')
        group.addParam('fitting', params.EnumParam, label="Fitting type: ", default=0, condition='doFit',
                       choices=['rot+trans', 'rotxy+transxy', 'translation', 'transxy', 'progressive'],
                       help='Fitting technique to the initial structure. '
                            'https://amberhub.chpc.utah.edu/cpptraj/running-cpptraj/')
        group = form.addGroup('Cutting')
        group.addParam('doDrop', params.BooleanParam, label="Cut trajectory?: ", default=False,
                       help='Cut a trajectory, saving only from first to last frames ')
        line = group.addLine('Frames: ', condition='doDrop',
                             help='First and last frames to save from trajectory (0 from first, 0 until last)')
        line.addParam('firstFrame', params.FloatParam, label="First: ", default=0)
        line.addParam('lastFrame', params.FloatParam, label="Last: ", default=0)

        group = form.addGroup('Subsample')
        group.addParam('doSubsample', params.BooleanParam, label="Subsample trajectory?: ", default=False,
                       help='Subsample trajectory frames')
        group.addParam('subsampleF', params.IntParam, label="Subsample factor: ", default=10, condition='doSubsample',
                       help='Subsample factor. Take a frame for each x original frames')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('modifySystem')
        self._insertFunctionStep('createOutputStep')

    def modifySystem(self):
        inSystem = self.amberSystem.get()
        inputStructure = os.path.abspath(inSystem.getSystemFile())
        inputTopology = os.path.abspath(inSystem.getTopologyFile())
        inputTrajectory = os.path.abspath(inSystem.getTrajectoryFile())
        print('Inputs: ', inputStructure, inputTopology, inputTrajectory)

        cleanTop, cleanTraj, cleanStruct = \
            self.getCleanTopologyFile(), self.getCleanTrajectoryFile(), self.getCleanStructureFile()

        cleaningStr = ''
        if self.cleaningW:
            cleaningStr += 'strip :WAT\n'
        if self.cleaningI:
            ionResStr = ','.join(inSystem.getIonResNames())
            cleaningStr += 'strip :{}\n'.format(ionResStr)

        if cleaningStr:
            cppStr = 'parm {}\ntrajin {}\n{}\ntrajout {}\ntrajout {} start 1 stop 1\ngo'.\
              format(inputTopology, inputTrajectory, cleaningStr, cleanTraj, cleanStruct)
            amberPlugin.runCPPTRAJ(self, self._getExtraPath('cleaning.txt'), cppStr)

            cppStr = 'parm {}\nparmstrip {}\nparmwrite out {}\ngo'. \
                format(inputTopology, cleaningStr, cleanTop)
            amberPlugin.runCPPTRAJ(self, self._getExtraPath('cleaning.txt'), cppStr)
        else:
            cleanTop, cleanTraj, cleanStruct = inputTopology, inputTrajectory, inputStructure

        if inputTrajectory:
            if self.doFit or self.doDrop or self.doSubsample:
                outTrj = os.path.abspath(self._getPath('modified.nc'))

                trjArgs = "parm {}\ntrajin {}\n".format(cleanTop, cleanTraj)
                if self.doFit:
                    trjArgs += 'autoimage\n'

                trjArgs += 'trajout {} '.format(outTrj)
                if self.doDrop:
                    trjArgs += 'start {} stop {}'.format(self.firstFrame.get(), self.lastFrame.get())

                if self.doSubsample:
                    trjArgs += 'offset {}'.format(self.subsampleF.get())
                trjArgs += '\n'

                amberPlugin.runCPPTRAJ(self, self._getExtraPath('modified.txt'), trjArgs, cwd=self._getPath())
            else:
                outTrj = cleanTraj
          
        self.saveFiles(cleanStruct, cleanTop, outTrj)


    def createOutputStep(self):
      outStruct, outTop, outTraj = self.getOutFiles()
        
      outSystem = AmberSystem()
      outSystem.setSystemFile(outStruct)
      outSystem.setTopologyFile(outTop)
      if self.amberSystem.get().getTrajectoryFile():
          outSystem.setTrajectoryFile(outTraj)

      self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def getCleanStructureFile(self):
        inputStructure = self.amberSystem.get().getFileName()
        name, ext = os.path.splitext(inputStructure)
        return os.path.abspath(self._getPath(os.path.basename(name)) + '.rst')

    def getCleanTrajectoryFile(self):
        inputTrajectory = self.amberSystem.get().getTrajectoryFile()
        return os.path.abspath(self._getPath(os.path.basename(inputTrajectory.split(".")[0])) + '.nc')

    def getCleanTopologyFile(self):
        inputTopology = self.amberSystem.get().getTopologyFile()
        return os.path.abspath(self._getPath(os.path.basename(inputTopology.split(".")[0])) + '.top')

    def saveFiles(self, finalStruct, finalTop, finalTraj):
        inSystem = self.amberSystem.get()
        inputStructure = os.path.abspath(inSystem.getSystemFile())
        inputTopology = os.path.abspath(inSystem.getTopologyFile())
        inputTrajectory = os.path.abspath(inSystem.getTrajectoryFile())
        
        outStruct, outTop, outTraj = self.getOutFiles()
        if finalStruct == inputStructure:
            shutil.copy(finalStruct, outStruct)
        else:
            os.rename(finalStruct, outStruct)
            
        if finalTop == inputTopology:
            shutil.copy(finalTop, outTop)
        else:
            os.rename(finalTop, outTop)
            
        if finalTraj == inputTrajectory:
            shutil.copy(finalTraj, outTraj)
        else:
            os.rename(finalTraj, outTraj)
        
    def getOutFiles(self):
        return os.path.abspath(self._getPath('outputSystem.rst')), os.path.abspath(self._getPath('outputSystem.top')), \
               os.path.abspath(self._getPath('outputSystem.nc'))
    