# **************************************************************************
# *
# * Authors:     Aida Pinacho Pérez
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
import os, glob, subprocess
import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
from pwchem.viewers import PyMolViewer, PyMolView, VmdViewPopen, MDSystemPViewer
# from pwchem.viewers.viewers_data import PML_MD_STR

from pwchem.utils import natural_sort
# from pwchem.constants import TCL_MD_STR, PML_MD_STR

from amber import Plugin
from ..objects import AmberSystem
from ..protocols import AmberMDSimulation
from ..constants import *

PML_MD_STR = '''load {}
load_traj {}, format=trj
hide everything, not br. all within 3 of (byres polymer & name CA)
set movie_fps, 15
'''

class AmberSystemViewer(pwviewer.Viewer):
  _label = 'Viewer Molecular Dynamics system'
  _environments = [pwviewer.DESKTOP_TKINTER]
  _targets = []

  def _visualize(self, obj, onlySystem=False, trjFile=None, **kwargs):
    systemFile = os.path.abspath(obj.getSystemFile())
    topoFile = os.path.abspath(obj.getTopologyFile())
    if not trjFile:
        trjFile = obj.hasTrajectory()

    if not trjFile or onlySystem:
        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(topoFile, cwd=os.path.dirname(topoFile))

    else:
        trjFile = os.path.abspath(obj.getTrajectoryFile())
        outPml = os.path.join(os.path.dirname(trjFile), 'pymolSimulation.pml')
        with open(outPml, 'w') as f:
          f.write(PML_MD_STR.format(os.path.abspath(topoFile),
                                    os.path.abspath(trjFile)))

        return [PyMolView(os.path.abspath(outPml), cwd=os.path.dirname(trjFile))]

class AmberSystemPViewer(MDSystemPViewer):
    """ Visualize the output of Molecular Dynamics simulation """
    _label = 'Viewer Molecular Dynamics System'
    _targets = [AmberSystem]

    def __init__(self, **args):
      super().__init__(**args)

    def _defineParams(self, form):
      super()._defineParams(form)

    def _getAnalysisTopFile(self):
      return self.getMDSystem().getTopologyFile()

    def _showMdPymol(self, paramName=None):
      system = self.getMDSystem()
      return AmberSystemViewer(project=self.getProject())._visualize(system)

    def _showMdVMD(self, paramName=None):
      system = self.getMDSystem()

      outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
      sysExt = os.path.splitext(system.getTopologyFile())[1][1:]
      trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
      self.writeTCL(outTcl, system.getTopologyFile(), sysExt, system.getTrajectoryFile(), 'netcdf')

      args = '-e {}'.format(outTcl)
      return [VmdViewPopen(args)]


class AmberSimulationViewer(AmberSystemPViewer):
    """ Visualize an Amber MD simulation, allowing the analysis to be restricted to the
    trajectory of a single simulation stage (minimization / heating / production / custom). """
    _label = 'Viewer Amber Simulation'
    _targets = [AmberMDSimulation]

    def __init__(self, **args):
      super().__init__(**args)

    def _defineSimParams(self, form):
      '''Mirror the base "Open MD simulation" group but prepend a stage selector so the
      trajectory shown with PyMol/VMD can be limited to one simulation stage.'''
      group = form.addGroup('Open MD simulation')
      group.addParam('chooseStage', params.EnumParam,
                     choices=self._getStagesWTrj(), default=0,
                     label='Choose the stage to display: ',
                     help='Restrict trajectory visualization to a single simulation stage.\n'
                          '"All" (default) uses the full concatenated production trajectory, '
                          'preserving the previous behavior.')
      group.addParam('displayMdPymol', params.LabelParam,
                     label='Display trajectory with PyMol: ',
                     help='Display the selected stage trajectory with PyMol.')
      group.addParam('displayMdVMD', params.LabelParam,
                     label='Display trajectory with VMD: ',
                     help='Display the selected stage trajectory with VMD.')

    def _defineMDTrajParams(self, form):
      '''Add the inherited MDTraj analysis section and place an independent stage selector.'''
      super()._defineMDTrajParams(form)
      group = form.getParam('MDTraj_analysis') or form.addGroup('MDTraj analysis')
      group.addParam('chooseStageAnalysis', params.EnumParam,
                     choices=self._getStagesWTrj(), default=0,
                     label='Choose the stage to analyze: ',
                     help='Restrict the MDTraj analysis to a single simulation stage.\n'
                          '"All" (default) analyzes the full concatenated production trajectory.')

    def _getStagesWTrj(self):
      '''Stages with a saved trajectory, plus the "All" (full trajectory) default.'''
      return ['All'] + self.protocol.getStagesWithTrj()

    def getMDSystem(self, objType=AmberSystem):
      system = super().getMDSystem(objType)
      stage = getattr(self, '_activeStage', 'All')
      if system and stage and stage != 'All':
        trjFile = self.protocol.getStageTrjFile(stage)
        if trjFile and os.path.exists(trjFile):
          system = system.clone()
          system.setTrajectoryFile(trjFile)
      return system

    def _withStage(self, stageParam, func, paramName):
      '''Run an inherited display/analysis method with getMDSystem scoped to the stage chosen
      in stageParam, then restore the default scope.'''
      self._activeStage = self.getEnumText(stageParam)
      try:
        return func(paramName)
      finally:
        self._activeStage = 'All'

    # Display (PyMol / VMD) -> scoped by 'chooseStage'
    def _showMdPymol(self, paramName=None):
      return self._withStage('chooseStage', super()._showMdPymol, paramName)

    def _showMdVMD(self, paramName=None):
      return self._withStage('chooseStage', super()._showMdVMD, paramName)

    # MDTraj analysis -> scoped by 'chooseStageAnalysis'
    def _showMDTrajAnalysis(self, paramName=None):
      return self._withStage('chooseStageAnalysis', super()._showMDTrajAnalysis, paramName)