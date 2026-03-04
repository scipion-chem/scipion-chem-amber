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