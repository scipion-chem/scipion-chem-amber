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
from pwchem.viewers import VmdViewPopen
from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import ChimeraViewer

from pwchem.viewers import PyMolViewer, PyMolView
from pwchem.utils import natural_sort

from amber import Plugin
from ..objects import AmberSystem
from ..protocols import AmberMDSimulation
from ..constants import *


class AmberSystemViewer(pwviewer.Viewer):
    _label = 'Viewer AMBER system'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [AmberSystem]

    def _visualize(self, obj, **kwargs):
        amberFile = os.path.abspath(obj.getCheckFile())

        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(amberFile, cwd=os.path.dirname(amberFile))


class AmberMDSimulationViewer(pwviewer.ProtocolViewer):
    ''' Visualize the trajectory output'''
    _label = 'Viewer AMBER Simulation'
    _targets = [AmberMDSimulation]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        form.addSection(label='Visualization of AMBER Simulation')
        group = form.addGroup('Open MD simulation')
        group.addParam('chooseStage', params.EnumParam,
                       choices=self._getStagesWTrj(), default=0,
                       label='Choose the stage to analyze: ',
                       help='Choose the simulation stage to analyze')
        group.addParam('displayMdVMD', params.LabelParam,
                       label='Display trajectory with VMD: ',
                       help='Display trajectory with VMD. \n'
                            'Protein represented as NewCartoon and waters as dots'
                       )

    def _getVisualizeDict(self):
        return {
            'displayMdVMD': self._showMdVMD,
        }

    def _showMdVMD(self, paramName=None):
        stage = self.getEnumText('chooseStage')
        _, trjFile = self.getStageFiles(stage)
        system = self.protocol.outputSystem
        topFile = self.AmberSystem.get().getTopologyFile()

        params = 'mol addrep 0' \
                 'mol new s%ds type {parm7} first 0 last -1 step 1 waitfor 1' \
                 'mol addfile s%ds type {netcdf} first 0 last -1 step 1 waitfor 1 0' % (topFile, trjFile)

        outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
        with open(outTcl, 'w') as f:
            f.write(params)
        args = '-e {}'.format(outTcl)

        return [VmdViewPopen(args)]

    ################################# UTILS #########################

    def _getStagesWTrj(self):
        '''Return stages with a saved trajectory'''
        stages = ['All']
        for stDir in natural_sort(glob.glob(self.protocol._getExtraPath('stage_*'))):
            stage = os.path.basename(stDir)
            trjFile = '{}/{}.trr'.format(stDir, stage)
            if os.path.exists(trjFile):
                stages.append(stage)
        return stages

    def getStageFiles(self, stage):
        if stage == 'All':
            _, _, tprFile = self.protocol.getPrevFinishedStageFiles(reverse=True)
            trjFile = self.protocol._getPath('outputTrajectory.nc')
        else:
            _, _, tprFile = self.protocol.getPrevFinishedStageFiles(stage)
            trjFile = self.protocol._getExtraPath('{}/{}_corrected.nc'.format(stage, stage))

        return tprFile, trjFile
