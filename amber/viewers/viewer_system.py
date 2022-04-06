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
import pyworkflow.protocol.params as params
from pwchem.viewers import VmdViewPopen
from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import ChimeraViewer

from pwchem.viewers import PyMolViewer, PyMolView
from pwchem.utils import natural_sort

from amber import Plugin
from ..objects import AmberSystem
from ..protocols import AmberSystemPrep
from ..constants import *


class AmberSystemViewer(pwviewer.Viewer):
    _label = 'Viewer AMBER system'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [AmberSystem]

    def _visualize(self, obj, **kwargs):
        amberFile = os.path.abspath(obj.getCheckFile())

        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(amberFile, cwd=os.path.dirname(amberFile))
