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
This module will perform energy minimizations and equilibrium for the system befor MD simultion
"""
import os, glob, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message, runJob, createLink
from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from amber.objects import *
from amber.constants import *
from amber import Plugin as amberPlugin


class AmberMDSimulation(EMProtocol):
    """
        This protocol will perform energy minimization and equilibrium on the system previously prepared by the protocol
         "system prepartion". This step is necessary to energy minimize the system in order to avoid unwanted conformations.
    """

    _label = 'Molecular dynamics simulation'

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection('Minimization')
        form.addParam('AmberSystem', params.PointerParam, label="Input Amber System: ",
                      pointerClass='AmberSystem',
                      allowsNull=True,
                      help='Amber solvated system to be simulated')
        form.addSection('Equilibrium')

        form.addSection('Simulation')

