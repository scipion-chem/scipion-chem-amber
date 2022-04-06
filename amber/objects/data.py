# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Aida Pinacho Pérez
# *
# *
# * your institution
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
import os, shutil
from subprocess import check_call
import pwem.objects.data as data
import pyworkflow.object as pwobj
from pwchem.objects import MDSystem

class AmberSystem(MDSystem):
    """A system atom structure (prepared for MD) in the file format of AMBER
   crd : cordinate file .crd
   top : topology file .prmtop
   check : PDB file to visualize the structure
   """

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        self._checkFile = pwobj.String(kwargs.get('checkFile', None))

    def __str__(self):
        return '{} ({}, hasTrj={})'.format(self.getClassName(), os.path.basename(self.getSystemFile()),
                                           self.hasTrajectory())

    def getCheckFile(self):
        return self._checkFile.get()

    def setCheckFile(self, value):
        self._checkFile.set(value)
