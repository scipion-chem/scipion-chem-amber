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
import os
from subprocess import check_call
import pyworkflow.object as pwobj
from pwchem.objects import MDSystem

from amber.constants import _cations, _anions, PROT_RESNAMES

class AmberSystem(MDSystem):
    """A system atom structure (prepared for MD) in the file format of AMBER
   crd : cordinate file .crd
   top : topology file .prmtop
   check : PDB file to visualize the structure
   """

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        self._checkFile = pwobj.String(kwargs.get('checkFile', None))
        self._libFile = pwobj.String(kwargs.get('libFile', None))
        self._originFile = pwobj.String(kwargs.get('originFile', None))
        self._missingFile = pwobj.String(kwargs.get('missingFile', None))
        self._resNames = pwobj.String(kwargs.get('resNames', None))

    def __str__(self):
        cn = self.getClassName()
        try:
            sf = os.path.basename(self.getSystemFile())
        except:
            sf = self.getSystemFile()

        ht = self.hasTrajectory()
        return '{} ({}, hasTrj={})'.format(self.getClassName(), sf,
                                           self.hasTrajectory())

    def getCheckFile(self):
        return self._checkFile.get()

    def setCheckFile(self, value):
        self._checkFile.set(value)

    def getLibFile(self):
        return self._libFile.get()

    def setLibFile(self, value):
        self._libFile.set(value)

    def getOriginFile(self):
        return self._originFile.get()

    def setOriginFile(self, value):
        self._originFile.set(value)

    def getMissingFile(self):
        return self._missingFile.get()

    def setMissingFile(self, value):
        self._missingFile.set(value)

    def getResNames(self):
        if self._resNames.get():
            return self._resNames.get().split(',')

    def setResNames(self, resNames=[], parse=False):
        if parse:
            resNames = self.parseResidueNames()
        self._resNames.set(','.join(resNames))

    def parseResidueNames(self):
        if self.getTopologyFile():
            com = 'cat {} | grep -zoP "%FLAG RESIDUE_LABEL.*\\n%[^%]*" > residueStr.txt'.\
                format(os.path.abspath(self.getTopologyFile()))
            check_call(com, cwd='/tmp', shell=True)

            resNames = set()
            with open('/tmp/residueStr.txt') as fIn:
                for line in fIn:
                    if not '%' in line:
                        for ele in line.split():
                            resNames.add(ele)
            if '\x00' in resNames:
                resNames.remove('\x00')
            return resNames
        else:
            return []

    def getIonResNames(self):
        ionResNames = []
        _ions = list(_cations.values()) + list(_anions.values())
        for rn in self.getResNames():
            if rn in _ions:
                ionResNames.append(rn)
        return ionResNames

    def getLigResNames(self):
        ligResNames = []
        nonLig = list(_cations.values()) + list(_anions.values()) + PROT_RESNAMES + ['WAT']
        for rn in self.getResNames():
            if rn not in nonLig:
                ligResNames.append(rn)
        return ligResNames

