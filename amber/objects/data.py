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
import os, shutil, re
from subprocess import check_call
import pwem.objects.data as data
import pyworkflow.object as pwobj
from pwchem.objects import MDSystem
from pyworkflow.object import Float, Integer, String, Object
from amber import Plugin as amberPlugin

class AmberSystem(MDSystem):
    """A system atom structure (prepared for MD) in the file format of AMBER
   crd : cordinate file .crd
   top : topology file .top
   """

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        self._libFile = pwobj.String(kwargs.get('libFile', None))
        self._crdFile = pwobj.String(kwargs.get('crdFile', None))
        self._nFrames = pwobj.Integer(kwargs.get('nFrames', None))
        self._nTime = pwobj.Float(kwargs.get('nTime', None))

    def __str__(self):
        strStr = '{} ({}'.format(self.getClassName(), os.path.basename(self.getSystemFile()))

        if self.hasTrajectory():
            strStr += ', frames: {}, time(ps): {:.1f}'.format(
                self.getNFrames(),
                self.getNTime()
            )
        strStr += ')'
        return strStr

    def getCrdFile(self):
        return self._crdFile.get()

    def setCrdFile(self, value):
        self._crdFile.set(value)

    def getNFrames(self):
        return self._nFrames.get()

    def setNFrames(self, value):
        self._nFrames.set(value)

    def getNTime(self):
        return self._nTime.get()

    def setNTime(self, value):
        self._nTime.set(value)

    def getFrameIdxs(self):
        """Return [firstFrame, lastFrame] as Python ints.
        Returns [0, 0] if not yet populated.
        """
        return [self._firstFrame.get(), self._lastFrame.get()]

    def hasTrjInfo(self):
        """True only when readTrjInfo has been called and found valid data."""
        return self._lastFrame.get() > 0


    def readTrjInfo(self, protocol, nTime, outDir=None):
        topFile = os.path.abspath(self.getTopologyFile())
        trjFile = os.path.abspath(self.getTrajectoryFile())

        nFrames = self._cpptrajGetNFrames(protocol, topFile, trjFile)

        self.setNFrames(nFrames)
        self.setNTime(nTime)


    # ── helpers ─────────────────────────────

    def _cpptrajGetNFrames(self, protocol, topFile, trjFile):
        """
        Run ``cpptraj -p <top> -y <traj> -tl`` and parse ``Frames: <N>``
        from STDOUT.  Returns None if parsing fails.
        """
        args = '-p {} -y {} -tl'.format(topFile, trjFile)
        amberPlugin.runAmbertools(protocol, program='cpptraj',
                                    args=args)

        for logName in ('run.stdout', 'run.stderr'):
            candidate = protocol._getPath('logs', logName)
            if os.path.exists(candidate):
                with open(candidate) as fh:
                    for line in fh:
                        m = re.search(r'Frames:\s*(\d+)', line)
                        if m:
                            nFrames = int(m.group(1))
                            break
                break

        return nFrames