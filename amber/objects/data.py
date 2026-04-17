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

        self._firstFrame = Integer(0)
        self._lastFrame = Integer(0)

    def __str__(self):
        strStr = '{} ({}'.format(self.getClassName(), os.path.basename(self.getSystemFile()))

        if self.hasTrajectory():
            strStr += ', frames: {} - {}, time(ps): {:.1f}'.format(
                *self.getFrameIdxs(),
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

    def getNTime(self):
        return self._nTime.get()

    def setNTime(self, value):
        self._nTime.set(value)

    def setFrameIdxs(self, idxs):
        """Store first and last frame indices (1-based integers).

        Parameters
        ----------
        idxs : list | tuple
            [firstFrame, lastFrame]
        """
        self._firstFrame.set(int(idxs[0]))
        self._lastFrame.set(int(idxs[1]))

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
        outDir = os.path.dirname(self.getTrajectoryFile()) if not outDir else outDir

        nFrames = self._cpptrajGetNFrames(protocol, topFile, trjFile, outDir)

        firstFrame = 1
        lastFrame = nFrames if nFrames is not None else 1
        self.setFrameIdxs([firstFrame, lastFrame])

        self.setNTime(nTime)


    # ── helpers (add these as methods of AmberSystem) ─────────────────────────────

    def _cpptrajGetNFrames(self, protocol, topFile, trjFile, outDir):
        """
        Run ``cpptraj -p <top> -y <traj> -tl`` and parse ``Frames: <N>``
        from STDOUT.  Returns None if parsing fails.
        """

        logFile = os.path.abspath(os.path.join(outDir, 'trjinfo_tl.log'))

        # runAmberProgram must redirect stdout; we capture via -o flag if available,
        # otherwise read the scipion stdout log produced by the protocol runner.
        args = '-p {} -y {} -tl -o {}'.format(topFile, trjFile, logFile)
        amberPlugin.runAmbertools(protocol, program='cpptraj',
                                    args=args, cwd=outDir)

        nFrames = None
        if os.path.exists(logFile):
            with open(logFile) as fh:
                for line in fh:
                    m = re.search(r'Frames:\s*(\d+)', line)
                    if m:
                        nFrames = int(m.group(1))
                        break

        # Fallback: parse from the protocol run.stdout / run.stderr logs
        if nFrames is None:
            for logName in ('run.stdout', 'run.stderr'):
                candidate = protocol._getPath('logs', logName)
                if os.path.exists(candidate):
                    with open(candidate) as fh:
                        for line in fh:
                            m = re.search(r'Frames:\s*(\d+)', line)
                            if m:
                                nFrames = int(m.group(1))
                                break
                if nFrames is not None:
                    break

        return nFrames