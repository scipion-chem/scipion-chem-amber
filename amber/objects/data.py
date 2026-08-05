# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
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
import os, shutil, re
import numpy as np

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
        self._nTime = pwobj.Float(kwargs.get('nTimeNs', None))  # total trajectory time, ns

    def __str__(self):
        strStr = '{} ({}'.format(self.getClassName(), os.path.basename(self.getSystemFile()))

        if self.hasTrajectory():
            strStr += ', frames: {}'.format(self.getNFrames())
            nTimeNs = self.getNTimeNs()
            if nTimeNs is not None:
                strStr += ', time: {:.3f} ns'.format(nTimeNs)
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

    def getNTimeNs(self):
        """Total trajectory time in nanoseconds, or None if unknown."""
        return self._nTime.get()

    def setNTimeNs(self, valueNs):
        """Store the total trajectory time, in nanoseconds."""
        self._nTime.set(valueNs)

    def getTimeStepNs(self):
        """Time between consecutive saved frames (ns), or None when unknown.
        Derived from the stored total time and frame count"""
        nFrames, nTimeNs = self.getNFrames(), self.getNTimeNs()
        if nFrames and nTimeNs:
            return nTimeNs / nFrames
        return None

    def hasTrjInfo(self):
        """True once readTrjInfo has populated a positive frame count."""
        return bool(self.getNFrames())

    def readTrjInfo(self, protocol, nTimeNs=None):
        """Create trajectory metadata (number of frames and total time in ns)
        by reading it from the trajectory file."""
        topFile = os.path.abspath(self.getTopologyFile())
        trjFile = os.path.abspath(self.getTrajectoryFile())

        self.setNFrames(self._cpptrajGetNFrames(protocol, topFile, trjFile))

        trjTimeNs = self._readTrajTimeNs(trjFile)
        self.setNTimeNs(trjTimeNs if trjTimeNs is not None else nTimeNs)

    # ── helpers ─────────────────────────────

    def _cpptrajGetNFrames(self, protocol, topFile, trjFile):
        """
        Run ``cpptraj -p <top> -y <traj> -tl`` and parse ``Frames: <N>``
        from STDOUT.  Returns None if parsing fails.
        """
        args = '-p {} -y {} -tl'.format(topFile, trjFile)
        amberPlugin.runAmbertools(protocol, program='cpptraj',
                                    args=args)

        nFrames = None
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

    def _readTrajTimeNs(trjFile):
        """Return the total elapsed time (ns) actually stored in the trajectory:
        the time stamp of its last frame. Returns None when no usable time is present"""
        if not trjFile or not str(trjFile).lower().endswith(('.nc', '.netcdf', '.ncdf')):
            return None
        try:
            from scipy.io import netcdf_file
            nc = netcdf_file(trjFile, 'r', mmap=False)
            try:
                if 'time' not in nc.variables:
                    return None
                times = np.array(nc.variables['time'][:], dtype=float)
            finally:
                nc.close()
        except Exception:
            return None

        if times.size == 0:
            return None
        lastPs = float(times[-1])
        if lastPs <= 0:
            mxPs = float(times.max())
            lastPs = mxPs if mxPs > 0 else 0.0
        return lastPs / 1000.0 if lastPs > 0 else None