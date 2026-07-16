# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
This module modifies an Amber MD system trajectory and/or coordinates using cpptraj.
"""

import os
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from amber import Plugin as amberPlugin
from amber.objects import AmberSystem
from amber.constants import ENV_RES, WATER_RES, ION_RES


# ── Enum index constants (keep in sync with the choices lists in _defineParams) ──
FIT_MASK_BACKBONE = 0
FIT_MASK_CA       = 1
FIT_MASK_ALL      = 2
FIT_MASK_CUSTOM   = 3

STRIP_SOLVENT     = 0   # waters + ions
STRIP_ENV         = 1   # full environment (waters, ions, lipids) → ENV_RES
STRIP_CUSTOM      = 2

OUT_FMT_NC        = 0
OUT_FMT_DCD       = 1
OUT_FMT_XTC       = 2
OUT_FMT_CRD       = 3    # Amber ASCII trajectory (cpptraj keyword 'crd')
OUT_FMT_PDB       = 4

TIME_UNIT_FS = 0
TIME_UNIT_PS = 1
TIME_UNIT_NS = 2


class AmberModifySystem(EMProtocol):
    """
    Modifies an Amber MD system trajectory and/or coordinates using **cpptraj**.

    The requested operations are combined into a single cpptraj script and applied
    in the order cpptraj processes them:

      1. Cutting / subsampling  – frame selection at ``trajin`` (start / stop / offset)
      2. Imaging                – fix PBC artefacts with ``autoimage``
      3. Stripping              – remove solvent / ions; writes a new stripped topology
      4. Fitting                – RMS-fit every frame onto the first frame
      5. Running average        – coordinate running average with ``runavg``
      6. Output                 – write the processed trajectory in the chosen format
    """

    _label = 'System modification'

    # ── param definitions ────────────────────────────────────────────────────
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('amberSystem', params.PointerParam,
                      label='Input Amber System: ',
                      pointerClass='AmberSystem',
                      help='Amber solvated system (topology + trajectory) to process.')

        # ── 1. Cutting ────────────────────────────────────────────────────
        group = form.addGroup('Cutting')
        group.addParam('doDrop', params.BooleanParam,
                       label='Cut trajectory?: ', default=False,
                       help='Keep only frames within a time (or frame-number) window.')
        group.addParam('cutByTime', params.BooleanParam,
                       label='Define window by time?: ', default=True,
                       condition='doDrop',
                       help='If Yes, the window is given in simulation time and converted '
                            'to frame indices using the input trajectory metadata '
                            '(frames / time). If No, frame indices are used directly.')
        line = group.addLine('Time window: ',
                             condition='doDrop and cutByTime',
                             help='Start and end times (0 = use the trajectory limits).')
        line.addParam('firstTime', params.FloatParam, label='Start: ', default=0.0)
        line.addParam('lastTime',  params.FloatParam, label='End: ',   default=0.0)
        line.addParam('timeUnit',  params.EnumParam,
                      label='Units: ', default=TIME_UNIT_NS,
                      choices=['fs', 'ps', 'ns'])
        line2 = group.addLine('Frame window: ',
                              condition='doDrop and not cutByTime',
                              help='First and last frame numbers (1-based; '
                                   '0 = use the trajectory limits).')
        line2.addParam('firstFrame', params.IntParam, label='Start: ', default=1)
        line2.addParam('lastFrame',  params.IntParam, label='End: ',   default=0)

        # ── 2. Subsampling ────────────────────────────────────────────────
        group = form.addGroup('Subsampling')
        group.addParam('doSubsample', params.BooleanParam,
                       label='Subsample trajectory?: ', default=False,
                       help='Keep only every N-th frame from the (possibly cut) trajectory.')
        group.addParam('subsampleF', params.IntParam,
                       label='Keep every N frames: ', default=10,
                       condition='doSubsample',
                       help='Frame stride. E.g. 10 = keep 1 frame out of 10.')

        # ── 3. Imaging ────────────────────────────────────────────────────
        group = form.addGroup('Imaging (PBC fix)')
        group.addParam('doAutoimage', params.BooleanParam,
                       label='Autoimage trajectory?: ', default=True,
                       help='Re-image molecules back into the primary unit cell using the '
                            'cpptraj *autoimage* command. Recommended as the first action '
                            'for solvated systems. Requires box information in the trajectory.')
        group.addParam('autoimageAnchor', params.StringParam,
                       label='Anchor mask: ', default='',
                       condition='doAutoimage',
                       help='Amber mask for the anchor molecule (the component kept centred). '
                            'Leave blank to let cpptraj choose automatically. Example: ``:1-300``')

        # ── 4. Stripping ──────────────────────────────────────────────────
        group = form.addGroup('Stripping')
        group.addParam('doStrip', params.BooleanParam,
                       label='Strip atoms?: ', default=False,
                       help='Remove unwanted atoms (solvent, ions, ...). '
                            'A new stripped topology (.prmtop) is written automatically so '
                            'the output system stays self-consistent.')
        group.addParam('stripSelection', params.EnumParam,
                       label='Strip selection: ', default=STRIP_SOLVENT,
                       condition='doStrip',
                       choices=['Solvent + ions',
                                'Full environment (waters, ions, lipids)',
                                'Custom mask'],
                       help='Predefined or custom Amber mask of atoms to remove.\n'
                            '- Solvent + ions: waters and monatomic ions.\n'
                            '- Full environment: everything in the plugin ENV_RES list '
                            '(waters, ions and lipids).')
        group.addParam('stripMaskCustom', params.StringParam,
                       label='Custom strip mask: ', default=':WAT',
                       condition='doStrip and stripSelection=={}'.format(STRIP_CUSTOM),
                       help='Amber mask of residues/atoms to strip, e.g. ``:WAT,Na+,Cl-``')

        # ── 5. Fitting ────────────────────────────────────────────────────
        group = form.addGroup('Fitting')
        group.addParam('doFit', params.BooleanParam,
                       label='RMS-fit trajectory?: ', default=False,
                       help='Superpose every frame onto the first frame of the (processed) '
                            'trajectory using the cpptraj *rms* command. Removes global '
                            'rotation/translation. Applied after stripping, so the fit is '
                            'computed on the retained atoms only.')
        group.addParam('fitMaskType', params.EnumParam,
                       label='Fit atom selection: ', default=FIT_MASK_BACKBONE,
                       condition='doFit',
                       choices=['Backbone (@N,CA,C,O)',
                                'C-alpha only (@CA)',
                                'All heavy atoms (!@H=)',
                                'Custom mask'],
                       help='Atom selection used for the least-squares fit.')
        group.addParam('fitMaskCustom', params.StringParam,
                       label='Custom fit mask: ', default='@CA',
                       condition='doFit and fitMaskType=={}'.format(FIT_MASK_CUSTOM),
                       help='Amber atom mask string, e.g. ``:1-250@CA``')

        # ── 6. Running average ─────────────────────────────────────────────
        group = form.addGroup('Running average')
        group.addParam('doRunAvg', params.BooleanParam,
                       label='Running average?: ', default=False,
                       help='Apply a coordinate running average with the cpptraj *runavg* '
                            'action. Smooths high-frequency motion; useful before visual '
                            'inspection. Reduces the number of output frames.')
        group.addParam('runAvgWindow', params.IntParam,
                       label='Window (frames): ', default=5,
                       condition='doRunAvg',
                       help='Number of frames averaged in each window (>= 2). '
                            'Larger values = stronger smoothing.')

        # ── Output format ─────────────────────────────────────────────────
        group = form.addGroup('Output')
        group.addParam('outputFormat', params.EnumParam,
                       label='Output trajectory format: ', default=OUT_FMT_NC,
                       choices=['NetCDF (.nc)', 'Charmm DCD (.dcd)', 'Gromacs XTC (.xtc)',
                                'Amber ASCII (.crd)', 'Multi-model PDB (.pdb)'],
                       help='Format for the output trajectory file.')

    # ── STEPS ────────────────────────────────────────────────────────────────
    def _insertAllSteps(self):
        self._insertFunctionStep('modifySystemStep')
        self._insertFunctionStep('createOutputStep')

    # ── main processing step ─────────────────────────────────────────────────
    def modifySystemStep(self):
        inSystem = self.amberSystem.get()
        topFile = os.path.abspath(inSystem.getTopologyFile())
        trajSource, hasTraj = self._getTrajSource()

        lines = ['parm {}'.format(topFile)]

        # --- input trajectory / coordinates (with cut + subsample) ----------
        trajinLine = 'trajin {}'.format(os.path.abspath(trajSource))
        if hasTraj:
            frameSpec = self._getTrajinFrameSpec()
            if frameSpec:
                trajinLine += ' ' + frameSpec
        lines.append(trajinLine)

        # --- actions (order matters in cpptraj) -----------------------------
        # 1. autoimage first, so subsequent actions see whole molecules
        if self.doAutoimage and hasTraj:
            aiLine = 'autoimage'
            anchor = self.autoimageAnchor.get().strip()
            if anchor:
                aiLine += ' anchor {}'.format(anchor)
            lines.append(aiLine)

        # 2. strip before fitting (do not fit onto solvent) and write new topology
        if self.doStrip:
            lines.append('strip {} parmout {}'.format(
                self._getStripMask(), os.path.abspath(self.getCleanTopologyFile())))

        # 3. RMS fit onto the first (processed) frame
        if self.doFit and hasTraj:
            lines.append('rms first {}'.format(self._getFitMask()))

        # 4. running average of coordinates
        if self.doRunAvg and hasTraj:
            lines.append('runavg window {}'.format(self.runAvgWindow.get()))

        # --- output --------------------------------------------------------
        if hasTraj:
            lines.append('trajout {} {}'.format(
                os.path.abspath(self.getCleanTrajectoryFile()), self._getOutputFormatKw()))
        else:
            # no trajectory: only a processed structure is produced
            lines.append('trajout {} pdb'.format(
                os.path.abspath(self.getCleanStructureFile())))

        lines += ['run', 'quit']

        cpptrajIn = os.path.abspath(self._getPath('cpptraj_modify.in'))
        with open(cpptrajIn, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

        amberPlugin.runAmbertools(self, program='cpptraj',
                                  args='-i {}'.format(cpptrajIn),
                                  cwd=self._getPath())

        # If no stripping was done, the output topology equals the input one.
        if not self.doStrip:
            shutil.copy(topFile, self.getCleanTopologyFile())

        # Extract a single-frame structure (consistent with the output topology) from
        # the processed trajectory. Done in a separate pass so it does not interact
        # with frame-buffering actions such as runavg.
        if hasTraj:
            self.writeSystemStructure()

    def writeSystemStructure(self):
        """Write the first frame of the processed trajectory to a PDB structure file
        that matches the (possibly stripped) output topology."""
        lines = ['parm {}'.format(os.path.abspath(self.getCleanTopologyFile())),
                 'trajin {} 1 1'.format(os.path.abspath(self.getCleanTrajectoryFile())),
                 'trajout {} pdb'.format(os.path.abspath(self.getCleanStructureFile())),
                 'run', 'quit']
        structIn = os.path.abspath(self._getPath('cpptraj_structure.in'))
        with open(structIn, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        amberPlugin.runAmbertools(self, program='cpptraj',
                                  args='-i {}'.format(structIn),
                                  cwd=self._getPath())

    # ── output step ──────────────────────────────────────────────────────────
    def createOutputStep(self):
        inSystem = self.amberSystem.get()
        _, hasTraj = self._getTrajSource()

        outSystem = AmberSystem()
        outSystem.setSystemFile(os.path.abspath(self.getCleanStructureFile()))
        outSystem.setTopologyFile(os.path.abspath(self.getCleanTopologyFile()))

        # carry over force-field / ligand metadata (unchanged by this protocol)
        if inSystem.getForceField():
            outSystem.setForceField(inSystem.getForceField())
        if inSystem.getWaterForceField():
            outSystem.setWaterForceField(inSystem.getWaterForceField())
        if inSystem.getLigandID():
            outSystem.setLigandID(inSystem.getLigandID())
        if inSystem.getLigTopologyFile():
            outSystem.setLigTopologyFile(inSystem.getLigTopologyFile())
        # the original coordinates only stay valid if no atoms were stripped
        if not self.doStrip and inSystem.getCrdFile():
            outSystem.setCrdFile(inSystem.getCrdFile())

        if hasTraj:
            outSystem.setTrajectoryFile(os.path.abspath(self.getCleanTrajectoryFile()))
            outSystem.readTrjInfo(protocol=self, nTime=self._estimateOutTime(),
                                  outDir=self._getExtraPath())

        self._defineOutputs(outputSystem=outSystem)

    # ── validation / info ──────────────────────────────────────────────────────
    def _validate(self):
        errs = []
        inSystem = self.amberSystem.get()

        if self.doRunAvg and self.runAvgWindow.get() < 2:
            errs.append('The running-average window must be >= 2 frames.')
        if self.doSubsample and self.subsampleF.get() < 2:
            errs.append('The subsampling stride must be >= 2.')
        if self.doDrop and not self.cutByTime:
            f0, f1 = self.firstFrame.get(), self.lastFrame.get()
            if f1 != 0 and f1 <= f0:
                errs.append('The last frame must be greater than the first frame '
                            '(set Last = 0 to use the end of the trajectory).')
        if self.doDrop and self.cutByTime and inSystem is not None:
            if not inSystem.getNFrames() or not inSystem.getNTime():
                errs.append('Time-based cutting needs the input trajectory metadata '
                            '(frames / time), which is not available. Use frame-based '
                            'cutting instead.')
        if self.doFit and self.fitMaskType.get() == FIT_MASK_CUSTOM \
                and not self.fitMaskCustom.get().strip():
            errs.append('A custom fit mask must be provided.')
        if self.doStrip and self.stripSelection.get() == STRIP_CUSTOM \
                and not self.stripMaskCustom.get().strip():
            errs.append('A custom strip mask must be provided.')
        return errs

    def _warnings(self):
        warns = []
        if not self.amberSystem.get().hasTrajectory():
            for flag, label in [(self.doFit,       'fitting'),
                                 (self.doDrop,      'cutting'),
                                 (self.doSubsample, 'subsampling'),
                                 (self.doRunAvg,    'running average')]:
                if flag:
                    warns.append('The input system has no trajectory - {} will be '
                                 'skipped (only the coordinates are processed).'.format(label))
        return warns

    def _summary(self):
        summary = []
        ops = []
        if self.doDrop:      ops.append('cut trajectory')
        if self.doSubsample: ops.append('subsample (1/{})'.format(self.subsampleF.get()))
        if self.doAutoimage: ops.append('autoimage')
        if self.doStrip:     ops.append('strip ({})'.format(self._getStripMask()))
        if self.doFit:       ops.append('RMS fit ({})'.format(self._getFitMask()))
        if self.doRunAvg:    ops.append('running average (window={})'.format(self.runAvgWindow.get()))
        if ops:
            summary.append('Applied operations: ' + ', '.join(ops))
        return summary

    def _methods(self):
        return ['The MD system trajectory was processed with cpptraj (AmberTools).']

    # ── helpers ──────────────────────────────────────────────────────────────
    def _getTrajSource(self):
        """Return (path, hasTrajectory). Falls back to coordinates/structure when the
        input system has no trajectory."""
        inSystem = self.amberSystem.get()
        trj = inSystem.getTrajectoryFile()
        if trj:
            return trj, True
        return inSystem.getCrdFile() or inSystem.getSystemFile(), False

    def _getTrajinFrameSpec(self):
        """Build the ``[start [stop [offset]]]`` part of the trajin command from the
        cutting and subsampling options. Returns '' when the full trajectory is used."""
        start, stop, offset = 1, None, None

        if self.doDrop:
            if self.cutByTime:
                psPerFrame = self._psPerFrame()
                t0 = self._timeToPs(self.firstTime.get(), self.getEnumText('timeUnit'))
                t1 = self._timeToPs(self.lastTime.get(),  self.getEnumText('timeUnit'))
                if self.firstTime.get() > 0 and psPerFrame:
                    start = max(1, int(round(t0 / psPerFrame)))
                if self.lastTime.get() > 0 and psPerFrame:
                    stop = max(start, int(round(t1 / psPerFrame)))
            else:
                start = self.firstFrame.get() if self.firstFrame.get() > 0 else 1
                stop = self.lastFrame.get() if self.lastFrame.get() > 0 else None

        if self.doSubsample:
            offset = self.subsampleF.get()

        # trajin requires start & stop to be present before an offset can be given
        if stop is None and offset is None and start == 1:
            return ''
        parts = [str(start), 'last' if stop is None else str(stop)]
        if offset is not None:
            parts.append(str(offset))
        return ' '.join(parts)

    def _psPerFrame(self):
        inSystem = self.amberSystem.get()
        nFrames, nTime = inSystem.getNFrames(), inSystem.getNTime()
        if nFrames and nTime:
            return nTime / nFrames
        return None

    def _estimateOutTime(self):
        """Best-effort total time (ps) of the output trajectory, derived from the input
        per-frame spacing, the applied stride and the resulting frame count."""
        psPerFrame = self._psPerFrame()
        if not psPerFrame:
            return self.amberSystem.get().getNTime() or 0.0
        stride = self.subsampleF.get() if self.doSubsample else 1
        outFrames = self._countFrames(os.path.abspath(self.getCleanTrajectoryFile()))
        if outFrames:
            return psPerFrame * stride * outFrames
        return self.amberSystem.get().getNTime() or 0.0

    def _countFrames(self, trjFile):
        """Return the number of frames in trjFile using ``cpptraj -tl``, or None."""
        import re
        topFile = os.path.abspath(self.getCleanTopologyFile())
        amberPlugin.runAmbertools(self, program='cpptraj',
                                  args='-p {} -y {} -tl'.format(topFile, trjFile))
        for logName in ('run.stdout', 'run.stderr'):
            candidate = self._getPath('logs', logName)
            if os.path.exists(candidate):
                with open(candidate) as fh:
                    for line in fh:
                        m = re.search(r'Frames:\s*(\d+)', line)
                        if m:
                            return int(m.group(1))
        return None

    def getCleanTopologyFile(self):
        top = self.amberSystem.get().getTopologyFile()
        base = os.path.splitext(os.path.basename(top))[0]
        return os.path.abspath(self._getPath(base + '_modified.prmtop'))

    def getCleanTrajectoryFile(self):
        extMap = {OUT_FMT_NC:  '.nc',  OUT_FMT_DCD: '.dcd', OUT_FMT_XTC: '.xtc',
                  OUT_FMT_CRD: '.crd', OUT_FMT_PDB: '.pdb'}
        trj = self.amberSystem.get().getTrajectoryFile()
        base = os.path.splitext(os.path.basename(trj))[0] if trj else 'trajectory'
        ext = extMap.get(self.outputFormat.get(), '.nc')
        return os.path.abspath(self._getPath(base + '_modified' + ext))

    def getCleanStructureFile(self):
        sysFile = self.amberSystem.get().getSystemFile()
        base = os.path.splitext(os.path.basename(sysFile))[0] if sysFile else 'system'
        return os.path.abspath(self._getPath(base + '_modified.pdb'))

    def _getStripMask(self):
        sel = self.stripSelection.get()
        if sel == STRIP_SOLVENT:
            return ':{},{}'.format(WATER_RES, ION_RES)
        elif sel == STRIP_ENV:
            return ':{}'.format(ENV_RES)
        return self.stripMaskCustom.get().strip()

    def _getFitMask(self):
        sel = self.fitMaskType.get()
        if sel == FIT_MASK_BACKBONE:
            return '@N,CA,C,O'
        elif sel == FIT_MASK_CA:
            return '@CA'
        elif sel == FIT_MASK_ALL:
            return '!@H='
        return self.fitMaskCustom.get().strip()

    def _getOutputFormatKw(self):
        kws = {OUT_FMT_NC:  'netcdf', OUT_FMT_DCD: 'dcd', OUT_FMT_XTC: 'xtc',
               OUT_FMT_CRD: 'crd',    OUT_FMT_PDB: 'pdb'}
        return kws.get(self.outputFormat.get(), 'netcdf')

    @staticmethod
    def _timeToPs(value, unitStr):
        """Convert a time value to picoseconds."""
        factors = {'fs': 1e-3, 'ps': 1.0, 'ns': 1e3}
        return value * factors.get(unitStr, 1.0)
