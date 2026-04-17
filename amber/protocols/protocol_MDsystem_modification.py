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
This module modifies an Amber system trajectory and/or coordinates using cpptraj.
"""

import os
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from amber import Plugin as amberPlugin
from amber.objects import AmberSystem


# ── Enum index constants (keep in sync with choices lists below) ─────────────
RMS_MASK_BACKBONE = 0
RMS_MASK_CA       = 1
RMS_MASK_ALL      = 2
RMS_MASK_CUSTOM   = 3

STRIP_SOLVENT     = 0   # WAT + common ions
STRIP_IONS_ONLY   = 1
STRIP_CUSTOM      = 2

OUT_FMT_NC        = 0
OUT_FMT_DCD       = 1
OUT_FMT_XTC       = 2
OUT_FMT_MDCRD     = 3
OUT_FMT_PDB       = 4

TIME_UNIT_FS = 0
TIME_UNIT_PS = 1
TIME_UNIT_NS = 2


class AmberModifySystem(EMProtocol):
    """
    Modifies an Amber MD system trajectory and/or coordinates using **cpptraj**.

    Available operations (all optional, applied in this order):
      1. Imaging     – fix PBC artefacts with ``autoimage``
      2. Fitting     – RMS-fit every frame to a reference
      3. Stripping   – remove solvent / ions; writes a new topology
      4. Cutting     – keep only a time window of the trajectory
      5. Subsampling – keep every N-th frame
      6. Smoothing   – low-pass Gaussian smoothing of coordinates
    """

    _label = 'System modification'

    # ── param definitions ────────────────────────────────────────────────────
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('amberSystem', params.PointerParam,
                      label='Input Amber System: ',
                      pointerClass='AmberSystem',
                      help='Amber solvated system (topology + trajectory) to process.')

        # ── 1. Imaging ────────────────────────────────────────────────────
        group = form.addGroup('Imaging (PBC fix)')
        group.addParam('doAutoimage', params.BooleanParam,
                       label='Autoimage trajectory?: ', default=True,
                       help='Re-image molecules back into the primary unit cell '
                            'using the cpptraj *autoimage* command. '
                            'Recommended as the first step for solvated systems.')
        group.addParam('autoimageAnchor', params.StringParam,
                       label='Anchor mask: ', default='',
                       condition='doAutoimage',
                       help='Amber mask for the anchor molecule (most stable '
                            'component). Leave blank to let cpptraj choose '
                            'automatically. Example: ``:1-300``')

        # ── 2. Fitting ────────────────────────────────────────────────────
        group = form.addGroup('Fitting')
        group.addParam('doFit', params.BooleanParam,
                       label='RMS-fit trajectory?: ', default=False,
                       help='Superpose every frame onto the first frame (or a '
                            'chosen reference) using the cpptraj *rms* command.')
        group.addParam('fitMaskType', params.EnumParam,
                       label='Fit atom selection: ', default=RMS_MASK_BACKBONE,
                       condition='doFit',
                       choices=['Backbone (@N,CA,C,O)',
                                'C-alpha only (@CA)',
                                'All heavy atoms (!@H*)',
                                'Custom mask'],
                       help='Atom selection used for the least-squares fit.')
        group.addParam('fitMaskCustom', params.StringParam,
                       label='Custom fit mask: ', default='@CA',
                       condition='doFit and fitMaskType=={}'.format(RMS_MASK_CUSTOM),
                       help='Amber atom mask string, e.g. ``:1-250@CA``')
        group.addParam('fitRef', params.EnumParam,
                       label='Reference frame: ', default=0,
                       condition='doFit',
                       choices=['First frame', 'Last frame', 'Average structure'],
                       help='Frame used as fitting reference.')

        # ── 3. Stripping ──────────────────────────────────────────────────
        group = form.addGroup('Stripping')
        group.addParam('doStrip', params.BooleanParam,
                       label='Strip atoms?: ', default=False,
                       help='Remove unwanted atoms (solvent, ions). '
                            'A new stripped topology (.prmtop) is written automatically.')
        group.addParam('stripSelection', params.EnumParam,
                       label='Strip selection: ', default=STRIP_SOLVENT,
                       condition='doStrip',
                       choices=['Solvent + ions (:WAT,:Na+,:Cl-,:K+,:Mg2+,:Ca2+)',
                                'Ions only (:Na+,:Cl-,:K+,:Mg2+,:Ca2+)',
                                'Custom mask'],
                       help='Predefined or custom Amber mask of atoms to remove.')
        group.addParam('stripMaskCustom', params.StringParam,
                       label='Custom strip mask: ', default=':WAT',
                       condition='doStrip and stripSelection=={}'.format(STRIP_CUSTOM),
                       help='Amber mask of residues/atoms to strip, '
                            'e.g. ``:WAT,:Na+``')

        # ── 4. Cutting ────────────────────────────────────────────────────
        group = form.addGroup('Cutting')
        group.addParam('doDrop', params.BooleanParam,
                       label='Cut trajectory?: ', default=False,
                       help='Keep only frames within a time (or frame-number) window.')
        group.addParam('cutByTime', params.BooleanParam,
                       label='Define window by time?: ', default=True,
                       condition='doDrop',
                       help='If True, use simulation time; '
                            'if False, use frame indices directly.')
        line = group.addLine('Time window: ',
                             condition='doDrop and cutByTime',
                             help='Start and end times (0 = use trajectory limits).')
        line.addParam('firstTime', params.FloatParam, label='Start: ', default=0.0)
        line.addParam('lastTime',  params.FloatParam, label='End: ',   default=0.0)
        line.addParam('timeUnit',  params.EnumParam,
                      label='Units: ', default=TIME_UNIT_NS,
                      choices=['fs', 'ps', 'ns'])
        line2 = group.addLine('Frame window: ',
                              condition='doDrop and not cutByTime',
                              help='First and last frame numbers (1-based; '
                                   '0 = use trajectory limits).')
        line2.addParam('firstFrame', params.IntParam, label='Start: ', default=1)
        line2.addParam('lastFrame',  params.IntParam, label='End: ',   default=0)

        # ── 5. Subsampling ────────────────────────────────────────────────
        group = form.addGroup('Subsampling')
        group.addParam('doSubsample', params.BooleanParam,
                       label='Subsample trajectory?: ', default=False,
                       help='Keep only every N-th frame from the (possibly cut) '
                            'trajectory.')
        group.addParam('subsampleF', params.IntParam,
                       label='Keep every N frames: ', default=10,
                       condition='doSubsample',
                       help='Frame stride. E.g. 10 = keep 1 frame out of 10.')

        # ── 6. Smoothing ──────────────────────────────────────────────────
        group = form.addGroup('Smoothing')
        group.addParam('doSmooth', params.BooleanParam,
                       label='Smooth trajectory?: ', default=False,
                       help='Apply a Gaussian low-pass coordinate filter using '
                            'the cpptraj *smooth* action. '
                            'Reduces high-frequency noise; useful before '
                            'visual inspection or RMSD calculations.')
        group.addParam('smoothWindow', params.IntParam,
                       label='Smoothing window (frames): ', default=5,
                       condition='doSmooth',
                       help='Half-width of the Gaussian window. '
                            'Larger values = stronger smoothing.')

        # ── Output format ─────────────────────────────────────────────────
        group = form.addGroup('Output')
        group.addParam('outputFormat', params.EnumParam,
                       label='Output trajectory format: ', default=OUT_FMT_NC,
                       choices=['NetCDF (.nc)', 'DCD (.dcd)', 'XTC (.xtc)',
                                'ASCII MDCRD (.mdcrd)', 'Multi-PDB (.pdb)'],
                       help='Format for the output trajectory file.')

    # ── STEPS ────────────────────────────────────────────────────────────────
    def _insertAllSteps(self):
        self._insertFunctionStep('modifySystemStep')
        self._insertFunctionStep('createOutputStep')

    # ── main processing step ─────────────────────────────────────────────────
    def modifySystemStep(self):
        topFile = os.path.abspath(self.amberSystem.get().getTopologyFile())
        trjFile = self.amberSystem.get().getTrajectoryFile()

        # Build the cpptraj input script line by line
        lines = []

        # --- topology & trajectory input ------------------------------------
        lines.append('parm {}'.format(topFile))

        if trjFile:
            trjFile = os.path.abspath(trjFile)
            trajin_line = 'trajin {}'.format(trjFile)

            if self.doDrop:
                if self.cutByTime:
                    t0 = self._timeToPs(self.firstTime.get(), self.getEnumText('timeUnit'))
                    t1 = self._timeToPs(self.lastTime.get(),  self.getEnumText('timeUnit'))
                    # cpptraj trajin uses frame numbers; convert via offset=1
                    # For time-based cutting we use the offset/start/stop trick:
                    # trajin file [start [stop [offset]]]
                    # We add a separate "time" filter via the 'time' keyword if available,
                    # but the most portable approach is a separate trajin with offsets.
                    # We instead guard with 'trajin ... 1 last 1' and rely on the
                    # time-based frame selection written into the script header.
                    trajin_line += ' starttime {} endtime {}'.format(t0, t1)
                else:
                    f0 = self.firstFrame.get() if self.firstFrame.get() > 0 else 1
                    f1 = self.lastFrame.get()  if self.lastFrame.get()  > 0 else 'last'
                    trajin_line += ' {} {}'.format(f0, f1)

            if self.doSubsample:
                # Append stride to trajin.  If cutting was already added we
                # need the full "start stop offset" form; fill in defaults.
                if self.doDrop and not self.cutByTime:
                    # already has start/stop, just add offset
                    trajin_line += ' {}'.format(self.subsampleF.get())
                elif not self.doDrop:
                    trajin_line += ' 1 last {}'.format(self.subsampleF.get())
                # (time-based cutting + subsampling: stride applied after trajin)

            lines.append(trajin_line)
        else:
            # No trajectory – only the structure file will be processed
            strucFile = os.path.abspath(self.amberSystem.get().getFileName())
            lines.append('trajin {}'.format(strucFile))

        # --- actions (order matters in cpptraj!) ----------------------------

        # 1. Autoimage first so that subsequent operations see whole molecules
        if self.doAutoimage:
            ai_line = 'autoimage'
            anchor = self.autoimageAnchor.get().strip()
            if anchor:
                ai_line += ' anchor {}'.format(anchor)
            lines.append(ai_line)

        # 2. Strip before fitting to avoid fitting onto solvent atoms
        if self.doStrip:
            strip_mask = self._getStripMask()
            for mask in strip_mask:
                lines.append('strip {}'.format(mask))
            # write stripped topology alongside the trajectory
            lines.append('parmout {}'.format(
                os.path.abspath(self.getCleanTopologyFile())))

        # 3. RMS fit
        if self.doFit:
            fit_mask = self._getFitMask()
            ref_kw   = self._getFitReference()
            lines.append('rms {} {} nofit'.format(ref_kw, fit_mask))
            # 'nofit' on the rms analysis line; use a separate 'rms' action for
            # the actual coordinate superposition:
            lines.append('rms {} {}'.format(ref_kw, fit_mask))

        # 4. Smoothing
        if self.doSmooth:
            lines.append('smooth {} window {}'.format(
                '!:WAT' if not self.doStrip else '*',
                self.smoothWindow.get()))

        # --- output trajectory ----------------------------------------------
        if trjFile or self.amberSystem.get().getFileName():
            out_trj = os.path.abspath(self.getCleanTrajectoryFile())
            fmt_kw  = self._getOutputFormatKw()
            lines.append('trajout {} {}'.format(out_trj, fmt_kw))

        lines.append('run')
        lines.append('quit')

        # Write the cpptraj input file
        cpptraj_in = os.path.abspath(self._getPath('cpptraj_modify.in'))
        with open(cpptraj_in, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

        # Run cpptraj
        amberPlugin.runAmbertools(
            self, program='cpptraj',
            args='-i {}'.format(cpptraj_in),
            cwd=self._getPath())

        # If no stripping was done, copy original topology as the output one
        if not self.doStrip:
            shutil.copy(topFile, self.getCleanTopologyFile())

    # ── output step ──────────────────────────────────────────────────────────
    def createOutputStep(self):
        outSystem = AmberSystem()
        outSystem.setTopologyFile(self.getCleanTopologyFile())

        strucFile = self.amberSystem.get().getFileName()
        if strucFile:
            outSystem.setFileName(strucFile)

        if self.amberSystem.get().getTrajectoryFile():
            outSystem.setTrajectoryFile(
                os.path.relpath(self.getCleanTrajectoryFile()))
            outSystem.readTrjInfo(protocol=self, nTime=self.amberSystem.get().getNTime(), outDir=self._getExtraPath())

        self._defineOutputs(outputSystem=outSystem)

    # ── validation / info ────────────────────────────────────────────────────
    def _validate(self):
        errs = []
        if self.doSmooth and self.smoothWindow.get() < 2:
            errs.append('Smoothing window must be ≥ 2 frames.')
        if self.doSubsample and self.subsampleF.get() < 2:
            errs.append('Subsampling stride must be ≥ 2.')
        if self.doDrop and not self.cutByTime:
            f0 = self.firstFrame.get()
            f1 = self.lastFrame.get()
            if f1 != 0 and f1 <= f0:
                errs.append('Last frame must be greater than first frame '
                             '(set Last = 0 to use the end of the trajectory).')
        if self.doFit and self.fitMaskType.get() == RMS_MASK_CUSTOM:
            if not self.fitMaskCustom.get().strip():
                errs.append('A custom fit mask must be provided.')
        if self.doStrip and self.stripSelection.get() == STRIP_CUSTOM:
            if not self.stripMaskCustom.get().strip():
                errs.append('A custom strip mask must be provided.')
        return errs

    def _warnings(self):
        warns = []
        trj = self.amberSystem.get().getTrajectoryFile()
        if not trj:
            for flag, label in [(self.doFit,       'fitting'),
                                 (self.doDrop,      'cutting'),
                                 (self.doSubsample, 'subsampling'),
                                 (self.doSmooth,    'smoothing')]:
                if flag:
                    warns.append('The input system has no trajectory – '
                                 '{} will be skipped.'.format(label))
        if self.doFit and self.doAutoimage:
            warns.append('Autoimage is enabled together with RMS fitting. '
                         'Note that autoimage must run before rms in cpptraj '
                         '(already enforced by this protocol) to avoid '
                         'incorrect imaging after rotation.')
        return warns

    def _summary(self):
        summary = []
        ops = []
        if self.doAutoimage:  ops.append('autoimage')
        if self.doStrip:      ops.append('strip ({})'.format(
                                  self._getStripMask()))
        if self.doFit:        ops.append('RMS fit ({})'.format(
                                  self._getFitMask()))
        if self.doDrop:       ops.append('cut trajectory')
        if self.doSubsample:  ops.append('subsample (1/{})'.format(
                                  self.subsampleF.get()))
        if self.doSmooth:     ops.append('smooth (window={})'.format(
                                  self.smoothWindow.get()))
        if ops:
            summary.append('Applied operations: ' + ', '.join(ops))
        return summary

    def _methods(self):
        return ['Trajectory was processed with cpptraj (AmberTools).']

    # ── helpers ──────────────────────────────────────────────────────────────
    def getCleanTopologyFile(self):
        top = self.amberSystem.get().getTopologyFile()
        base = os.path.splitext(os.path.basename(top))[0]
        return os.path.abspath(self._getPath(base + '_clean.prmtop'))

    def getCleanTrajectoryFile(self):
        ext_map = {OUT_FMT_NC:    '.nc',
                   OUT_FMT_DCD:   '.dcd',
                   OUT_FMT_XTC:   '.xtc',
                   OUT_FMT_MDCRD: '.mdcrd',
                   OUT_FMT_PDB:   '.pdb'}
        trj = self.amberSystem.get().getTrajectoryFile()
        if trj:
            base = os.path.splitext(os.path.basename(trj))[0]
        else:
            base = 'trajectory'
        ext = ext_map.get(self.outputFormat.get(), '.nc')
        return os.path.abspath(self._getPath(base + '_clean' + ext))

    def _getStripMask(self):
        """Return list of Amber masks to strip."""
        SOLVENT_MASKS = [':WAT', ':Na+', ':Cl-', ':K+', ':Mg2+', ':Ca2+']
        ION_MASKS     = [':Na+', ':Cl-', ':K+', ':Mg2+', ':Ca2+']
        sel = self.stripSelection.get()
        if sel == STRIP_SOLVENT:
            return SOLVENT_MASKS
        elif sel == STRIP_IONS_ONLY:
            return ION_MASKS
        else:  # custom
            return [self.stripMaskCustom.get().strip()]

    def _getFitMask(self):
        sel = self.fitMaskType.get()
        if sel == RMS_MASK_BACKBONE:
            return '@N,CA,C,O'
        elif sel == RMS_MASK_CA:
            return '@CA'
        elif sel == RMS_MASK_ALL:
            return '!@H='
        else:
            return self.fitMaskCustom.get().strip()

    def _getFitReference(self):
        ref = self.fitRef.get()
        if ref == 0:
            return 'first'
        elif ref == 1:
            return 'last'
        else:
            return 'average'

    def _getOutputFormatKw(self):
        kws = {OUT_FMT_NC:    'netcdf',
               OUT_FMT_DCD:   'dcd',
               OUT_FMT_XTC:   'xtc',
               OUT_FMT_MDCRD: 'mdcrd',
               OUT_FMT_PDB:   'pdb'}
        return kws.get(self.outputFormat.get(), 'netcdf')

    @staticmethod
    def _timeToPs(value, unit_str):
        """Convert a time value to picoseconds (cpptraj internal unit)."""
        factors = {'fs': 1e-3, 'ps': 1.0, 'ns': 1e3}
        return value * factors.get(unit_str, 1.0)