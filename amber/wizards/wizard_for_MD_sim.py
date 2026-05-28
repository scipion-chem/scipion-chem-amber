# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
import tkinter as tk
from tkinter import messagebox
from pyworkflow.gui import ListTreeProviderString, dialog
from pwem.objects import Pointer, String

from pwchem.wizards import DeleteElementWizard, VariableWizard, WatchElementWizard
from pwchem.utils import cifFromASFile

from ..protocols import AmberMDSimulation, AmberSystemPrep
from ..constants import *

class AmberAddElementSummaryWizard(VariableWizard):
    """Add a step of the workflow in the defined position"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        numSteps = protocol.countSteps()

        if 'min' in form.wizParamName:
            stageType = 'Minimization'
        elif 'heat' in form.wizParamName:
            stageType = 'Heating'
        elif 'sim' in form.wizParamName:
            stageType = 'Simulation'
        elif 'custom' in form.wizParamName:
            stageType = 'Custom'

        if getattr(protocol, inputParam[0]).get().strip() != '':
            index = int(getattr(protocol, inputParam[0]).get())
        else:
            index = numSteps + 1

        msjDic = protocol.getStageParamsDic(stageType)

        if index > numSteps:
            prevStr = getattr(protocol, outputParam[0]).get() \
                if getattr(protocol, outputParam[0]).get() is not None else ''
            form.setVar(outputParam[0], prevStr + str(msjDic) + '\n')

            newSum = protocol.createSummary()
            form.setVar(outputParam[1], newSum)

        elif numSteps >= index > 0:
            workSteps = getattr(protocol, outputParam[0]).get().split('\n')
            workSteps.insert(index - 1, str(msjDic))
            form.setVar(outputParam[0], '\n'.join(workSteps))

            newSum = protocol.createSummary()
            form.setVar(outputParam[1], newSum)

AmberAddElementSummaryWizard().addTarget(protocol=AmberMDSimulation,
                                         targets=['minInsertStep'],
                                         inputs=['minInsertStep'],
                                         outputs=['workFlowSteps', 'summarySteps'])

AmberAddElementSummaryWizard().addTarget(protocol=AmberMDSimulation,
                                         targets=['heatInsertStep'],
                                         inputs=['heatInsertStep'],
                                         outputs=['workFlowSteps', 'summarySteps'])

AmberAddElementSummaryWizard().addTarget(protocol=AmberMDSimulation,
                                         targets=['simInsertStep'],
                                         inputs=['simInsertStep'],
                                         outputs=['workFlowSteps', 'summarySteps'])

AmberAddElementSummaryWizard().addTarget(protocol=AmberMDSimulation,
                                         targets=['customInsertStep'],
                                         inputs=['customInsertStep'],
                                         outputs=['workFlowSteps', 'summarySteps'])

DeleteElementWizard().addTarget(protocol=AmberMDSimulation,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])



class AmberAddDefaultWorkflow(VariableWizard):
    """Add a default workflow based on system type"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        # Determine workflow type from wizard parameter name
        if 'protein' in form.wizParamName:
            workflowSteps = PROTWORK
        elif 'protLig' in form.wizParamName:
            workflowSteps = PROTLIGWORK
        elif 'memLig' in form.wizParamName:
            workflowSteps = MEMPROTLIGWORK
        elif 'membrane' in form.wizParamName:
            workflowSteps = MEMPROTWORK
        else:
            workflowSteps = PROTWORK  # Default fallback

        # Set the workflow steps
        form.setVar(outputParam[0], workflowSteps)

        # Generate and set summary using the unified function
        newSum = protocol.createSummary(workflowSteps)
        form.setVar(outputParam[1], newSum)

AmberAddDefaultWorkflow().addTarget(protocol=AmberMDSimulation,
                                    targets=['proteinDefault'],
                                    inputs=[''],
                                    outputs=['workFlowSteps', 'summarySteps'])

AmberAddDefaultWorkflow().addTarget(protocol=AmberMDSimulation,
                                    targets=['protLigDefault'],
                                    inputs=[''],
                                    outputs=['workFlowSteps', 'summarySteps'])

AmberAddDefaultWorkflow().addTarget(protocol=AmberMDSimulation,
                                    targets=['membraneDefault'],
                                    inputs=[''],
                                    outputs=['workFlowSteps', 'summarySteps'])

AmberAddDefaultWorkflow().addTarget(protocol=AmberMDSimulation,
                                    targets=['memLigDefault'],
                                    inputs=[''],
                                    outputs=['workFlowSteps', 'summarySteps'])


class DisulfideBondWizard(VariableWizard):
    """Wizard to select multiple pairs of CYS residues for disulfide bond creation"""
    _targets, _inputs, _outputs = [], {}, {}

    def parseCysPairs(self, inputFile, protocol):
        """Parse PDB/mmCIF file and return CYS pairs whose SG atoms are 1.5–2.5 Å apart
        (typical disulfide bond distance ~2.05 Å).

        Returns: [((chainA, resIdxA, resNameA), (chainB, resIdxB, resNameB)), ...]
        Standardizes any input to mmCIF format using pwem utilities before parsing.
        """
        # key: (chainId, resIdx) -> (x, y, z)
        sgCoords = {}

        base, _ = os.path.splitext(os.path.basename(inputFile))
        tmpCifPath = os.path.abspath(protocol.getProject().getTmpPath(f'{base}_temp.cif'))

        cifFile = cifFromASFile(inputFile, tmpCifPath)
        if not cifFile or not os.path.exists(cifFile):
            print(f"ERROR: Conversion failed or mmCIF file missing: {cifFile}")
            return []

        with open(cifFile, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                try:
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    if parts[5].strip() != 'CYS':
                        continue
                    if parts[3].strip() != 'SG':  # only the sulfur atom
                        continue

                    chainId = parts[6].strip()
                    if not chainId or chainId == '.':
                        chainId = '_'

                    resNumStr = parts[8].strip()
                    if not resNumStr or not resNumStr.replace('-', '').isdigit():
                        continue

                    resIdx = int(resNumStr)
                    x, y, z = float(parts[10]), float(parts[11]), float(parts[12])

                    sgCoords[(chainId, resIdx)] = (x, y, z)

                except (ValueError, IndexError) as e:
                    print(f"WARNING: Could not parse line: {line.rstrip()}\n  Error: {e}")
                    continue

        if cifFile == tmpCifPath and os.path.exists(tmpCifPath):
            try:
                os.remove(tmpCifPath)
            except OSError:
                pass

        # Pairwise SG–SG distance filter (1.5–2.5 Å = disulfide bond range)
        pairs = []
        keys = list(sgCoords.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                kA, kB = keys[i], keys[j]
                xA, yA, zA = sgCoords[kA]
                xB, yB, zB = sgCoords[kB]
                dist = ((xA - xB) ** 2 + (yA - yB) ** 2 + (zA - zB) ** 2) ** 0.5
                if 1.5 <= dist <= 2.5:
                    chainA, resIdxA = kA
                    chainB, resIdxB = kB
                    pairs.append((
                        (chainA, resIdxA, 'CYS'),
                        (chainB, resIdxB, 'CYS'),
                    ))

        return pairs

    def _pairLabel(self, pair):
        (chainA, idxA, nameA), (chainB, idxB, nameB) = pair
        return f"{nameA} {chainA}:{idxA} \u2194 {nameB} {chainB}:{idxB}"

    def show(self, form, *params):
        """Main wizard entry point"""
        inputParams, outputParams = self.getInputOutput(form)
        protocol = form.protocol

        inputObj = getattr(protocol, inputParams[0]).get()

        if hasattr(inputObj, 'getFileName'):
            pdbFile = inputObj.getFileName()
        else:
            dialog.showError("Missing Input",
                         "Please select a valid input structure first.",
                               form.root)

        cysPairs = self.parseCysPairs(pdbFile, protocol)

        if not cysPairs:
            dialog.showInfo("Disulfide Bonds",
                            "No CYS pairs with SG distance in [1.5, 2.5] Å found.",
                            form.root)
            return

        pairLabels = [self._pairLabel(p) for p in cysPairs]
        labelStrings = [String(lbl) for lbl in pairLabels]

        provider = ListTreeProviderString(labelStrings)
        dlg = dialog.ListDialog(
            form.root, "Select Disulfide Bonds", provider,
            "Select which disulfide bonds to form\n"
            "(Ctrl+Click or Shift+Click for multiple):"
        )

        if dlg.values:
            selectedLabels = {val.get() for val in dlg.values}

            # Encode each selected pair as "chainA_resIdxA-chain_resIdxB"
            selectedTokens = []
            for pair, label in zip(cysPairs, pairLabels):
                if label in selectedLabels:
                    (chainA, idxA, _), (chainB, idxB, _) = pair
                    selectedTokens.append(f"{chainA}_{idxA}-{chainB}_{idxB}")

            if selectedTokens:
                form.setVar(outputParams[0], '/'.join(selectedTokens))

DisulfideBondWizard().addTarget(
    protocol=AmberSystemPrep,
    targets=['disulfideBridgesNumber'],
    inputs=['inputStructure'],
    outputs=['disulfideBridgesNumber']
)