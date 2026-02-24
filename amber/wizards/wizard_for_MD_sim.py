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

# Imports
from pwchem.wizards import DeleteElementWizard, VariableWizard, WatchElementWizard

from ..protocols import AmberMDSimulation
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

# WatchElementWizard().addTarget(protocol=AmberMDSimulation,
#                                 targets=['watchStep'],
#                                 inputs=['watchStep'],
#                                 outputs=['workFlowSteps', 'summarySteps'])


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