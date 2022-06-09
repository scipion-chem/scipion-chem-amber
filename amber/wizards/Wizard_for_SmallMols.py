# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Alba Lomas
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

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pyworkflow.wizard import Wizard

import amber
from pwem.wizards import GetStructureChainsWizard, EmWizard
from pwem.convert import AtomicStructHandler
import pyworkflow.object as pwobj
import pyworkflow.wizard as pwizard
import os, requests
import pwchem.protocols as emprot
from amber.protocols import AmberSystemPrep


class SelectFileNameWizard(pwizard.Wizard):
    """Base wizard for selecting an attribute from those contained in the items of a given input
  inputParam: Name of the input parameter where the items are stored
  outputParam:
  """
    _targets, _inputs, _outputs = [], {}, {}

    def addTarget(self, protocol, targets, inputs, outputs):
        self._targets += [(protocol, targets)]
        self._inputs.update({protocol: inputs})
        self._outputs.update({protocol: outputs})

    def getInputOutput(self, form):
        '''Retrieving input and output corresponding to the protocol where the wizard is used'''
        outParam = ''
        for target in self._targets:
            if form.protocol.__class__ == target[0]:
                inParam, outParam = self._inputs[target[0]], self._outputs[target[0]]
        return inParam, outParam

    def getInputFiles(self, form, inputParam):

        inputPointer = getattr(form.protocol, inputParam)

        # fileNames = []

        # for mol in inputPointer:
        # mol = mol.get()
        # fileNames.append(os.path.abspath(mol.getFileName()))

        return inputPointer.get()

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        molList = self.getInputFiles(form, inputParam[0])
        finalFilesList = []

        for i in molList:
            finalFilesList.append(pwobj.String(i.getUniqueName()))

        provider = ListTreeProviderString(finalFilesList)

        dlg = dialog.ListDialog(form.root, "Filter set", provider,
                                "Select one of the molecules")

        form.setVar(outputParam[0], dlg.values[0].get())


# Defining target for the SelectAttributeWizard
SelectFileNameWizard().addTarget(protocol=amber.protocols.AmberSystemPrep,
                                 targets=['inputLigandSelect'],
                                 inputs=['inputLigand'],
                                 outputs=['inputLigandSelect'])


