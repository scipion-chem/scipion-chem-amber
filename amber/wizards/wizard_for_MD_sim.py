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
import os
import tkinter as tk
from tkinter import messagebox
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

import tkinter as tk
import tkinter.ttk as ttk
from pyworkflow.gui.tree import Tree

import tkinter as tk
from tkinter import ttk, messagebox


class DisulfideBondWizard(VariableWizard):
    """Wizard to select multiple pairs of CYS residues for disulfide bond creation"""
    _targets, _inputs, _outputs = [], {}, {}

    def parseCysResidues(self, inputFile, protocol):
        """Parse PDB/mmCIF file and return CYS residues grouped by chain
        Returns: {chainName: [(resIdx, resName), ...], ...}
        Standardizes any input to mmCIF format using pwem utilities before parsing.
        """
        cysDict = {}

        # 1. FIX: Use .cif extension for the temporary file instead of .pdb
        base, _ = os.path.splitext(os.path.basename(inputFile))
        tmpCifPath = os.path.abspath(protocol.getProject().getTmpPath(f'{base}_temp.cif'))

        # Convert/copy input to mmCIF
        cifFile = cifFromASFile(inputFile, tmpCifPath)

        if not cifFile or not os.path.exists(cifFile):
            print(f"ERROR: Conversion failed or mmCIF file missing: {cifFile}")
            return {}

        with open(cifFile, 'r') as f:
            for line in f:
                # Standard mmCIF coordinate rows start with ATOM or HETATM tokens
                if not line.startswith(('ATOM', 'HETATM')):
                    continue

                try:
                    parts = line.split()
                    if len(parts) < 10:
                        continue

                    # 2. FIX: Extract the variables from standard mmCIF columns
                    # parts[5] = Residue Name (e.g., CYS)
                    # parts[6] = Chain ID (e.g., A)
                    # parts[8] = Sequence Index / Residue Number (e.g., 14)

                    resName = parts[5].strip()
                    if resName != 'CYS':
                        continue

                    chainId = parts[6].strip()
                    if not chainId or chainId == '.':
                        chainId = '_'

                    resNumStr = parts[8].strip()
                    if not resNumStr or not resNumStr.replace('-', '').isdigit():
                        continue
                    resIdx = int(resNumStr)

                    # Add to dictionary
                    if chainId not in cysDict:
                        cysDict[chainId] = set()

                    cysDict[chainId].add((resIdx, resName))

                except (ValueError, IndexError) as e:
                    print(f"WARNING: Could not parse line: {line.rstrip()}\n  Error: {e}")
                    continue

        # 3. FIX: Safe cleanup. ONLY delete if it's the temporary file we created
        if cifFile == tmpCifPath and os.path.exists(tmpCifPath):
            try:
                os.remove(tmpCifPath)
            except OSError:
                pass

        # Sort values sequentially by residue index per chain
        for chain in cysDict:
            cysDict[chain] = sorted(list(cysDict[chain]), key=lambda x: x[0])

        return cysDict

    def createSelectionDialog(self, cysDict, protocol):
        """Create dialog with two tables for CYS selection and a list for multiple bonds"""

        dialog = tk.Toplevel()
        dialog.title("Select Disulfide Bond Residues")
        dialog.geometry("850x750")  # Slightly taller to accommodate the new listbox

        mainFrame = tk.Frame(dialog)
        mainFrame.pack(fill='both', expand=True, padx=10, pady=10)

        instrLabel = tk.Label(mainFrame,
                              text="1. Select one CYS from each table.\n2. Click 'Add Bond' to save the pair.\n3. Click OK when finished.",
                              font=('Arial', 10, 'bold'), justify='left')
        instrLabel.pack(pady=(0, 10), anchor='w')

        # --- TABLES FRAME ---
        tablesFrame = tk.Frame(mainFrame)
        tablesFrame.pack(fill='both', expand=True)

        # Left panel
        leftFrame = tk.Frame(tablesFrame)
        leftFrame.pack(side='left', fill='both', expand=True, padx=(0, 5))
        tk.Label(leftFrame, text="First CYS Residue", font=('Arial', 9, 'bold')).pack()

        self.leftTree = ttk.Treeview(leftFrame, columns=('Chain', 'Residue', 'Name'), show='headings',
                                     selectmode='browse')
        for col in ('Chain', 'Residue', 'Name'):
            self.leftTree.heading(col, text=col)
            self.leftTree.column(col, width=80, anchor='center')

        leftScrollbar = ttk.Scrollbar(leftFrame, orient='vertical', command=self.leftTree.yview)
        self.leftTree.configure(yscrollcommand=leftScrollbar.set)
        self.leftTree.pack(side='left', fill='both', expand=True)
        leftScrollbar.pack(side='right', fill='y')

        # Right panel
        rightFrame = tk.Frame(tablesFrame)
        rightFrame.pack(side='right', fill='both', expand=True, padx=(5, 0))
        tk.Label(rightFrame, text="Second CYS Residue", font=('Arial', 9, 'bold')).pack()

        self.rightTree = ttk.Treeview(rightFrame, columns=('Chain', 'Residue', 'Name'), show='headings',
                                      selectmode='browse')
        for col in ('Chain', 'Residue', 'Name'):
            self.rightTree.heading(col, text=col)
            self.rightTree.column(col, width=80, anchor='center')

        rightScrollbar = ttk.Scrollbar(rightFrame, orient='vertical', command=self.rightTree.yview)
        self.rightTree.configure(yscrollcommand=rightScrollbar.set)
        self.rightTree.pack(side='left', fill='both', expand=True)
        rightScrollbar.pack(side='right', fill='y')

        # Populate trees
        for chain in sorted(cysDict.keys()):
            for resIdx, resName in cysDict[chain]:
                values = (chain, resIdx, resName)
                self.leftTree.insert('', 'end', values=values)
                self.rightTree.insert('', 'end', values=values)

        # --- ADD BUTTON & SELECTION DISPLAY ---
        selectionFrame = tk.Frame(mainFrame)
        selectionFrame.pack(fill='x', pady=10)

        self.selectionLabel = tk.Label(selectionFrame, text="Current Selection: None - None", font=('Arial', 10))
        self.selectionLabel.pack()

        self.leftTree.bind('<<TreeviewSelect>>', self.updateSelection)
        self.rightTree.bind('<<TreeviewSelect>>', self.updateSelection)

        # The list to store formatted bond strings
        self.addedBonds = []

        def onAddBond():
            leftSelection = self.leftTree.selection()
            rightSelection = self.rightTree.selection()

            if leftSelection and rightSelection:
                leftItem = self.leftTree.item(leftSelection[0])
                rightItem = self.rightTree.item(rightSelection[0])

                bond = f"{leftItem['values'][0]}_{leftItem['values'][1]}-{rightItem['values'][0]}_{rightItem['values'][1]}"
                reverse_bond = f"{rightItem['values'][0]}_{rightItem['values'][1]}-{leftItem['values'][0]}_{leftItem['values'][1]}"

                # Prevent duplicates (checking both directions A-B and B-A)
                if bond in self.addedBonds or reverse_bond in self.addedBonds:
                    messagebox.showinfo("Duplicate", "This disulfide bond has already been added.")
                elif leftItem['values'] == rightItem['values']:
                    messagebox.showwarning("Invalid Bond", "Cannot bond a residue to itself.")
                else:
                    self.addedBonds.append(bond)
                    self.bondsListbox.insert(tk.END, bond)
            else:
                messagebox.showwarning("Incomplete Selection", "Please select one residue from each table first.")

        addBtn = tk.Button(selectionFrame, text="↓ Add Selected Bond ↓", command=onAddBond, width=25, bg='#e0e0e0',
                           font=('Arial', 9, 'bold'))
        addBtn.pack(pady=5)

        # --- LISTBOX FOR ADDED BONDS ---
        listFrame = tk.Frame(mainFrame)
        listFrame.pack(fill='both', expand=True, pady=5)

        tk.Label(listFrame, text="Added Disulfide Bonds:", font=('Arial', 9, 'bold')).pack(
            anchor='w')

        self.bondsListbox = tk.Listbox(listFrame, height=6, font=('Courier', 10))
        self.bondsListbox.pack(side='left', fill='both', expand=True)

        listScroll = ttk.Scrollbar(listFrame, orient='vertical', command=self.bondsListbox.yview)
        self.bondsListbox.configure(yscrollcommand=listScroll.set)
        listScroll.pack(side='right', fill='y')

        def onRemoveBond():
            selection = self.bondsListbox.curselection()
            if selection:
                idx = selection[0]
                self.bondsListbox.delete(idx)
                self.addedBonds.pop(idx)

        removeBtn = tk.Button(mainFrame, text="Remove Selected from List", command=onRemoveBond)
        removeBtn.pack(anchor='e')

        # --- OK / CANCEL BUTTONS ---
        buttonFrame = tk.Frame(mainFrame)
        buttonFrame.pack(pady=15)

        self.finalBondString = None

        def onOk():
            if not self.addedBonds:
                if not messagebox.askyesno("No Bonds", "No bonds have been added. Proceed with empty selection?"):
                    return
                self.finalBondString = ""
            else:
                # This joins the list with a forward slash: A_25-A_29/B_25-B_54
                self.finalBondString = "/".join(self.addedBonds)
            dialog.destroy()

        def onCancel():
            self.finalBondString = None
            dialog.destroy()

        okButton = tk.Button(buttonFrame, text="OK", command=onOk, width=12, font=('Arial', 9, 'bold'))
        okButton.pack(side='left', padx=10)

        cancelButton = tk.Button(buttonFrame, text="Cancel", command=onCancel, width=12)
        cancelButton.pack(side='left', padx=10)

        dialog.transient()
        dialog.grab_set()
        dialog.wait_window()

        return self.finalBondString

    def updateSelection(self, event=None):
        """Update the label showing what is currently highlighted in the trees"""
        leftSelection = self.leftTree.selection()
        rightSelection = self.rightTree.selection()

        left_text = f"{self.leftTree.item(leftSelection[0])['values'][0]}_{self.leftTree.item(leftSelection[0])['values'][1]}" if leftSelection else "None"
        right_text = f"{self.rightTree.item(rightSelection[0])['values'][0]}_{self.rightTree.item(rightSelection[0])['values'][1]}" if rightSelection else "None"

        self.selectionLabel.config(text=f"Current Selection: {left_text} - {right_text}")

    def show(self, form, *params):
        """Main wizard entry point"""
        inputParams, outputParams = self.getInputOutput(form)
        protocol = form.protocol

        inputObj = getattr(protocol, inputParams[0]).get()

        if hasattr(inputObj, 'getFileName'):
            pdbFile = inputObj.getFileName()
        else:
            print("ERROR: Could not get PDB file from input object")
            return

        cysDict = self.parseCysResidues(pdbFile, protocol)

        if not cysDict:
            messagebox.showwarning("No CYS Found", "No cysteine residues found in the structure")
            return

        # Capture the joined string
        finalBondString = self.createSelectionDialog(cysDict, protocol)

        if finalBondString is not None:
            print(f"Selected disulfide bonds: {finalBondString}")
            if outputParams:
                form.setVar(outputParams[0], finalBondString)
                param = getattr(protocol, outputParams[0])
                param.set(finalBondString)
        else:
            print("Selection cancelled")

# Register the Wizard (Update these parameters to match your specific plugin structure)
DisulfideBondWizard().addTarget(
    protocol=AmberSystemPrep,  # Replace with your actual protocol class
    targets=['disulfideBridgesNumber'],  # The UI element this wizard is attached to
    inputs=['inputStructure'],  # The parameter holding the PDB file
    outputs=['disulfideBridgesNumber']  # Where the string result will be saved
)