# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Aida Pinacho Pérez
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
This module will perform energy minimizations and equilibrium for the system befor MD simultion
"""
import os, glob, shutil
from os.path import relpath

from pyworkflow.mapper.sqlite import SELF
from pyworkflow.protocol import params
from pyworkflow.utils import Message, runJob, createLink

import amber
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.utils import natural_sort

from amber.objects import *
from amber.constants import *
from amber import Plugin as amberPlugin


class AmberMDSimulation(EMProtocol):
    """
    This protocol will perform energy minimization heating and similation on the system previously prepared by the protocol
    "system prepartion".
    """
    _amberEngines = ['sander', 'pmemd']
    _label = 'run MD simulation'
    _ensemTypes = ['NVT', 'NPT']
    _thermostats = ['Andersen', 'Langevin', 'Nose-Hoover', 'Nose-Hoover RESPA', 'Berendsen']
    _barostats = ['Berendsen', 'Monte Carlo']

    _coupleStyle = ['No pressure scaling', 'isotropic', 'anisotropic', 'semiisotropic']

    _shakeAlgorithm = ['Shake not performed', 'Bonds involving hydrogens are constrains', 'all bonds are constrained']

    _omitParamNames = ['amberSystem', 'runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                       'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi', 'minInsertStep', 'heatInsertStep',
                       'simInsertStep', 'customInsertStep']

    _key_map = {'Minimization': 'min', 'Heating': 'heat', 'Simulation': 'sim', 'Custom': 'custom'}

    _restrained_groups = list(RESTRAINS_DIC.keys())

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="This protocol has both CPU and GPU implementation.\
                                                         Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addSection('Input')
        form.addParam('amberSystem', params.PointerParam, label="Input Amber System: ",
                      pointerClass='AmberSystem',
                      allowsNull=True,
                      help='Amber solvated system to be simulated')
        group = form.addGroup('Minimization')
        group.addParam('energyMin', params.BooleanParam,
                       label='Energy Minimization: ', default=True,
                       help='Energy minimization before the simulation (recommended)')
        line = group.addLine('Minimization settings: ', condition='energyMin',
                             help='The first 20 cycles will utilize the steepest descent'
                                  'algorithm before shifting to the conjugate gradient '
                                  'algorithm for the remaining cycles\nThe first x cycles will utilize the steepest descent'
                                  'algorithm before shifting to the conjugate gradient '
                                  'algorithm for the remaining cycles\n'
                                  'Sphere of influence for each atom during energy minimization')
        line.addParam('minMaxCycles', params.IntParam, default=1000,
                      label='Maximum cycles:', condition='energyMin')
        line.addParam('minSdCycles', params.IntParam, default=500,
                      label='Steepest Descent cycles:', condition='energyMin')
        line.addParam('minIntCutoff', params.FloatParam, default=8.0,
                      label='Interaction cutoff', condition='energyMin')

        group.addParam('minRestraint', params.BooleanParam, default=False,
                       label='Add restrains',
                       help='Restraining specified atoms in Cartesian space using a harmonic potential')
        lineMin = group.addLine('Restrains in minimization: ', condition='minRestraint',
                                help="Specify the components of the system to be restraint and the associated force constant.")
        lineMin.addParam('minRestrAtoms', params.EnumParam, choices=self._restrained_groups, default=0,
                         label='Atoms to restrain')
        lineMin.addParam('minRestrForce', params.FloatParam, default=50.0,
                         label='Force (kcal·mol-1·Å-2)')
        group.addParam('minInsertStep', params.StringParam, default='',
                       label='Insert Minimization step number: ',
                       help='Insert the defined Minimization step into the workflow on the defined position (number).\n'
                            'The default (when empty) is the last position')
        group = form.addGroup('Heating - NVT')
        line = group.addLine('Heating simulation time: ',
                             help='Time settings\n'
                                  'Number of MD steps to run '
                                  'Time step in ps. (Number of MD steps * Time step = run length in ps)'
                                  'Trajectory step size: The trajectory coordinates are written to a traj file every x steps.')
        line.addParam('heatMDSteps', params.IntParam, default=10000,
                      label='Number of MD steps:',
                      help='Number of MD steps in run (x * time step = run length in ps)')
        line.addParam('heatTimeStep', params.FloatParam, default=0.002,
                      label='Time step (ps)')
        group.addParam('heatSaveTrj', params.BooleanParam, default=False,
                       label='Save trajectory: ',
                       help='Save trajectory coordinates during the heating stage.')
        group.addParam('heatTraj', params.IntParam, default=1000,
                       label='Trajectory step size', condition='heatSaveTrj',
                       help='The coordinates are written to a mdcrd file every x steps.')

        line = group.addLine('Temperature increase: ',
                             help='Initial and final temperature (K)')
        line.addParam('heatInTemp', params.FloatParam, default=0,
                      label='Initial temperature (K)')
        line.addParam('heatFiTemp', params.FloatParam, default=300,
                      label='Final temperature (K)')

        line = group.addLine('Heating temperature control: ',
                             help='Thermostat type and associated parameters')
        line.addParam('heatThermostat', params.EnumParam, default=1,
                      label='Thermostat: ', choices=self._thermostats)
        line.addParam('heatCollisFreq', params.FloatParam, default=2.0,
                      label='Collision frequency (1/ps): ', condition='heatThermostat==1')
        line.addParam('heatCoupConst', params.FloatParam, default=2.0,
                      label='Coupling constant (1/ps): ', condition='heatThermostat==2')
        line.addParam('heatFricConst', params.FloatParam, default=2.0,
                      label='Friction constant (1/ps): ', condition='heatThermostat==3')

        group.addParam('heatRestraint', params.BooleanParam, default=False,
                       label='Add restrains',
                       help='Restraining specified atoms in Cartesian space using a harmonic potential')
        line = group.addLine('Restrains in heating: ', condition='heatRestraint',
                             help="Specify the components of the system to be restraint and the associated force constant.")
        line.addParam('heatRestrAtoms', params.EnumParam, default=3, choices=self._restrained_groups,
                      label='Atoms to restrain')
        line.addParam('heatRestrForce', params.FloatParam, default=50.0,
                      label='Force (kcal·mol-1·Å-2)')

        group.addParam('heatInsertStep', params.StringParam, default='',
                       label='Insert Heating step number: ',
                       help='Insert the defined Heating step into the workflow on the defined position (number).\n'
                            'The default (when empty) is the last position')

        group = form.addGroup('Simulation - NVT or NPT')
        line = group.addLine('Simulation time: ',
                             help='Time setting\n'
                                  'Number of MD steps to run '
                                  'Time step in ps. (Number of MD steps * Time step = run length in ps)'
                                  'Trajectory step size: The trajectory coordinates are written to a traj file every x steps.')
        line.addParam('simMDSteps', params.IntParam, default=100000,
                      label='Number of MD steps:',
                      help='Number of MD steps in run (nstlim * dt = run length in ps)')
        line.addParam('simTimeStep', params.FloatParam, default=0.002,
                      label='Time step (ps)')
        group.addParam('simSaveTrj', params.BooleanParam, default=True,
                       label='Save trajectory: ',
                       help='Save trajectory coordinates during the simulation stage.')
        group.addParam('simTrajStep', params.IntParam, default=1000,
                       label='Trajectory step size', condition='simSaveTrj',
                       help='The trajectory coordinates are written to a traj file every x steps.')

        group.addParam('simEnsemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=1,
                       help='Type of simulation to perform: NVT or NPT\n')
        line = group.addLine('Simulation temperature control: ',
                             help='Thermostat type and associated params')
        line.addParam('simThermostat', params.EnumParam, default=1,
                      label='Thermostat: ', choices=self._thermostats)
        line.addParam('simCollisFreq', params.FloatParam, default=2.0,
                      label='Collision frequency (1/ps): ', condition='simThermostat==1')
        line.addParam('simCoupConst', params.FloatParam, default=2.0,
                      label='Coupling constant (1/ps): ', condition='simThermostat==2')
        line.addParam('simFricConst', params.FloatParam, default=2.0,
                      label='Friction constant (1/ps): ', condition='simThermostat==3')

        line = group.addLine('Pressure control: ', condition='simEnsemType==1',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Pressure scaling: Recomemded Isotropic in all cases but membrane'
                                  'systems, where semiisotropic is recommended.')
        line.addParam('simPressure', params.FloatParam, default=1.0, condition='simEnsemType==1',
                      label='Pressure (bar): ')
        line.addParam('simBarostat', params.EnumParam, default=1, condition='simEnsemType==1',
                      label='Barostat type: ', choices=self._barostats)
        line.addParam('simPressureScaling', params.EnumParam, default=1, condition='simEnsemType==1',
                      label='Pressure scaling: ', choices=self._coupleStyle)

        group.addParam('simRestraint', params.BooleanParam, default=False,
                       label='Add restrains',
                       help='Restraining specified atoms in Cartesian space using a harmonic potential')
        lineSim = group.addLine('Restrains in simulation: ', condition='simRestraint',
                                help="Specify the components of the system to be restraint and the associated force constant.")
        lineSim.addParam('simRestrAtoms', params.EnumParam, default=3, choices=self._restrained_groups,
                         label='Atoms to restrain')
        lineSim.addParam('simRestrForce', params.FloatParam, default=50.0,
                         label='Force (kcal·mol-1·Å-2)')

        group.addParam('simInsertStep', params.StringParam, default='',
                       label='Insert Simulation step number: ',
                       help='Insert the defined Simulation step into the workflow on the defined position (number).\n'
                            'The default (when empty) is the last position')

        group = form.addGroup('Custom input', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('customIn', params.TextParam, width=60, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED,
                       help=f'Upload a custom configuration file a sander/pmemd MD step.\n'
                            'For detailed syntax and options, refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')
        group.addParam('customInsertStep', params.StringParam, default='',
                       label='Insert Simulation step number: ',
                       help='Insert the custom step into the workflow on the defined position (number).\n'
                            'The default (when empty) is the last position')

        group = form.addGroup('Summary')
        group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
                       label='Summary of steps',
                       help='Summary of the defined steps. \nManual modification will have no '
                            'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete step number: ',
                       help='Delete the step of the specified index from the workflow.')
        # group.addParam('watchStep', params.StringParam, default='',
        #                label='Watch relaxation step number: ',
        #                help='''Watch the parameters step of the specified index from the workflow..\n
        #                                This might be useful if you want to change some parameters of a predefined step.\n
        #                                However, the parameters are not changed until you add the new step (and probably\n
        #                                you may want to delete the previous unchanged step)''')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

        form.addSection('Example workflows')
        group = form.addGroup('Protein')
        group.addParam('proteinDefault', params.LabelParam, label='Protein default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system'
                            'Summary of steps is updated')
        group = form.addGroup('Protein+Ligand')
        group.addParam('protLigDefault', params.LabelParam, label='Protein+Ligand default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system'
                            'Summary of steps is updated')
        group = form.addGroup('Transmembrane protein')
        group.addParam('membraneDefault', params.LabelParam, label='Transmembrane default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system'
                            'Summary of steps is updated')
        group = form.addGroup('Transmembrane protein+Ligand')
        group.addParam('memLigDefault', params.LabelParam, label='Transmembrane+Ligand default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system'
                            'Summary of steps is updated')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        print('Printing each step specific params:')
        print(self.workFlowSteps.get())
        self.createGUISummary()
        i = 1
        for wStep in self.workFlowSteps.get().strip().split('\n'):
            self._insertFunctionStep(self.simulateStageStep, wStep, i)
            i += 1

        self._insertFunctionStep(self.createOutputStep)

    def simulateStageStep(self, wStep, i):
        msjDic = eval(wStep)
        mdpFile = self.generateMDPFile(msjDic, str(i))
        self.callAmber(mdpFile, saveTrj=self.shouldSaveTrj(msjDic))

    def createOutputStep(self):
        lastCrdFile, lastTopoFile, lastOutFile = self.getPrevFinishedStageFiles()
        oriSystemFile = self.amberSystem.get().getSystemFile()
        ligTopFile = self.amberSystem.get().getLigTopologyFile()
        ligID = self.amberSystem.get().getLigandID()

        localCrdFile, localTopFile = self._getPath('outputSystem.rst7'), self._getPath('systemTopology.parm7')
        shutil.copyfile(lastCrdFile, localCrdFile), shutil.copyfile(lastTopoFile, localTopFile)

        mFF, wFF = self.getFFFiles()

        outSystem = AmberSystem(filename=relpath(oriSystemFile), ff=mFF, wff=wFF, nTime=self.calculateTotalSimTime(),
                                ligTopFile=ligTopFile, ligName=ligID)

        outSystem.setTopologyFile(localTopFile)
        outSystem.setCrdFile(localCrdFile)

        concatTrjFile = self.prepareSimTrj()
        if concatTrjFile is not None:
            outputTrajectory = self._getPath('outputTrajectory.nc')
            shutil.copyfile(concatTrjFile, outputTrajectory)
            outSystem.setTrajectoryFile(outputTrajectory)
            outSystem.readTrjInfo(protocol=self, nTime=self.calculateSavedTrjTime(),
                                  outDir=self._getExtraPath())

        finalPdbFile = self.crdToPDB(localCrdFile, localTopFile)
        finalAtomStruct = AtomStruct(filename=relpath(finalPdbFile))
        self._defineOutputs(outputSystem=outSystem, lastFrameStruct=finalAtomStruct)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        fnSummary = self._getExtraPath("summary.txt")
        if not os.path.exists(fnSummary):
            summary = ["No summary information yet."]
        else:
            fhSummary = open(fnSummary, "r")
            summary = []
            for line in fhSummary.readlines():
                summary.append(line.rstrip())
            fhSummary.close()
        return summary

    def createSummary(self, workSteps=None):
        """Creates the displayed summary from workflow steps string."""
        if workSteps is None:
            workSteps = self.workFlowSteps.get()
        if not workSteps or not workSteps.strip():
            return ''

        sumStr = ''
        lastTemp = 300
        for i, dicLine in enumerate(workSteps.split('\n')):
            if not dicLine.strip():
                continue
            msjDic = self.addDefaultForMissing(eval(dicLine))
            stepType = msjDic.get('stepType', 'Step')
            lineText = f'{i + 1}) {stepType} - '

            if stepType == 'Minimization':
                lineText += f"Max Cycles: {msjDic.get('MaxCycles', 0)}"
                if msjDic.get('Restraint'):
                    lineText += f", restraint on {msjDic.get('RestrAtoms')}"

            elif stepType == 'Heating':
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lastTemp = msjDic.get('FiTemp', 300)
                lineText += f"Sim. time: {nTime} ps, NVT ensemble, {msjDic.get('InTemp', 0)} K to {lastTemp} K"
                if not self.shouldSaveTrj(msjDic):
                    lineText += ', trajectory not saved'
                if msjDic.get('Restraint'):
                    lineText += f", restraint on {msjDic.get('RestrAtoms')}"

            elif stepType == 'Simulation':  # Simulation
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lineText += f"Sim. time: {nTime} ps, {msjDic.get('EnsemType', 'NPT')} ensemble, {lastTemp} K"
                if not self.shouldSaveTrj(msjDic):
                    lineText += ', trajectory not saved'
                if msjDic.get('Restraint'):
                    lineText += f", restraint on {msjDic.get('RestrAtoms')}"

            elif stepType == 'Custom':
                # Show the first non-empty line of the .in file as a hint
                firstLine = next(
                    (l.strip() for l in msjDic.get('In', '').splitlines() if l.strip()), 'custom input'
                )
                lineText += f'Custom input: "{firstLine}"'

            sumStr += lineText + '\n'
        return sumStr

    def createGUISummary(self):
        with open(self._getExtraPath("summary.txt"), 'w') as f:
            f.write(self.createSummary())

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append('The methods used to perform the ')

        return methods

    ######################## UTILS ##################################
    def countSteps(self):
        """Count steps from workFlowSteps"""
        workStepsStr = self.workFlowSteps.get() if self.workFlowSteps.get() is not None else ''
        if workStepsStr.strip() == '':
            return 0
        steps = [step for step in workStepsStr.split('\n') if step.strip() != '']
        return len(steps)

    def getStageParamsDic(self, stageType):
        """Return {cleanParamName: value} for parameters belonging to the given stage type."""
        prefix = self._key_map.get(stageType, '')
        paramsDic = {}
        paramsDic['stepType'] = stageType
        for paramName, param in self._definition.iterAllParams():
            if paramName in self._omitParamNames:
                continue
            if isinstance(param, (params.Group, params.Line)):
                continue
            if prefix and paramName.startswith(prefix):
                cleanName = paramName[len(prefix):]
                if isinstance(param, params.EnumParam):
                    paramsDic[cleanName] = self.getEnumText(paramName)
                else:
                    paramsDic[cleanName] = getattr(self, paramName).get()
        return paramsDic

    def addDefaultForMissing(self, msjDic):
        """Add default values for any parameters missing from msjDic."""
        if msjDic.get('stepType') == 'Custom':
            return msjDic
        for paramName, param in self._definition.iterAllParams():
            if paramName in self._omitParamNames:
                continue
            if isinstance(param, (params.Group, params.Line)):
                continue
            if paramName not in msjDic:
                msjDic[paramName] = param.default
        return msjDic

    def customMDPFile(self, msjDic, type):
        stageDir = self._getExtraPath(type)
        os.makedirs(stageDir, exist_ok=True)

        key = f"{self._key_map.get(type)}CustomIn"

        # Get content, clean spaces, and ensure exactly one blank line at the end
        params = msjDic.get(key, "").replace('\xa0', ' ').rstrip() + '\n\n'

        mdpFile = os.path.join(stageDir, f"{type}.in")
        with open(mdpFile, 'w') as f:
            f.write(params)

        return mdpFile

    def shouldSaveTrj(self, msjDic):
        """Return whether a workflow step should write a trajectory file.

        Old workflow dictionaries do not contain SaveTrj. Treat them as True
        to preserve the previous Amber protocol behavior.
        """
        return msjDic.get('SaveTrj', True)

    def generateMDPFile(self, msjDic, i):
        '''Generate .in file'''
        stepType = msjDic['stepType']
        stageDir = self._getExtraPath('{}_{}'.format(i, stepType))
        if os.path.exists(stageDir):
            shutil.rmtree(stageDir)
        os.makedirs(stageDir)

        mdpFile = os.path.join(stageDir, '{}_{}.in'.format(i, stepType))

        params = ''
        if stepType == 'Custom':
            params = msjDic['In'].rstrip() + '\n'

        if stepType == 'Minimization':
            params = 'MINIMIZATION\n&cntrl \n' \
                     'imin=1, ntx=1, irest=0, maxcyc={}, ncyc={}, ntpr=100,' \
                     ' ntwx=0, cut={}'.format(msjDic['MaxCycles'], msjDic['SdCycles'], msjDic['IntCutoff'])
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' /".format(msjDic['RestrForce'],
                                                                                  RESTRAINS_DIC[msjDic['RestrAtoms']])

        elif stepType == 'Heating':
            params = 'HEATING\n&cntrl \n' \
                     'imin=0, nstlim={}, dt={}, ntf=2, ntc=2, tempi={}, ' \
                     'temp0={}, ntpr={} , ntwx={}, ntb=1, ntp=0, ig=-1, ' \
                     'cut=8.0 '.format(msjDic['MDSteps'],
                                       msjDic['TimeStep'],
                                       msjDic['InTemp'],
                                       msjDic['FiTemp'],
                                       msjDic['Traj'],
                                       msjDic['Traj'] if self.shouldSaveTrj(msjDic) else 0)
            params += self.addThermostatParams(msjDic)
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' ".format(msjDic['RestrForce'],
                                                                                 RESTRAINS_DIC[msjDic['RestrAtoms']])

            params += '/\n&wt type=\'TEMP0\', istep1=0, istep2={}, value1={}, value2={} /\n'.format(
                msjDic['MDSteps'],
                msjDic['InTemp'],
                msjDic['FiTemp'],
                msjDic['FiTemp'])

        elif stepType == 'Simulation':
            params = 'MD SIMULATION\n&cntrl \n' \
                     'imin=0, ntx=5, irest=1, nstlim={}, dt={}, ntf=2, ntc=2, ' \
                     'temp0={}, ntpr={} , ntwx={}, ig=-1, ' \
                     'cut=8.0'.format(msjDic['MDSteps'],
                                      msjDic['TimeStep'],
                                      self.getLastHeatingTemp(),
                                      msjDic['TrajStep'],
                                      msjDic['TrajStep'] if self.shouldSaveTrj(msjDic) else 0)
            params += self.addThermostatParams(msjDic)
            if msjDic['EnsemType'] == 'NPT':
                params += self.addBarostatParams(msjDic)
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' /".format(msjDic['RestrForce'],
                                                                                  RESTRAINS_DIC[msjDic['RestrAtoms']])

        if stepType != 'Custom':
            params += '\n&end \nEND'
        with open(mdpFile, 'w') as f:
            f.write(params)

        print('sander/pmemd.cuda input params are:\n {}'.format(params))

        return mdpFile

    def addThermostatParams(self, msjDic):
        ntt = None
        gamma_ln = None
        thermostat = msjDic['Thermostat']
        if thermostat == 'Andersen':
            ntt = 2
        elif thermostat == 'Langevin':
            ntt = 3
            gamma_ln = msjDic.get('CollisFreq')
        elif thermostat == 'Nose-Hoove':
            ntt = 9
            gamma_ln = msjDic.get('CoupConst')
        elif thermostat == 'Nose-Hoover RESPA':
            ntt = 10
            gamma_ln = msjDic.get('FricConst')
        elif thermostat == 'Berendsen':
            ntt = 11

        params = []
        if ntt is not None:
            params.append(f"\nntt={ntt}")
        if gamma_ln is not None:
            params.append(f"gamma_ln={gamma_ln}")

        return ", ".join(params)

    def addBarostatParams(self, msjDic):
        params = []
        extraMembParams = ""

        if msjDic['EnsemType'] == 'NPT':
            pres0 = msjDic['Pressure']
            ntb = 2

            pressScaParam = msjDic['PressureScaling']
            barosParam = msjDic['Barostat']

            if pressScaParam == 'isotropic':
                ntp = 1
            elif pressScaParam == 'anisotropic':
                ntp = 2
            elif pressScaParam == 'semiisotropic':
                # used for membranes
                ntp = 3
                extraMembParams = ", csurften=3, gamma_ten=0.0"
            else:
                ntp = 0

            if barosParam == 'Berendsen':
                barostat = 1
            elif barosParam == 'Monte Carlo':
                barostat = 2
            else:
                barostat = 0

            params.append(
                f"\nntb={ntb}, ntp={ntp}, pres0={pres0}, barostat={barostat}{extraMembParams}"
            )

        else:
            ntp = 0
            ntb = 1
            params.append(f"\nntb={ntb}, ntp={ntp}")

        return ", ".join(params)

    def callAmber(self, mdpFile, saveTrj=True):
        inputFile = os.path.abspath(mdpFile)
        stageDir = os.path.dirname(mdpFile)
        stage = os.path.split(stageDir)[-1]
        stageType = stage.split('_')[1]
        stageNum = stage.split('_')[0]
        outFile = os.path.join(stage, '{}.o'.format(stage, type))
        crdFile, topFile, _ = self.getPrevFinishedStageFiles(stageNum)

        if stageType == 'Minimization':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst -o {}.o' \
                      ' -inf {}.inf'.format(inputFile, crdFile, topFile, crdFile, *[stage] * 4)
        elif stageType == 'Heating':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst' \
                      ' -o {}.o '.format(inputFile, crdFile, topFile, crdFile, *[stage] * 2)
            if saveTrj:
                command += ' -x {}.netcdf'.format(stage)
            command += ' -inf {}.inf'.format(stage)

        elif stageType == 'Simulation':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst' \
                      ' -o {}.o'.format(inputFile, crdFile, topFile, crdFile, *[stage] * 2)
            if saveTrj:
                command += ' -x {}.netcdf'.format(stage)
            command += ' -inf {}.inf'.format(stage)

        elif stageType == 'Custom':
            # The .in content determines what sander/pmemd actually runs;
            command = f'-i {inputFile} -c {crdFile} -p {topFile} -ref {crdFile} -r {stage}.ncrst -o {stage}.o -inf {stage}.inf' \
                      f' -x {stage}.netcdf'

        if self.useGpu.get():
            os.environ["CUDA_VISIBLE_DEVICES"] = self.gpuList.get()
            amberPlugin.runPmemd(self, ' -O ', args=command, cwd=stageDir)

        else:
            # Use sander for CPU execution
            engine = 'sander'
            amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)

        return os.path.join(stageDir, outFile)

    def getPrevFinishedStageFiles(self, stageNum=None):
        '''Return the previous .crd and topology files if number stage is provided.
        If not, returns the ones of the lastest stage'''
        topFile = self.amberSystem.get().getTopologyFile()
        if stageNum:
            if stageNum == '1':
                crdFile = os.path.abspath(self.amberSystem.get().getCrdFile())
                outFile = None

            else:
                prevDir = self.getStageDir((int(stageNum) - 1))
                for file in os.listdir(prevDir):
                    if '.ncrst' in file:
                        crdFile = os.path.join(prevDir, file)
                    elif '.o' in file:
                        outFile = os.path.join(prevDir, file)
        else:
            stageDir = self.getLastStageDir()
            for file in os.listdir(stageDir):
                if '.ncrst' in file:
                    crdFile = os.path.join(stageDir, file)
                elif '.o' in file:
                    outFile = os.path.join(stageDir, file)

        return os.path.abspath(crdFile), os.path.abspath(topFile), outFile

    def crdToPDB(self, crdFile, topFile):
        """Run cpptraj to strip waters and ions and save the coordinates to a PDB file.
        """
        systemName = os.path.splitext(os.path.basename(self.amberSystem.get().getSystemFile()))[0]
        finalPdbFile = self._getPath(f'{systemName}_final.pdb')

        cpptrajCmds = f"""parm {os.path.abspath(topFile)}
            trajin {os.path.abspath(crdFile)}
            strip :{ENV_RES}
            trajout {os.path.abspath(finalPdbFile)} pdb
            run
            quit
            """
        scriptPath = self._getExtraPath('strip_solvent_cpptraj.in')
        with open(scriptPath, 'w') as f:
            f.write(cpptrajCmds)

        cmd = f'-i {os.path.abspath(scriptPath)}'
        amberPlugin.runAmbertools(self, 'cpptraj', args=cmd, cwd=self._getExtraPath())

        return finalPdbFile

    def getFFFiles(self):
        system = self.amberSystem.get()
        return system.getForceField(), system.getWaterForceField()

    def getStageDir(self, stage):
        """Returns the directory path for a given stage number"""
        pattern = self._getExtraPath('{}*'.format(stage))
        matchingDirs = glob.glob(pattern)

        if matchingDirs: return os.path.abspath(matchingDirs[0])

        # 3. Si no existe, retornamos None o podrías lanzar un error
        print("Warning: Directory for stage {} not found.".format(stage))
        return None

    def getLastStageDir(self):
        """Returns the directory path of the stage with the highest number"""
        pattern = self._getExtraPath('*_*')
        allDirs = glob.glob(pattern)

        if not allDirs:
            return None
        stageDirs = natural_sort(allDirs, rev=True)

        return os.path.abspath(stageDirs[0])

    def prepareSimTrj(self):
        """Concatenate all .netcdf trajectory files from Simulation/Custom stages using cpptraj and wrap atoms
        for visualization"""
        simDirs = self.getSavedTrjStageDirs()
        if not simDirs:
            print("Warning: No saved Simulation/Custom trajectories found")
            return None

        trjFiles = [
            os.path.abspath(os.path.join(d, f))
            for d in simDirs
            for f in os.listdir(d) if f.endswith('.netcdf')
        ]
        if not trjFiles:
            print("Warning: No .netcdf files found")
            return None

        outputTrj = os.path.abspath(self._getExtraPath('prepSimulation.nc'))
        topFile = os.path.abspath(self.amberSystem.get().getTopologyFile())
        cpptrajInParams = ['autoimage']
        cpptrajInParams.append(f"trajout {outputTrj}")
        cpptrajInParams.append("run")
        cpptrajInPath = os.path.abspath(self._getExtraPath('cpptraj.in'))
        with open(cpptrajInPath, 'w') as f:
            f.write("\n".join(cpptrajInParams))

        cmd = f'-p {topFile} -y {" ".join(trjFiles)} -i {cpptrajInPath}'
        amberPlugin.runAmbertools(self, 'cpptraj', cmd, cwd=self._getExtraPath())

        if os.path.exists(outputTrj):
            print(f"Successfully concatenated {len(trjFiles)} trajectory files: {' '.join(trjFiles)}")
            return outputTrj
        print("Error: Concatenation failed")
        return None

    def getSavedTrjStageDirs(self):
        """Return the final contiguous block of Simulation/Custom dirs with saved trajectories."""
        stageDirs = natural_sort(glob.glob(self._getExtraPath('*_Simulation')) +
                                 glob.glob(self._getExtraPath('*_Custom')), rev=True)
        savedStageDirs = []
        for stageDir in stageDirs:
            trjFiles = [f for f in os.listdir(stageDir) if f.endswith('.netcdf')]
            if not trjFiles:
                break
            savedStageDirs.append(stageDir)
        savedStageDirs.reverse()
        return savedStageDirs

    def calculateSavedTrjTime(self):
        """Calculate the simulation time represented by the final saved trajectory block."""
        savedStageNames = {os.path.basename(stageDir) for stageDir in self.getSavedTrjStageDirs()}
        totalPs = 0.0
        workSteps = self.workFlowSteps.get()

        for i, dicLine in enumerate(workSteps.split('\n'), start=1):
            if not dicLine.strip():
                continue

            msjDic = eval(dicLine)
            stageName = '{}_{}'.format(i, msjDic.get('stepType'))
            if stageName in savedStageNames and msjDic.get('stepType') in ['Simulation', 'Custom']:
                totalPs += msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.0)

        return totalPs

    def calculateTotalSimTime(self):
        """
        Calculates the total simulation time (in ps) by summing
        MDSteps * TimeStep for all production and heating stages.
        """
        totalPs = 0.0
        workSteps = self.workFlowSteps.get()

        for dicLine in workSteps.split('\n'):
            if not dicLine.strip(): continue

            msjDic = eval(dicLine)

            if msjDic.get('stepType') in ['Simulation', 'Custom']:
                totalPs += msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.0)

        return totalPs

    def getLastHeatingTemp(self):
        """Get the final temperature from the last Heating step in the workflow"""
        workSteps = self.workFlowSteps.get()
        lines = workSteps.split('\n')
        lastHeatingTemp = 300  # Default

        for dicLine in reversed(lines):
            if dicLine.strip() == '':
                continue
            try:
                msjDic = eval(dicLine)
                if msjDic.get('stepType') == 'Heating':
                    lastHeatingTemp = msjDic.get('FiTemp', 300)
                    break
            except:
                continue

        return lastHeatingTemp
