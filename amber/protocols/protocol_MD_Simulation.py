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
from cProfile import label
from email.policy import default
from os.path import join

from pyworkflow.protocol import params
from pyworkflow.utils import Message, runJob, createLink

import amber
from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from amber.objects import *
from amber.constants import *
from amber import Plugin as amberPlugin


class AmberMDSimulation(EMProtocol):
    """
        This protocol will perform energy minimization and equilibrium on the system previously prepared by the protocol
         "system prepartion". This step is necessary to energy minimize the system in order to avoid unwanted conformations.
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
                       'simInsertStep']

    _key_map = {'Minimization': 'min', 'Heating': 'heat', 'Simulation': 'sim'}

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

        group.addParam('minCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED, condition='energyMin',
                       help='Upload a custom configuration file for the Minimization step.'
                            'Providing a file here will override all other Minimization parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')
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
        line.addParam('heatTraj', params.IntParam, default=1000,
                       label='Trajectory step size', help='The coordinates are written to a mdcrd file x times.')

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

        group.addParam('heatCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED,
                       help='Upload a custom configuration file for the Heating step.'
                            'Providing a file here will override all other Heating parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')

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
        line.addParam('simTrajStep', params.IntParam, default=1000,
                      label='Trajectory step size',
                      help='The trajectory coordinates are written to a traj file every x steps.')

        group.addParam('simEnsemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
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
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('simPressure', params.FloatParam, default=1.0, condition='simEnsemType==1',
                      label='Pressure (bar): ')
        line.addParam('simBarostat', params.EnumParam, default=1, condition='simEnsemType==1',
                      label='Barostat type: ', choices=self._barostats)
        line.addParam('simPressureScaling', params.EnumParam, default=1, condition='simEnsemType==1',
                      label='Pressure scaling: ', choices=self._coupleStyle)

        group.addParam('simCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED,
                       help='Upload a custom configuration file for the MD simulation.'
                            'Providing a file here will override all other simulation parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')

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

        group = form.addGroup('Summary')
        group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
                       label='Summary of steps',
                       help='Summary of the defined steps. \nManual modification will have no '
                            'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete relaxation step number: ',
                       help='Delete the step of the specified index from the workflow.')
        # group.addParam('watchStep', params.StringParam, default='',
        #                label='Watch relaxation step number: ',
        #                help='''Watch the parameters step of the specified index from the workflow..\n
        #                                This might be useful if you want to change some parameters of a predefined step.\n
        #                                However, the parameters are not changed until you add the new step (and probably\n
        #                                you may want to delete the previous unchanged step)''')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

        form.addSection('Recommended workflows')
        group = form.addGroup('Protein')
        group.addParam('proteinDefault', params.LabelParam, label='Protein default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system')
        group = form.addGroup('Protein+Ligand')
        group.addParam('protLigDefault', params.LabelParam, label='Protein+Ligand default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system')
        group = form.addGroup('Transmembrane protein')
        group.addParam('membraneDefault', params.LabelParam, label='Transmembrane default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system')
        group = form.addGroup('Transmembrane protein+Ligand')
        group.addParam('memLigDefault', params.LabelParam, label='Transmembrane+Ligand default workflow',
                       help='Click the wizard to set a simluation with default params for a Protein system')


    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        print(self.workFlowSteps.get())
        self.createGUISummary()
        i = 1
        for wStep in self.workFlowSteps.get().strip().split('\n'):
            self._insertFunctionStep('simulateStageStep', wStep, i)
            i += 1

        self._insertFunctionStep('createOutputStep')

    def simulateStageStep(self, wStep, i):
        if wStep in ['', None]:
            msjDic = self.createMSJDic()
        else:
            msjDic = eval(wStep)
            mdpFile = self.generateMDPFile(msjDic, str(i))

            tprFile = self.callAmber(mdpFile)
            # self.callMDRun(tprFile, saveTrj=msjDic['saveTrj'])

    def minimizationStep(self):
        msjDic = self.getStageParamsDicNew('Minimization')
        if self.minCustomIn.get():
            mdpFile = self.customMDPFile(msjDic, 'Minimization')
        else:
            mdpFile = self.generateMDPFile(msjDic, 'Minimization')

        outFile = self.callAmber(mdpFile, 'Minimization')

    def heatingStep(self):
        msjDic = self.getStageParamsDicNew('Heating')
        if self.heatCustomIn.get():
            mdpFile = self.customMDPFile(msjDic, 'Heating')
        else:
            mdpFile = self.generateMDPFile(msjDic, 'Heating')

        outFile = self.callAmber(mdpFile, 'Heating')

    def simulationStep(self):
        msjDic = self.getStageParamsDicNew('Simulation')
        if self.simCustomIn.get():
            mdpFile = self.customMDPFile(msjDic, 'Simulation')
        else:
            mdpFile = self.generateMDPFile(msjDic, 'Simulation')

        outFile = self.callAmber(mdpFile, 'Simulation')

    # def createOutputStep(self):
    #     CrdAmberFile, localTopFile = self._getPath('crdFile.crd'), self._getPath('systemTopology.parm7')
    #     shutil.copyfile(self.amberSystem.get().getCrdFile(), CrdAmberFile)
    #     shutil.copyfile(self.amberSystem.get().getTopologyFile(), localTopFile)
    #
    #     outTrj = self.getSimTrajStepFile()
    #     outputTrajectory = self._getPath('outputTrajectory.netcdf')
    #     shutil.copyfile(outTrj, outputTrajectory)
    #
    #     system_visualization = self._getPath(os.path.basename(self.amberSystem.get().getFileName()))
    #     shutil.copyfile(self.amberSystem.get().getFileName(), system_visualization)
    #
    #     mFF, wFF = self.getFFFiles()
    #     nFrames = self.getNFrames()
    #     nTime = nFrames * self.simTimeStep.get()
    #
    #     outSystem = AmberSystem(filename=system_visualization, ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)
    #
    #     outSystem.setTopologyFile(localTopFile)
    #     if outTrj:
    #         outSystem.setTrajectoryFile(outputTrajectory)
    #
    #     self._defineOutputs(outputSystem=outSystem)

    def createOutputStep(self):
        lastCrdFile, lastTopoFile, lastOutFile = self.getPrevFinishedStageFiles()
        oriSystemFile = self.amberSystem.get().getSystemFile()

        localCrdFile, localTopFile = self._getPath('outputSystem.rst7'), self._getPath('systemTopology.parm7')
        shutil.copyfile(lastCrdFile, localCrdFile), shutil.copyfile(lastTopoFile, localTopFile)
        # outTrj = self.concatTrjFiles(outTrj='outputTrajectory.xtc', tprFile=lastTprFile)

        mFF, wFF = self.getFFFiles()

        outSystem = AmberSystem(filename=oriSystemFile, ff=mFF, wff=wFF)

        outSystem.setTopologyFile(localTopFile)

        outputTrajectory = self._getPath('outputTrajectory.netcdf')

        concatTrjFile = self.concatSimTrj()
        ############### gestionar las trayectorias - concatenar???
        shutil.copyfile(concatTrjFile, outputTrajectory)
        outSystem.setTrajectoryFile(outputTrajectory)

        self._defineOutputs(outputSystem=outSystem)

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
        '''Creates the displayed summary from workflow steps'''
        sumStr = ''
        lastTemp = 300  # default

        # If no workSteps provided, use the protocol's workFlowSteps
        if workSteps is None:
            workSteps = self.workFlowSteps.get()

        if not workSteps or workSteps.strip() == '':
            return ''

        lines = workSteps.split('\n')
        for i, dicLine in enumerate(lines):
            if dicLine.strip() == '':
                continue

            msjDic = eval(dicLine)
            msjDic = self.addDefaultForMissing(msjDic)
            stepType = msjDic.get('stepType', 'Step')

            lineText = '{}) {} - '.format(i + 1, stepType)

            if stepType == 'Minimization':
                lineText += 'Max Cycles: {}'.format(msjDic.get('MaxCycles', 0))
                if msjDic['Restraint']:
                    lineText += ', restraint on {}'.format(msjDic.get('RestrAtoms'))

            elif stepType == 'Heating':
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lastTemp = msjDic.get('FiTemp', 300)
                lineText += 'Sim. time: {} ps, NVT ensemble, {} K to {} K'.format(
                    nTime, msjDic.get('InTemp', 0), lastTemp)
                if msjDic['Restraint']:
                    lineText += ', restraint on {}'.format(msjDic.get('RestrAtoms'))

            else:  # Simulation
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lineText += 'Sim. time: {} ps, {} ensemble, {} K'.format(
                    nTime, msjDic.get('EnsemType', 'NPT'), lastTemp)
                if msjDic['Restraint']:
                    lineText += ', restraint on {}'.format(msjDic.get('RestrAtoms'))

            sumStr += lineText + '\n'
        return sumStr

    def createDefaultSummary(self, workSteps):
        '''Creates the default summary from the internal state of the steps'''
        sumStr = ''
        lastTemp = 300 # default
        lines = workSteps.split('\n')
        for i, dicLine in enumerate(lines):
            if dicLine.strip() == '':
                continue

            msjDic = eval(dicLine)
            msjDic = self.addDefaultForMissing(msjDic)
            stepType = msjDic.get('stepType', 'Step')

            lineText = '{}) {} - '.format(i + 1, stepType)

            if stepType == 'Minimization':
                lineText += 'Max Cycles: {}'.format(msjDic.get('maxCycles', 0))

            elif stepType == 'Heating':
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lastTemp = msjDic.get('FiTemp')
                lineText += 'Sim. time: {} ps, NVT ensemble, {} K to {} K'.format(
                    nTime, msjDic.get('InTemp'), lastTemp)

            else:
                nTime = msjDic.get('MDSteps', 0) * msjDic.get('TimeStep', 0.002)
                lineText += 'Sim. time: {} ps, {} ensemble, {} K'.format(
                    nTime, msjDic.get('EnsemType'), lastTemp)

            sumStr += lineText + '\n'
        return sumStr

    def createGUISummary(self):
        with open(self._getExtraPath("summary.txt"), 'w') as f:
            if self.workFlowSteps.get():
                f.write(self.createSummary())
            else:
                f.write(self.createSummary(self.createMSJDic()))

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

    def getStageParamsDic(self, type='All'):
        '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
        Type'''
        paramsDic = {}
        for paramName, param in self._definition.iterAllParams():
            if not paramName in self._omitParamNames and not isinstance(param, params.Group) and not isinstance(param,
                                                                                                                params.Line):
                if type == 'All':
                    paramsDic[paramName] = param
                elif type == 'Enum' and isinstance(param, params.EnumParam):
                    paramsDic[paramName] = param
                elif type == 'Normal' and not isinstance(param, params.EnumParam):
                    paramsDic[paramName] = param
        return paramsDic

    def createMSJDic(self, stageType):
        msjDic = {}
        for pName in self.getStageParamsDicNew(type='Normal').keys():
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self.getStageParamsDicNew(type='Enum').keys():
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def getStageParamsDicNew(self, type):
        '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
        Type'''
        paramsDic = {}
        if type == 'Minimization':
            prefix = "min"
        elif type == 'Heating':
            prefix = "heat"
        elif type == 'Simulation':
            prefix = "sim"

        for paramName, param in self._definition.iterAllParams():
            if not paramName in self._omitParamNames and not isinstance(param, params.Group) and not isinstance(param,
                                                                                                                params.Line):
                # if type == 'All':
                #     paramsDic[paramName] = param

                if prefix and paramName.startswith(prefix):
                    cleanName = paramName[len(prefix):]

                    if isinstance(param, params.EnumParam):
                        paramsDic[cleanName] = self.getEnumText(paramName)
                    else:
                        paramsDic[cleanName] = getattr(self, paramName).get()

        paramsDic['stepType'] = type

        return paramsDic

    def addDefaultForMissing(self, msjDic):
        '''Add default values for missing parameters in the msjDic'''
        paramDic = self.getStageParamsDic()
        for pName in paramDic.keys():
            if not pName in msjDic:
                msjDic[pName] = paramDic[pName].default
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

    def generateMDPFile(self, msjDic, i):
        '''Generate .in file'''
        stepType = msjDic['stepType']
        stageDir = self._getExtraPath('{}_{}'.format(i, stepType))
        if os.path.exists(stageDir):
            shutil.rmtree(stageDir)
        os.makedirs(stageDir)

        mdpFile = os.path.join(stageDir, '{}_{}.in'.format(i, stepType))

        params = ''
        if stepType == 'Minimization':
            params = 'MINIMIZATION\n&cntrl \n' \
                     'imin=1, ntx=1, irest=0, maxcyc={}, ncyc={}, ntpr=100,' \
                     ' ntwx=0, cut={}'.format(msjDic['MaxCycles'], msjDic['SdCycles'], msjDic['IntCutoff'])
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' /".format(msjDic['RestrForce'],
                                                                                  RESTRAINS_DIC[msjDic['RestrAtoms']])

        if stepType == 'Heating':
            params = 'HEATING\n&cntrl \n' \
                     'imin=0, nstlim={}, dt={}, ntf=2, ntc=2, tempi={}, ' \
                     'temp0={}, ntpr={} , ntwx={}, ntb=1, ntp=0, ig=-1, ' \
                     'cut=8.0 /'.format(msjDic['MDSteps'],
                                        msjDic['TimeStep'],
                                        msjDic['InTemp'],
                                        msjDic['FiTemp'],
                                        msjDic['Traj'],
                                        msjDic['Traj'])
            params += self.addThermostatParams(msjDic)
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' /".format(msjDic['RestrForce'],
                                                                                  RESTRAINS_DIC[msjDic['RestrAtoms']])

            params += '\n&wt type=\'TEMP0\', istep1=0, istep2={}, value1={}, value2={} /\n'.format(
                msjDic['MDSteps'],
                msjDic['InTemp'],
                msjDic['FiTemp'],
                msjDic['FiTemp'])

        if stepType == 'Simulation':
            params = 'MD SIMULATION\n&cntrl \n' \
                     'imin=0, ntx=5, irest=1, nstlim={}, dt={}, ntf=2, ntc=2, ' \
                     'temp0={}, ntpr={} , ntwx={}, ig=-1, ' \
                     'cut=8.0'.format(msjDic['MDSteps'],
                                      msjDic['TimeStep'],
                                      self.getLastHeatingTemp(),
                                      msjDic['TrajStep'],
                                      msjDic['TrajStep'])
            params += self.addThermostatParams(msjDic)
            if msjDic['EnsemType']== 'NPT':
                params += self.addBarostatParams(msjDic)
            if msjDic['Restraint']:
                params += ", ntr=1, restraint_wt={}, restraintmask='{}' /".format(msjDic['simRestrForce'],
                                                                                  RESTRAINS_DIC[msjDic['RestrAtoms']])

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
                ntp = 3
            else:
                ntp = 0

            if barosParam == 'Berendsen':
                barostat = 1
            elif barosParam == 'Monte Carlo':
                barostat = 2
            else:
                barostat = 0

            params.append(
                f"\nntb={ntb}, ntp={ntp}, pres0={pres0}, barostat={barostat}"
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

        if os.path.exists(outFile):
            return topFile
        crdFile, topFile, _ = self.getPrevFinishedStageFiles(stageNum)

        if stageType == 'Minimization':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst -o {}.o' \
                      ' -inf {}.inf'.format(inputFile, crdFile, topFile, crdFile, *[stage] * 4)
        elif stageType == 'Heating':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst' \
                      ' -o {}.o ' \
                      ' -x {}.netcdf -inf {}.inf'.format(inputFile, crdFile, topFile, crdFile, *[stage] * 5)

        elif stageType == 'Simulation':
            command = '-i {} -c {} -p {} -ref {} -r {}.ncrst' \
                      ' -o {}.o' \
                      ' -x {}.netcdf -inf {}.inf'.format(inputFile, crdFile, topFile, crdFile, *[stage] * 5)

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

    def checkIfPrevTrj(self, stageNum):
        if stageNum == '1':
            return False
        else:
            prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
            for file in os.listdir(prevDir):
                if '.netcdf' in file:
                    return os.path.join(prevDir, file)
        return False

    def getTrjFiles(self):
        trjFiles = []
        stagesDirs = natural_sort(glob.glob(self._getExtraPath()), rev=True)
        for sDir in stagesDirs:
            cont = False
            for file in os.listdir(sDir):
                if '.netcdf' in file:
                    trjFiles.append(os.path.abspath(os.path.join(sDir, file)))
                    cont = True
            if not cont:
                break
        trjFiles.reverse()
        return trjFiles

    def getNFrames(self):
        nFrames = self.simMDSteps.get() // self.simTrajStep.get()
        return nFrames

    def getFFFiles(self):
        system = self.amberSystem.get()
        return system.getForceField(), system.getWaterForceField()

    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def getStageDir(self, stage):
        """ Returns the directory path for a given stage number.
        """
        pattern = self._getExtraPath('{}*'.format(stage))
        matchingDirs = glob.glob(pattern)

        if matchingDirs: return os.path.abspath(matchingDirs[0])

        # 3. Si no existe, retornamos None o podrías lanzar un error
        print("Warning: Directory for stage {} not found.".format(stage))
        return None

    def getLastStageDir(self):
        """ Returns the directory path of the stage with the highest number. """
        pattern = self._getExtraPath('*_*')
        allDirs = glob.glob(pattern)

        if not allDirs:
            return None
        stageDirs = natural_sort(allDirs, rev=True)

        return os.path.abspath(stageDirs[0])

    def concatSimTrj(self):
        """Concatenate all .netcdf trajectory files from Simulation stages using cpptraj"""
        pattern = self._getExtraPath('*_Simulation')
        simDirs = glob.glob(pattern)

        if not simDirs:
            print("Warning: No Simulation directories found")
            return None

        simDirs = natural_sort(simDirs)

        # Collect all .netcdf files from Simulation stages
        trjFiles = []
        for sDir in simDirs:
            for file in os.listdir(sDir):
                if file.endswith('.netcdf'):
                    trjFiles.append(os.path.abspath(os.path.join(sDir, file)))

        if not trjFiles:
            print("Warning: No .netcdf trajectory files found in Simulation stages")
            return None

        if len(trjFiles) == 1:
            # Only one trajectory, no need to concatenate
            print("Only one Simulation trajectory found, skipping concatenation")
            return trjFiles[0]

        outputTrj = os.path.abspath(self._getExtraPath('concatSimulation.nc'))
        topFile = self.amberSystem.get().getTopologyFile()

        # command for cpptraj
        command = '-p {} -y {} -x {}'.format(topFile, ' '.join(trjFiles), outputTrj)

        amberPlugin.runAmbertools(self, 'cpptraj', command, cwd=self._getExtraPath())

        if os.path.exists(outputTrj):
            print(f"Successfully concatenated {len(trjFiles)} trajectory files")
            return outputTrj
        else:
            print("Error: Concatenation failed, output file not created")
            return None


    def getLastHeatingTemp(self):
        """Get the final temperature from the last Heating step in the workflow"""
        workSteps = self.workFlowSteps.get()
        lines = workSteps.split('\n')
        lastHeatingTemp = 300  # Default

        # Search for the last Heating step and get its FiTemp
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