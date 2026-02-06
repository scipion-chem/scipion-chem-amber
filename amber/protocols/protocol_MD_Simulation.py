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
    _ensemTypes = [ 'NVT', 'NPT']
    _thermostats = ['Andersen', 'Langevin', 'Nose-Hoover', 'Nose-Hoover RESPA', 'Berendsen']
    _barostats = ['Berendsen', 'Monte Carlo']

    _coupleStyle = ['No pressure scaling', 'isotropic', 'anisotropic', 'semiisotropic']

    _shakeAlgorithm = ['Shake not performed', 'Bonds involving hydrogens are constrains', 'all bonds are constrained']

    _omitParamNames = ['runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                       'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi']

    _key_map = {'Minimization': 'min', 'Heating': 'heat', 'Simulation': 'sim'}

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
        group.addParam('minCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED, condition='energyMin',
                       help='Upload a custom configuration file for the Minimization step.'
                            'Providing a file here will override all other Minimization parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')

        group = form.addGroup('Heating')
        group.addParam('heatMDSteps', params.IntParam, default=5000,
                       label='Number of MD steps:',
                       help='Number of MD steps in run (x * time step = run length in ps)')
        group.addParam('heatTimeStep', params.FloatParam, default=0.002,
                       label='Time step (ps)')
        group.addParam('heatInTemp', params.FloatParam, default=0,
                       label='Initial temperature (K)')
        group.addParam('heatFiTemp', params.FloatParam, default=300,
                       label='Final temperature (K)')
        group.addParam('heatTraj', params.IntParam, default=1000,
                       label='Trajectory step size', help='The coordinates are written to a mdcrd file x times.')
        group.addParam('heatCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED,
                       help='Upload a custom configuration file for the Heating step.'
                            'Providing a file here will override all other Heating parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')

        group = form.addGroup('Simulation')
        group.addParam('simMDSteps', params.IntParam, default=10000,
                       label='Number of MD steps:',
                       help='Number of MD steps in run (nstlim * dt = run length in ps)')
        group.addParam('simTimeStep', params.FloatParam, default=0.002,
                       label='Time step (ps)')
        group.addParam('simTrajStep', params.IntParam, default=1000,
                       label='Trajectory step size', help='The trajectory coordinates are written to a traj file every x steps.')
        group.addParam('simCustomIn', params.TextParam, width=60, readOnly=True, default=None,
                       label='Input sander/pmemd', expertLevel=params.LEVEL_ADVANCED,
                       help='Upload a custom configuration file for the MD simulation.'
                            'Providing a file here will override all other simulation parameters defined in the interface. For detailed syntax and options, '
                            'refer to the Amber Manual https://ambermd.org/doc12/Amber25.pdf.')
        group = form.addGroup('Ensemble')

        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform: no periodicity, NVT or NPT\n')
        line = group.addLine('Temperature control: ',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Relaxation time constant for thermostat (ps)')
        line.addParam('thermostat', params.EnumParam, default=1,
                      label='Thermostat: ', choices=self._thermostats)
        line.addParam('collisFreq', params.FloatParam, default=2.0,
                      label='Collision frequency (1/ps): ', condition='thermostat==1')
        line.addParam('coupConst', params.FloatParam, default=2.0,
                      label='Coupling constant (1/ps): ', condition='thermostat==2')
        line.addParam('fricConst', params.FloatParam, default=2.0,
                      label='Friction constant (1/ps): ', condition='thermostat==3')

        line = group.addLine('Pressure control: ', condition='ensemType==1',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', params.FloatParam, default=1.0, condition='ensemType==1',
                      label='Pressure (bar): ')
        line.addParam('barostat', params.EnumParam, default=1, condition='ensemType==1',
                      label='Barostat type: ', choices=self._barostats)
        line.addParam('pressureScaling', params.EnumParam, default=2, condition='ensemType==1',
                      label='Pressure scaling: ', choices=self._coupleStyle)
        # line = group.addLine('SHAKE algorithm : ', help='In SHAKE algorithm, the system of non-linear constraint '
        #                                                 'equations is solved using the Gauss–Seidel method which '
        #                                                 'approximates the solution of the linear system of equations '
        #                                                 'using the Newton–Raphson method ')
        # line.addParam('Shake', params.EnumParam, default=0,
        #               label='SHAKE algortihm: ', choices=self._shakeAlgorithm)
        #
        # group = form.addGroup('Summary')
        # group.addParam('insertStep', params.StringParam, default='',
        #                label='Insert relaxation step number: ',
        #                help='Insert the defined relaxation step into the workflow on the defined position.\n'
        #                     'The default (when empty) is the last position')
        # group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
        #                label='Summary of steps',
        #                help='Summary of the defined steps. \nManual modification will have no '
        #                     'effect, use the wizards to add / delete the steps')
        # group.addParam('deleteStep', params.StringParam, default='',
        #                label='Delete relaxation step number: ',
        #                help='Delete the step of the specified index from the workflow.')
        # group.addParam('watchStep', params.StringParam, default='',
        #                label='Watch relaxation step number: ',
        #                help='''Watch the parameters step of the specified index from the workflow..\n
        #                                This might be useful if you want to change some parameters of a predefined step.\n
        #                                However, the parameters are not changed until you add the new step (and probably\n
        #                                you may want to delete the previous unchanged step)''')
        # group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        if self.energyMin.get():
            self._insertFunctionStep('minimizationStep')
        self._insertFunctionStep('heatingStep')
        self._insertFunctionStep('simulationStep')
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
            mdpFile = self.generateMDPFileNew(msjDic, 'Minimization')

        outFile = self.callAmber(mdpFile, 'Minimization')

    def heatingStep(self):
        msjDic = self.getStageParamsDicNew('Heating')
        if self.heatCustomIn.get():
            mdpFile = self.customMDPFile(msjDic, 'Heating')
        else:
            mdpFile = self.generateMDPFileNew(msjDic, 'Heating')

        outFile = self.callAmber(mdpFile, 'Heating')

    def simulationStep(self):
        msjDic = self.getStageParamsDicNew('Simulation')
        if self.simCustomIn.get():
            mdpFile = self.customMDPFile(msjDic, 'Simulation')
        else:
            mdpFile = self.generateMDPFileNew(msjDic, 'Simulation')

        outFile = self.callAmber(mdpFile, 'Simulation')

    def createOutputStep(self):
        CrdAmberFile, localTopFile = self._getPath('crdFile.crd'), self._getPath('systemTopology.parm7')
        shutil.copyfile(self.amberSystem.get().getCrdFile(), CrdAmberFile)
        shutil.copyfile(self.amberSystem.get().getTopologyFile(), localTopFile)

        outTrj = self.getSimTrajStepFile()
        outputTrajectory = self._getPath('outputTrajectory.netcdf')
        shutil.copyfile(outTrj, outputTrajectory)

        system_visualization = self._getPath(os.path.basename(self.amberSystem.get().getFileName()))
        shutil.copyfile(self.amberSystem.get().getFileName(), system_visualization)

        mFF, wFF = self.getFFFiles()
        nFrames = self.getNFrames()
        nTime = nFrames * self.simTimeStep.get()

        outSystem = AmberSystem(filename=system_visualization, ff=mFF, wff=wFF, nFrames=nFrames, nTime=nTime)

        outSystem.setTopologyFile(localTopFile)
        if outTrj:
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

    def createSummary(self, msjDic=None):
        '''Creates the displayed summary from the internal state of the steps'''
        sumStr = ''
        if not msjDic:
            for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
                if dicLine != '':
                    msjDic = eval(dicLine)
                    msjDic = self.addDefaultForMissing(msjDic)
                    method, ensemType = msjDic['thermostat'], msjDic['ensemType']
                    sumStr += '{}) Sim. time ({}): {} ns, {} ensemble'. \
                        format(i + 1, msjDic['integrator'], msjDic['simTime'], ensemType)
                    sumStr += ', {} K\n'.format(msjDic['temperature'])
        else:
            msjDic = self.addDefaultForMissing(msjDic)
            method, ensemType = msjDic['thermostat'], msjDic['ensemType']
            sumStr += 'Sim. time ({}): {} ns, {} ensemble'. \
                format(msjDic['simTime'], msjDic['integrator'], ensemType, method)
            sumStr += ', {} K\n'.format(msjDic['temperature'])
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
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

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

    def getStageParamsDicNew(self, type):
        '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
        Type'''
        paramsDic = {}
        for paramName, param in self._definition.iterAllParams():
            if not paramName in self._omitParamNames and not isinstance(param, params.Group) and not isinstance(param,
                                                                                                                params.Line):
                if type == 'All':
                    paramsDic[paramName] = param
                elif type == 'Minimization' and paramName.startswith("min"):
                    paramsDic[paramName] = getattr(self, paramName).get()
                elif type == 'Heating' and paramName.startswith("heat"):
                    paramsDic[paramName] = getattr(self, paramName).get()
                elif type == 'Simulation' and paramName.startswith("sim"):
                    paramsDic[paramName] = getattr(self, paramName).get()
        return paramsDic

    def createMSJDic(self):
        msjDic = {}
        for pName in self.getStageParamsDic(type='Normal').keys():
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self.getStageParamsDic(type='Enum').keys():
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def addDefaultForMissing(self, msjDic):
        '''Add default values for missing parameters in the msjDic'''
        paramDic = self.getStageParamsDic()
        for pName in paramDic.keys():
            if not pName in msjDic:
                msjDic[pName] = paramDic[pName].default
        return msjDic

    def generateMDPFile(self, msjDic, mdpStage):
        '''Generate .in file'''
        stageDir = self._getExtraPath('stage_{}'.format(mdpStage))
        os.mkdir(stageDir)
        mdpFile = os.path.join(stageDir, 'stage_{}.in'.format(mdpStage))

        params = '\n&cntrl \n' \
                 '      dt={},' \
                 ' nstlim={}, ntwr=50, ntwx=50, ntwe=50, '.format(msjDic['timeStep'],
                                                                  int(msjDic['simTime'] / msjDic['timeStep']))

        if msjDic['EnergyMin']:
            params += ' imin=1,'
        else:
            params += ' imin=0,'

        if msjDic['ensemType'] == 'NVT':
            params += ' ntb=1,'
        elif msjDic['ensemType'] == 'NPT':
            params += ' ntb=2,'

        if msjDic['thermostat'] == 'no':
            params += ' ntt=0,'
        if msjDic['thermostat'] == 'Andersen':
            params += ' ntt=2,'
        if msjDic['thermostat'] == 'Langevin':
            params += ' ntt=3, gamma ln=2.0'
        if msjDic['thermostat'] == 'Nose-Hoove':
            params += ' ntt=9,'
        if msjDic['thermostat'] == 'Nose-Hoover RESPA':
            params += ' ntt=10,'
        if msjDic['thermostat'] == 'Berendsen':
            params += ' ntt=11,'

        params += ' temp0={}, pres0={},'.format(msjDic['temperature'], msjDic['pressure'])

        if msjDic['barostat'] == 'Berendsen':
            params += ' barostat=1,'
        elif msjDic['barostat'] == 'Monte Carlo':
            params += ' barostat=2,'
        else:
            params += ''

        if msjDic['ensemType'] == 'NPT':

            if msjDic['pressureDynamics'] == 'isotropic':
                params += ' ntp=1,'
            if msjDic['pressureDynamics'] == 'anisotropic':
                params += ' ntp=2,'
            if msjDic['pressureDynamics'] == 'semiisotropic':
                params += ' ntp=3,'
        else:
            params += ' ntp=0,'

        if msjDic['Shake'] == 'Shake not performed':
            params += ' ntc=1,'
        if msjDic['Shake'] == 'Bonds involving hydrogens are constrains':
            params += ' ntc=2,'
        if msjDic['Shake'] == 'all bonds are constrained':
            params += ' ntc=3,'

        params += '\n&end \nEND'
        print(msjDic)
        with open(mdpFile, 'w') as f:
            f.write(params)

        return mdpFile

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

    def generateMDPFileNew(self, msjDic, type):
        '''Generate .in file'''
        stageDir = self._getExtraPath('{}'.format(type))
        os.mkdir(stageDir)
        mdpFile = os.path.join(stageDir, '{}.in'.format(type))

        params = ''
        if type == 'Minimization':
            params = 'MINIMIZATION\n&cntrl \n' \
                     'imin=1, ntx=1, irest=0, maxcyc={}, ncyc={},ntpr=100,' \
                     ' ntwx=0, cut={}'.format(msjDic['minMaxCycles'], msjDic['minSdCycles'], msjDic['minIntCutoff'])

        if type == 'Heating':
            params = 'HEATING\n&cntrl \n' \
                     'imin=0, nstlim={}, dt={}, ntf=2, ntc=2, tempi={}, ' \
                     'temp0={}, ntpr={} , ntwx={}, ntb=1, ntp=0, ig=-1, ' \
                     'cut=8.0 /'.format(msjDic['heatMDSteps'],
                                           msjDic['heatTimeStep'],
                                           msjDic['heatInTemp'],
                                           msjDic['heatFiTemp'],
                                           msjDic['heatTraj'],
                                           msjDic['heatTraj'])
            params += self.addThermostatParams(msjDic)
            params += '\n&wt type=\'TEMP0\', istep1=0, istep2={}, value1={}, value2={} /\n'.format(
                msjDic['heatMDSteps'],
                msjDic['heatInTemp'],
                msjDic['heatFiTemp'],
                msjDic['heatFiTemp'])
        if type == 'Simulation':
            params = 'MD SIMULATION\n&cntrl \n' \
                     'imin=0, ntx=5, irest=1, nstlim={}, dt={}, ntf=2, ntc=2, ' \
                     'temp0={}, ntpr={} , ntwx={}, ig=-1, ' \
                     'cut=8.0'.format(msjDic['simMDSteps'],
                                       msjDic['simTimeStep'],
                                       self.heatFiTemp.get(),
                                       msjDic['simTrajStep'],
                                       msjDic['simTrajStep'])
            params += self.addThermostatParams(msjDic)
            if self.getEnumText('ensemType') == 'NPT':
                params += self.addBarostatParams(msjDic)

        params += '\n&end \nEND'
        with open(mdpFile, 'w') as f:
            f.write(params)

        print('sander/pmemd.cuda input params are:\n {}'.format(params))

        return mdpFile

    def addThermostatParams(self, msjDic):
        ntt = None
        gamma_ln = None
        thermostat = self.getEnumText('thermostat')
        if thermostat == 'Andersen':
            ntt = 2
        elif thermostat == 'Langevin':
            ntt = 3
            gamma_ln = self.collisFreq.get()
        elif thermostat == 'Nose-Hoove':
            ntt = 9
            gamma_ln = self.coupConst.get()
        elif thermostat == 'Nose-Hoover RESPA':
            ntt = 10
            gamma_ln = self.fricConst.get()
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

        if self.getEnumText('ensemType') == 'NPT':
            pres0 = self.pressure.get()
            ntb = 2

            pressScaParam = self.getEnumText('pressureScaling')
            barosParam = self.getEnumText('barostat')

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

    # def callAmber(self, mdpFile, saveTrj=True):
    #     inputStructure = os.path.abspath(self.amberSystem.get().getFileName())
    #     systemBasename = os.path.basename(inputStructure.split(".")[0])
    #
    #     stageDir = os.path.dirname(mdpFile)
    #     stage = os.path.split(stageDir)[-1]
    #     stageNum = stage.replace('stage_', '').strip()
    #     amberFile = self.getPrevFinishedStageFiles(stage)
    #     outFile = '{}.in'.format(stage)
    #     topFile = self.amberSystem.get().getTopologyFile()
    #     crdFile = self.amberSystem.get().get
    #     print(stageDir)
    #
    #     if self.checkIfPrevTrj(stageNum):
    #         prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
    #     else:
    #         prevTrjStr = ''
    #
    #     command = '-i {} -c {} -p {} -r {}.r \
    #                 -o {}.o -x {}.netcdf -e {}.e -ref {}.crd \
    #                 -inf min.inf'.format(outFile, amberFile, topFile, *[stage] * 5)
    #
    #     # Manage warnings
    #     nWarns = self.countWarns(stageNum)
    #     print('{} warnings in stage {}'.format(nWarns, stageNum))
    #     if nWarns >= 1:
    #         command += ' -maxwarn {}'.format(nWarns)
    #
    #     amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)
    #     if not saveTrj:
    #         trjFile = os.path.join(stageDir, '{}.trr'.format(stage))
    #         os.remove(trjFile)
    #
    #     return os.path.join(stageDir, outFile)

    def callAmber(self, mdpFile, type, saveTrj=True):
        inputStructure = os.path.abspath(self.amberSystem.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])
        stageDir = os.path.dirname(mdpFile)

        inputFile = '{}.in'.format(type)
        topFile = self.amberSystem.get().getTopologyFile()
        crdFile = self.amberSystem.get().getCrdFile()
        outFile = '{}.o'.format(type)

        # if self.checkIfPrevTrj(stageNum):
        #     prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        # else:
        #     prevTrjStr = ''
        if type == 'Minimization':
            command = '-i {} -c {} -p {} -r {}.ncrst -o {}.o' \
                       ' -ref {}.crd' \
                       ' -inf {}.inf'.format(inputFile, crdFile, topFile, *[type] * 4)
        elif type == 'Heating':
            if self.energyMin.get():
                crdFile = os.path.abspath(os.path.join(self._getExtraPath(), 'Minimization', 'Minimization.ncrst'))
            command = '-i {} -c {} -p {} -r {}.ncrst' \
                      ' -o {}.o -ref {}.crd' \
                      ' -x {}.netcdf -inf {}.inf'.format(inputFile, crdFile, topFile, *[type] * 5)

        elif type == 'Simulation':
            crdFile = os.path.abspath(os.path.join(self._getExtraPath(), 'Heating', 'Heating.ncrst'))
            command = '-i {} -c {} -p {} -r {}.ncrst' \
                      ' -o {}.o -ref {}.crd' \
                      ' -x {}.netcdf -inf {}.inf'.format(inputFile, crdFile, topFile, *[type] * 5)

        if self.useGpu.get():
            os.environ["CUDA_VISIBLE_DEVICES"] = self.gpuList.get()
            amberPlugin.runPmemd(self, ' -O ', args=command, cwd=stageDir)

        else:
            # Use sander for CPU execution
            engine = 'sander'
            amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)

        if not saveTrj:
            trjFile = os.path.join(type, '{}.trr'.format(type))
            os.remove(trjFile)

        return os.path.join(stageDir, outFile)

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

    def getSimTrajStepFile(self):
        trjFile = os.path.join(self._getExtraPath(), 'Simulation', 'Simulation.netcdf')
        return trjFile

    def getNFrames(self):
      nFrames = self.simMDSteps.get() // self.simTrajStep.get()
      return nFrames

    def getFFFiles(self):
      system = self.amberSystem.get()
      return system.getForceField(), system.getWaterForceField()
