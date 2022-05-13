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

    _label = 'Molecular dynamics simulation'
    _ensemTypes = ['no periodicity', 'NVT', 'NPT']

    _thermostats = ['no', 'Andersen', 'Langevin', 'Nose-Hoover', 'Nose-Hoover RESPA', 'Berendsen']
    _barostats = ['no', 'Berendsen', 'Monte Carlo']
    _coupleStyle = ['No pressure scaling', 'isotropic', 'anisotropic', 'semiisotropic']

    _shakeAlgorithm = ['Shake not performed', 'Bonds involving hydrogens are constrains', 'all bonds are constrained']

    _paramNames = ['simTime', 'timeStep', 'timeNeigh', 'saveTrj', 'trajInterval', 'temperature', 'tempRelaxCons',
                   'tempCouple', 'pressure', 'presRelaxCons', 'presCouple', 'EnergyMin']
    _enumParamNames = ['integrator', 'ensemType', 'thermostat', 'barostat', 'pressureDynamics', 'Shake']
    _defParams = {'simTime': 100, 'timeStep': 0.002, 'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0,
                  'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'integrator': 'md',
                  'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restrainForce': 50.0,
                  'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman',
                  'restrains': 'None', 'pressureDynamics': 'anisotropic', 'Shake': 'Shake not performed'}

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection('Minimization')
        form.addParam('AmberSystem', params.PointerParam, label="Input Amber System: ",
                      pointerClass='AmberSystem',
                      allowsNull=True,
                      help='Amber solvated system to be simulated')
        group = form.addGroup('Trajectory')
        group.addParam('saveTrj', params.BooleanParam, default=self._defParams['saveTrj'],
                       label="Save trajectory: ",
                       help='Save trajectory of the atoms during stage simulation.'
                            'The output will concatenate those trajectories which appear after the last stage '
                            'where the trajectory was not saved.')
        group.addParam('trajInterval', params.FloatParam, default=self._defParams['trajInterval'],
                       label='Interval time (ps):', condition='saveTrj',
                       help='Time between each frame recorded in the simulation (ps)')

        group = form.addGroup('Simulation time')
        group.addParam('simTime', params.FloatParam, default=100,
                       label='Simulation time (ps):',
                       help='Total time of the simulation stage (ps)')
        group.addParam('timeStep', params.FloatParam, default=0.002,
                       label='Simulation time steps (ps)[dt]:',
                       help='Time of the steps for simulation (ps)[dt] \n 0.002ps is recommended if SHAKE '
                            'algorithm is chosen \n'
                            'if not, 0.001ps is recommended ')

        group = form.addGroup('Ensemble')
        group.addParam('EnergyMin', params.BooleanParam,
                       label='Energy Minimization: ', default=True,
                       help='Whether this step should perform energy minimization')

        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform in the step: Energy minimization, NVT or NPT\n')

        line = group.addLine('Temperature settings: ', condition='ensemType!=0',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Relaxation time constant for thermostat (ps)')
        line.addParam('temperature', params.FloatParam, default=300, condition='ensemType!=0',
                      label='Temperature: ')
        line.addParam('thermostat', params.EnumParam, default=5, condition='ensemType!=0',
                      label='Thermostat: ', choices=self._thermostats)

        line = group.addLine('Pressure settings: ', condition='ensemType==2',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', params.FloatParam, default=1.0,
                      label='   Pressure (bar):   ')
        line.addParam('pressureDynamics', params.EnumParam, default=0,
                      label='  Pressure dynamics type:   ', choices=self._coupleStyle)
        line.addParam('barostat', params.EnumParam, default=1,
                      label='  Barostat type:   ', choices=self._barostats)
        line = group.addLine('SHAKE algorithm : ', help='In SHAKE algorithm, the system of non-linear constraint '
                                                        'equations is solved using the Gauss–Seidel method which '
                                                        'approximates the solution of the linear system of equations '
                                                        'using the Newton–Raphson method ')
        line.addParam('Shake', params.EnumParam, default=0,
                      label='SHAKE algortihm: ', choices=self._shakeAlgorithm)

        group = form.addGroup('Summary')
        group.addParam('insertStep', params.StringParam, default='',
                       label='Insert relaxation step number: ',
                       help='Insert the defined relaxation step into the workflow on the defined position.\n'
                            'The default (when empty) is the last position')
        group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
                       label='Summary of steps',
                       help='Summary of the defined steps. \nManual modification will have no '
                            'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete relaxation step number: ',
                       help='Delete the step of the specified index from the workflow.')
        group.addParam('watchStep', params.StringParam, default='',
                       label='Watch relaxation step number: ',
                       help='''Watch the parameters step of the specified index from the workflow..\n
                                       This might be useful if you want to change some parameters of a predefined step.\n
                                       However, the parameters are not changed until you add the new step (and probably\n
                                       you may want to delete the previous unchanged step)''')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
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

    def createOutputStep(self):
        CrdAmberFile, localTopFile = self._getPath('CrdFile.crd'), self._getPath('systemTopology.top')

        shutil.copyfile(self.AmberSystem.get().getSystemFile(), CrdAmberFile)
        shutil.copyfile(self.AmberSystem.get().getTopologyFile(), localTopFile)

        outTrj = self.getTrjFiles()
        outputTrajectory = self._getPath('outputTrajectory.nc')
        shutil.copyfile(outTrj[-1], outputTrajectory)

        outSystem = AmberSystem(filename=CrdAmberFile)
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
                    sumStr += '{}) Sim. time ({}): {} ps, {} ensemble'. \
                        format(i + 1, msjDic['integrator'], msjDic['simTime'], ensemType)
                    sumStr += ', {} K\n'.format(msjDic['temperature'])
        else:
            msjDic = self.addDefaultForMissing(msjDic)
            method, ensemType = msjDic['thermostat'], msjDic['ensemType']
            sumStr += 'Sim. time ({}): {} ps, {} ensemble'. \
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
            methods.append('The methods used to perform the energy minimization protocol have been *"gmx grompp"* to '
                           'create the tpr files and *"gmx mdrun"* to run the minimization simulation.')

        return methods

    ######################## UTILS ##################################
    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def createMSJDic(self):
        msjDic = {}
        for pName in self._paramNames:
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self._enumParamNames:
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def addDefaultForMissing(self, msjDic):
        '''Add default values for missing parameters in the msjDic'''
        for pName in [*self._paramNames, *self._enumParamNames]:
            if not pName in msjDic:
                msjDic[pName] = self._defParams[pName]
        return msjDic

    def generateMDPFile(self, msjDic, mdpStage):
        stageDir = self._getExtraPath('stage_{}'.format(mdpStage))
        os.mkdir(stageDir)
        mdpFile = os.path.join(stageDir, 'stage_{}.in'.format(mdpStage))

        params = '\n &cntrl \n' \
                 '      dt={},' \
                 ' nstlim={}, ntwr=50, ntwx=50, ntwe=50, '.format(msjDic['timeStep'],
                                                                  int(msjDic['simTime'] / msjDic['timeStep']))

        if msjDic['EnergyMin']:
            params += ' imin=1,'
        else:
            params += ' imin=0,'

        if msjDic['ensemType'] == 'no periodicity':
            params += ' ntb=0, cut=99'
        elif msjDic['ensemType'] == 'NVT':
            params += ' ntb=1,'
        elif msjDic['ensemType'] == 'NPT':
            params += ' ntb=2,'

        if msjDic['thermostat'] == 'no':
            params += ' ntt=0,'
        if msjDic['thermostat'] == 'Andersen':
            params += ' ntt=2,'
        if msjDic['thermostat'] == 'Langevin':
            params += ' ntt=3,'
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

        params += '\n &end \n END'
        print(msjDic)
        with open(mdpFile, 'w') as f:
            f.write(params)

        return mdpFile

    def callAmber(self, mdpFile, saveTrj=True):
        inputStructure = os.path.abspath(self.AmberSystem.get().getFileName())
        systemBasename = os.path.basename(inputStructure.split(".")[0])

        stageDir = os.path.dirname(mdpFile)
        stage = os.path.split(stageDir)[-1]
        stageNum = stage.replace('stage_', '').strip()
        amberFile = self.getPrevFinishedStageFiles(stage)
        outFile = '{}.in'.format(stage)
        topFile = self.AmberSystem.get().getTopologyFile()
        crdFile = self.AmberSystem.get().getSystemFile()

        if self.checkIfPrevTrj(stageNum):
            prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        else:
            prevTrjStr = ''

        command = '-i {} -c {} -p {} -r {}.r \
                                       -o {}.o \
                                       -x {}.1.nc \
                                       -e {}.1.e \
                                       -ref {}.crd \
                                       -inf min.inf'.format(outFile, amberFile, topFile, *[stage] * 5)

        # Manage warnings
        nWarns = self.countWarns(stageNum)
        print('{} warnings in stage {}'.format(nWarns, stageNum))
        if nWarns >= 1:
            command += ' -maxwarn {}'.format(nWarns)

        amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)
        if not saveTrj:
            trjFile = os.path.join(stageDir, '{}.trr'.format(stage))
            os.remove(trjFile)

        return os.path.join(stageDir, outFile)

    def getPrevFinishedStageFiles(self, stage=None, reverse=False):
        '''Return the previous .gro and topology files if number stage is provided.
        If not, returns the ones of the lastest stage'''
        print(stage)
        if stage:
            stageNum = stage.replace('stage_', '').strip()
            if stageNum == '1':
                amberFile = os.path.abspath(self.AmberSystem.get().getSystemFile())

            else:
                prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
                for file in os.listdir(prevDir):
                    if '.r' in file:
                        amberFile = os.path.join(prevDir, file)

        else:
            stageDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=reverse)
            for file in os.listdir(stageDirs[-1]):
                if '.r' in file:
                    amberFile = os.path.join(stageDirs[-1], file)

        return os.path.abspath(amberFile)

    def checkIfPrevTrj(self, stageNum):
        if stageNum == '1':
            return False
        else:
            prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
            for file in os.listdir(prevDir):
                if '.nc' in file:
                    return os.path.join(prevDir, file)
        return False

    def getTrjFiles(self):
        trjFiles = []
        stagesDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=True)
        for sDir in stagesDirs:
            cont = False
            for file in os.listdir(sDir):
                if '.nc' in file:
                    trjFiles.append(os.path.abspath(os.path.join(sDir, file)))
                    cont = True
            if not cont:
                break

        trjFiles.reverse()
        return trjFiles

    def countWarns(self, stageNum):
        nWarns = 0
        for warn in self._warnings():
            if warn.split()[1] in ['all', str(stageNum)]:
                nWarns += 1
        return nWarns
