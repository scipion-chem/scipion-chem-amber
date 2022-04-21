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
    _ensemTypes = ['Energy min', 'NVT', 'NPT']

    _thermostats = ['no', 'Andersen', 'Langevin', 'Nose-Hoover', 'Nose-Hoover RESPA', 'Berendsen']
    _barostats = ['no', 'Berendsen', 'Monte Carlo']
    _coupleStyle = ['No pressure scaling', 'isotropic', 'anisotropic', 'semiisotropic']

    _shakeAlgorithm = ['Shake not performed', 'Bonds involving hydrogens are constrains', 'all bonds are constrained']

    _paramNames = ['simTime', 'timeStep', 'timeNeigh', 'saveTrj', 'trajInterval', 'temperature', 'tempRelaxCons',
                   'tempCouple', 'pressure', 'presRelaxCons', 'presCouple']
    _enumParamNames = ['integrator', 'ensemType', 'thermostat', 'barostat']
    _defParams = {'simTime': 100, 'timeStep': 0.002, 'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0,
                  'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'integrator': 'md',
                  'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restrainForce': 50.0,
                  'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman',
                  'restrains': 'None'}

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
        line.addParam('pressureDynamics', params.EnumParam, default=2,
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

            tprFile = self.callGROMPP(mdpFile)
            self.callMDRun(tprFile, saveTrj=msjDic['saveTrj'])

    def createOutputStep(self):
        lastAmberFile, lastTopoFile, lastTprFile = self.getPrevFinishedStageFiles()
        localAmberFile, localTopFile = self._getPath('outputSystem.nc'), self._getPath('systemTopology.top')
        shutil.copyfile(lastAmberFile, localAmberFile)
        shutil.copyfile(lastTopoFile, localTopFile)

        outTrj = self.concatTrjFiles(outTrj='outputTrajectory.nc', tprFile=lastTprFile)

        outSystem = AmberSystem(filename=localAmberFile)
        outSystem.setTopologyFile(localTopFile)
        if outTrj:
            outSystem.setTrajectoryFile(outTrj)

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

    def _warnings(self):
        warns = []
        # Global warnings
        inpSystem = self.gromacsSystem.get()
        if str(inpSystem.getForceField()).startswith('gromos'):
            warns.append('\nStep all : GROMOS force field is not recommended by GROMACS: '
                         'https://chemrxiv.org/engage/chemrxiv/article-details/60c74701bdbb895afaa38ce2')
        elif str(inpSystem.getForceField()).startswith('charmm36') and 'CA' in inpSystem.getIons() and \
                inpSystem.getIons()['CA'] > 20:
            warns.append('\nStep all : More than 20 non-matching atom names because of (Cal - CA) in charmm36. '
                         'Cal from topology file will be used')

        prevTrj = False
        for step, wStep in enumerate(self.workFlowSteps.get().strip().split('\n')):
            if wStep in ['', None]:
                msjDic = self.createMSJDic()
            else:
                msjDic = eval(wStep)

            if msjDic['ensemType'] == 'NPT' and msjDic['barostat'] != 'Berendsen' and not prevTrj:
                warns.append(
                    '\nStep {} : Berendsen is the barostat recommended for system equilibration, {} might not be'
                    ' the best option for the first trajectory saved'.format(step + 1, msjDic['barostat']))
            if msjDic['saveTrj']:
                prevTrj = True

            if msjDic['ensemType'] != 'Energy min':
                tCoup = msjDic['timeNeigh'] if msjDic['tempCouple'] == -1 else msjDic['tempCouple']
                if msjDic['thermostat'] == 'Nose-Hoover' and 20 * tCoup * msjDic['timeStep'] > msjDic['tempRelaxCons']:
                    warns.append('\nStep {} : For proper integration of the Nose-Hoover thermostat, tau-t ({}) should '
                                 'be at least 20 times larger than nsttcouple*dt ({}*{})'.format(step + 1,
                                                                                                 msjDic[
                                                                                                     'tempRelaxCons'],
                                                                                                 tCoup,
                                                                                                 msjDic['timeStep']))

        return warns

    def _validate(self):
        vals = []
        for step, wStep in enumerate(self.workFlowSteps.get().strip().split('\n')):
            if wStep in ['', None]:
                msjDic = self.createMSJDic()
            else:
                msjDic = eval(wStep)
            if msjDic['ensemType'] != 'Energy min':
                if 'Andersen' in msjDic['thermostat'] and msjDic['integrator'] == 'md':
                    vals.append(
                        'Step {} : Andersen temperature control not supported for integrator md.'.format(step + 1))
        return vals

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
        mdpFile = os.path.join(stageDir, 'stage_{}.mdp'.format(mdpStage))

        params = '&cntrl \n' \
                 'dt = {},' \
                 'nstlim = {},'.format(msjDic['timeStep'], msjDic['simTime'])

        if msjDic['ensemType'] == 'Energy min':
            params += 'imin = 1,ntb = 2,'
        if msjDic['ensemType'] == 'NVT':
            params += 'imin = 0, ntb = 1,'
        if msjDic['ensemType'] == 'NPT':
            params += 'imin = 0, ntb = 2,'

        if msjDic['thermostat'] == 'no':
            params += 'ntt = 0,'
        if msjDic['thermostat'] == 'Andersen':
            params += 'ntt = 2,'
        if msjDic['thermostat'] == 'Langevin':
            params += 'ntt = 3,'
        if msjDic['thermostat'] == 'Nose-Hoove':
            params += 'ntt = 9,'
        if msjDic['thermostat'] == 'Nose-Hoover RESPA':
            params += 'ntt = 10,'
        if msjDic['thermostat'] == 'Berendsen':
            params += 'ntt = 11,'

        params += 'temp0 = {}, pres0 = {},'.format(msjDic['temperature'], msjDic['pressure'])

        if msjDic['barostat'] == 'Berendsen':
            params += 'barostat = 1,'
        elif msjDic['barostat'] == 'Monte Carlo':
            params += 'barostat = 2,'
        else:
            params += ''

        if msjDic['pressureDynamic'] == 'No pressure scaling':
            params += 'ntp = 0,'
        if msjDic['pressureDynamic'] == 'isotropic':
            params += 'ntp = 1,'
        if msjDic['pressureDynamic'] == 'anisotropic':
            params += 'ntp = 2,'
        if msjDic['pressureDynamic'] == 'semiisotropic':
            params += 'ntp = 3,'

        if msjDic['Shake'] == 'Shake not performed':
            params += 'ntp = 1,'
        if msjDic['Shake'] == 'Bonds involving hydrogens are constrains':
            params += 'ntp = 2,'
        if msjDic['Shake'] == 'all bonds are constrained':
            params += 'ntp = 3,'

        params += '\n &end'

        with open(mdpFile, 'w') as f:
            f.write(params)

        return mdpFile

    def callAmber(self, mdpFile):
        stageDir = os.path.dirname(mdpFile)
        stage = os.path.split(stageDir)[-1]
        stageNum = stage.replace('stage_', '').strip()
        amberFile, topFile, _ = self.getPrevFinishedStageFiles(stage)
        outFile = '{}.in'.format(stage)
        topFile = self.AmberSystem.get().getTopologyFile()
        crdFile = self.AmberSystem.get().getfilename()

        if self.checkIfPrevTrj(stageNum):
            prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        else:
            prevTrjStr = ''

        command = '-i {} -c {} -p {}'.format(outFile, crdFile, topFile)
        # Manage warnings
        nWarns = self.countWarns(stageNum)
        print('{} warnings in stage {}'.format(nWarns, stageNum))
        if nWarns >= 1:
            command += ' -maxwarn {}'.format(nWarns)

        amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)
        return os.path.join(stageDir, outFile)


    def getPrevFinishedStageFiles(self, stage=None, reverse=False):
        '''Return the previous .gro and topology files if number stage is provided.
        If not, returns the ones of the lastest stage'''
        topFile = self.AmberSystem.get().getTopologyFile()
        if stage:
            stageNum = stage.replace('stage_', '').strip()
            if stageNum == '1':
                groFile = os.path.abspath(self.AmberSystem.get().getSystemFile())
                tprFile = None

            else:
                prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
                for file in os.listdir(prevDir):
                    if '.gro' in file:
                        groFile = os.path.join(prevDir, file)
                    elif '.tpr' in file:
                        tprFile = os.path.join(prevDir, file)
        else:
            stageDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=reverse)
            for file in os.listdir(stageDirs[-1]):
                if '.gro' in file:
                    groFile = os.path.join(stageDirs[-1], file)
                elif '.tpr' in file:
                    tprFile = os.path.join(stageDirs[-1], file)

        return os.path.abspath(groFile), os.path.abspath(topFile), tprFile

    def checkIfPrevTrj(self, stageNum):
        if stageNum == '1':
            return False
        else:
            prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
            for file in os.listdir(prevDir):
                if '.cpt' in file:
                    return os.path.join(prevDir, file)
        return False

    def getTrjFiles(self):
        trjFiles = []
        stagesDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=True)
        for sDir in stagesDirs:
            cont = False
            for file in os.listdir(sDir):
                if '.trr' in file:
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
