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

from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from amber.objects import *
from amber import Plugin as amberPlugin

EMIN, NVE, NVT, NPT = 0, 1, 2, 3                                    # Ensembles
CG, ST_CG, ST, XMIN, LMOD = 0, 1, 2, 3, 4                           # Integrators
WCOUP, ANDER, LANG, NOSE, NOSE_RES, BEREN = 1, 2, 3, 9, 10, 11      # Thermostats
BEREN_P, MONCAR = 1, 2                                              # Barostats
NO_SHAKE, SHAKE_H, SHAKE_ALL = 0, 1, 2                              # SHAKE options
SOL, PROT, BB, CA, CUST = 0, 1, 2, 3, 4                     # Restraints options


class AmberMDSimulation(EMProtocol):
    """
        This protocol will perform energy minimization and equilibrium on the system previously prepared by the protocol
         "system preparation". This step is necessary to energy minimize the system in order to avoid unwanted conformations.
    """

    _label = 'Molecular dynamics simulation'
    _ensemTypes = ['Energy Minimization', 'NVE', 'NVT', 'NPT']
    _integrators = ['Conjugate Gradient', 'Steep and CG', 'Steepest Descent', 'XMIN', 'LMOD']

    _thermostats = ['Weak-Coupling', 'Andersen', 'Langevin', 'Nose-Hoover', 'Nose-Hoover RESPA', 'Berendsen']
    _map_therm = {0: WCOUP, 1: ANDER, 2: LANG, 3: NOSE, 4: NOSE_RES, 5: BEREN}

    _barostats = ['Berendsen', 'Monte Carlo']
    _map_bar = {0: BEREN_P, 1: MONCAR}

    _shakeAlgorithm = ['SHAKE not performed', 'Bonds involving hydrogen are constrained', 'All bonds are constrained']
    _restraints = ['Solute', 'Protein', 'Ligand', 'Backbone', 'Alpha carbons', 'Custom']

    _paramNames = ['temperature', 'thermostat', 'tautp', 'gamma_ln', 'vlimit',                  # Temperature settings
                   'pressure', 'barostat', 'taup',                                              # Pressure settings
                   'saveTrj', 'trajInterval',                                                   # Trajectory settings
                   'integrator', 'ensemType',                                                   # Ensemble settings
                   'maxcyc', 'ncyc', 'dx0', 'drms',                                             # Minimization settings
                   'simTime', 'timeStep',                                                       # MD settings
                   'shake', 'restraint', 'restraint_enum', 'restraintmask', 'restraint_wt', 'restraint_heavy',     # Restraints settings
                   'extraParams']

    _enumParamNames = []
    _defParams = {'temperature': 300.0, 'thermostat': 0, 'tautp': 1.0, 'gamma_ln': 1.0, 'vlimit': 20.0,
                  'pressure': 1.0, 'barostat': 0, 'taup': 1.0,
                  'saveTrj': False, 'trajInterval': 1.0,
                  'integrator': 1, 'ensemType': 0,
                  'maxcyc': 100, 'ncyc': 20, 'dx0': 0.01, 'drms': 0.0001,
                  'simTime': 100, 'timeStep': 0.002,
                  'shake': 0, 'restraint': False, 'restraint_enum': 0, 'restraintmask': '', 'restraint_wt': 50,
                  'restraint_heavy': True, 'extraParams': ''
    }

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection('Minimization')
        form.addParam('AmberSystem', params.PointerParam, label="Input Amber System: ",
                      pointerClass='AmberSystem', help='Amber solvated system to be simulated')

        group = form.addGroup('Ensemble')
        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform in the step: Energy minimization, NVT or NPT\n')

        line = group.addLine('Temperature settings: ', condition='ensemType in [{}, {}]'.format(NVT, NPT),
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Time constant (ps) for heat bath coupling for the system\n'
                                  'The collision frequency γ, in ps −1 , when ntt = 3\n')
        line.addParam('temperature', params.FloatParam, default=300.0,
                      condition='ensemType in [{}, {}]'.format(NVT, NPT), label='Temperature: ')
        line.addParam('thermostat', params.EnumParam, default=5, condition='ensemType in [{}, {}]'.format(NVT, NPT),
                      label='Thermostat: ', choices=self._thermostats)
        line.addParam('tautp', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                      condition='ensemType in [{}, {}] and thermostat==0'.format(NVT, NPT))
        line.addParam('gamma_ln', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED,
                      condition='ensemType in [{}, {}] and thermostat==2'.format(NVT, NPT))
        group.addParam('vlimit', params.FloatParam, default=20.0, label='Velocity limit:',
                       condition='ensemType in [{}, {}]'.format(NVT, NPT), expertLevel=params.LEVEL_ADVANCED,
                       help='If not equal to 0.0, then any component of the velocity that is greater than abs(VLIMIT) '
                            'will be reduced to VLIMIT (preserving the sign)')

        line = group.addLine('Pressure settings: ', condition='ensemType=={}'.format(NPT),
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', params.FloatParam, default=1.0,
                      label='   Pressure (bar):   ')
        line.addParam('barostat', params.EnumParam, default=1,
                      label='  Barostat type:   ', choices=self._barostats)
        line.addParam('taup', params.FloatParam, default=1.0, expertLevel=params.LEVEL_ADVANCED)

        group = form.addGroup('Trajectory', condition='ensemType!=0',)
        group.addParam('saveTrj', params.BooleanParam, default=self._defParams['saveTrj'],
                       label="Save trajectory: ",
                       help='Save trajectory of the atoms during stage simulation.'
                            'The output will concatenate those trajectories which appear after the last stage '
                            'where the trajectory was not saved.')
        group.addParam('trajInterval', params.FloatParam, default=self._defParams['trajInterval'],
                       label='Interval time (ps):', condition='saveTrj',
                       help='Time between each frame recorded in the simulation (ps)')

        group = form.addGroup('Simulation time')
        # ENERGY MINIMIZATION
        group.addParam('integrator', params.EnumParam,
                       label='Simulation integrator: ', condition='ensemType=={}'.format(EMIN),
                       choices=self._integrators, default=1,
                       help='Type of integrator to use in simulation.')

        group.addParam('maxcyc', params.IntParam, default=100,
                       label='Minimization cycles:', condition='ensemType=={}'.format(EMIN),
                       help='Maximum number of minimization cycles')
        group.addParam('ncyc', params.IntParam, default=10, label='Steep cycles before CG:',
                       condition='ensemType=={} and integrator=={}'.format(EMIN, ST_CG),
                       help='For NCYC cycles the steepest descent method is used then conjugate gradient '
                            'is switched on')
        group.addParam('dx0', params.FloatParam, default=0.01, expertLevel=params.LEVEL_ADVANCED,
                       label='Length of minimization step [dx]:', condition='ensemType=={}'.format(EMIN),
                       help='The initial step length. If the initial step length is too big then will give a huge '
                            'energy; however the minimizer is smart enough to adjust itself')
        group.addParam('drms', params.FloatParam, default=0.0001, expertLevel=params.LEVEL_ADVANCED,
                       label='Convergence energy derivative (kcal/molkcal·Å):', condition='ensemType=={}'.format(EMIN),
                       help='The convergence criterion for the energy Derivative: minimization will halt when the '
                            'Root-Mean-Square of the Cartesian elements of the gradient of the energy is less '
                            'than this')

        # MOLECULAR DYNAMICS
        group.addParam('simTime', params.FloatParam, default=100,
                       label='Simulation time (ps):', condition='ensemType!={}'.format(EMIN),
                       help='Total time of the simulation stage (ps)')
        group.addParam('timeStep', params.FloatParam, default=0.002,
                       label='Simulation time steps (ps/step)[dt]:', condition='ensemType!={}'.format(EMIN),
                       help='Time of the steps for simulation (ps/step)[dt] \n 0.002ps is recommended if SHAKE '
                            'algorithm is chosen; if not, 0.001ps is recommended ')

        group = form.addGroup('Restraints')
        group.addParam('shake', params.EnumParam, default=0,
                       label='Use SHAKE constraints: ', choices=self._shakeAlgorithm,
                       help='Perform bond length constraints. should be used for most MD calculations. '
                            'The size of the MD timestep is determined by the fastest motions in the system. '
                            'SHAKE removes the bond stretching freedom, which is the fastest motion, and consequently '
                            'allows a larger timestep to be used')
        group.addParam('restraint', params.BooleanParam, default=False, label="Perform restraining: ",
                       help='Restraining specified atoms in Cartesian space using a harmonic potential')
        group.addParam('restraint_enum', params.EnumParam, default=0,
                       label='Restraint on: ', choices=self._restraints, condition='restraint',
                       help='Perform the restraints on the selected group of atoms')
        group.addParam('restraintmask', params.StringParam, default='',
                       label="Restrained atoms: ", condition='restraint and restraint_enum=={}'.format(CUST),
                       help='Restraining specified atoms in Cartesian space using a harmonic potential in Chimera '
                            'specifier: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/select.html')
        group.addParam('restraint_heavy', params.BooleanParam, default=True, condition='restraint',
                       label='Only heavy atoms: ', help='Restraint only heavy atoms (not hydrogens)')
        group.addParam('restraint_wt', params.FloatParam, default=50,
                       label="Restraint weight: ", condition='restraint',
                       help='(The weight for the positional restraints (kcal/mol/Å2)')

        group = form.addGroup('Extra parameters')
        group.addParam('extraParams', params.TextParam, width=120,
                       label='Extra params: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Included any extra parameters you want the input file to have for sander.\n'
                            'Be aware that you need to specify them with the proper syntax Sander expects')

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
        rstAmberFile, localTopFile = self._getPath('rstFile.rst7'), self._getPath('systemTopology.prmtop')

        shutil.copyfile(self.AmberSystem.get().getSystemFile(), rstAmberFile)
        shutil.copyfile(self.AmberSystem.get().getTopologyFile(), localTopFile)

        outTrjs = self.getTrjFiles()
        outputTrajectory = self._getPath('outputTrajectory.nc')
        self.combineTrajectories(outTrjs, outputTrajectory)

        outSystem = AmberSystem(filename=rstAmberFile)
        outSystem.setTopologyFile(localTopFile)
        if outTrjs:
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
                    sumStr += '{}) {}\n'.format(i+1, self.createSummaryLine(eval(dicLine)))
        else:
            sumStr += self.createSummaryLine(msjDic)
        return sumStr

    def createSummaryLine(self, msjDic):
        msjDic = self.addDefaultForMissing(msjDic)
        if msjDic['ensemType'] == EMIN:
            inte = self._integrators[msjDic['integrator']]
            steps = '{}-{}'.format(msjDic['ncyc'], msjDic['maxcyc']) if msjDic['integrator'] == ST_CG \
              else msjDic['maxcyc']
            lineStr = 'Minimization ({}): {} steps, {} convergence energy/dx'.format(inte, steps, msjDic['drms'])

        else:
            lineStr = 'MD simulation ({}): {} ps'.format(self._ensemTypes[msjDic['ensemType']], msjDic['simTime'])

        if msjDic['restraint']:
          lineStr += ', restraint on {}'.format(self._restraints[msjDic['restraint_enum']])
        lineStr += ', {} K'.format(msjDic['temperature'])
        if msjDic['saveTrj'] and not msjDic['ensemType'] == EMIN:
            lineStr += ', save traj'
        return lineStr

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

        params = f'Stage {mdpStage}\n &cntrl'

        if msjDic['ensemType'] == EMIN:
            # Energy minimization
            params += '\nimin=1, ntmin={}, '.format(msjDic['integrator'])
            if msjDic['integrator'] == 'Steep and CG':
                params += 'ncyc={}, '.format(msjDic['ncyc'])
            params += 'maxcyc={}, dx0={}, drms={}, '.format(msjDic['maxcyc'], msjDic['dx0'], msjDic['drms'])
        else:
            # MD simulation
            params += '\nimin=0, '
            params += 'dt={}, nstlim={}, '.format(msjDic['timeStep'], int(msjDic['simTime'] / msjDic['timeStep']))

            if msjDic['ensemType'] == NVE:
                params += '\nntt=0, '
            else:
                thermostat = self._map_therm[msjDic['thermostat']]
                params += '\nntt={}, temp0={}, '.format(thermostat, msjDic['temperature'])
                if thermostat == WCOUP:
                    params += 'tautp={}, '.format(msjDic['tautp'])
                elif thermostat == LANG:
                    params += 'gamma_ln={}, '.format(msjDic['gamma_ln'])

            if msjDic['ensemType'] == NPT:
                barostat = self._map_bar[msjDic['barostat']]
                params += '\nntp=1, barostat={}, pres0={}, taup={},'.\
                    format(barostat, msjDic['pressure'], msjDic['taup'])

        # OUTPUT SETTINGS
        if msjDic['saveTrj'] and not msjDic['ensemType'] == EMIN:
            stepsSave = int(msjDic['trajInterval'] / msjDic['timeStep'])
            params += '\nntwr={}, ntwx={}, ntwe={}, '.format(*[stepsSave]*3)
            pass

        params += '\nntc={},'.format(msjDic['shake']+1)
        if msjDic['restraint']:
            params += "\nntr=1, restraint_wt={}, restraintmask='{}', ".\
                format(msjDic['restraint_wt'], self.getRestraintMask(msjDic))

        if msjDic['extraParams']:
            params += '\n{}'.format(msjDic['extraParams'])
        params += '\n &end \n END'
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
        inFile = '{}.in'.format(stage)
        topFile = self.AmberSystem.get().getTopologyFile()
        rstFile = self.AmberSystem.get().getSystemFile()
        print(stageDir)

        if self.checkIfPrevTrj(stageNum):
            prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        else:
            prevTrjStr = ''

        command = '-i {} -c {} -p {} -r {}.r -o {}.o ' \
                  '-x {}.nc -e {}.e -ref {} -inf min.inf'.\
            format(inFile, amberFile, topFile, *[stage] * 4, rstFile)

        # Manage warnings
        nWarns = self.countWarns(stageNum)
        print('{} warnings in stage {}'.format(nWarns, stageNum))
        if nWarns >= 1:
            command += ' -maxwarn {}'.format(nWarns)

        amberPlugin.runAmbertools(self, 'sander -O ', command, cwd=stageDir)
        if not saveTrj:
            trjFile = os.path.join(stageDir, '{}.trr'.format(stage))
            os.remove(trjFile)

        return os.path.join(stageDir, inFile)

    def getRestraintMask(self, msjDic):
        inputSystem = self.AmberSystem.get()
        ionsResNames = inputSystem.getIonResNames()
        ligResName = inputSystem.getLigResNames()

        ionStr, nonProtStr = '', ''
        if ionsResNames:
          ionStr = ',' + ','.join(ionsResNames)
        if ionsResNames or ligResName:
          nonProtStr = ',' + ','.join(ionsResNames + ligResName)

        if msjDic['restraint_enum'] == SOL:
            restrStr = '!:WAT{}'.format(ionStr)

        elif msjDic['restraint_enum'] == PROT:
            restrStr = '!:WAT{}'.format(nonProtStr)

        elif msjDic['restraint_enum'] == BB:
            restrStr = '!:WAT{} & @n,ca,c,o'.format(nonProtStr)
        elif msjDic['restraint_enum'] == CA:
            restrStr = '!:WAT{} & @ca'.format(nonProtStr)
        elif msjDic['restraint_enum'] == CUST:
            restrStr = msjDic['restraintmask']

        if msjDic['restraint_heavy']:
            restrStr += ' & !@H='
        return restrStr


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

    def combineTrajectories(self, trajectoryFiles, outTrajectory):
        inputTop = self.AmberSystem.get().getTopologyFile()
        readTrajs = ['trajin {}'.format(trjFile) for trjFile in trajectoryFiles]
        combineStr = 'parm {}\n{}\nautoimage\ntrajout {}\ngo'.format(inputTop, '\n'.join(readTrajs), outTrajectory)
        cppTrajFile = self._getExtraPath('combineTraj.cpptraj')
        with open(cppTrajFile, 'w') as f:
            f.write(combineStr)

        amberPlugin.runAmbertools(self, 'cpptraj', '-i {}'.format(cppTrajFile))
        return outTrajectory
