# **************************************************************************
# *
# * Authors:     Aida Pinacho
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

import os, subprocess

import pwem

from pyworkflow import join

from amber.constants import AMBER_HOME, V2020, AMBER, AMBER_DEFAULT_VERSION

_logo = "icon.png"
_references = ['Salomon-Ferrer2013']


class Plugin(pwem.Plugin):
    _homeVar = AMBER_HOME
    _pathVars = [AMBER_HOME]
    _supportedVersions = [V2020]
    _amberName = 'ambertools-21'
    _pluginHome = join(pwem.Config.EM_ROOT, _amberName)
    _Ambertools21Env = "Ambertools21"

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(AMBER_HOME, cls._amberName)
        cls._defineVar("AMBERTOOLS_ENV_ACTIVATION", 'conda activate %s' % cls._Ambertools21Env)

    @classmethod
    def getAmbertoolsEnvActivation(cls):
        activation = cls.getVar("AMBERTOOLS_ENV_ACTIVATION")
        return activation

    @classmethod
    def defineBinaries(cls, env, default=False):
        # Creating a new conda enviroment for Ambertools21
        AMBER_INSTALLED = '%s_installed' % AMBER
        ambertools_commands = 'conda create -y -n %s && ' % cls._Ambertools21Env
        ambertools_commands += '%s %s && ' % (cls.getCondaActivationCmd(), cls.getAmbertoolsEnvActivation())
        ambertools_commands += 'conda install -y -c conda-forge ambertools=21 compilers && '
        ambertools_commands += 'touch {}'.format(AMBER_INSTALLED)  # Flag installation finished

        ambertools_commands = [(ambertools_commands, AMBER_INSTALLED)]
        env.addPackage('ambertools', version='21',
                       tar='void.tgz',
                       commands=ambertools_commands,
                       default=True)

    @classmethod
    def getPluginHome(cls, path=""):
        import amber
        fnDir = os.path.split(amber.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def getScriptsDir(cls, scriptName=''):
        return cls.getPluginHome('scripts/%s' % scriptName)

    @classmethod
    def getEnvActivation(cls, env):
        activation = cls.getVar("{}_ENV_ACTIVATION".format(env.upper()))
        return activation

    @classmethod
    def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False):
        """ Run modeller command from a given protocol. """
        scriptName = cls.getScriptsDir(scriptName)
        fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getEnvActivation(env), 'python', scriptName)
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(fullProgram + args, cwd=cwd, shell=True)

    @classmethod
    def runAmbertools(cls, protocol, program, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getAmbertoolsEnvActivation(), program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron_amber(), cwd=cwd)

    @classmethod
    def runCPPTRAJ(cls, protocol, cppTrajFile, cppTrajStr='', cwd=None):
        """ Run Ambertools command from a given protocol. """
        if cppTrajStr:
            with open(os.path.abspath(cppTrajFile), 'w') as f:
                f.write(cppTrajStr)

        args = '-i {}'.format(os.path.abspath(cppTrajFile))
        cls.runAmbertools(protocol, 'cpptraj', args, cwd)

    @classmethod
    def runAmberPrintf(cls, protocol, program, printfValues, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        AmberPath = join(cls._pluginHome, 'bin/{}'.format(program))
        program = 'printf "{}\n" | {}'.format('\n'.join(printfValues), AmberPath)
        protocol.runJob(program, args, cwd=cwd)



    @classmethod  # Test that
    def getEnviron_amber(cls):
        pass

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def runAmber(cls, self, param, params, cwd):
        pass
