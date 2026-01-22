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

import pwem

from pyworkflow import join

from amber.constants import *
from pwem.convert.atom_struct import getEnviron

import pwchem

from scipion.install.funcs import InstallHelper

_logo = "icon.png"
_references = ['Salomon-Ferrer2013']

AMBER_DIC = {'name': 'amber', 'version': AMBER_DEFAULT_VERSION, 'home': 'AMBER_HOME', 'pmemd_home': 'PMEMD_HOME'}

class Plugin(pwchem.Plugin):
    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(AMBER_DIC['home'], cls.getEnvName(AMBER_DIC))
        cls._defineVar("AMBERTOOLS_ENV_ACTIVATION", 'conda activate %s' % cls.getEnvName(AMBER_DIC))
        cls._defineEmVar(AMBER_DIC['pmemd_home'], 'pmemd24')


    @classmethod
    def getAmbertoolsEnvActivation(cls):
        activation = cls.getVar("AMBERTOOLS_ENV_ACTIVATION")
        return activation

    @classmethod
    def defineBinaries(cls, env, default=False):
        # Creating a new conda enviroment for Ambertools25
        AMBER_INSTALLED = '%s_%s_installed' % (AMBER, V2025)
        ambertools_commands = 'conda create -y -n %s && ' % cls.getEnvName(AMBER_DIC)
        ambertools_commands += '%s %s && ' % (cls.getCondaActivationCmd(), cls.getAmbertoolsEnvActivation())
        ambertools_commands += 'conda install -y -c dacase -c conda-forge ambertools-dac=25 compilers && '
        ambertools_commands += 'touch {}'.format(AMBER_INSTALLED)  # Flag installation finished

        ambertools_commands = [(ambertools_commands, AMBER_INSTALLED)]
        env.addPackage(AMBER_DIC['name'], version=AMBER_DIC['version'],
                       tar='void.tgz',
                       commands=ambertools_commands,
                       default=True)

        installer = InstallHelper(AMBER_DIC['name'], packageHome=cls.getVar(AMBER_DIC['home']),
                                  packageVersion=AMBER_DIC['version'])

        AMBER_INSTALLED = f'{AMBER}_{V2025}_installed'

        installer.addCommand(f'conda create -y -n {cls.getEnvName(AMBER_DIC)}','AMBER_ENV_CREATED'
            ).addCommand(f'{cls.getCondaActivationCmd()} {cls.getAmbertoolsEnvActivation()} && '
            f'conda install -y -c dacase -c conda-forge ambertools-dac=25 compilers',
            'AMBERTOOLS_INSTALLED'
            ).addCommand(f'touch {AMBER_INSTALLED}',AMBER_INSTALLED
            ).addPackage(env, dependencies=['conda'], default=True)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def runAmbertools(cls, protocol, program, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        fullProgram = f' {cls.getEnvActivationCommand(AMBER_DIC)} && {program}'
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runPmemd(cls, protocol, program, args, gpu=False, cwd=None):
        fullProgram = f' {cls.getEnvActivationCommand(AMBER_DIC)} && {cls.getPmemdBin()}'
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def getPmemdBin(cls, prog='pmemd.cuda '):
        """ Retorna la ruta al binario pmemd configurado. """
        return join(cls.getVar(AMBER_DIC['pmemd_home']), 'bin', prog)
