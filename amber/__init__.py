# **************************************************************************
# *
# * Authors:     Aida Pinacho (you@yourinstitution.email)
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

from amber.constants import AMBER_HOME, V2020, AMBER, AMBER_DEFAULT_VERSION

_logo = "icon.png"
_references = ['you2019']


class Plugin(pwem.Plugin):
    _homeVar = AMBER_HOME
    _pathVars = [AMBER_HOME]
    _supportedVersions = [V2020]
    _amberName = AMBER + '-' + AMBER_DEFAULT_VERSION
    _pluginHome = join(pwem.Config.EM_ROOT, _amberName)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(AMBER_HOME, cls._amberName)

    @classmethod
    def defineBinaries(cls, env):
        cMakeCmd = 'mkdir build && cd build && '
        cMakeCmd += 'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON ' \
                           '-DCMAKE_INSTALL_PREFIX={} > cMake.log'.format(cls._pluginHome)
        makeCmd = 'cd build && make -j {} > make.log && make check'.format(env.getProcessors())
        makeInstallCmd = 'cd build && make install'

        # Creating validation file
        AMBER_INSTALLED = '%s_installed' % AMBER
        installationCmd = 'touch %s' % AMBER_INSTALLED  # Flag installation finished

        env.addPackage(AMBER,
                       version=AMBER_DEFAULT_VERSION,
                       url=cls._getAmberDownloadUrl(),
                       commands=[(cMakeCmd, []),
                                 (makeCmd, []),
                                 (makeInstallCmd, []),
                                 (installationCmd, AMBER_INSTALLED)],
                       default=True)

    @classmethod
    def runAmbertools(cls, protocol, program, args, cwd=None):
        """ Run Ambertools command from a given protocol. """
        protocol.runJob(join(cls._pluginHome, 'bin/{}'.format(program)), args, cwd=cwd)

    @classmethod
    def runAmberPrintf(cls, protocol, program, printfValues, args, cwd=None):
      """ Run Ambertools command from a given protocol. """
      AmberPath = join(cls._pluginHome, 'bin/{}'.format(program))
      program = 'printf "{}\n" | {}'.format('\n'.join(printfValues), AmberPath)
      protocol.runJob(program, args, cwd=cwd)

    @classmethod  # Test that
    def getEnviron(cls):
        pass

    @classmethod
    def _getAmberDownloadUrl(cls):
        return 'https://ambermd.org/GetAmber.php{}.tar.gz'.format(AMBER_DEFAULT_VERSION)

    @classmethod
    def _getAmberTar(cls):
        return cls._pluginHome + '/' + cls._amberName + '.tar.gz'

    # ---------------------------------- Utils functions  -----------------------





