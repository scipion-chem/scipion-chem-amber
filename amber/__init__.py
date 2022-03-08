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
import pwchem
import pyworkflow

import pyworkflow


_logo = "icon.png"
_references = ['you2019']


class Plugin(pwchem.Plugin):
    _homeVar = AMBER_HOME
    _pathVars = [AMBER_HOME]
    _supportedVersions = [V2020]
    _gromacsName = AMBER + '-' + AMBER_DEFAULT_VERSION
    _pluginHome = join(pwem.Config.EM_ROOT, _amberName)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(AMBER_HOME, cls._amberName)




