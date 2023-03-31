# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

'''Script to describe a SetOfSmallMolecules docked using the Open Drug Discovery Toolkit
(ODDT, https://github.com/oddt/oddt). It must be launch with the conda environment with rdkit and oddt'''

import sys, os
import BioSimSpace as BSS

def parseParams(paramsFile):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split(':')
      paramsDic[key] = value.strip()
  return paramsDic

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])

    ligFile = paramsDic['ligandFile']
    ligName = os.path.basename(os.path.splitext(ligFile)[0])
    molecule = BSS.IO.readPDB(ligFile, pdb4amber=True).getMolecule(0)
    molecule = BSS.Parameters.parameterise(molecule, paramsDic['ff']).getMolecule()

    topFile, crdFile = BSS.IO.saveMolecules(ligName, molecule, ['PRM7', 'RST7'])


    with open(paramsDic['outputPath'], 'w') as f:
        f.write("Parametrised molecules:\n")
        f.write("Topology file:\t{}\nCoordinates file:\t{}\n".format(topFile, crdFile))
