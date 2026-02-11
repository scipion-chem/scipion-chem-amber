# **************************************************************************
# *
# * Authors: Aida Pinacho Pérez
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


AMBER_HOME = 'AMBER_HOME'

AMBER = 'amber'
V2020 = '2020.1'
V2025 = '2025.1'

AMBER_DEFAULT_VERSION = V2025

VMD_HOME = 'VMD_HOME'

AMBER_DIC = {'name': 'amber', 'version': AMBER_DEFAULT_VERSION, 'home': 'AMBER_HOME', 'pmemd_home': 'PMEMD_HOME'}

PROTEIN_RES = "ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL," \
              "HID,HIE,HIP,CYX,ASH,GLH,LYN,ARN,ACE,NME,NHE"
ENV_RES = "WAT,HOH,TIP3,SPC,SPCE,Na+,Cl-,K+,Cs+,Rb+,Li+,Mg+,Ca2+,Zn2+"
RESTRAINS_DIC = {'Protein + Ligand': f':{PROTEIN_RES},LIG & !@H=', 'Protein only': f':{PROTEIN_RES} & !@H=', 'Ligand': ':LIG',
                 'Backbone': '@CA,C,N,O', 'CA': '@CA', 'Everything except wat+ions': '!(:WAT,HOH,TIP3,Na+,Cl-,K+)'}