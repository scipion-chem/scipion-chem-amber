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
AMBER_DEFAULT_VERSION = V2020


#  Ions from Ambertools21/dat/leap/lib/atomic_ions.lib
_cations = {'Ag+': 'AG', 'Ag2+': 'Ag', 'Al3+': 'AL', 'Ba2+': 'BA', 'Be2+': 'Be', 'Ca2+': 'CA', 'Cd2+': 'CD',
            'Ce3+': 'CE', 'Ce4+': 'Ce', 'Co2+': 'CO', 'Cr2+': 'Cr', 'Cr3+': 'CR', 'Cs+': 'CS', 'Cu+': 'CU',
            'Cu2+': 'CU', 'Dy3+': 'Dy', 'Er3+': 'Er', 'Eu2+': 'EU', 'Eu3+': 'EU3', 'Fe2+': 'FE2', 'Fe3+': 'FE',
            'Gd3+': 'GD', 'H3O+': 'O', 'HE+': 'H', 'HZ+': 'H', 'Hf4+': 'Hf', 'Hg2+': 'HG', 'In3+': 'IN', 'K+': 'K',
            'La3+': 'LA', 'Li+': 'LI', 'Lu3+': 'LU', 'Mg2+': 'MG', 'Mn2+': 'MN', 'NH4+': 'N', 'Na+': 'Na+',
            'Nd3+': 'Nd', 'Ni2+': 'NI', 'Pb2+': 'PB', 'Pd2+': 'PD', 'Pr3+': 'PR', 'Pt2+': 'PT', 'Pu4+': 'Pu',
            'Ra2+': 'Ra', 'Rb+': 'RB', 'Sm2+': 'Sm', 'Sm3+': 'SM', 'Sn2+': 'Sn', 'Sr2+': 'SR', 'Tb3+': 'TB',
            'Th4+': 'Th', 'Tl+': 'TL', 'Tl3+': 'Tl', 'Tm3+': 'Tm', 'U4+': 'U', 'V2+': 'V2+', 'Y3+': 'Y', 'Yb2+': 'YB2',
            'Zn2+': 'ZN', 'Zr4+': 'Zr'}
_anions = {'Br-': 'BR', 'Cl-': 'Cl-', 'F-': 'F', 'I-': 'I'}

PROT_RESNAMES = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'CYX',  'MET',
                 'PHE', 'TYR', 'TRP', 'ASP', 'GLU', 'LYS', 'ARG', 'HID', 'HIE', 'HIP']

