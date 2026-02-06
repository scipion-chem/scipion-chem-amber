#!/usr/bin/env python
import sys
import math
import os.path
import argparse
import numpy as np

# * Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# Adapted script from https://github.com/callumjd/AMBER-Membrane_protein_tutorial
################################################################################
#
# alignMembrane.py - Align membrane systems using protein structure
# Run as:
# ./alignMembrane.py -i original_protein.pdb -m bilayer.pdb -o output_name.pdb
#
################################################################################

class Coord(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def atom_mass(str):
    weight = {'C': '12.01', 'H': '1.008', 'P': '30.97', 'N': '14.01',
              'O': '16.0', 'Cl': '35.45', 'F': '19.0', 'S': '32.065'}
    if str[0].isalpha() == True:
        if str[1].isalpha() == True:
            try:
                return weight[str[0:2]]
            except:
                return weight[str[0]]
        else:
            return weight[str[0]]
    elif str[0].isalpha() == False:
        if str[1].isalpha() == True and str[2].isalpha() == True:
            try:
                return weight[str[1:3]]
            except:
                return weight[str[1]]
        else:
            return weight[str[1]]


def get_protein_coords(file_in):
    """Extract protein atom coordinates (CA atoms for alignment)"""
    coords = []
    disallow = ['MEMB', 'TIP3', 'SOD', 'POT', 'CLA']

    with open(file_in, 'r') as f:
        for line in f:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                if line.split()[-1] not in disallow:
                    # Use CA atoms for alignment (or all atoms if you prefer)
                    atom_name = line[12:16].strip()
                    if atom_name == 'CA':  # Change to 'if True:' to use all protein atoms
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])

    return np.array(coords)


def kabsch_alignment(P, Q):
    """
    Calculate optimal rotation matrix using Kabsch algorithm
    P: reference structure (Nx3)
    Q: structure to align (Nx3)
    Returns: rotation matrix (3x3) and translation vector
    """
    # Center both structures
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Calculate covariance matrix
    H = Q_centered.T @ P_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Calculate rotation matrix
    # Check for reflection
    d = np.linalg.det(Vt.T @ U.T)
    if d < 0:
        Vt[-1, :] *= -1

    R = Vt.T @ U.T

    # Calculate translation
    t = centroid_P - R @ centroid_Q

    return R, t, centroid_Q


def apply_transformation(coords, R, t):
    """Apply rotation R and translation t to coordinates"""
    return (R @ coords.T).T + t


def calculate_rmsd(P, Q):
    """Calculate RMSD between two coordinate sets"""
    return np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1)))


def get_wat_size(file_in):
    """Get water box dimensions"""
    wat_xyz = []
    with open(file_in, 'r') as f_in:
        for line in f_in:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                if line.split()[-1] == 'TIP3':
                    x = float(line.split()[-6])
                    y = float(line.split()[-5])
                    z = float(line.split()[-4])
                    wat_xyz.append((x, y, z))

    if len(wat_xyz) > 0:
        wat_xyz_np = np.array(wat_xyz)
        box_x = abs(np.min(wat_xyz_np[:, 0])) + abs(np.max(wat_xyz_np[:, 0]))
        box_y = abs(np.min(wat_xyz_np[:, 1])) + abs(np.max(wat_xyz_np[:, 1]))
        box_z = abs(np.min(wat_xyz_np[:, 2])) + abs(np.max(wat_xyz_np[:, 2]))
        return Coord(box_x, box_y, box_z)
    else:
        return Coord(0, 0, 0)


################################################################################
#### Main program ####
################################################################################

parser = argparse.ArgumentParser(description='Align membrane system to reference protein structure')
parser.add_argument("-i", type=str, help="Reference protein PDB file", required=True)
parser.add_argument("-m", type=str, help="Protein+membrane PDB to align", required=True)
parser.add_argument("-o", type=str, help="Output aligned PDB", required=True)
parser.add_argument("--use-all-atoms", action='store_true',
                    help="Use all protein atoms instead of just CA (default: CA only)")

args = parser.parse_args()

# Check files exist
if not (os.path.isfile(args.i) and os.path.isfile(args.m)):
    print("Error: Input files not found")
    sys.exit()

prot_file = args.i
membrane_file = args.m
output_file = args.o

################################################################################
#### Calculate alignment transformation ####
################################################################################

print("Reading reference protein from:", prot_file)
print("Reading membrane system from:", membrane_file)
print()

# Get protein coordinates
ref_coords = get_protein_coords(prot_file)
mobile_coords = get_protein_coords(membrane_file)

print(f"Found {len(ref_coords)} alignment atoms in reference")
print(f"Found {len(mobile_coords)} alignment atoms in membrane system")
print()

if len(ref_coords) != len(mobile_coords):
    print("WARNING: Different number of atoms! Alignment may fail.")
    print("Make sure both files contain the same protein structure.")
    min_len = min(len(ref_coords), len(mobile_coords))
    ref_coords = ref_coords[:min_len]
    mobile_coords = mobile_coords[:min_len]

# Calculate optimal rotation and translation
R, t, centroid_mobile = kabsch_alignment(ref_coords, mobile_coords)

print("Rotation matrix:")
print(R)
print()
print("Translation vector:", t)
print()

# Calculate RMSD before and after
rmsd_before = calculate_rmsd(ref_coords, mobile_coords)
aligned_coords = apply_transformation(mobile_coords, R, t)
rmsd_after = calculate_rmsd(ref_coords, aligned_coords)

print(f"RMSD before alignment: {rmsd_before:.3f} Å")
print(f"RMSD after alignment:  {rmsd_after:.3f} Å")
print()

################################################################################
#### Apply transformation and write output ####
################################################################################

print(f"Writing aligned system to: {output_file}")

disallow = ['MEMB', 'TIP3', 'SOD', 'POT', 'CLA']

with open(output_file, 'w') as f_out:
    with open(membrane_file, 'r') as f_in:
        for line in f_in:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                # Only process membrane, water, and ion atoms (not protein)
                if line.split()[-1] in disallow:
                    # Extract coordinates
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # Apply transformation
                    coord = np.array([x, y, z])
                    new_coord = apply_transformation(coord.reshape(1, -1), R, t)[0]

                    # Write transformed coordinates
                    new_line = (f"{line[0:30]}{new_coord[0]:8.3f}{new_coord[1]:8.3f}"
                                f"{new_coord[2]:8.3f}{line[54:]}")
                    f_out.write(new_line)
                # Skip protein atoms - don't write them
            elif line.split()[0] == 'TER':
                f_out.write(line)
            elif line.split()[0] == 'END':
                f_out.write(line)

box_dimensions = get_wat_size(membrane_file)
if box_dimensions.x > 0:
    print(f'\nOriginal box X, Y, Z: {box_dimensions.x:.3f} {box_dimensions.y:.3f} {box_dimensions.z:.3f}')

print("\nAlignment complete!")