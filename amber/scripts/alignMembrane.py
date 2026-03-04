#!/usr/bin/env python
import sys
import os.path
import argparse
import numpy as np


# Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# Adapted script from https://github.com/callumjd/AMBER-Membrane_protein_tutorial
################################################################################
# Run as:
# ./align_protein_to_membrane.py -i protein.pdb -m membrane_system.pdb -o aligned_protein.pdb
#
# This aligns your protein to match the orientation of the protein in the
# packmol-memgen generated membrane system (where Z-axis is transmembrane).
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


def get_protein_coords(file_in, use_all=False):
    """Extract protein atom coordinates (CA atoms by default for alignment)"""
    coords = []
    disallow = ['MEMB', 'TIP3', 'SOD', 'POT', 'CLA']

    with open(file_in, 'r') as f:
        for line in f:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                if line.split()[-1] not in disallow:
                    atom_name = line[12:16].strip()
                    if use_all or atom_name == 'CA':
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])

    return np.array(coords)


def read_mol2_coords(file_in):
    """Read coordinates from mol2 file"""
    coords = []
    in_atom_section = False

    with open(file_in, 'r') as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atom_section = True
                continue
            elif '@<TRIPOS>' in line and in_atom_section:
                break
            elif in_atom_section:
                parts = line.split()
                if len(parts) >= 5:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    coords.append([x, y, z])

    return np.array(coords)


def write_mol2_transformed(input_file, output_file, R, t):
    """Write mol2 file with transformed coordinates"""
    with open(output_file, 'w') as f_out:
        with open(input_file, 'r') as f_in:
            in_atom_section = False

            for line in f_in:
                if '@<TRIPOS>ATOM' in line:
                    in_atom_section = True
                    f_out.write(line)
                    continue
                elif '@<TRIPOS>' in line and in_atom_section:
                    in_atom_section = False
                    f_out.write(line)
                    continue

                if in_atom_section:
                    parts = line.split()
                    if len(parts) >= 5:
                        # Extract and transform coordinates
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        coord = np.array([x, y, z])
                        new_coord = apply_transformation(coord.reshape(1, -1), R, t)[0]

                        # Reconstruct line with new coordinates
                        new_line = f"{parts[0]:>7} {parts[1]:<8} {new_coord[0]:>9.4f} {new_coord[1]:>9.4f} {new_coord[2]:>9.4f}"
                        if len(parts) > 5:
                            new_line += " " + " ".join(parts[5:])
                        new_line += "\n"
                        f_out.write(new_line)
                    else:
                        f_out.write(line)
                else:
                    f_out.write(line)


def read_sdf_coords(file_in):
    """Read coordinates from SDF file"""
    coords = []

    with open(file_in, 'r') as f:
        lines = f.readlines()

        # Find counts line (4th line typically)
        if len(lines) < 4:
            return np.array(coords)

        counts_line = lines[3].split()
        if len(counts_line) < 2:
            return np.array(coords)

        try:
            n_atoms = int(counts_line[0])
        except:
            return np.array(coords)

        # Read atom coordinates (start at line 4)
        for i in range(4, min(4 + n_atoms, len(lines))):
            parts = lines[i].split()
            if len(parts) >= 3:
                try:
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                    coords.append([x, y, z])
                except:
                    continue

    return np.array(coords)


def write_sdf_transformed(input_file, output_file, R, t):
    """Write SDF file with transformed coordinates"""
    with open(output_file, 'w') as f_out:
        with open(input_file, 'r') as f_in:
            lines = f_in.readlines()

            if len(lines) < 4:
                f_out.writelines(lines)
                return

            # Write header lines
            for i in range(4):
                f_out.write(lines[i])

            # Get number of atoms
            counts_line = lines[3].split()
            try:
                n_atoms = int(counts_line[0])
            except:
                f_out.writelines(lines[4:])
                return

            # Transform and write atom coordinates
            for i in range(4, min(4 + n_atoms, len(lines))):
                parts = lines[i].split()
                if len(parts) >= 3:
                    try:
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        coord = np.array([x, y, z])
                        new_coord = apply_transformation(coord.reshape(1, -1), R, t)[0]

                        # Reconstruct line with new coordinates
                        new_line = f"{new_coord[0]:10.4f}{new_coord[1]:10.4f}{new_coord[2]:10.4f}"
                        if len(parts) > 3:
                            new_line += " " + " ".join(parts[3:])
                        new_line += "\n"
                        f_out.write(new_line)
                    except:
                        f_out.write(lines[i])
                else:
                    f_out.write(lines[i])

            # Write remaining lines
            f_out.writelines(lines[4 + n_atoms:])


def kabsch_alignment(P, Q):
    """
    Calculate optimal rotation matrix using Kabsch algorithm
    P: target/reference structure (Nx3) - membrane system protein
    Q: structure to align (Nx3) - your protein to move
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

parser = argparse.ArgumentParser(description='Align protein to membrane system orientation')
parser.add_argument("-i", type=str, help="Input protein PDB file to align", required=True)
parser.add_argument("-m", type=str, help="Reference membrane system PDB (with correct orientation)", required=True)
parser.add_argument("-o", type=str, help="Output aligned protein PDB", required=True)
parser.add_argument("-l", type=str, help="Optional ligand file (.mol2 or .sdf) to align", default=None)
parser.add_argument("--ligand-out", type=str, help="Output aligned ligand file (same format as input)", default=None)
parser.add_argument("--use-all-atoms", action='store_true',
                    help="Use all protein atoms instead of just CA (default: CA only)")
parser.add_argument("--save-all", action='store_true',
                    help="Save all non-membrane atoms (protein + ligands/cofactors)")

args = parser.parse_args()

# Check files exist
if not (os.path.isfile(args.i) and os.path.isfile(args.m)):
    print("Error: Input files not found")
    sys.exit()

protein_file = args.i
membrane_file = args.m
output_file = args.o

################################################################################
#### Calculate alignment transformation ####
################################################################################

print("Reading protein to align from:", protein_file)
print("Reading reference membrane system from:", membrane_file)
print()

# Get protein coordinates
target_coords = get_protein_coords(membrane_file, args.use_all_atoms)
mobile_coords = get_protein_coords(protein_file, args.use_all_atoms)

print(f"Found {len(target_coords)} alignment atoms in membrane system")
print(f"Found {len(mobile_coords)} alignment atoms in protein to align")
print()

if len(target_coords) != len(mobile_coords):
    print("WARNING: Different number of atoms! Alignment may fail.")
    print("Make sure both files contain the same protein structure.")
    min_len = min(len(target_coords), len(mobile_coords))
    target_coords = target_coords[:min_len]
    mobile_coords = mobile_coords[:min_len]

# Calculate optimal rotation and translation
R, t, centroid_mobile = kabsch_alignment(target_coords, mobile_coords)

print("Rotation matrix:")
print(R)
print()
print("Translation vector:", t)
print()

# Calculate RMSD before and after
rmsd_before = calculate_rmsd(target_coords, mobile_coords)
aligned_coords = apply_transformation(mobile_coords, R, t)
rmsd_after = calculate_rmsd(target_coords, aligned_coords)

print(f"RMSD before alignment: {rmsd_before:.3f} Å")
print(f"RMSD after alignment:  {rmsd_after:.3f} Å")
print()

################################################################################
#### Apply transformation and write output ####
################################################################################

print(f"Writing aligned protein to: {output_file}")

disallow = ['MEMB', 'TIP3', 'SOD', 'POT', 'CLA']

with open(output_file, 'w') as f_out:
    with open(protein_file, 'r') as f_in:
        for line in f_in:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                # Check if this is a protein/ligand atom (not membrane/water/ion)
                if line.split()[-1] not in disallow:
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
                # Skip membrane/water/ion atoms if present
            elif line.split()[0] == 'TER':
                f_out.write(line)
            elif line.split()[0] == 'END':
                f_out.write(line)

# Always write aligned membrane+waters to membrane.pdb
output_dir = os.path.dirname(output_file) or '.'
membrane_output = os.path.join(output_dir, 'membrane.pdb')

print(f"Writing membrane+waters (original orientation) to: {membrane_output}")

with open(membrane_output, 'w') as f_out:
    with open(membrane_file, 'r') as f_in:
        for line in f_in:
            if line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM':
                # Only save membrane, water, and ion atoms (no transformation)
                if line.split()[-1] in disallow:
                    f_out.write(line)
                # Skip protein atoms
            elif line.split()[0] == 'TER':
                f_out.write(line)
            elif line.split()[0] == 'END':
                f_out.write(line)

box_dimensions = get_wat_size(membrane_file)
if box_dimensions.x > 0:
    print(f'\nMembrane system box X, Y, Z: {box_dimensions.x:.3f} {box_dimensions.y:.3f} {box_dimensions.z:.3f}')

print("\nAlignment complete!")
print(f"Output files created:")
print(f"  - Aligned protein: {output_file}")
print(f"  - Aligned membrane+waters: {membrane_output}")

# Process ligand if provided
if args.l and os.path.isfile(args.l):
    ligand_file = args.l
    ligand_ext = os.path.splitext(ligand_file)[1].lower()

    # Determine output file name
    if args.ligand_out:
        ligand_output = args.ligand_out
    else:
        # Auto-generate output name
        base_name = os.path.splitext(ligand_file)[0]
        ligand_output = f"{base_name}_aligned{ligand_ext}"

    print(f"\nProcessing ligand file: {ligand_file}")

    try:
        if ligand_ext == '.mol2':
            write_mol2_transformed(ligand_file, ligand_output, R, t)
            print(f"Aligned ligand saved to: {ligand_output}")
        elif ligand_ext == '.sdf':
            write_sdf_transformed(ligand_file, ligand_output, R, t)
            print(f"Aligned ligand saved to: {ligand_output}")
        else:
            print(f"WARNING: Unsupported ligand format '{ligand_ext}'. Only .mol2 and .sdf are supported.")
    except Exception as e:
        print(f"ERROR processing ligand file: {e}")
elif args.l and not os.path.isfile(args.l):
    print(f"\nWARNING: Ligand file specified but not found: {args.l}")
