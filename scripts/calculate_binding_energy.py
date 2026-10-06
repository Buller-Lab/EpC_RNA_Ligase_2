import pyrosetta
import mdtraj as mdt
import csv
import os
from pyrosetta import rosetta
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.select import *

# Define a function to calculate the total binding energy
def calculate_binding_energy(protein_pose, ligand_pose, complex_pose):
    # Score the entire complex pose
    scorefxn(complex_pose)
    total_energy = complex_pose.energies().total_energy()

    # Score the protein and ligand individually
    scorefxn(protein_pose)
    protein_energy = protein_pose.energies().total_energy()

    scorefxn(ligand_pose)
    ligand_energy = ligand_pose.energies().total_energy()

    # Calculate the binding energy
    binding_energy = total_energy - (protein_energy + ligand_energy)

    return binding_energy

# List of input PDB files
input_files = [f"productive_frames_replicate_{i}.pdb" for i in range(6)]  # Adjust range as needed

# Open a CSV file for appending binding energies
with open('binding_energies.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Input File", "Frame", "Binding Energy (kcal/mol)"])  # Write header

    for input_file in input_files:
        if not os.path.exists(input_file):
            print(f"Warning: File {input_file} does not exist. Skipping.")
            continue

        # Load the HDF5 file (trajectory) using MDTraj
        traj = mdt.load(input_file)  # Ensure the correct topology file if needed

        # Separate protein and ligand chains using MDTraj
        protein_chains = ['A', 'B']  # Adjust these chain IDs based on your structure
        ligand_chains = ['C', 'D', 'E']  # Adjust these chain IDs based on your structure

        # Extract protein and ligand based on chain IDs
        protein_residues = traj.topology.select('chainid 0 or chainid 1')  # Protein chains 'A' and 'B'
        ligand_residues = traj.topology.select('chainid 2 or chainid 3 or chainid 4')  # Ligand chains 'C', 'D', 'E'

        # Convert selection indices to integers
        protein_residues = protein_residues.astype(int)
        ligand_residues = ligand_residues.astype(int)

        # Iterate over each frame in the trajectory
        for frame_idx in range(traj.n_frames):
            # Extract protein and ligand from the frame
            complex_frame = traj.slice([frame_idx])
            protein_frame = traj.slice([frame_idx]).atom_slice(protein_residues)
            ligand_frame = traj.slice([frame_idx]).atom_slice(ligand_residues)

            # Save the protein and ligand frames as temporary PDB files
            complex_pdb = f'complex_frame_{frame_idx}.pdb'
            protein_pdb = f'protein_frame_{frame_idx}.pdb'
            ligand_pdb = f'ligand_frame_{frame_idx}.pdb'

            complex_frame.save(complex_pdb)
            protein_frame.save(protein_pdb)
            ligand_frame.save(ligand_pdb)

            # Initialize PyRosetta without the '-use_rna' option
            pyrosetta.init(extra_options=f'-ex1 -ex2aro -use_input_sc -native {complex_pdb} -add_orbitals')

            # Create score function (adjust this based on the exact scoring function you want)
            scorefxn = pyrosetta.create_score_function('rna_res_level_energy4.wts')

            # Load the frames into PyRosetta poses
            protein_pose = pyrosetta.pose_from_pdb(protein_pdb)
            ligand_pose = pyrosetta.pose_from_pdb(ligand_pdb)

            # Combine the protein and ligand poses into a complex pose
            complex_pose = protein_pose.clone()  # Copy the protein pose
            complex_pose.append_pose_by_jump(ligand_pose, 1)  # Append ligand pose to protein pose

            # Calculate the binding energy for this frame
            binding_energy = calculate_binding_energy(protein_pose, ligand_pose, complex_pose)

            # Print the binding energy
            print(f"File: {input_file}, Frame {frame_idx}: {binding_energy} kcal/mol")

            # Append the binding energy to the CSV file
            writer.writerow([input_file, frame_idx, binding_energy])  # Save the file name, frame, and binding energy

            # Delete the temporary PDB files
            os.remove(complex_pdb)
            os.remove(protein_pdb)
            os.remove(ligand_pdb)
