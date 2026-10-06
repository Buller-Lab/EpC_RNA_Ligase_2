import mdtraj as md
import numpy as np
import pandas as pd
import pyemma

# List to store all productive frames
data_records = []

# Loop through the trajectory files
for i in range(5):
    traj = md.load(f'production_{i}.h5')
    topology = traj.topology
    
    # Initialize selections
    o5_selection = []
    o3_selection = []
    
    # Select O5' atom from residue A13 (resSeq 13)
    for atom in topology.atoms:
        if atom.residue.resSeq == 13 and atom.name == "O5'":
            o5_selection = [atom.index]
    
    # Select O3' atom from residue U22 (resSeq 22)
    for atom in topology.atoms:
        if atom.residue.resSeq == 22 and atom.name == "O3'":
            o3_selection = [atom.index]
    
    # Define the Mg atom selections
    mg_selection = topology.select("resname MG")
    
    # Print selections for verification
    print(f"\nProcessing trajectory production_{i}.h5")
    print("Mg atoms (resname MG):")
    for idx in mg_selection:
        atom = topology.atom(idx)
        print(f"  - {atom.name} (index: {idx}, residue: {atom.residue})")
    
    print("O5 atom (resname A13, name O5'):")
    if not o5_selection:
        print(f"  - Warning: No O5' atom found in resname A13")
    for idx in o5_selection:
        atom = topology.atom(idx)
        print(f"  - {atom.name} (index: {idx}, residue: {atom.residue})")
    
    print("O3 atom (resname U22, name O3'):")
    if not o3_selection:
        print(f"  - Warning: No O3' atom found in resname U22")
    for idx in o3_selection:
        atom = topology.atom(idx)
        print(f"  - {atom.name} (index: {idx}, residue: {atom.residue})")
    
    # Check if selections are valid
    if not o5_selection:
        print(f"Warning: Skipping trajectory {i} due to missing O5' atom in A13")
        continue
    if not o3_selection:
        print(f"Warning: Skipping trajectory {i} due to missing O3' atom in U22")
        continue
    if len(o5_selection) > 1:
        print(f"Warning: Multiple O5' atoms found in A13: {o5_selection}. Using first one.")
        o5_selection = [o5_selection[0]]
    if len(o3_selection) > 1:
        print(f"Warning: Multiple O3' atoms found in U22: {o3_selection}. Using first one.")
        o3_selection = [o3_selection[0]]
    
    # Define atom pairs for distance calculations
    atom_pairs = [[mg, o5] for mg in mg_selection for o5 in o5_selection] + \
                 [[mg, o3] for mg in mg_selection for o3 in o3_selection]
    o5_o3_pair = [[o5, o3] for o5 in o5_selection for o3 in o3_selection][0]
    
    if not atom_pairs:
        print(f"Skipping trajectory {i} due to empty atom pairs")
        continue
    
    # Compute distances with periodic boundary conditions
    distances = md.compute_distances(traj, atom_pairs, periodic=True)
    o5_o3_distances = md.compute_distances(traj, [o5_o3_pair], periodic=True)[:, 0]
    
    # Split distances into Mg-O5' and Mg-O3' pairs
    num_mg_o5_pairs = len(mg_selection) * len(o5_selection)
    distances_mg_o5 = distances[:, :num_mg_o5_pairs]
    distances_mg_o3 = distances[:, num_mg_o5_pairs:]
    
    # Compute mean distance per frame
    mean_distances = (np.mean(distances_mg_o5, axis=1) + np.mean(distances_mg_o3, axis=1)) / 2
    
    # Select top 100 frames with lowest mean distances
    top_indices = np.argsort(mean_distances)[:100]
    
    # Extract the top 100 frames
    top_frames = traj[top_indices]
    molecules = top_frames.topology.find_molecules()
    anchor_molecules = [molecules[0]]
    top_frames = top_frames.image_molecules(inplace=False, anchor_molecules=anchor_molecules, make_whole=True)
    # Select atoms for alignment (using all protein heavy atoms - modify as needed)
    align_indices = top_frames.topology.select('protein and not type H')
    
    # Superpose all frames to the first frame
    top_frames.superpose(top_frames, frame=0, atom_indices=align_indices)
  
    # Save the aligned trajectory, overwriting the original file
    top_frames.save_pdb(f'productive_frames_replicate_{i + 1}.pdb')

    # Reload the PDB file
    aligned_traj = md.load(f'productive_frames_replicate_{i + 1}.pdb')

    # Perform clustering using pyemma
    cluster_model = pyemma.coordinates.cluster_kmeans(aligned_traj.xyz.reshape(len(aligned_traj), -1), k=1)
    cluster_centroid_idx = np.argmin(np.sum(cluster_model.dtrajs[0][:, None] - cluster_model.dtrajs[0], axis=1))
    centroid_traj = aligned_traj[cluster_centroid_idx]

    # Save the centroid PDB
    centroid_traj.save_pdb(f'productive_cluster_replicate_{i + 1}_centroid.pdb')
    
    # Store frame data
    for frame_idx in top_indices:
        data_records.append([f'production_{i}', frame_idx, mean_distances[frame_idx] * 10, o5_o3_distances[frame_idx] * 10])

# Create and save DataFrame
df = pd.DataFrame(data_records, columns=['trajectory', 'frame', 'mean_distance (Å)', 'distance_O5_O3 (Å)'])
df.to_excel('top_100_productive_data.xlsx', index=False)