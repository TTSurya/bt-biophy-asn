
def analyze_protein_trajectory(dna_coords, protein_coords, target_site):
    coord_to_indices = {}
    for idx, coord in enumerate(dna_coords):
        coord_tuple = tuple(coord)
        if coord_tuple not in coord_to_indices:
            coord_to_indices[coord_tuple] = set()
        coord_to_indices[coord_tuple].add(idx)
    target_coord = tuple(dna_coords[target_site])
    associations = 0
    dissociations = 0
    one_d_steps = 0
    three_d_steps = 0
    previous_dna_indices = set()
    previous_on_dna = False
    for protein_pos in protein_coords:
        pos_tuple = tuple(protein_pos)
        current_dna_indices = coord_to_indices.get(pos_tuple, set())
        current_on_dna = bool(current_dna_indices)
        if current_on_dna:
            if not previous_on_dna:
                associations += 1
            else:
                found_1d_slide = False
                for curr_idx in current_dna_indices:
                    for prev_idx in previous_dna_indices:
                        if abs(curr_idx - prev_idx) == 1:
                            one_d_steps += 1
                            found_1d_slide = True
                            break
                    if found_1d_slide:
                        break
                if not found_1d_slide and previous_dna_indices != current_dna_indices:
                    three_d_steps += 1
        else:
            if previous_on_dna:
                dissociations += 1
            three_d_steps += 1
        previous_dna_indices = current_dna_indices
        previous_on_dna = current_on_dna
        if pos_tuple == target_coord:
            print(f"Target reached at site with DNA coordinate {target_coord}!")
            break
    results = {
        "Number of Associations": associations,
        "Number of Dissociations": dissociations,
        "Number of 1D Steps (Sliding)": one_d_steps,
        "Number of 3D Steps (Jumping/Hopping)": three_d_steps
    }
    return results


# ---------------------------------------------------------
# ---------------------------------------------------------

# Only make changes in this section
# Put your dataset here, simply change the values in dna_coords and protein_coords
# In case of any confusion on how to put the right values, look for the dataset for EP22B050 at page number: 336
# For any other confusion, please ask on the group instead of personally messaging Lakshman.

dna_coords = [
    [0,3,2], [1,2,3], [2,1,4], [1,0,5], [0,1,4],
    [1,2,5], [2,3,4], [3,4,5], [2,5,4], [1,4,3],
    [0,5,2], [1,4,1], [2,5,0], [3,4,1], [4,3,2],
    [3,4,3], [4,5,4], [5,4,5], [4,3,4], [3,2,5]
]
protein_coords = [
    [4,5,2], [5,4,1], [4,3,2], [3,4,3], [4,5,4], 
    [3,4,3], [4,5,4], [2,3,4], [1,2,5], [2,3,4],
    [1,2,5], [1,2,5], [0,1,4], [1,2,5], [2,3,4],
    [2,1,4], [1,2,3], [2,1,4], [1,2,3], [0,3,2],
    [1,2,1], [0,3,0], [1,2,1], [2,1,0], [1,2,1], 
    [0,3,0], [0,5,2]
]
TARGET_SITE = 10 # This looks same for everyone

# ----------------------------------------------------------
# ----------------------------------------------------------

final_counts = analyze_protein_trajectory(dna_coords, protein_coords, TARGET_SITE)
print("\n--- Final Counts ---")
for event, count in final_counts.items():
    print(f"{event}: {count}")