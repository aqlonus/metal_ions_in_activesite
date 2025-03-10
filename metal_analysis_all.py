"""
Complete Metal Ion Analysis Pipeline
----------------------------------
This script performs a comprehensive analysis of metal ions in protein structures:
1. Identifies metal ions and their neighbors in PDB files
2. Categorizes results based on residue numbers
3. Performs statistical analysis of distances
"""

from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import pandas as pd
import os
import shutil

# Configuration
METALS = ["MG", "CO", "NI", "ZN", "MN", "CA", "FE", "CU", "CD"]
DISTANCE_THRESHOLD = 2.8
BASE_PATH = "C:/Users/Agnieszka/Desktop/komputer_work/analiza_SAXS/XI/XI_metal_analysis"
PDB_PATH = os.path.join(BASE_PATH, "pdb_files")
RESULTS_PATH = os.path.join(PDB_PATH, "results")
METAL1_PATH = os.path.join(RESULTS_PATH, "metal_1")
METAL2_PATH = os.path.join(RESULTS_PATH, "metal_2")

# Residue numbers for categorization
RES_NUMBERS_ATOM_1 = [257, 255, 220]
RES_NUMBERS_ATOM_2 = [245, 244, 180, 181]


def create_directories():
    """Create necessary directories if they don't exist."""
    for path in [RESULTS_PATH, METAL1_PATH, METAL2_PATH]:
        os.makedirs(path, exist_ok=True)


def Metal_Ion(analysed_metal, protein, parser, path):
    """Find metal ions and their coordinates in PDB file."""
    protein_path = os.path.join(path, protein)
    structure = parser.get_structure(protein, protein_path)
    metal_occurrences = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == analysed_metal:
                        metal_occurrences.append(
                            {
                                "Atom ID": atom.id,
                                "Coordinates": atom.coord,
                                "Chain ID": chain.id,
                                "Residue Number": residue.id[1],
                                "Residue ID": residue,
                                "Structure Name": os.path.basename(protein),
                            }
                        )
    return metal_occurrences


def Calculate_distance(coord1, coord2):
    """Calculate distance between two 3D coordinates."""
    return np.linalg.norm(coord1 - coord2)


def Find_Closest_Neighbors(metal_occurrence, structure, distance_threshold):
    """Find neighboring atoms within distance threshold."""
    metal_coord = metal_occurrence["Coordinates"]
    ns = NeighborSearch(list(structure.get_atoms()))

    closest_neighbors = []
    for atom in ns.search(metal_coord, distance_threshold):
        if atom.element != metal_occurrence["Atom ID"]:
            distance = Calculate_distance(metal_coord, atom.coord)
            closest_neighbors.append(
                {
                    "Atom ID": atom.id,
                    "Coordinates": atom.coord,
                    "Distance": distance,
                    "Chain ID": atom.parent.parent.id,
                    "Residue Number": atom.parent.id[1],
                    "Residue ID": atom.parent,
                    "Structure Name": metal_occurrence["Structure Name"],
                }
            )

    metal_occurrence["Neighbors"] = closest_neighbors
    return closest_neighbors


def categorize_and_copy_files(input_dir):
    """Categorize files based on residue numbers and copy to appropriate directories."""
    metal_1 = []
    metal_2 = []
    other = []

    for file_name in os.listdir(input_dir):
        if not file_name.endswith(".csv"):
            continue

        file_path = os.path.join(input_dir, file_name)
        try:
            df = pd.read_csv(file_path)
            added = False

            for _, row in df.iterrows():
                res_number = row["Residue Number"]
                if res_number in RES_NUMBERS_ATOM_1:
                    metal_1.append(file_name)
                    added = True
                    break
                elif res_number in RES_NUMBERS_ATOM_2:
                    metal_2.append(file_name)
                    added = True
                    break

            if not added:
                other.append(file_name)

        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}")

    # Copy files to respective directories
    for file in set(metal_1):
        shutil.copy(os.path.join(input_dir, file), METAL1_PATH)
    for file in set(metal_2):
        shutil.copy(os.path.join(input_dir, file), METAL2_PATH)

    return other


def analyze_distances(metal_files_dir):
    """Analyze distances for metal coordination."""
    for metal in METALS:
        print(f"\nAnalyzing {metal} coordination...")

        # Read and combine all CSV files
        all_data = []
        for file in os.listdir(metal_files_dir):
            if file.endswith(".csv"):
                df = pd.read_csv(os.path.join(metal_files_dir, file))
                all_data.append(df)

        if not all_data:
            continue

        combined_df = pd.concat(all_data)
        filtered_df = combined_df[combined_df["Atom ID"] != metal]

        # Group and analyze
        output_text = ""
        for (atom_id, res_num), group in filtered_df.groupby(
            ["Atom ID", "Residue Number"]
        ):
            if len(group) == 0:
                continue

            distances = group["Distance"]
            stats = {
                "mean": np.mean(distances),
                "std": np.std(distances),
                "median": np.median(distances),
                "max": distances.max(),
                "min": distances.min(),
                "max_structure": group.loc[distances.idxmax(), "Structure Name"],
                "min_structure": group.loc[distances.idxmin(), "Structure Name"],
            }

            output_text += (
                f"\n{atom_id}, residue number={res_num} len={len(group)}, \n"
                f"max dist={stats['max']} from {stats['max_structure']} "
                f"and min dist = {stats['min']} from {stats['min_structure']}\n"
                f"mean = {stats['mean']:.3f} with std {stats['std']:.3f}, "
                f"median = {stats['median']:.3f}\n"
            )

        # Save results
        if output_text:
            output_file = os.path.join(metal_files_dir, f"output_{metal}_analysis.txt")
            with open(output_file, "w") as f:
                f.write(output_text)
            print(f"Results saved to {output_file}")


def main():
    try:
        print("Starting metal ion analysis pipeline...")
        create_directories()

        # Part 1: Initial metal ion analysis
        parser = PDBParser(QUIET=True)
        pdb_files = [f for f in os.listdir(PDB_PATH) if f.endswith(".pdb")]

        print(f"Found {len(pdb_files)} PDB files to analyze")

        # Process each metal
        for metal in METALS:
            print(f"\nAnalyzing {metal}...")
            all_metal_occurrences = []

            for protein in pdb_files:
                metal_occurrences = Metal_Ion(metal, protein, parser, PDB_PATH)
                all_metal_occurrences.extend(metal_occurrences)

                # Find neighbors for each occurrence
                for occurrence in metal_occurrences:
                    structure = parser.get_structure(
                        occurrence["Structure Name"],
                        os.path.join(PDB_PATH, occurrence["Structure Name"]),
                    )
                    Find_Closest_Neighbors(occurrence, structure, DISTANCE_THRESHOLD)

            # Save results for this metal
            for i, occurrence in enumerate(all_metal_occurrences):
                csv_filename = (
                    f"{occurrence['Structure Name']}_metal_occurrence_{i}.csv"
                )
                csv_path = os.path.join(RESULTS_PATH, csv_filename)

                with open(csv_path, "w") as f:
                    f.write(
                        "Atom ID,Coordinates,Distance,Chain ID,Residue Number,Residue ID,Structure Name\n"
                    )
                    f.write(
                        f"{occurrence['Atom ID']},{occurrence['Coordinates']},0,"
                        f"{occurrence['Chain ID']},{occurrence['Residue Number']},"
                        f"{occurrence['Residue ID']},{occurrence['Structure Name']}\n"
                    )
                    for neighbor in occurrence["Neighbors"]:
                        f.write(
                            f"{neighbor['Atom ID']},{neighbor['Coordinates']},"
                            f"{neighbor['Distance']},{neighbor['Chain ID']},"
                            f"{neighbor['Residue Number']},{neighbor['Residue ID']},"
                            f"{neighbor['Structure Name']}\n"
                        )

        # Part 2: Categorize results
        print("\nCategorizing results...")
        other_files = categorize_and_copy_files(RESULTS_PATH)
        if other_files:
            print("\nFiles requiring manual check:", other_files)

        # Part 3: Statistical analysis
        print("\nPerforming statistical analysis...")
        analyze_distances(METAL1_PATH)
        analyze_distances(METAL2_PATH)

        print("\nAnalysis pipeline completed successfully!")

    except Exception as e:
        print(f"An error occurred: {str(e)}")


if __name__ == "__main__":
    main()
