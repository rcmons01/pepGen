import os
import itertools
import random
from pathlib import Path
import argparse
import stk
from rdkit import Chem
from rdkit.Chem import AllChem

def read_smiles_file(smiles_file):
    """Read SMILES strings and residue names from a file and return a list of tuples (name, smiles)."""
    residues = []
    with open(smiles_file, 'r') as file:
        for line in file:
            name, smiles = line.strip().split()
            residues.append((name, smiles))
    return residues

def create_building_blocks(residues):
    """Create BuildingBlock objects from SMILES strings with defined functional groups."""
    building_blocks = []
    for name, smiles in residues:
        building_block = stk.BuildingBlock(
            smiles,
            functional_groups=[
                stk.BromoFactory(),
                stk.AmideFactory()
            ]
        )
        building_blocks.append((name, building_block))
    return building_blocks

def create_polymer(sequence, building_blocks):
    """Create a polymer with the specified sequence."""
    selected_blocks = [building_blocks[idx][1] for idx in sequence]
    repeating_unit = [i for i in range(len(sequence))]

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=selected_blocks,
            repeating_unit=repeating_unit,
            num_repeating_units=1,
            optimizer=stk.Collapser(scale_steps=False),
        )
    )
    return polymer

def optimize_molecule(molecule):
    """Optimize the molecule using RDKit."""
    try:
        Chem.SanitizeMol(molecule)  # Ensure valences are calculated
        
        if molecule.GetNumAtoms() > 200:  # Skip large molecules
            print("Skipping large molecule with more than 200 atoms.")
            return None
        
        rdkit_mol = Chem.AddHs(molecule)
        
        # Attempt to embed the molecule in 3D
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d
        params.numThreads = 0  # Use all available threads
        
        if AllChem.EmbedMolecule(rdkit_mol, params=params) != 0:
            print("3D embedding failed, falling back to 2D coordinates.")
            AllChem.Compute2DCoords(rdkit_mol)
        else:
            # Optimize using MMFF94 force field
            if AllChem.MMFFOptimizeMolecule(rdkit_mol) != 0:
                raise ValueError("MMFF optimization failed.")
        
    except Exception as e:
        print(f"Error optimizing molecule: {e}")
        return None  # Return None if optimization fails
    
    return rdkit_mol

def save_molecule_to_file(molecule, output_dir, filename):
    """Save the optimized molecule to an SDF file."""
    optimized_mol = optimize_molecule(molecule)
    
    if optimized_mol is None:
        print(f"Skipping molecule {filename} due to optimization failure.")
        return
    
    # Write the molecule to .sdf file using RDKit
    sdf_file = os.path.join(output_dir, filename + '.sdf')
    writer = Chem.SDWriter(sdf_file)
    writer.write(optimized_mol)
    writer.close()
    
    print(f"Generated {sdf_file}")

def remove_first_bromine(molecule):
    """Remove the bromine atom from the first residue in the molecule."""
    rw_mol = Chem.RWMol(molecule)
    bromine_idx = next((atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetSymbol() == 'Br'), None)
    if bromine_idx is not None:
        rw_mol.RemoveAtom(bromine_idx)
    Chem.SanitizeMol(rw_mol)  # Sanitize after modification
    return rw_mol

def generate_random_sequences(num_blocks, length, number):
    """Generate a specified number of random sequences."""
    sequences = []
    for _ in range(number):
        sequence = [random.randint(0, num_blocks-1) for _ in range(length)]
        sequences.append(sequence)
    return sequences

def main(args):
    residues = read_smiles_file(args.smiles_file)
    print(f"Residues: {residues}")
    building_blocks = create_building_blocks(residues)

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    if args.mode == "combinations":
        sequences = itertools.product(range(len(building_blocks)), repeat=args.length)
    elif args.mode == "single":
        residue_names = args.sequence.split('_')
        sequences = [tuple(building_blocks.index((name, bb)) for name, bb in building_blocks if name in residue_names)]
    elif args.mode == "random":
        sequences = generate_random_sequences(len(building_blocks), args.length, args.number)

    for seq in sequences:
        polymer = create_polymer(seq, building_blocks)
        sequence_names = '_'.join(building_blocks[idx][0] for idx in seq)
        polymer_rdkit_mol = polymer.to_rdkit_mol()
        polymer_rdkit_mol = Chem.RWMol(polymer_rdkit_mol)
        polymer_rdkit_mol = remove_first_bromine(polymer_rdkit_mol)
        save_molecule_to_file(polymer_rdkit_mol, args.output_dir, f'peptoid_{sequence_names}')

    print("All molecules have been generated and saved.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate peptoid molecules.")
    parser.add_argument("--smiles_file", required=True, help="Path to the SMILES file.")
    parser.add_argument("--output_dir", required=True, help="Name of the output directory.")
    parser.add_argument("--mode", choices=["combinations", "single", "random"], required=True, help="Mode of generation: 'combinations' for all n-length sequences, 'single' for a specified sequence, 'random' for random sequences.")
    parser.add_argument("--length", type=int, help="Length of the peptoid sequence (required if mode is 'combinations' or 'random').")
    parser.add_argument("--number", type=int, help="Number of random sequences to generate (required if mode is 'random').")
    parser.add_argument("--sequence", help="Residue names sequence (required if mode is 'single').")

    args = parser.parse_args()
    main(args)
