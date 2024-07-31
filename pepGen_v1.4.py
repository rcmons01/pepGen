import os
import itertools
from pathlib import Path
import argparse
import stk
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel

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
    Chem.SanitizeMol(molecule)  # Ensure valences are calculated
    rdkit_mol = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(rdkit_mol, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(rdkit_mol)
    rdkit_mol = Chem.RemoveHs(rdkit_mol)
    return rdkit_mol

def save_molecule_to_file(molecule, output_dir, filename, file_format="mol2"):
    """Convert and save molecule to specified file format."""
    # Optimize the molecule
    optimized_mol = optimize_molecule(molecule)
    
    if file_format == "mol2":
        # Write the molecule to a temporary .mol file using RDKit
        temp_mol_file = os.path.join(output_dir, filename + '.mol')
        Chem.MolToMolFile(optimized_mol, temp_mol_file)
        
        # Convert from .mol to .mol2 with OpenBabel
        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats("mol", "mol2")
        ob_mol = openbabel.OBMol()
        ob_conversion.ReadFile(ob_mol, temp_mol_file)
        
        # Set the force field to MMFF94 and optimize
        ff = openbabel.OBForceField.FindForceField("MMFF94")
        ff.Setup(ob_mol)
        ff.ConjugateGradients(250)  # Perform a force field optimization
        ff.GetCoordinates(ob_mol)
        
        mol2_file = os.path.join(output_dir, filename + '.mol2')
        ob_conversion.WriteFile(ob_mol, mol2_file)
        
        # Remove the temporary .mol file
        os.remove(temp_mol_file)
        
        print(f"Generated {mol2_file}")
    
    elif file_format == "sdf":
        # Generate 2D coordinates
        AllChem.Compute2DCoords(optimized_mol)
        
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

def swap_terminal_hydroxyl_with_nitrogen(molecule):
    """Swap the terminal hydroxyl group with a nitrogen atom in the last residue."""
    rw_mol = Chem.RWMol(molecule)
    carbon_idx = next(
        (atom.GetIdx() for atom in rw_mol.GetAtoms()
         if atom.GetSymbol() == 'C' and any(
             neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1
             for neighbor in atom.GetNeighbors()
         )),
        None
    )
    if carbon_idx is not None:
        oxygen_idx = next(
            (atom.GetIdx() for atom in rw_mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1),
            None
        )
        if oxygen_idx is not None:
            rw_mol.ReplaceAtom(oxygen_idx, Chem.Atom('N'))
    Chem.SanitizeMol(rw_mol)  # Sanitize after modification
    return rw_mol

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

    for seq in sequences:
        polymer = create_polymer(seq, building_blocks)
        sequence_names = '_'.join(building_blocks[idx][0] for idx in seq)
        polymer_rdkit_mol = polymer.to_rdkit_mol()
        polymer_rdkit_mol = Chem.RWMol(polymer_rdkit_mol)
        polymer_rdkit_mol = remove_first_bromine(polymer_rdkit_mol)
        polymer_rdkit_mol = swap_terminal_hydroxyl_with_nitrogen(polymer_rdkit_mol)
        save_molecule_to_file(polymer_rdkit_mol, args.output_dir, f'peptoid_{sequence_names}', file_format=args.format)

    print("All molecules have been generated and saved.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate peptoid molecules.")
    parser.add_argument("--smiles_file", required=True, help="Path to the SMILES file.")
    parser.add_argument("--output_dir", required=True, help="Name of the output directory.")
    parser.add_argument("--mode", choices=["combinations", "single"], required=True, help="Mode of generation: 'combinations' for all n-length sequences, 'single' for a specified sequence.")
    parser.add_argument("--length", type=int, help="Length of the peptoid sequence (required if mode is 'combinations').")
    parser.add_argument("--sequence", help="Residue names sequence (required if mode is 'single').")
    parser.add_argument("--format", choices=["mol2", "sdf"], default="mol2", help="Output file format: 'mol2' or 'sdf'.")

    args = parser.parse_args()
    main(args)
