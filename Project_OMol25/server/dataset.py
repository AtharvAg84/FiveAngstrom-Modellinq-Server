import os
import logging
from fairchem.core.datasets import AseDBDataset
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Configuration
DATASET_PATH = "./train_4M/"  # Update with actual path

def get_molecule_details(dataset_path, molecule_index=0, smiles=None):
    """
    Retrieve detailed information about a molecule from the OMol25 dataset.
    
    Args:
        dataset_path (str): Path to the dataset directory.
        molecule_index (int): Index of the molecule in the dataset (used if smiles is None).
        smiles (str): SMILES string to find a specific molecule (optional).
    
    Returns:
        dict: Dictionary containing molecular details.
    """
    try:
        # Load dataset
        logger.info(f"Loading dataset from {dataset_path}...")
        dataset = AseDBDataset(config=dict(src=dataset_path))
        logger.info("Dataset loaded successfully.")

        # Select molecule
        if smiles:
            # Search for molecule by SMILES
            for i in range(len(dataset)):
                atoms = dataset.get_atoms(i)
                if atoms.info.get("smiles") == smiles:
                    molecule_index = i
                    break
            else:
                logger.error(f"No molecule found with SMILES: {smiles}")
                return None
        else:
            if molecule_index >= len(dataset):
                logger.error(f"Index {molecule_index} out of range for dataset size {len(dataset)}")
                return None

        # Get atoms object
        atoms = dataset.get_atoms(molecule_index)
        smiles = atoms.info.get("smiles", "C")  # Fallback SMILES
        logger.info(f"Processing molecule with SMILES: {smiles}")

        # Initialize RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None

        # Add hydrogens and compute 3D conformation (optional)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates if needed

        # Extract properties
        details = {
            "smiles": smiles,
            "index": molecule_index,
            "energy": atoms.info.get("energy", 0.0) / 1000.0,  # Normalize energy (eV)
            "num_atoms": mol.GetNumAtoms(),
            "molecular_weight": Descriptors.MolWt(mol),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": Chem.GetSSSR(mol),  # Number of rings
            "logP": Descriptors.MolLogP(mol),  # Partition coefficient
            "positions": atoms.positions.tolist(),  # 3D coordinates
            "atomic_numbers": atoms.get_atomic_numbers().tolist(),
            "atom_symbols": atoms.get_chemical_symbols(),
            "bonds": [],
            "formal_charges": [],
            "degrees": []  # Number of bonds per atom
        }

        # Extract bond information
        for bond in mol.GetBonds():
            details["bonds"].append({
                "begin_atom": bond.GetBeginAtomIdx(),
                "end_atom": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType())
            })

        # Extract atom properties
        for atom in mol.GetAtoms():
            details["formal_charges"].append(atom.GetFormalCharge())
            details["degrees"].append(atom.GetDegree())

        # Additional descriptors
        details["tpsa"] = Descriptors.TPSA(mol)  # Topological polar surface area
        details["num_h_donors"] = Descriptors.NumHDonors(mol)
        details["num_h_acceptors"] = Descriptors.NumHAcceptors(mol)

        logger.info(f"Successfully retrieved details for molecule at index {molecule_index}")
        return details

    except Exception as e:
        logger.error(f"Error processing molecule: {str(e)}")
        return None

def print_molecule_details(details):
    """Pretty-print molecule details."""
    if not details:
        print("No molecule details available.")
        return

    print("\nMolecule Details:")
    print(f"Index: {details['index']}")
    print(f"SMILES: {details['smiles']}")
    print(f"Energy (eV): {details['energy']:.4f}")
    print(f"Number of Atoms: {details['num_atoms']}")
    print(f"Molecular Weight: {details['molecular_weight']:.2f} g/mol")
    print(f"Number of Bonds: {details['num_bonds']}")
    print(f"Number of Rings: {details['num_rings']}")
    print(f"LogP: {details['logP']:.4f}")
    print(f"TPSA (Å²): {details['tpsa']:.2f}")
    print(f"H-bond Donors: {details['num_h_donors']}")
    print(f"H-bond Acceptors: {details['num_h_acceptors']}")
    print("\nAtomic Details:")
    for i, (symbol, charge, degree) in enumerate(zip(details['atom_symbols'], details['formal_charges'], details['degrees'])):
        print(f"Atom {i}: {symbol}, Formal Charge: {charge}, Degree: {degree}")
    print("\n3D Coordinates:")
    for i, pos in enumerate(details['positions']):
        print(f"Atom {i}: {pos}")
    print("\nBonds:")
    for bond in details['bonds']:
        print(f"Bond {bond['begin_atom']} - {bond['end_atom']}: {bond['bond_type']}")

if __name__ == "__main__":
    # Example usage: Get details for molecule at index 0 or by SMILES
    molecule_index = 0
    # smiles = "CCO"  # Uncomment to search by SMILES
    details = get_molecule_details(DATASET_PATH, molecule_index=molecule_index)  # Add smiles=smiles if searching by SMILES
    print_molecule_details(details)