import os
import logging
import torch
import numpy as np
from fairchem.core.datasets import AseDBDataset
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pickle
from pathlib import Path
import time

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Configuration
DATASET_PATH = "./train_4M/"  # Update with actual path
INDEX_FILE = "smiles_index.pkl"  # Cache file for SMILES-to-index mapping
BATCH_SIZE = 1000  # Batch size for dataset loading

def create_smiles_index(dataset_path, index_file):
    """Create and cache a SMILES-to-index mapping for fast lookup."""
    index_path = Path(index_file)
    if index_path.exists():
        logger.info(f"Loading SMILES index from {index_file}")
        with open(index_path, "rb") as f:
            return pickle.load(f)

    logger.info("Creating SMILES index...")
    dataset = AseDBDataset(config=dict(src=dataset_path))
    smiles_to_index = {}
    for i in range(len(dataset)):
        atoms = dataset.get_atoms(i)
        smiles = atoms.info.get("smiles", None)
        if smiles:
            smiles_to_index[smiles] = i

    with open(index_path, "wb") as f:
        pickle.dump(smiles_to_index, f)
    logger.info(f"SMILES index saved to {index_file}")
    return smiles_to_index

def smiles_to_tensor(smiles, device):
    """Convert SMILES string to a tensor of ASCII values for GPU processing."""
    ascii_vals = [ord(c) for c in smiles]
    return torch.tensor(ascii_vals, dtype=torch.int32, device=device)

@torch.no_grad()
def get_molecule_details_by_smiles(dataset_path, smiles, index_file=INDEX_FILE):
    """
    Retrieve detailed information about a molecule by SMILES using GPU acceleration.
    
    Args:
        dataset_path (str): Path to the dataset directory.
        smiles (str): SMILES string of the molecule.
        index_file (str): Path to SMILES index file.
    
    Returns:
        dict: Molecular details, or None if not found/invalid.
    """
    start_time = time.time()
    try:
        # Check for GPU availability
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info(f"Using device: {device}")

        # Validate SMILES
        if not smiles or not isinstance(smiles, str):
            logger.error("Invalid SMILES: Must be a non-empty string")
            return None

        # Load SMILES index
        smiles_to_index = create_smiles_index(dataset_path, index_file)
        molecule_index = smiles_to_index.get(smiles)
        if molecule_index is None:
            logger.error(f"No molecule found with SMILES: {smiles}")
            return None

        # Load dataset
        dataset = AseDBDataset(config=dict(src=dataset_path))
        if molecule_index >= len(dataset):
            logger.error(f"Index {molecule_index} out of range")
            return None

        # Get atoms object
        atoms = dataset.get_atoms(molecule_index)
        logger.info(f"Processing molecule with SMILES: {smiles} at index {molecule_index}")

        # Initialize RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None

        # Add hydrogens (CPU-based, as RDKit is not GPU-accelerated)
        mol = Chem.AddHs(mol)

        # Extract properties
        positions = torch.tensor(atoms.positions, dtype=torch.float32, device=device)
        atomic_numbers = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.int32, device=device)

        details = {
            "smiles": smiles,
            "index": molecule_index,
            "energy": atoms.info.get("energy", 0.0) / 1000.0,  # Normalize energy (eV)
            "num_atoms": mol.GetNumAtoms(),
            "molecular_weight": Descriptors.MolWt(mol),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": Chem.GetSSSR(mol),
            "logP": Descriptors.MolLogP(mol),
            "positions": positions.cpu().tolist(),  # Move back to CPU for output
            "atomic_numbers": atomic_numbers.cpu().tolist(),
            "atom_symbols": atoms.get_chemical_symbols(),
            "bonds": [],
            "formal_charges": [],
            "degrees": [],
            "tpsa": Descriptors.TPSA(mol),
            "num_h_donors": Descriptors.NumHDonors(mol),
            "num_h_acceptors": Descriptors.NumHAcceptors(mol)
        }

        # Extract bond and atom properties
        for bond in mol.GetBonds():
            details["bonds"].append({
                "begin_atom": bond.GetBeginAtomIdx(),
                "end_atom": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType())
            })
        for atom in mol.GetAtoms():
            details["formal_charges"].append(atom.GetFormalCharge())
            details["degrees"].append(atom.GetDegree())

        logger.info(f"Processing completed in {time.time() - start_time:.2f} seconds")
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
    print(f"Dataset Index: {details['index']}")
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
    # Prompt user for SMILES input
    smiles = input("Enter the SMILES string of the molecule: ").strip()
    if not smiles:
        logger.error("No SMILES provided")
        print("Error: Please provide a valid SMILES string.")
    else:
        details = get_molecule_details_by_smiles(DATASET_PATH, smiles)
        print_molecule_details(details)