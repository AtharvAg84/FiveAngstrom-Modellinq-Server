import os
import logging
import torch
import numpy as np
from fairchem.core.datasets import AseDBDataset
from rdkit import Chem
import requests
import pickle
from pathlib import Path
import time
import json

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Configuration
DATASET_PATH = "./train_4M/"  # Update with actual path
INDEX_FILE = "smiles_API_index.pkl"  # Cache file for SMILES-to-index mapping
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"

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

def fetch_pubchem_properties(smiles):
    """
    Fetch molecular properties from PubChem PUG REST API.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        dict: Molecular properties, or None if request fails.
    """
    try:
        # Canonicalize SMILES for consistency
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.warning(f"Invalid SMILES for PubChem: {smiles}")
            return None
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

        # API request for properties
        url = f"{PUBCHEM_API_BASE}/{canonical_smiles}/property/MolecularWeight,LogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract properties
        props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
        return {
            "molecular_weight": props.get("MolecularWeight"),
            "logP": props.get("LogP"),
            "tpsa": props.get("TPSA"),
            "num_h_donors": props.get("HBondDonorCount"),
            "num_h_acceptors": props.get("HBondAcceptorCount"),
            "cid": props.get("CID")  # PubChem Compound ID
        }
    except requests.RequestException as e:
        logger.error(f"PubChem API request failed: {str(e)}")
        return None
    except Exception as e:
        logger.error(f"Error processing PubChem response: {str(e)}")
        return None

@torch.no_grad()
def get_molecule_details_by_smiles_api(dataset_path, smiles, index_file=INDEX_FILE):
    """
    Retrieve molecular details using PubChem API and local OMol25 dataset with GPU acceleration.
    
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

        # Fetch properties from PubChem API
        logger.info(f"Fetching properties from PubChem for SMILES: {smiles}")
        pubchem_props = fetch_pubchem_properties(smiles)
        if not pubchem_props:
            logger.warning("Falling back to RDKit for molecular properties")

        # Load SMILES index and check local dataset
        smiles_to_index = create_smiles_index(dataset_path, index_file)
        molecule_index = smiles_to_index.get(smiles)
        atoms = None
        if molecule_index is not None:
            dataset = AseDBDataset(config=dict(src=dataset_path))
            if molecule_index < len(dataset):
                atoms = dataset.get_atoms(molecule_index)
                logger.info(f"Found molecule in dataset at index {molecule_index}")
            else:
                logger.warning(f"Index {molecule_index} out of range")
                molecule_index = None

        # Initialize RDKit molecule for local computations
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logger.error(f"Invalid SMILES: {smiles}")
            return None
        mol = Chem.AddHs(mol)

        # Extract local dataset properties if available
        details = {
            "smiles": smiles,
            "index": molecule_index if molecule_index is not None else -1,
            "energy": atoms.info.get("energy", 0.0) / 1000.0 if atoms else None,
            "positions": atoms.positions.tolist() if atoms else None,
            "atomic_numbers": atoms.get_atomic_numbers().tolist() if atoms else None,
            "atom_symbols": atoms.get_chemical_symbols() if atoms else None
        }

        # Use GPU for coordinate processing if available
        if details["positions"]:
            positions = torch.tensor(details["positions"], dtype=torch.float32, device=device)
            details["positions"] = positions.cpu().tolist()  # Move back to CPU for output
        if details["atomic_numbers"]:
            atomic_numbers = torch.tensor(details["atomic_numbers"], dtype=torch.int32, device=device)
            details["atomic_numbers"] = atomic_numbers.cpu().tolist()

        # Populate molecular properties (API or RDKit fallback)
        if pubchem_props:
            details.update({
                "molecular_weight": pubchem_props["molecular_weight"],
                "logP": pubchem_props["logP"],
                "tpsa": pubchem_props["tpsa"],
                "num_h_donors": pubchem_props["num_h_donors"],
                "num_h_acceptors": pubchem_props["num_h_acceptors"],
                "pubchem_cid": pubchem_props["cid"]
            })
        else:
            details.update({
                "molecular_weight": Descriptors.MolWt(mol),
                "logP": Descriptors.MolLogP(mol),
                "tpsa": Descriptors.TPSA(mol),
                "num_h_donors": Descriptors.NumHDonors(mol),
                "num_h_acceptors": Descriptors.NumHAcceptors(mol),
                "pubchem_cid": None
            })

        # Extract RDKit-specific properties
        details.update({
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": Chem.GetSSSR(mol),
            "bonds": [],
            "formal_charges": [],
            "degrees": []
        })

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
    print(f"SMILES: {details['smiles']}")
    print(f"Dataset Index: {details['index'] if details['index'] >= 0 else 'Not in local dataset'}")
    print(f"PubChem CID: {details['pubchem_cid'] if details['pubchem_cid'] else 'Not available'}")
    print(f"Energy (eV): {details['energy']:.4f}" if details['energy'] is not None else "Energy: Not available")
    print(f"Number of Atoms: {details['num_atoms']}")
    print(f"Molecular Weight: {details['molecular_weight']:.2f} g/mol")
    print(f"Number of Bonds: {details['num_bonds']}")
    print(f"Number of Rings: {details['num_rings']}")
    print(f"LogP: {details['logP']:.4f}")
    print(f"TPSA (Å²): {details['tpsa']:.2f}")
    print(f"H-bond Donors: {details['num_h_donors']}")
    print(f"H-bond Acceptors: {details['num_h_acceptors']}")
    if details["atom_symbols"]:
        print("\nAtomic Details:")
        for i, (symbol, charge, degree) in enumerate(zip(details['atom_symbols'], details['formal_charges'], details['degrees'])):
            print(f"Atom {i}: {symbol}, Formal Charge: {charge}, Degree: {degree}")
    if details["positions"]:
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
        details = get_molecule_details_by_smiles_api(DATASET_PATH, smiles)
        print_molecule_details(details)