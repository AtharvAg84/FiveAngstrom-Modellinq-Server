import requests
import urllib.parse
import re
from rdkit import Chem
from rdkit.Chem import AllChem

def validate_smiles(smiles):
    """
    Basic validation for SMILES string.
    Returns True if the string seems valid, False otherwise.
    """
    if not smiles or not isinstance(smiles, str):
        return False
    allowed_chars = r'^[A-Za-z0-9@#%\[\]\(\)=\-+/*.]+$'
    return bool(re.match(allowed_chars, smiles))

def name_to_smiles(name):
    """
    Convert a chemical name to SMILES using PubChem.

    Args:
        name (str): Chemical name.

    Returns:
        str: SMILES string, or None if failed.
    """
    try:
        encoded_name = urllib.parse.quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/JSON"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        return None
    except Exception as e:
        print(f"Name to SMILES failed: {e}")
        return None

def get_pubchem_coordinates_with_ids(smiles):
    """
    Fetch 2D coordinates from PubChem using SMILES, returning atom IDs, symbols, and coordinates.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        tuple: (mol, coords) where mol is RDKit molecule and coords is list of dictionaries with coordinate_id, atom_symbol, x, y, or None if failed.
    """
    if not validate_smiles(smiles):
        raise ValueError("Invalid SMILES string provided")

    try:
        encoded_smiles = urllib.parse.quote(smiles)

        # Get CID
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/cids/JSON"
        cid_resp = requests.get(cid_url, timeout=10)
        cid_resp.raise_for_status()
        cid_data = cid_resp.json()
        if "IdentifierList" not in cid_data or not cid_data["IdentifierList"].get("CID"):
            raise ValueError("No CID found for the given SMILES")
        cid = cid_data["IdentifierList"]["CID"][0]

        # Get 2D SDF from PubChem
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=2d"
        sdf_resp = requests.get(sdf_url, timeout=10)
        sdf_resp.raise_for_status()

        # Parse SDF with RDKit
        mol = Chem.MolFromMolBlock(sdf_resp.text, removeHs=False)
        if mol is None:
            raise ValueError("Failed to parse molecule from PubChem SDF")

        conformer = mol.GetConformer()
        coords = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conformer.GetAtomPosition(idx)
            coords.append({
                "coordinate_id": idx,
                "atom_symbol": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y
            })

        return (mol, coords)

    except requests.exceptions.HTTPError as e:
        raise ValueError(f"PubChem HTTP Error: {e} (Status code: {cid_resp.status_code if 'cid_resp' in locals() else sdf_resp.status_code})")
    except requests.exceptions.ConnectionError:
        raise ValueError("PubChem: Failed to connect. Check your internet connection.")
    except requests.exceptions.Timeout:
        raise ValueError("PubChem: Request timed out.")
    except requests.exceptions.RequestException as e:
        raise ValueError(f"PubChem: Error fetching data: {e}")
    except ValueError as e:
        raise
    except Exception as e:
        raise ValueError(f"PubChem: Unexpected error: {e}")

def get_rdkit_coordinates_with_ids(smiles):
    """
    Generate 2D coordinates using RDKit, returning atom IDs, symbols, and coordinates.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        tuple: (mol, coords) where mol is RDKit molecule and coords is list of dictionaries with coordinate_id, atom_symbol, x, y, or None if failed.
    """
    if not validate_smiles(smiles):
        raise ValueError("Invalid SMILES string provided")

    try:
        # Parse SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string provided")

        # Generate 2D coordinates
        if AllChem.Compute2DCoords(mol) != 0:
            raise ValueError("Failed to compute 2D coordinates")

        conformer = mol.GetConformer()
        coords = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conformer.GetAtomPosition(idx)
            coords.append({
                "coordinate_id": idx,
                "atom_symbol": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y
            })

        return (mol, coords)

    except Exception as e:
        raise ValueError(f"RDKit: Error processing SMILES: {e}")

def get_3d_coordinates(smiles):
    """
    Generate 3D coordinates using RDKit, returning atom IDs, symbols, and coordinates.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        tuple: (mol, coords) where mol is RDKit molecule and coords is list of dictionaries with coordinate_id, atom_symbol, x, y, z, or None if failed.
    """
    if not validate_smiles(smiles):
        raise ValueError("Invalid SMILES string provided")

    try:
        # Parse SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string provided")

        # Add hydrogens for 3D embedding
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            raise ValueError("Failed to embed molecule in 3D")

        # Optimize 3D coordinates
        AllChem.MMFFOptimizeMolecule(mol)

        conformer = mol.GetConformer()
        coords = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conformer.GetAtomPosition(idx)
            coords.append({
                "coordinate_id": idx,
                "atom_symbol": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z
            })

        return (mol, coords)

    except Exception as e:
        raise ValueError(f"RDKit 3D: Error processing SMILES: {e}")

def coords_to_pdb(smiles, coords, mol, name="MOLECULE"):
    """
    Convert coordinates to PDB format with connectivity.

    Args:
        smiles (str): SMILES string of the compound.
        coords (list): List of dictionaries with coordinate_id, atom_symbol, x, y, and optional z.
        mol (RDKit Mol): RDKit molecule object for bond information.
        name (str): Name of the molecule for the PDB header.

    Returns:
        str: PDB-formatted string.
    """
    if not coords or not mol:
        return ""

    # Initialize PDB content
    pdb_lines = []
    pdb_lines.append(f"HEADER    {name.upper()}")
    dimension = "3D" if any("z" in coord for coord in coords) else "2D"
    pdb_lines.append(f"TITLE     {dimension} Structure of {name}")
    pdb_lines.append("COMPND    MOL")

    # Generate ATOM records
    atom_counts = {}  # Track atom counts for unique names (e.g., C1, C2)
    for coord in coords:
        idx = coord["coordinate_id"]
        symbol = coord["atom_symbol"]
        x = coord["x"]
        y = coord["y"]
        z = coord.get("z", 0.0)  # Default to 0.0 for 2D
        # Generate unique atom name (e.g., C1, O1)
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        atom_name = f"{symbol}{atom_counts[symbol]}"
        # PDB ATOM format: ATOM serial name resName resSeq x y z occ temp element
        pdb_lines.append(
            f"ATOM  {idx + 1:>5} {atom_name:<4} MOL     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {symbol}"
        )

    # Generate CONECT records
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx() + 1  # 1-based indexing for PDB
        atom2 = bond.GetEndAtomIdx() + 1
        pdb_lines.append(f"CONECT{atom1:>5}{atom2:>5}")
        # Add reverse direction for symmetry
        pdb_lines.append(f"CONECT{atom2:>5}{atom1:>5}")

    pdb_lines.append("END")
    return "\n".join(pdb_lines)

def get_2d_coordinates(query):
    """
    Attempt to fetch 2D coordinates from PubChem, fall back to RDKit if failed.
    Handles SMILES or chemical names.

    Args:
        query (str): SMILES string or chemical name.

    Returns:
        tuple: (source, coords, mol) where source is 'PubChem' or 'RDKit', coords is a list of dictionaries, and mol is RDKit molecule, or (None, None, None) if failed.
    """
    # Try as SMILES first
    smiles = query
    if validate_smiles(smiles):
        try:
            mol, coords = get_pubchem_coordinates_with_ids(smiles)
            return ("PubChem", coords, mol)
        except ValueError as e:
            print(f"PubChem failed: {e}")

        try:
            mol, coords = get_rdkit_coordinates_with_ids(smiles)
            return ("RDKit", coords, mol)
        except ValueError as e:
            print(f"RDKit failed: {e}")

    # Try as chemical name
    smiles = name_to_smiles(query)
    if smiles:
        try:
            mol, coords = get_pubchem_coordinates_with_ids(smiles)
            return ("PubChem", coords, mol)
        except ValueError as e:
            print(f"PubChem failed: {e}")

        try:
            mol, coords = get_rdkit_coordinates_with_ids(smiles)
            return ("RDKit", coords, mol)
        except ValueError as e:
            print(f"RDKit failed: {e}")

    return (None, None, None)