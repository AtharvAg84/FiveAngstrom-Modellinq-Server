from Bio.PDB import PDBList
import os
import requests
import re

class PDBFetcher:
    """
    A class to fetch PDB files and FASTA sequences from the RCSB PDB database using Biopython and REST API.
    """
    def __init__(self, temp_dir: str = "temp_pdb"):
        """
        Initialize the PDBFetcher with a temporary directory for file storage.
        
        Args:
            temp_dir (str): Directory to store temporary PDB files. Defaults to 'temp_pdb'.
        """
        self.temp_dir = temp_dir
        os.makedirs(self.temp_dir, exist_ok=True)
        self.pdb_list = PDBList()

    def fetch_pdb_file(self, pdb_id: str) -> str:
        """
        Retrieve the PDB file content for a given PDB ID.
        
        Args:
            pdb_id (str): The PDB ID (e.g., '1A3N').
        
        Returns:
            str: The PDB file content as a string, or empty string if failed.
        """
        try:
            # Validate PDB ID (basic check for 4-character alphanumeric)
            if not isinstance(pdb_id, str) or not len(pdb_id) == 4 or not pdb_id.isalnum():
                print(f"Invalid PDB ID: {pdb_id}")
                return ""

            # Download PDB file to the temporary directory
            file_path = self.pdb_list.retrieve_pdb_file(
                pdb_id,
                pdir=self.temp_dir,
                file_format="pdb",
                overwrite=True
            )

            # Read the PDB file content
            if os.path.exists(file_path):
                with open(file_path, 'r') as f:
                    pdb_content = f.read()
                # Clean up the temporary file
                os.remove(file_path)
                return pdb_content
            else:
                print(f"PDB file not found for ID: {pdb_id}")
                return ""
        except Exception as e:
            print(f"Error fetching PDB file for {pdb_id}: {str(e)}")
            return ""

    def fetch_fasta_sequence(self, pdb_id: str) -> str:
        """
        Retrieve the FASTA sequence for a given PDB ID using the RCSB PDB REST API.
        
        Args:
            pdb_id (str): The PDB ID (e.g., '1A3N').
        
        Returns:
            str: The FASTA sequence as a string (amino acids only), or empty string if failed.
        """
        try:
            # Validate PDB ID
            if not isinstance(pdb_id, str) or not len(pdb_id) == 4 or not pdb_id.isalnum():
                print(f"Invalid PDB ID: {pdb_id}")
                return ""

            # Fetch FASTA from RCSB PDB REST API
            url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
            response = requests.get(url)
            response.raise_for_status()

            # Parse FASTA content
            fasta_content = response.text
            # Extract sequence (remove header lines starting with '>')
            sequence = "".join(line.strip() for line in fasta_content.split("\n") if not line.startswith(">"))
            # Remove any non-amino acid characters (basic cleanup)
            sequence = re.sub(r"[^A-Za-z]", "", sequence)
            
            if not sequence:
                print(f"No valid FASTA sequence found for ID: {pdb_id}")
                return ""
            return sequence
        except requests.RequestException as e:
            print(f"Error fetching FASTA sequence for {pdb_id}: {str(e)}")
            return ""
        except Exception as e:
            print(f"Unexpected error for {pdb_id}: {str(e)}")
            return ""