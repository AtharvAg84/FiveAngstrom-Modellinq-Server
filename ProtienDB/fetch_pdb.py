from Bio.PDB import PDBList
import os

class PDBFetcher:
    """
    A class to fetch PDB files from the RCSB PDB database using Biopython.
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