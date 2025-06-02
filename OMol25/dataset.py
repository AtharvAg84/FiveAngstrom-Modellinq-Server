import os
from datasets import load_dataset
from rdkit import Chem
import torch
from torch.utils.data import DataLoader
from memory_profiler import profile
from itertools import islice

# Configuration
HF_TOKEN = os.getenv("HF_TOKEN")  # Ensure HF_TOKEN is set in your environment
BATCH_SIZE = 32
MAX_EXAMPLES = 1000  # Limit for demonstration; remove for full dataset

@profile  # Monitor memory usage
def process_omol25():
    # Load dataset in streaming mode
    dataset = load_dataset(
        "facebook/OMol25",
        split="train",
        streaming=True,
        use_auth_token=HF_TOKEN
    )

    # Define processing function for molecular data
    def process_molecule(example):
        try:
            # Convert SMILES to RDKit molecule and compute number of atoms
            mol = Chem.MolFromSmiles(example["smiles"])
            num_atoms = mol.GetNumAtoms() if mol else 0
            # Normalize energy (hypothetical example)
            energy = example["energy"] / 1000.0 if "energy" in example else 0.0
            return {
                "smiles": example["smiles"],
                "num_atoms": num_atoms,
                "energy": energy
            }
        except:
            # Handle invalid SMILES
            return {
                "smiles": example["smiles"],
                "num_atoms": 0,
                "energy": 0.0
            }

    # Apply preprocessing
    processed_dataset = dataset.map(process_molecule)

    # Create DataLoader for batching
    def collate_fn(examples):
        smiles = [ex["smiles"] for ex in examples]
        num_atoms = torch.tensor([ex["num_atoms"] for ex in examples], dtype=torch.int32)
        energies = torch.tensor([ex["energy"] for ex in examples], dtype=torch.float32)
        return {"smiles": smiles, "num_atoms": num_atoms, "energies": energies}

    # Wrap dataset in DataLoader
    dataloader = DataLoader(
        processed_dataset,
        batch_size=BATCH_SIZE,
        collate_fn=collate_fn
    )

    # Process batches (limited to MAX_EXAMPLES for demo)
    example_count = 0
    for batch in dataloader:
        # Simulate model inference
        print(f"Batch SMILES (first 2): {batch['smiles'][:2]}")
        print(f"Batch Num Atoms (first 2): {batch['num_atoms'][:2]}")
        print(f"Batch Energies (first 2): {batch['energies'][:2]}")
        
        example_count += len(batch["smiles"])
        if example_count >= MAX_EXAMPLES:
            break

    print(f"Processed {example_count} examples.")

# Run the processing
if __name__ == "__main__":
    process_omol25()