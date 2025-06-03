from fairchem.core.datasets import AseDBDataset

dataset_path = "./train_4M/"
dataset = AseDBDataset({"src": dataset_path})
print(dataset)
for i in range(0,100):
# index the dataset to get a structure
    atoms = dataset.get_atoms(i)
    print(atoms.get_())
    print("\n --------------- \n")

