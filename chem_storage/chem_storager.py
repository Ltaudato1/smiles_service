from rdkit import Chem


class SmilesStorage:
    def __init__(self):
        self.storage = {}

    def add_molecule(self, id: str, smiles: str) -> bool:
        self.storage[id] = smiles
        return True

    def get_molecule(self, id: str) -> str:
        return self.storage.get(id, "")

    def update_molecule(self, id: str, smiles: str) -> bool:
        if id in self.storage:
            self.storage[id] = smiles
            return True
        return False

    def delete_molecule(self, id: str) -> bool:
        if id in self.storage:
            del self.storage[id]
            return True
        return False

    def substructure_search(self, molecule: str) -> dict:
        target_molecule = Chem.MolFromSmiles(molecule)
        answer = {}
        for id, smiles in self.storage.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol.HasSubstructMatch(target_molecule):
                answer[id] = smiles
        return answer
