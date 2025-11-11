from rdkit import Chem
from src.chems.chem import Molecule


class SmilesStorage:
    def __init__(self, db_session):
        self.db_session = db_session

    def add_molecule(self, id: str, smiles: str) -> bool:
        db = self.db_session
        if Chem.MolFromSmiles(smiles) is None:
            return False # Invalid smiles
        if db.query(Molecule).filter(Molecule.id == id).first():
            return False # DB contains this molecule
        
        db.add(Molecule(id=id, smiles=smiles))
        db.commit()
        return True

    def get_molecule(self, id: str) -> str:
        db = self.db_session
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        return molecule.smiles if molecule else ""

    def get_all_molecules(self) -> list:
        db = self.db_session
        molecules = db.query(Molecule).all()
        return [{'id': mol.id, 'smiles': mol.smiles} for mol in molecules]

    def update_molecule(self, id: str, smiles: str) -> bool:
        db = self.db_session
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        if molecule:
            molecule.smiles = smiles
            self.db_session.commit()
            self.db_session.refresh(molecule)
            return True
        return False

    def delete_molecule(self, id: str) -> bool:
        db = self.db_session
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        if molecule:
            db.delete(molecule)
            db.commit()
            return True
        return False

    def substructure_search(self, molecule: str) -> dict:
        target_molecule = Chem.MolFromSmiles(molecule)
        molecules = self.get_all_molecules()
        answer = {}
        for mol in molecules:
            mol_struct = Chem.MolFromSmiles(mol['smiles'])
            if mol_struct and mol_struct.HasSubstructMatch(target_molecule):
                answer[mol['id']] = mol['smiles']
        return answer
