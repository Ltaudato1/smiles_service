from rdkit import Chem
from chems.chem import Molecule
from db import Base, engine, SessionLocal

Base.metadata.create_all(bind=engine)


class SmilesStorage:
    def __init__(self):
        pass

    def add_molecule(self, id: str, smiles: str) -> bool:
        db = SessionLocal()
        db.add(Molecule(id=id, smiles=smiles))
        db.commit()
        db.close()
        return True

    def get_molecule(self, id: str) -> str:
        db = SessionLocal()
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        db.close()
        return molecule.smiles if molecule else ""

    def get_all_molecules(self) -> list:
        db = SessionLocal()
        molecules = db.query(Molecule).all()
        db.close()
        return [{'id': mol.id, 'smiles': mol.smiles} for mol in molecules]

    def update_molecule(self, id: str, smiles: str) -> bool:
        db = SessionLocal()
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        if molecule:
            molecule.smiles = smiles
            db.commit()
            db.refresh(molecule)
            db.close()
            return True
        db.close()
        return False

    def delete_molecule(self, id: str) -> bool:
        db = SessionLocal()
        molecule = db.query(Molecule).filter(Molecule.id == id).first()
        if molecule:
            db.delete(molecule)
            db.commit()
            db.close()
            return True
        db.close()
        return False

    def substructure_search(self, molecule: str) -> dict:
        target_molecule = Chem.MolFromSmiles(molecule)
        molecules = self.get_all_molecules()
        answer = {}
        for mol in molecules:
            mol_struct = Chem.MolFromSmiles(mol['smiles'])
            if mol_struct.HasSubstructMatch(target_molecule):
                answer[mol['id']] = mol['smiles']
        return answer
