from fastapi import APIRouter
from rdkit import Chem
from chem_storage.chem_storager import SmilesStorage

router = APIRouter()
storage = SmilesStorage()

@router.post("/add")
def add_molecule(id: str, smiles: str):
    storage.add_molecule(id, smiles)
    return {"message": "Molecule added successfully."}

@router.get("/get")
def get_molecule(id: str):
    smiles = storage.get_molecule(id)
    if smiles:
        return {"id": id, "smiles": smiles}
    return {"message": "Molecule not found."}

@router.put("/update")
def update_molecule(id: str, smiles: str):
    if storage.update_molecule(id, smiles):
        return {"message": "Molecule updated successfully."}
    return {"message": "Molecule not found."}

@router.delete("/delete")
def delete_molecule(id: str):
    if storage.delete_molecule(id):
        return {"message": "Molecule deleted successfully."}
    return {"message": "Molecule not found."}

@router.get("/all")
def get_all_molecules():
    return {"molecules": storage.storage}

@router.get("/search")
def substructure_search(molecule: str):
    return {"matches": storage.substructure_search(molecule)}