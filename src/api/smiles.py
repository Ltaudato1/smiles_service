from fastapi import APIRouter, File, UploadFile
from rdkit import Chem
import csv
from chems.chem_storager import SmilesStorage
from db import Base, engine, SessionLocal

Base.metadata.create_all(bind=engine)

router = APIRouter()

# ======================== POST =====================================


@router.post("/add")
async def add_molecule(id: str, smiles: str):
    db = SessionLocal()
    storage = SmilesStorage(db)
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        storage.add_molecule(id, smiles)
        db.close()
        return {"message": "Molecule added successfully."}
    return {"message": "Invalid SMILES string."}


@router.post("/upload_csv")
async def upload_csv(file: UploadFile = File(...)):
    content = await file.read()
    decoded = content.decode("utf-8").splitlines()
    reader = csv.DictReader(decoded)
    db = SessionLocal()
    storage = SmilesStorage(db)

    added = []
    errors = []

    for row in reader:
        mol_id = row.get("id")
        smiles = row.get("smiles")
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES")
            storage.add_molecule(mol_id, smiles)
            added.append({'id': mol_id, 'smiles': smiles})
        except Exception as e:
            errors.append({"id": mol_id, "smiles": smiles, "error": str(e)})
    db.close()

    return {
        "message": f"Uploaded {len(added)} molecules",
        "added": added,
        "errors": errors
    }


# ======================== GET =====================================


@router.get("/get")
async def get_molecule(id: str):
    db = SessionLocal()
    storage = SmilesStorage(db)
    smiles = storage.get_molecule(id)
    db.close()
    if smiles:
        return {"id": id, "smiles": smiles}
    return {"message": "Molecule not found."}


@router.get("/all")
async def get_all_molecules():
    db = SessionLocal()
    storage = SmilesStorage(db)
    molecules = storage.get_all_molecules()
    db.close()
    return molecules


@router.get("/search")
async def substructure_search(molecule: str):
    mol = Chem.MolFromSmiles(molecule)
    db = SessionLocal()
    storage = SmilesStorage(db)
    if mol:
        matches = storage.substructure_search(molecule)
        db.close()
        return {"matches": matches}
    db.close()
    return {"message": "Invalid SMILES string."}


# ======================== PUT =====================================


@router.put("/update")
async def update_molecule(id: str, smiles: str):
    db = SessionLocal()
    storage = SmilesStorage(db)
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        flag = storage.update_molecule(id, smiles)
        db.close()
        if flag:
            return {"message": "Molecule updated successfully."}
        else:
            return {"message": "Molecule not found."}
    else:
        return {"message": "Invalid SMILES string."}


# ======================== DELETE =====================================


@router.delete("/delete")
async def delete_molecule(id: str):
    db = SessionLocal()
    storage = SmilesStorage(db)
    if storage.delete_molecule(id):
        db.close()
        return {"message": "Molecule deleted successfully."}
    db.close()
    return {"message": "Molecule not found."}
