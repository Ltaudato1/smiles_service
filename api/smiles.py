from fastapi import APIRouter, File, UploadFile
from rdkit import Chem
import csv
from chem_storage.chem_storager import SmilesStorage

router = APIRouter()
storage = SmilesStorage()

# ======================== POST =====================================

@router.post("/add")
async def add_molecule(id: str, smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        storage.add_molecule(id, smiles)
        return {"message": "Molecule added successfully."}
    return {"message": "Invalid SMILES string."}

@router.post("/upload_csv")
async def upload_csv(file: UploadFile = File(...)):
    content = await file.read()
    decoded = content.decode("utf-8").splitlines()
    reader = csv.DictReader(decoded)

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
            added.append(mol_id)
        except Exception as e:
            errors.append({"id": mol_id, "smiles": smiles, "error": str(e)})

    return {
        "message": f"Загружено {len(added)} молекул",
        "added": added,
        "errors": errors
    }

# ======================== PUT =====================================

@router.put("/update")
async def update_molecule(id: str, smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        if storage.update_molecule(id, smiles):
            return {"message": "Molecule updated successfully."}
        else:
            return {"message": "Molecule not found."}
    else:
        return {"message": "Invalid SMILES string."}

# ======================== DELETE =====================================

@router.delete("/delete")
async def delete_molecule(id: str):
    if storage.delete_molecule(id):
        return {"message": "Molecule deleted successfully."}
    return {"message": "Molecule not found."}

# ======================== GET =====================================

@router.get("/get")
async def get_molecule(id: str):
    smiles = storage.get_molecule(id)
    if smiles:
        return {"id": id, "smiles": smiles}
    return {"message": "Molecule not found."}

@router.get("/all")
async def get_all_molecules():
    return {"molecules": storage.storage}

@router.get("/search")
async def substructure_search(molecule: str):
    mol = Chem.MolFromSmiles(molecule)
    if mol:
        return {"matches": storage.substructure_search(molecule)}
    return {"message": "Invalid SMILES string."}