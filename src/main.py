from fastapi import FastAPI
from api import smiles_router
from os import getenv

app = FastAPI(
    title="Smiles Service API",
    description="API для поиска и работы с молекулами",
    version="1.0.0"
)

app.include_router(smiles_router, prefix="/smiles", tags=["smiles"])


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}
