from sqlalchemy import Column, Integer, String
from db import Base


class Molecule(Base):
    __tablename__ = "smiles"

    id = Column(
        Integer,
        primary_key=True,
        index=True
    )

    smiles = Column(
        String,
        unique=True,
        index=True
    )
