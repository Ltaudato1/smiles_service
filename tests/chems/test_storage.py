import pytest
from unittest.mock import MagicMock

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from src.chems.chem_storager import SmilesStorage
from src.chems.chem import Molecule


@pytest.fixture
def mock_db():
    db_mock = MagicMock()
    
    db_mock.query().filter().first.return_value = None
    db_mock.query().all.return_value = []
    
    return db_mock


def test_add_molecule(mock_db):
    storage = SmilesStorage(mock_db)
    result = storage.add_molecule("1", "CCO")

    assert result is True
    mock_db.add.assert_called_once()
    mock_db.commit.assert_called_once()


def test_get_molecule_found(mock_db):
    storage = SmilesStorage(mock_db)

    # Настраиваем возвращаемое значение
    fake_mol = Molecule(id="1", smiles="CCO")
    mock_db.query().filter().first.return_value = fake_mol

    result = storage.get_molecule("1")
    assert result == "CCO"


def test_get_molecule_not_found(mock_db):
    storage = SmilesStorage(mock_db)

    mock_db.query().filter().first.return_value = None
    result = storage.get_molecule("999")

    assert result == ""


def test_get_all_molecules(mock_db):
    storage = SmilesStorage(mock_db)

    fake_molecules = [Molecule(id="1", smiles="CCO"),
                      Molecule(id="2", smiles="c1ccccc1")]
    mock_db.query().all.return_value = fake_molecules

    result = storage.get_all_molecules()
    assert result == [
        {"id": "1", "smiles": "CCO"},
        {"id": "2", "smiles": "c1ccccc1"}
    ]


def test_update_molecule_success(mock_db):
    storage = SmilesStorage(mock_db)

    fake_mol = Molecule(id="1", smiles="CCO")
    mock_db.query().filter().first.return_value = fake_mol

    result = storage.update_molecule("1", "CCC")
    assert result is True
    assert fake_mol.smiles == "CCC"
    mock_db.commit.assert_called_once()
    mock_db.refresh.assert_called_once_with(fake_mol)


def test_update_molecule_not_found(mock_db):
    storage = SmilesStorage(mock_db)

    mock_db.query().filter().first.return_value = None
    result = storage.update_molecule("999", "CCC")

    assert result is False
    mock_db.commit.assert_not_called()


def test_delete_molecule_success(mock_db):
    storage = SmilesStorage(mock_db)

    fake_mol = Molecule(id="1", smiles="CCO")
    mock_db.query().filter().first.return_value = fake_mol

    result = storage.delete_molecule("1")
    assert result is True
    mock_db.delete.assert_called_once_with(fake_mol)
    mock_db.commit.assert_called_once()


def test_delete_molecule_not_found(mock_db):
    storage = SmilesStorage(mock_db)

    mock_db.query().filter().first.return_value = None
    result = storage.delete_molecule("999")

    assert result is False
    mock_db.delete.assert_not_called()


def test_substructure_search_single_match(monkeypatch, mock_db):
    storage = SmilesStorage(mock_db)
    monkeypatch.setattr(
        storage,
        "get_all_molecules",
        lambda: [
            {"id": 1, "smiles": "CCO"},
            {"id": 2, "smiles": "c1ccccc1"}
        ]
    )

    result = storage.substructure_search("CC")
    assert result == {1: "CCO"}


def test_substructure_search_multiple_matches(monkeypatch, mock_db):
    storage = SmilesStorage(mock_db)
    monkeypatch.setattr(
        storage,
        "get_all_molecules",
        lambda: [
            {"id": 1, "smiles": "CCO"},
            {"id": 2, "smiles": "CCN"},
            {"id": 3, "smiles": "c1ccccc1"}
        ]
    )
    result = storage.substructure_search("CC")
    assert result == {1: "CCO", 2: "CCN"}  # два совпадения


def test_substructure_search_no_match(monkeypatch, mock_db):
    storage = SmilesStorage(mock_db)
    monkeypatch.setattr(
        storage,
        "get_all_molecules",
        lambda: [
            {"id": 1, "smiles": "CCO"},
            {"id": 2, "smiles": "c1ccccc1"}
        ]
    )
    result = storage.substructure_search("NN")
    assert result == {}  # ничего не найдено


def test_substructure_search_invalid_smiles(monkeypatch, mock_db):
    storage = SmilesStorage(mock_db)
    monkeypatch.setattr(
        storage,
        "get_all_molecules",
        lambda: [
            {"id": 1, "smiles": "CCO"},
            {"id": 2, "smiles": "invalid_smiles"}
        ]
    )
    result = storage.substructure_search("CC")
    assert result == {1: "CCO"}


def test_add_duplicate_molecule(mock_db):
    storage = SmilesStorage(mock_db)
    mock_db.query().filter().first.return_value = Molecule(id="1", smiles="CCO")
    
    result = storage.add_molecule("1", "CCC")
    assert result is False


def test_add_invalid_smiles(mock_db):
    storage = SmilesStorage(mock_db)
    result = storage.add_molecule("1", "invalid_smiles")
    assert result is False