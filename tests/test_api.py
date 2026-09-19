"""接口与配平逻辑测试。"""
import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.parser import parse_formula

client = TestClient(app)


# ---------- 解析器 ----------

@pytest.mark.parametrize(
    "formula,expected,charge",
    [
        ("H2O", {"H": 2, "O": 1}, 0),
        ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}, 0),
        ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}, 0),
        ("((NH4)2)3[Fe(CN)6]2", {"N": 18, "H": 24, "Fe": 2, "C": 12}, 0),
        ("CuSO4.5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}, 0),
        ("Na2CO3·10H2O", {"Na": 2, "C": 1, "O": 13, "H": 20}, 0),
        ("SO4^2-", {"S": 1, "O": 4}, -2),
        ("Fe^3+", {"Fe": 1}, 3),
        ("NH4^+", {"N": 1, "H": 4}, 1),
        ("Na+", {"Na": 1}, 1),
        ("SO4--", {"S": 1, "O": 4}, -2),
        ("e-", {}, -1),
    ],
)
def test_parse_formula(formula, expected, charge):
    counts, q = parse_formula(formula)
    assert dict(counts) == expected
    assert q == charge


# ---------- 配平接口 ----------

def balance(reactants, products):
    resp = client.post("/api/v1/balance", json={"reactants": reactants, "products": products})
    assert resp.status_code == 200, resp.text
    return resp.json()


def test_balance_simple():
    data = balance(["H2", "O2"], ["H2O"])
    assert data["coefficients"] == [2, 1, 2]
    assert data["equation"] == "2H2 + O2 -> 2H2O"
    assert data["balanced"] is True
    assert all(e["balanced"] for e in data["elements"])
    assert data["charge"] == {"reactants": 0, "products": 0, "balanced": True}


def test_balance_kmno4_hcl():
    data = balance(["KMnO4", "HCl"], ["KCl", "MnCl2", "H2O", "Cl2"])
    assert data["coefficients"] == [2, 16, 2, 2, 8, 5]


def test_balance_nested_parens():
    data = balance(["Ca(OH)2", "H3PO4"], ["Ca3(PO4)2", "H2O"])
    assert data["coefficients"] == [3, 2, 1, 6]


def test_balance_hydrate():
    data = balance(["CuSO4.5H2O"], ["CuSO4", "H2O"])
    assert data["coefficients"] == [1, 1, 5]


def test_balance_ionic_with_electron():
    data = balance(["MnO4-", "H+", "e-"], ["Mn^2+", "H2O"])
    assert data["coefficients"] == [1, 8, 5, 1, 4]
    assert data["charge"]["balanced"] is True
    assert data["charge"]["reactants"] == 2


def test_balance_ionic_disproportionation():
    data = balance(["Cl2", "OH-"], ["Cl-", "ClO3-", "H2O"])
    assert data["coefficients"] == [3, 6, 5, 1, 3]


def test_balance_nullity_two_returns_422():
    resp = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O", "H2O2"]},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


def test_balance_not_all_positive_returns_422():
    # He 是旁观物质：零空间一维但 He 的系数被强制为 0，无法全正
    resp = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2", "O2", "He"], "products": ["H2O"]},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 1
    assert "全正" in resp.json()["detail"]["message"]


def test_balance_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2O$"], "products": ["H2", "O2"]},
    )
    assert resp.status_code == 400


def test_balance_unknown_element_returns_400():
    resp = client.post(
        "/api/v1/balance",
        json={"reactants": ["Xx2"], "products": ["Xx"]},
    )
    assert resp.status_code == 400


# ---------- 校验接口 ----------

def test_verify_correct():
    resp = client.post(
        "/api/v1/verify",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "reactant_coefficients": [2, 1],
            "product_coefficients": [2],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["balanced"] is True
    assert data["is_reduced"] is True


def test_verify_unreduced_coefficients():
    resp = client.post(
        "/api/v1/verify",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "reactant_coefficients": [4, 2],
            "product_coefficients": [4],
        },
    )
    data = resp.json()
    assert data["balanced"] is True
    assert data["gcd"] == 2
    assert data["is_reduced"] is False


def test_verify_wrong_coefficients():
    resp = client.post(
        "/api/v1/verify",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "reactant_coefficients": [1, 1],
            "product_coefficients": [1],
        },
    )
    data = resp.json()
    assert data["balanced"] is False
    assert any(not e["balanced"] for e in data["elements"])


def test_verify_charge_imbalance():
    resp = client.post(
        "/api/v1/verify",
        json={
            "reactants": ["Fe^2+"],
            "products": ["Fe^3+"],
            "reactant_coefficients": [1],
            "product_coefficients": [1],
        },
    )
    data = resp.json()
    assert data["balanced"] is False
    assert data["charge"]["balanced"] is False


def test_verify_length_mismatch_returns_400():
    resp = client.post(
        "/api/v1/verify",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "reactant_coefficients": [2],
            "product_coefficients": [2],
        },
    )
    assert resp.status_code == 400


def test_health():
    assert client.get("/health").json() == {"status": "ok"}
