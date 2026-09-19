"""化学计量接口测试。

内置原子量：H 1.008 / O 15.999 / Cl 35.45 / K 39.0983 / Mn 54.938
=> M(H2)=2.016, M(O2)=31.998, M(H2O)=18.015, M(HCl)=36.458,
   M(KMnO4)=158.0323, M(Cl2)=70.9
"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def stoich(reactants, products, feed, precision=None):
    payload = {"reactants": reactants, "products": products, "feed": feed}
    if precision is not None:
        payload["precision"] = precision
    resp = client.post("/api/v1/stoichiometry", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 正常计算 ----------

def test_basic_mol_feed():
    data = stoich(
        ["H2", "O2"],
        ["H2O"],
        [
            {"formula": "H2", "amount": 2, "unit": "mol"},
            {"formula": "O2", "amount": 2, "unit": "mol"},
        ],
    )
    assert data["equation"] == "2H2 + O2 -> 2H2O"
    assert data["coefficients"] == [2, 1, 2]
    assert data["precision"] == 4
    assert data["extent"] == "1.0000"
    assert data["limiting_reactants"] == ["H2"]

    (h2o,) = data["products"]
    assert h2o["formula"] == "H2O"
    assert h2o["coefficient"] == 2
    assert h2o["theoretical_moles"] == "2.0000"
    assert h2o["theoretical_grams"] == "36.0300"  # 2 * 18.015

    h2, o2 = data["reactants"]
    assert h2["limiting"] is True
    assert h2["feed_moles"] == "2.0000"
    assert h2["feed_grams"] == "4.0320"  # 2 * 2.016
    assert h2["consumed_moles"] == "2.0000"
    assert h2["remaining_moles"] == "0.0000"
    assert h2["remaining_grams"] == "0.0000"
    assert o2["limiting"] is False
    assert o2["feed_grams"] == "63.9960"  # 2 * 31.998
    assert o2["remaining_moles"] == "1.0000"
    assert o2["remaining_grams"] == "31.9980"


def test_gram_feed_with_tie_limiting():
    # 4.032 g H2 = 2 mol；31.998 g O2 = 1 mol —— 反应进度均为 1，并列限制
    data = stoich(
        ["H2", "O2"],
        ["H2O"],
        [
            {"formula": "H2", "amount": "4.032", "unit": "g"},
            {"formula": "O2", "amount": "31.998", "unit": "g"},
        ],
    )
    assert data["extent"] == "1.0000"
    assert data["limiting_reactants"] == ["H2", "O2"]
    assert all(r["remaining_moles"] == "0.0000" for r in data["reactants"])
    assert data["products"][0]["theoretical_moles"] == "2.0000"


def test_mixed_units_feed():
    # H2 10 mol；O2 100 g = 100/31.998 = 3.1251953... mol -> O2 限制
    data = stoich(
        ["H2", "O2"],
        ["H2O"],
        [
            {"formula": "H2", "amount": 10, "unit": "mol"},
            {"formula": "O2", "amount": 100, "unit": "g"},
        ],
    )
    assert data["limiting_reactants"] == ["O2"]
    assert data["extent"] == "3.1252"
    h2, o2 = data["reactants"]
    assert o2["feed_moles"] == "3.1252"
    assert o2["feed_grams"] == "100.0000"
    assert h2["remaining_moles"] == "3.7496"  # 10 - 2*3.1251953...
    assert h2["remaining_grams"] == "7.5592"  # 3.7496093... * 2.016
    (h2o,) = data["products"]
    assert h2o["theoretical_moles"] == "6.2504"
    assert h2o["theoretical_grams"] == "112.6008"


def test_rounding_is_half_up():
    # 0.12345 在第 5 位恰为 5：四舍五入（HALF_UP）应收尾为 0.1235
    data = stoich(
        ["H2", "O2"],
        ["H2O"],
        [
            {"formula": "H2", "amount": "0.12345", "unit": "mol"},
            {"formula": "O2", "amount": 10, "unit": "mol"},
        ],
    )
    h2 = data["reactants"][0]
    assert h2["feed_moles"] == "0.1235"
    assert data["extent"] == "0.0617"  # 0.061725
    assert data["products"][0]["theoretical_moles"] == "0.1235"


def test_precision_parameter():
    feed = [
        {"formula": "H2", "amount": 2, "unit": "mol"},
        {"formula": "O2", "amount": 2, "unit": "mol"},
    ]
    data = stoich(["H2", "O2"], ["H2O"], feed, precision=2)
    assert data["precision"] == 2
    assert data["products"][0]["theoretical_grams"] == "36.03"
    data = stoich(["H2", "O2"], ["H2O"], feed, precision=8)
    assert data["precision"] == 8
    assert data["products"][0]["theoretical_grams"] == "36.03000000"
    assert data["extent"] == "1.00000000"


def test_kmno4_hcl_stoichiometry():
    data = stoich(
        ["KMnO4", "HCl"],
        ["KCl", "MnCl2", "H2O", "Cl2"],
        [
            {"formula": "KMnO4", "amount": 2, "unit": "mol"},
            {"formula": "HCl", "amount": 20, "unit": "mol"},
        ],
    )
    assert data["coefficients"] == [2, 16, 2, 2, 8, 5]
    assert data["limiting_reactants"] == ["KMnO4"]
    assert data["extent"] == "1.0000"
    cl2 = data["products"][3]
    assert cl2["theoretical_moles"] == "5.0000"
    assert cl2["theoretical_grams"] == "354.5000"  # 5 * 70.9
    hcl = data["reactants"][1]
    assert hcl["remaining_moles"] == "4.0000"  # 20 - 16
    assert hcl["remaining_grams"] == "145.8320"  # 4 * 36.458


# ---------- 4xx 错误 ----------

def test_missing_feed_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [{"formula": "H2", "amount": 2, "unit": "mol"}],
        },
    )
    assert resp.status_code == 400
    assert "缺少" in resp.json()["detail"]
    assert "O2" in resp.json()["detail"]


def test_empty_feed_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={"reactants": ["H2", "O2"], "products": ["H2O"], "feed": []},
    )
    assert resp.status_code == 400


def test_duplicate_feed_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2", "amount": 2, "unit": "mol"},
                {"formula": "H2", "amount": 3, "unit": "g"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "重复投料" in resp.json()["detail"]


@pytest.mark.parametrize("amount", [0, -1, "-2.5"])
def test_non_positive_amount_returns_400(amount):
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2", "amount": amount, "unit": "mol"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "正数" in resp.json()["detail"]


def test_feed_of_non_reactant_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2", "amount": 2, "unit": "mol"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
                {"formula": "H2O", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "不在反应物" in resp.json()["detail"]


def test_invalid_unit_returns_422():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2", "amount": 2, "unit": "kg"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 422


def test_unsupported_atomic_weight_returns_400():
    # Zn 不在内置原子量表
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["Zn", "HCl"],
            "products": ["ZnCl2", "H2"],
            "feed": [
                {"formula": "Zn", "amount": 1, "unit": "mol"},
                {"formula": "HCl", "amount": 2, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "Zn" in resp.json()["detail"]
    assert "原子量" in resp.json()["detail"]


def test_unbalanceable_returns_422():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O", "H2O2"],
            "feed": [
                {"formula": "H2", "amount": 2, "unit": "mol"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


def test_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2O$", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2O$", "amount": 2, "unit": "mol"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
        },
    )
    assert resp.status_code == 400


@pytest.mark.parametrize("precision", [1, 9])
def test_precision_out_of_range_returns_422(precision):
    resp = client.post(
        "/api/v1/stoichiometry",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "feed": [
                {"formula": "H2", "amount": 2, "unit": "mol"},
                {"formula": "O2", "amount": 1, "unit": "mol"},
            ],
            "precision": precision,
        },
    )
    assert resp.status_code == 422
