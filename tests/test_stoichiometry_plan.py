"""批次投料规划接口测试。

内置原子量：H 1.008 / O 15.999 / Cl 35.45 / K 39.0983 / Mn 54.938
=> M(H2)=2.016, M(O2)=31.998, M(H2O)=18.015, M(HCl)=36.458,
   M(KMnO4)=158.0323, M(Cl2)=70.9
"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

H2O_PLAN = {
    "reactants": ["H2", "O2"],
    "products": ["H2O"],
    "target_product": "H2O",
    "target_mass": "36.03",
    "yield_fraction": 1,
    "materials": [
        {"formula": "H2", "purity": 1, "excess_percent": 0},
        {"formula": "O2", "purity": 1, "excess_percent": 0},
    ],
}


def plan(**overrides):
    payload = {**H2O_PLAN, **overrides}
    resp = client.post("/api/v1/stoichiometry/plan", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 正常计算 ----------

def test_basic_plan_full_yield_no_excess():
    # 36.03 g H2O = 2 mol；收率 1 -> 理论 2 mol -> ξ = 1
    data = plan()
    assert data["equation"] == "2H2 + O2 -> 2H2O"
    assert data["coefficients"] == [2, 1, 2]
    assert data["precision"] == 4
    assert data["target_product"] == "H2O"
    assert data["extent"] == "1.0000"

    h2, o2 = data["materials"]
    assert h2["formula"] == "H2"
    assert h2["coefficient"] == 2
    assert h2["stoichiometric_moles"] == "2.0000"
    assert h2["excess_moles"] == "2.0000"
    assert h2["actual_grams"] == "4.0320"  # 2 * 2.016 / 1
    assert o2["formula"] == "O2"
    assert o2["coefficient"] == 1
    assert o2["stoichiometric_moles"] == "1.0000"
    assert o2["excess_moles"] == "1.0000"
    assert o2["actual_grams"] == "31.9980"  # 1 * 31.998 / 1


def test_yield_excess_and_purity():
    # 100 g H2O / 18.015 = 5.5509... mol；收率 0.8 -> 理论 6.9387... mol -> ξ = 3.4693...
    data = plan(
        target_mass=100,
        yield_fraction="0.8",
        materials=[
            {"formula": "H2", "purity": "0.95", "excess_percent": 10},
            {"formula": "O2", "purity": 1, "excess_percent": 0},
        ],
    )
    assert data["extent"] == "3.4693"
    h2, o2 = data["materials"]
    assert h2["stoichiometric_moles"] == "6.9387"  # 2 * ξ
    assert h2["excess_moles"] == "7.6325"  # * 1.1
    assert h2["actual_grams"] == "16.1970"  # 7.6325... * 2.016 / 0.95
    assert o2["stoichiometric_moles"] == "3.4693"
    assert o2["excess_moles"] == "3.4693"
    assert o2["actual_grams"] == "111.0117"  # 3.4693... * 31.998 / 1


def test_kmno4_hcl_plan_with_precision():
    # 70.9 g Cl2 = 1 mol；收率 0.5 -> 理论 2 mol -> ξ = 0.4
    data = plan(
        reactants=["KMnO4", "HCl"],
        products=["KCl", "MnCl2", "H2O", "Cl2"],
        target_product="Cl2",
        target_mass="70.9",
        yield_fraction="0.5",
        materials=[
            {"formula": "KMnO4", "purity": "0.99", "excess_percent": 0},
            {"formula": "HCl", "purity": "0.36", "excess_percent": 20},
        ],
        precision=6,
    )
    assert data["coefficients"] == [2, 16, 2, 2, 8, 5]
    assert data["precision"] == 6
    assert data["extent"] == "0.400000"
    kmno4, hcl = data["materials"]
    assert kmno4["stoichiometric_moles"] == "0.800000"  # 2 * 0.4
    assert kmno4["excess_moles"] == "0.800000"
    assert kmno4["actual_grams"] == "127.702869"  # 0.8 * 158.0323 / 0.99
    assert hcl["stoichiometric_moles"] == "6.400000"  # 16 * 0.4
    assert hcl["excess_moles"] == "7.680000"  # * 1.2
    assert hcl["actual_grams"] == "777.770667"  # 7.68 * 36.458 / 0.36


def test_precision_parameter():
    data = plan(precision=2)
    assert data["precision"] == 2
    assert data["extent"] == "1.00"
    assert data["materials"][0]["actual_grams"] == "4.03"
    data = plan(precision=8)
    assert data["precision"] == 8
    assert data["extent"] == "1.00000000"
    assert data["materials"][1]["actual_grams"] == "31.99800000"


# ---------- 4xx 错误 ----------

def test_target_not_in_products_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "target_product": "O2"},
    )
    assert resp.status_code == 400
    assert "不在生成物" in resp.json()["detail"]
    assert "O2" in resp.json()["detail"]


@pytest.mark.parametrize("target_mass", [0, -1, "-2.5"])
def test_non_positive_target_mass_returns_400(target_mass):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "target_mass": target_mass},
    )
    assert resp.status_code == 400
    assert "正数" in resp.json()["detail"]


@pytest.mark.parametrize("yield_fraction", [0, -1, "1.5", 2])
def test_yield_out_of_range_returns_400(yield_fraction):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "yield_fraction": yield_fraction},
    )
    assert resp.status_code == 400
    assert "收率" in resp.json()["detail"]


@pytest.mark.parametrize("purity", [0, -1, "1.5"])
def test_purity_out_of_range_returns_400(purity):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "materials": [
                {"formula": "H2", "purity": purity, "excess_percent": 0},
                {"formula": "O2", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400
    assert "纯度" in resp.json()["detail"]


def test_negative_excess_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "materials": [
                {"formula": "H2", "purity": 1, "excess_percent": "-5"},
                {"formula": "O2", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400
    assert "过量" in resp.json()["detail"]


def test_missing_material_config_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "materials": [{"formula": "H2", "purity": 1, "excess_percent": 0}],
        },
    )
    assert resp.status_code == 400
    assert "缺少" in resp.json()["detail"]
    assert "O2" in resp.json()["detail"]


def test_empty_materials_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "materials": []},
    )
    assert resp.status_code == 400


def test_duplicate_material_config_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "materials": [
                {"formula": "H2", "purity": 1, "excess_percent": 0},
                {"formula": "H2", "purity": "0.9", "excess_percent": 5},
                {"formula": "O2", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400
    assert "重复配置" in resp.json()["detail"]


def test_material_of_non_reactant_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "materials": [
                {"formula": "H2", "purity": 1, "excess_percent": 0},
                {"formula": "O2", "purity": 1, "excess_percent": 0},
                {"formula": "H2O", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400
    assert "不在反应物" in resp.json()["detail"]


def test_unsupported_atomic_weight_returns_400():
    # Zn 不在内置原子量表
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["Zn", "HCl"],
            "products": ["ZnCl2", "H2"],
            "target_product": "H2",
            "target_mass": 2,
            "yield_fraction": 1,
            "materials": [
                {"formula": "Zn", "purity": 1, "excess_percent": 0},
                {"formula": "HCl", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400
    assert "Zn" in resp.json()["detail"]
    assert "原子量" in resp.json()["detail"]


def test_unbalanceable_returns_422():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "products": ["H2O", "H2O2"]},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


def test_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            **H2O_PLAN,
            "reactants": ["H2O$", "O2"],
            "materials": [
                {"formula": "H2O$", "purity": 1, "excess_percent": 0},
                {"formula": "O2", "purity": 1, "excess_percent": 0},
            ],
        },
    )
    assert resp.status_code == 400


@pytest.mark.parametrize("precision", [1, 9])
def test_precision_out_of_range_returns_422(precision):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={**H2O_PLAN, "precision": precision},
    )
    assert resp.status_code == 422
