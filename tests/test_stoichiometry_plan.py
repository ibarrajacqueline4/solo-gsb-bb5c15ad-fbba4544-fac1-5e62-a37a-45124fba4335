"""反应批次投料规划接口测试。

内置原子量：H 1.008 / O 15.999 / Cl 35.45 / K 39.0983 / Mn 54.938
=> M(H2)=2.016, M(O2)=31.998, M(H2O)=18.015, M(HCl)=36.458,
   M(KMnO4)=158.0323, M(Cl2)=70.9
"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def plan(reactants, products, target_product, target_mass, yield_fraction,
         materials, precision=None):
    payload = {
        "reactants": reactants,
        "products": products,
        "target_product": target_product,
        "target_mass": target_mass,
        "yield_fraction": yield_fraction,
        "materials": materials,
    }
    if precision is not None:
        payload["precision"] = precision
    resp = client.post("/api/v1/stoichiometry/plan", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 正常计算 ----------

def test_basic_plan_no_excess_full_purity():
    # 目标 36.03 g H2O（恰好 2 mol），收率 1、纯度 1、无过量
    # => 理论进度 1 mol，H2 2 mol / 4.032 g，O2 1 mol / 31.998 g
    data = plan(
        ["H2", "O2"],
        ["H2O"],
        "H2O",
        "36.03",
        "1",
        [
            {"formula": "H2", "purity": "1", "excess_percent": "0"},
            {"formula": "O2", "purity": "1", "excess_percent": "0"},
        ],
    )
    assert data["equation"] == "2H2 + O2 -> 2H2O"
    assert data["coefficients"] == [2, 1, 2]
    assert data["precision"] == 4
    assert data["target_product"] == "H2O"
    assert data["target_mass"] == "36.03"
    assert data["yield_fraction"] == "1"
    assert data["target_theoretical_moles"] == "2.0000"
    assert data["target_theoretical_grams"] == "36.0300"
    assert data["extent"] == "1.0000"

    h2, o2 = data["materials"]
    assert h2["formula"] == "H2"
    assert h2["coefficient"] == 2
    assert h2["purity"] == "1"
    assert h2["excess_percent"] == "0"
    assert h2["stoichiometric_moles"] == "2.0000"
    assert h2["moles_with_excess"] == "2.0000"
    assert h2["weigh_grams"] == "4.0320"  # 2 * 2.016
    assert o2["stoichiometric_moles"] == "1.0000"
    assert o2["weigh_grams"] == "31.9980"  # 1 * 31.998


def test_plan_with_yield_purity_and_excess():
    # 目标 36.03 g H2O，收率 0.9 => 理论 40.0333 g = 20/9 mol，进度 10/9 mol
    # H2：化学计量 20/9 mol，过量 10% => 22/9 mol，纯度 0.98
    #     => 称量 (22/9)*2.016/0.98 = 5.0285714... g
    data = plan(
        ["H2", "O2"],
        ["H2O"],
        "H2O",
        "36.03",
        "0.9",
        [
            {"formula": "H2", "purity": "0.98", "excess_percent": "10"},
            {"formula": "O2", "purity": "1", "excess_percent": "0"},
        ],
    )
    assert data["target_theoretical_grams"] == "40.0333"  # 36.03 / 0.9
    assert data["target_theoretical_moles"] == "2.2222"  # 20/9
    assert data["extent"] == "1.1111"  # 10/9

    h2, o2 = data["materials"]
    assert h2["purity"] == "0.98"
    assert h2["excess_percent"] == "10"
    assert h2["stoichiometric_moles"] == "2.2222"
    assert h2["moles_with_excess"] == "2.4444"  # 22/9
    assert h2["weigh_grams"] == "5.0286"
    assert o2["stoichiometric_moles"] == "1.1111"
    assert o2["moles_with_excess"] == "1.1111"
    assert o2["weigh_grams"] == "35.5533"  # (10/9) * 31.998


def test_plan_kmno4_hcl_target_cl2():
    # 2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 8H2O + 5Cl2
    # 目标 70.9 g Cl2（1 mol），收率 1 => 进度 0.2 mol
    data = plan(
        ["KMnO4", "HCl"],
        ["KCl", "MnCl2", "H2O", "Cl2"],
        "Cl2",
        "70.9",
        "1",
        [
            {"formula": "KMnO4", "purity": "1", "excess_percent": "0"},
            {"formula": "HCl", "purity": "0.36", "excess_percent": "20"},
        ],
    )
    assert data["coefficients"] == [2, 16, 2, 2, 8, 5]
    assert data["extent"] == "0.2000"
    kmno4, hcl = data["materials"]
    assert kmno4["stoichiometric_moles"] == "0.4000"
    assert kmno4["weigh_grams"] == "63.2129"  # 0.4 * 158.0323
    assert hcl["stoichiometric_moles"] == "3.2000"
    assert hcl["moles_with_excess"] == "3.8400"  # 3.2 * 1.2
    assert hcl["weigh_grams"] == "388.8853"  # 3.84 * 36.458 / 0.36


def test_plan_rounding_is_half_up():
    # 目标 0.0018015 g H2O = 0.0001 mol => 进度 0.00005 mol，
    # 第 5 位恰为 5：四舍五入（HALF_UP）应收尾为 0.0001
    data = plan(
        ["H2", "O2"],
        ["H2O"],
        "H2O",
        "0.0018015",
        "1",
        [
            {"formula": "H2", "purity": "1", "excess_percent": "0"},
            {"formula": "O2", "purity": "1", "excess_percent": "0"},
        ],
    )
    assert data["extent"] == "0.0001"
    assert data["materials"][0]["stoichiometric_moles"] == "0.0001"


def test_plan_precision_parameter():
    materials = [
        {"formula": "H2", "purity": "1", "excess_percent": "0"},
        {"formula": "O2", "purity": "1", "excess_percent": "0"},
    ]
    data = plan(["H2", "O2"], ["H2O"], "H2O", "36.03", "1", materials, precision=2)
    assert data["precision"] == 2
    assert data["extent"] == "1.00"
    assert data["materials"][0]["weigh_grams"] == "4.03"
    data = plan(["H2", "O2"], ["H2O"], "H2O", "36.03", "1", materials, precision=8)
    assert data["precision"] == 8
    assert data["extent"] == "1.00000000"
    assert data["materials"][0]["weigh_grams"] == "4.03200000"


# ---------- 4xx 错误 ----------

def test_target_not_in_products_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O2",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "不在生成物" in resp.json()["detail"]


@pytest.mark.parametrize("target_mass", [0, -1, "-2.5"])
def test_non_positive_target_mass_returns_400(target_mass):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": target_mass,
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "正数" in resp.json()["detail"]


@pytest.mark.parametrize("yield_fraction", [0, -0.5, "1.0001", 2])
def test_yield_out_of_range_returns_400(yield_fraction):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": yield_fraction,
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "收率" in resp.json()["detail"]


@pytest.mark.parametrize("purity", [0, -0.1, "1.0001", 2])
def test_purity_out_of_range_returns_400(purity):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": purity, "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "纯度" in resp.json()["detail"]


def test_negative_excess_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "-5"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "过量" in resp.json()["detail"]


def test_missing_material_config_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [{"formula": "H2", "purity": "1", "excess_percent": "0"}],
        },
    )
    assert resp.status_code == 400
    assert "缺少" in resp.json()["detail"]
    assert "O2" in resp.json()["detail"]


def test_empty_materials_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [],
        },
    )
    assert resp.status_code == 400


def test_duplicate_material_config_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "H2", "purity": "0.9", "excess_percent": "5"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_material_of_non_reactant_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
                {"formula": "H2O", "purity": "1", "excess_percent": "0"},
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
            "target_mass": "2",
            "yield_fraction": "1",
            "materials": [
                {"formula": "Zn", "purity": "1", "excess_percent": "0"},
                {"formula": "HCl", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400
    assert "Zn" in resp.json()["detail"]
    assert "原子量" in resp.json()["detail"]


def test_unbalanceable_returns_422():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O", "H2O2"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["null_space_dim"] == 2


def test_bad_formula_returns_400():
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2O$", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2O$", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
        },
    )
    assert resp.status_code == 400


@pytest.mark.parametrize("precision", [1, 9])
def test_precision_out_of_range_returns_422(precision):
    resp = client.post(
        "/api/v1/stoichiometry/plan",
        json={
            "reactants": ["H2", "O2"],
            "products": ["H2O"],
            "target_product": "H2O",
            "target_mass": "10",
            "yield_fraction": "1",
            "materials": [
                {"formula": "H2", "purity": "1", "excess_percent": "0"},
                {"formula": "O2", "purity": "1", "excess_percent": "0"},
            ],
            "precision": precision,
        },
    )
    assert resp.status_code == 422
