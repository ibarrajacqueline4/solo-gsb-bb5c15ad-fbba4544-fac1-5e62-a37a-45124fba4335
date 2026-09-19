"""经验式/分子式推断接口测试。

内置原子量：H 1.008 / C 12.011 / N 14.007 / O 15.999
=> M(CH2O)=30.026, M(H2O)=18.015, M(N2O3)=76.011, M(C6H12O6)=180.156
"""
import pytest
from fastapi.testclient import TestClient

from app.formula_infer import _candidate_subscripts, hill_order, hill_string
from app.main import app

client = TestClient(app)


# ---------- 内部工具 ----------

def test_hill_order():
    assert hill_order(["O", "C", "H"]) == ["C", "H", "O"]
    assert hill_order(["O", "H"]) == ["H", "O"]  # 无 C 时纯字母序
    assert hill_order(["S", "C", "Na"]) == ["C", "Na", "S"]
    assert hill_order(["Fe"]) == ["Fe"]


def test_hill_string():
    assert hill_string({"C": 1, "H": 2, "O": 1}) == "CH2O"
    assert hill_string({"H": 2, "O": 1}) == "H2O"
    assert hill_string({"C": 6, "H": 12, "O": 6}) == "C6H12O6"
    assert hill_string({"Na": 1, "Cl": 1}) == "ClNa"  # 字母序


def test_candidate_subscripts_coprime_and_deduped():
    from decimal import Decimal
    from math import gcd

    # 摩尔比 1:2:1：所有尺度都归约为同一个互质候选
    cands = _candidate_subscripts(
        {"C": Decimal(1), "H": Decimal(2), "O": Decimal(1)}, ["C", "H", "O"]
    )
    assert cands == [{"C": 1, "H": 2, "O": 1}]
    # 摩尔比 1:1.5：尺度 1 得 (1,2)，尺度 2 得互质的 (2,3)
    cands = _candidate_subscripts(
        {"N": Decimal(1), "O": Decimal("1.5")}, ["N", "O"]
    )
    assert {"N": 1, "O": 2} in cands
    assert {"N": 2, "O": 3} in cands
    # 全部候选互质（gcd == 1）且无重复
    assert all(gcd(*c.values()) == 1 for c in cands)
    keys = [tuple(c[e] for e in sorted(c)) for c in cands]
    assert len(keys) == len(set(keys))


def infer(percentages, tolerance=None, molar_mass=None, precision=None):
    payload = {"percentages": percentages}
    if tolerance is not None:
        payload["tolerance"] = tolerance
    if molar_mass is not None:
        payload["molar_mass"] = molar_mass
    if precision is not None:
        payload["precision"] = precision
    resp = client.post("/api/v1/formula/infer", json=payload)
    assert resp.status_code == 200, resp.text
    return resp.json()


# ---------- 正常推断 ----------

def test_glucose_empirical_only():
    data = infer({"C": "40.0", "H": "6.7", "O": "53.3"})
    assert data["empirical_formula"] == "CH2O"
    assert data["subscripts"] == {"C": 1, "H": 2, "O": 1}
    assert data["empirical_molar_mass"] == "30.0260"
    # 理论百分比（按内置原子量反算）
    assert data["theoretical_percentages"] == {
        "C": "40.0020",
        "H": "6.7142",
        "O": "53.2838",
    }
    # 逐元素偏差与最大偏差
    assert data["deviations"]["C"] == "0.0020"
    assert data["deviations"]["H"] == "0.0142"
    assert data["deviations"]["O"] == "0.0162"
    assert data["max_deviation"] == "0.0162"
    # 未提供 molar_mass：分子式相关字段为 null
    assert data["molecular_formula"] is None
    assert data["multiplier"] is None
    assert data["molecular_molar_mass"] is None
    assert data["molar_mass_relative_deviation"] is None


def test_glucose_with_molar_mass():
    data = infer({"C": "40.0", "H": "6.7", "O": "53.3"}, molar_mass="180.0")
    assert data["empirical_formula"] == "CH2O"
    assert data["molecular_formula"] == "C6H12O6"
    assert data["multiplier"] == 6
    assert data["molecular_molar_mass"] == "180.1560"  # 6 * 30.026
    # 相对偏差 |180.156 - 180| / 180 * 100 = 0.0867%
    assert data["molar_mass_relative_deviation"] == "0.0867"


def test_water():
    data = infer({"H": "11.19", "O": "88.81"})
    assert data["empirical_formula"] == "H2O"  # 无 C：Hill 顺序即字母序
    assert data["subscripts"] == {"H": 2, "O": 1}


def test_hill_order_with_carbon():
    # 含 C 时 C 优先、H 次之，其余按字母序
    data = infer({"C": "52.14", "H": "13.13", "O": "34.73"})
    assert data["empirical_formula"] == "C2H6O"
    assert list(data["subscripts"]) == ["C", "H", "O"]


def test_higher_scale_candidate_wins():
    # N2O3：摩尔比 N:O = 1 : 1.5，尺度 1 的 NO2 偏差过大，尺度 2 的 N2O3 胜出
    data = infer({"N": "36.86", "O": "63.14"})
    assert data["empirical_formula"] == "N2O3"
    assert data["subscripts"] == {"N": 2, "O": 3}
    assert data["empirical_molar_mass"] == "76.0110"


def test_single_element():
    data = infer({"Fe": "100"})
    assert data["empirical_formula"] == "Fe"
    assert data["subscripts"] == {"Fe": 1}
    assert data["theoretical_percentages"] == {"Fe": "100.0000"}
    assert data["deviations"] == {"Fe": "0.0000"}
    assert data["max_deviation"] == "0.0000"


def test_exact_theoretical_percentages_give_zero_deviation():
    # 用 CH2O 的理论百分比反推，偏差应为 0
    data = infer(
        {"C": "40.0020", "H": "6.7142", "O": "53.2838"}, tolerance="0.01"
    )
    assert data["empirical_formula"] == "CH2O"
    assert data["max_deviation"] == "0.0000"


def test_molar_mass_multiplier_one():
    data = infer({"H": "11.19", "O": "88.81"}, molar_mass="18.015")
    assert data["empirical_formula"] == "H2O"
    assert data["multiplier"] == 1
    assert data["molecular_formula"] == "H2O"
    assert data["molar_mass_relative_deviation"] == "0.0000"


def test_sum_within_tolerance_accepted():
    # 总和 99.7，与 100 的偏差 0.3 ≤ 默认 tolerance 0.5
    data = infer({"C": "39.9", "H": "6.6", "O": "53.2"})
    assert data["empirical_formula"] == "CH2O"


def test_precision_parameter():
    data = infer({"C": "40.0", "H": "6.7", "O": "53.3"}, precision=2)
    assert data["theoretical_percentages"]["C"] == "40.00"
    assert data["empirical_molar_mass"] == "30.03"
    data = infer({"C": "40.0", "H": "6.7", "O": "53.3"}, precision=8)
    assert data["empirical_molar_mass"] == "30.02600000"


# ---------- 400：输入数据不合法 ----------

def test_sum_deviation_exceeds_tolerance_returns_400():
    resp = client.post(
        "/api/v1/formula/infer",
        json={"percentages": {"C": "40.0", "H": "6.7", "O": "50.0"}},
    )
    assert resp.status_code == 400
    assert "总和" in resp.json()["detail"]


@pytest.mark.parametrize("bad", ["0", "-1.5"])
def test_non_positive_percentage_returns_400(bad):
    resp = client.post(
        "/api/v1/formula/infer",
        json={"percentages": {"C": bad, "O": "60"}},
    )
    assert resp.status_code == 400
    assert "正数" in resp.json()["detail"]


def test_unknown_element_returns_400():
    resp = client.post(
        "/api/v1/formula/infer",
        json={"percentages": {"Zn": "65.4", "S": "34.6"}},
    )
    assert resp.status_code == 400
    assert "Zn" in resp.json()["detail"]
    assert "原子量" in resp.json()["detail"]


def test_non_positive_molar_mass_returns_400():
    resp = client.post(
        "/api/v1/formula/infer",
        json={
            "percentages": {"C": "40.0", "H": "6.7", "O": "53.3"},
            "molar_mass": "0",
        },
    )
    assert resp.status_code == 400
    assert "摩尔质量" in resp.json()["detail"]


def test_empty_percentages_returns_422():
    resp = client.post("/api/v1/formula/infer", json={"percentages": {}})
    assert resp.status_code == 422


# ---------- 422：无解 ----------

def test_no_candidate_within_tolerance_returns_422():
    # C 42% / O 58%：1～12 倍尺度内最接近的候选 CO 偏差仍约 0.88 个百分点，
    # 超过默认 tolerance=0.5
    resp = client.post(
        "/api/v1/formula/infer",
        json={"percentages": {"C": "42.0", "O": "58.0"}},
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "无候选式" in detail["message"]
    assert detail["nearest_candidate"] == "CO"


def test_molar_mass_mismatch_returns_422():
    # CH2O 的整数倍摩尔质量中最接近 100 的是 3 倍 = 90.078，相对偏差约 9.9%
    resp = client.post(
        "/api/v1/formula/infer",
        json={
            "percentages": {"C": "40.0", "H": "6.7", "O": "53.3"},
            "molar_mass": "100",
        },
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "相对偏差" in detail["message"]
    assert detail["empirical_formula"] == "CH2O"
    assert detail["multiplier"] == 3


# ---------- 422：参数校验 ----------

def test_negative_tolerance_returns_422():
    resp = client.post(
        "/api/v1/formula/infer",
        json={"percentages": {"H": "11.19", "O": "88.81"}, "tolerance": "-0.1"},
    )
    assert resp.status_code == 422


@pytest.mark.parametrize("precision", [1, 9])
def test_precision_out_of_range_returns_422(precision):
    resp = client.post(
        "/api/v1/formula/infer",
        json={
            "percentages": {"H": "11.19", "O": "88.81"},
            "precision": precision,
        },
    )
    assert resp.status_code == 422
