"""守恒矩阵构建与有理数零空间求解。

约定：矩阵每列对应一种物质，反应物取正号、生成物取负号；
每行对应一种元素，最后一行为电荷行。零空间向量 x 满足 A·x = 0，
即各元素与电荷在反应前后守恒。仅当零空间为一维且基向量可化为
全正时，方程存在唯一（不计倍数）的全正整数配平。
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd, lcm

from .parser import parse_formula


class BalanceError(Exception):
    """方程无法配平（零空间维数不为 1，或无法化为全正系数）。"""

    def __init__(self, message: str, nullity: int | None = None):
        super().__init__(message)
        self.message = message
        self.nullity = nullity


def build_matrix(
    compositions: list[Counter], charges: list[int], n_reactants: int
) -> tuple[list[str], list[list[Fraction]]]:
    """构建元素守恒矩阵（末行为电荷行）。反应物列为正，生成物列为负。"""
    elements = sorted({el for comp in compositions for el in comp})
    n = len(compositions)
    matrix: list[list[Fraction]] = []
    for el in elements:
        row = []
        for j in range(n):
            sign = 1 if j < n_reactants else -1
            row.append(Fraction(sign * compositions[j].get(el, 0)))
        matrix.append(row)
    charge_row = []
    for j in range(n):
        sign = 1 if j < n_reactants else -1
        charge_row.append(Fraction(sign * charges[j]))
    matrix.append(charge_row)
    return elements, matrix


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    """有理数高斯消元，返回 (行最简形, 主元列下标列表)。"""
    mat = [row[:] for row in matrix]
    n_rows = len(mat)
    n_cols = len(mat[0]) if mat else 0
    pivots: list[int] = []
    r = 0
    for c in range(n_cols):
        if r >= n_rows:
            break
        piv = next((i for i in range(r, n_rows) if mat[i][c] != 0), None)
        if piv is None:
            continue
        mat[r], mat[piv] = mat[piv], mat[r]
        pv = mat[r][c]
        mat[r] = [x / pv for x in mat[r]]
        for i in range(n_rows):
            if i != r and mat[i][c] != 0:
                f = mat[i][c]
                mat[i] = [a - f * b for a, b in zip(mat[i], mat[r])]
        pivots.append(c)
        r += 1
    return mat, pivots


def null_space(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """返回零空间的一组基（每个自由变量对应一个基向量）。"""
    n_cols = len(matrix[0]) if matrix else 0
    mat, pivots = rref(matrix)
    pivot_set = set(pivots)
    free_cols = [c for c in range(n_cols) if c not in pivot_set]
    basis = []
    for f in free_cols:
        vec = [Fraction(0)] * n_cols
        vec[f] = Fraction(1)
        for row_idx, pcol in enumerate(pivots):
            vec[pcol] = -mat[row_idx][f]
        basis.append(vec)
    return basis


def to_integer_coeffs(vec: list[Fraction]) -> list[int]:
    """有理数向量 → 最大公约数约简后的最小正整数系数。"""
    denom_lcm = 1
    for x in vec:
        denom_lcm = lcm(denom_lcm, x.denominator)
    ints = [int(x * denom_lcm) for x in vec]
    g = 0
    for v in ints:
        g = gcd(g, v)
    g = abs(g) or 1
    return [v // g for v in ints]


def account(
    compositions: list[Counter],
    charges: list[int],
    coeffs: list[int],
    n_reactants: int,
    elements: list[str],
) -> tuple[list[dict], dict]:
    """逐元素与电荷核算：左右两侧各自的原子总数 / 总电荷。"""
    per_element = []
    for el in elements:
        left = sum(compositions[j].get(el, 0) * coeffs[j] for j in range(n_reactants))
        right = sum(
            compositions[j].get(el, 0) * coeffs[j]
            for j in range(n_reactants, len(compositions))
        )
        per_element.append(
            {"element": el, "reactants": left, "products": right, "balanced": left == right}
        )
    left_q = sum(charges[j] * coeffs[j] for j in range(n_reactants))
    right_q = sum(
        charges[j] * coeffs[j] for j in range(n_reactants, len(compositions))
    )
    charge = {"reactants": left_q, "products": right_q, "balanced": left_q == right_q}
    return per_element, charge


def format_equation(
    reactants: list[str], products: list[str], coeffs: list[int]
) -> str:
    def side(formulas: list[str], cs: list[int]) -> str:
        return " + ".join(
            f if c == 1 else f"{c}{f}" for f, c in zip(formulas, cs)
        )

    n_r = len(reactants)
    return f"{side(reactants, coeffs[:n_r])} -> {side(products, coeffs[n_r:])}"


def _parse_all(reactants: list[str], products: list[str]):
    formulas = list(reactants) + list(products)
    parsed = [parse_formula(f) for f in formulas]
    compositions = [c for c, _ in parsed]
    charges = [q for _, q in parsed]
    return compositions, charges


def balance_equation(reactants: list[str], products: list[str]) -> dict:
    """配平主流程：构建矩阵 → 零空间 → 一维且全正 → 最小整数系数。"""
    n_r = len(reactants)
    compositions, charges = _parse_all(reactants, products)
    elements, matrix = build_matrix(compositions, charges, n_r)
    basis = null_space(matrix)
    if len(basis) != 1:
        raise BalanceError(
            f"零空间维数为 {len(basis)}（必须恰好为 1 才能唯一配平）",
            nullity=len(basis),
        )
    vec = basis[0]
    if all(x < 0 for x in vec):
        vec = [-x for x in vec]
    if not all(x > 0 for x in vec):
        raise BalanceError(
            "零空间基向量无法化为全正系数（存在异号或零分量），"
            "该反应式不能用全部给定物质配平为全正整数系数",
            nullity=1,
        )
    coeffs = to_integer_coeffs(vec)
    per_element, charge = account(compositions, charges, coeffs, n_r, elements)
    return {
        "equation": format_equation(reactants, products, coeffs),
        "coefficients": coeffs,
        "reactant_coefficients": coeffs[:n_r],
        "product_coefficients": coeffs[n_r:],
        "elements": per_element,
        "charge": charge,
        "balanced": True,
    }


def verify_coefficients(
    reactants: list[str],
    products: list[str],
    reactant_coefficients: list[int],
    product_coefficients: list[int],
) -> dict:
    """校验给定系数是否满足逐元素与电荷守恒。"""
    n_r = len(reactants)
    compositions, charges = _parse_all(reactants, products)
    elements, _ = build_matrix(compositions, charges, n_r)
    coeffs = list(reactant_coefficients) + list(product_coefficients)
    per_element, charge = account(compositions, charges, coeffs, n_r, elements)
    balanced = charge["balanced"] and all(e["balanced"] for e in per_element)
    g = 0
    for v in coeffs:
        g = gcd(g, v)
    g = abs(g) or 1
    return {
        "balanced": balanced,
        "equation": format_equation(reactants, products, coeffs),
        "elements": per_element,
        "charge": charge,
        "gcd": g,
        "is_reduced": g == 1,
        "message": "配平正确" if balanced else "配平不正确：存在不守恒的元素或电荷",
    }
