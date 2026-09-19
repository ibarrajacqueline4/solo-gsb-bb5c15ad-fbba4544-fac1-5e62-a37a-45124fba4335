"""化学计量：限制反应物判定、理论产量与剩余量计算。

复用 parser 的分子式解析与 balancer 的配平结果，按内置原子量表
（H、C、N、O、Na、Mg、Al、S、Cl、K、Ca、Fe、Cu、Mn）换算物质的量。
全程使用 Decimal 计算，输出统一按请求 precision 四舍五入
（ROUND_HALF_UP）为定长小数字符串。
"""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP, localcontext
from typing import TYPE_CHECKING

from .balancer import balance_equation
from .parser import parse_formula

if TYPE_CHECKING:  # 避免在逻辑层硬依赖 API 模型
    from .schemas import FeedItem


class StoichiometryError(ValueError):
    """投料或原子量相关的输入错误（对应 HTTP 400）。"""


# 内置原子量表（g/mol，IUPAC 标准原子量常用值）
ATOMIC_WEIGHTS: dict[str, Decimal] = {
    symbol: Decimal(weight)
    for symbol, weight in {
        "H": "1.008",
        "C": "12.011",
        "N": "14.007",
        "O": "15.999",
        "Na": "22.98977",
        "Mg": "24.305",
        "Al": "26.982",
        "S": "32.06",
        "Cl": "35.45",
        "K": "39.0983",
        "Ca": "40.078",
        "Fe": "55.845",
        "Cu": "63.546",
        "Mn": "54.938",
    }.items()
}

_SUPPORTED_ELEMENTS = "、".join(ATOMIC_WEIGHTS)


def molar_mass(formula: str) -> Decimal:
    """按内置原子量表计算摩尔质量（g/mol）。

    离子电荷忽略不计（电子质量远小于原子质量）；电子 e 本身质量为 0。
    含表外元素时抛 StoichiometryError。
    """
    counts, _ = parse_formula(formula)
    mass = Decimal(0)
    for element, n in counts.items():
        weight = ATOMIC_WEIGHTS.get(element)
        if weight is None:
            raise StoichiometryError(
                f"物质 {formula!r} 含有元素 {element}，"
                f"不在内置原子量表（支持：{_SUPPORTED_ELEMENTS}）"
            )
        mass += weight * n
    return mass


def _round(value: Decimal, precision: int) -> str:
    """按 precision 四舍五入为定长小数字符串（负零归一化为正零）。"""
    quantum = Decimal(1).scaleb(-precision)
    rounded = value.quantize(quantum, rounding=ROUND_HALF_UP)
    if rounded == 0:
        rounded = rounded.copy_abs()
    return format(rounded, "f")


def compute_stoichiometry(
    reactants: list[str],
    products: list[str],
    feed: list[FeedItem],
    precision: int,
) -> dict:
    """主流程：配平 → 投料校验 → 物质的量换算 → 反应进度 → 产量与剩余量。

    FormulaError / BalanceError 继续向上抛，由路由层映射为 4xx。
    """
    balanced = balance_equation(reactants, products)
    r_coeffs = balanced["reactant_coefficients"]
    p_coeffs = balanced["product_coefficients"]

    # ---- 投料校验：重复投料 / 非反应物 / 非正数量 / 缺失投料 ----
    amounts: dict[str, tuple[Decimal, str]] = {}
    for item in feed:
        formula = item.formula.strip()
        if formula in amounts:
            raise StoichiometryError(f"同一物质重复投料: {formula!r}")
        if formula not in reactants:
            raise StoichiometryError(f"投料物质 {formula!r} 不在反应物列表中")
        if not item.amount.is_finite() or item.amount <= 0:
            raise StoichiometryError(
                f"投料量必须为正数: {formula!r} 的 amount 为 {item.amount}"
            )
        amounts[formula] = (item.amount, item.unit)
    missing = [r for r in reactants if r not in amounts]
    if missing:
        raise StoichiometryError(f"缺少反应物投料: {', '.join(missing)}")

    # ---- 摩尔质量：任一物质含表外元素即报错 ----
    masses = {f: molar_mass(f) for f in list(reactants) + list(products)}

    with localcontext() as ctx:
        # 内部精度远高于输出精度（2～8 位），避免中间舍入误差
        ctx.prec = 50

        feed_moles: dict[str, Decimal] = {}
        feed_grams: dict[str, Decimal] = {}
        for formula in reactants:
            amount, unit = amounts[formula]
            mass = masses[formula]
            if unit == "mol":
                feed_moles[formula] = amount
                feed_grams[formula] = amount * mass
            else:
                if mass == 0:
                    raise StoichiometryError(
                        f"物质 {formula!r} 摩尔质量为 0（如电子 e），不能以 g 投料"
                    )
                feed_moles[formula] = amount / mass
                feed_grams[formula] = amount

        # 反应进度 ξ_i = n_i / ν_i，最小者（可并列）为限制反应物
        extents = {
            formula: feed_moles[formula] / r_coeffs[i]
            for i, formula in enumerate(reactants)
        }
        extent = min(extents.values())
        limiting = [f for f in reactants if extents[f] == extent]

        product_rows = []
        for i, formula in enumerate(products):
            moles = p_coeffs[i] * extent
            product_rows.append(
                {
                    "formula": formula,
                    "coefficient": p_coeffs[i],
                    "theoretical_moles": _round(moles, precision),
                    "theoretical_grams": _round(moles * masses[formula], precision),
                }
            )

        reactant_rows = []
        for i, formula in enumerate(reactants):
            consumed = r_coeffs[i] * extent
            if formula in limiting:
                remaining = Decimal(0)
            else:
                # 除法舍入噪声下 remaining 不会真正为负，防御性截断
                remaining = max(feed_moles[formula] - consumed, Decimal(0))
            reactant_rows.append(
                {
                    "formula": formula,
                    "coefficient": r_coeffs[i],
                    "feed_moles": _round(feed_moles[formula], precision),
                    "feed_grams": _round(feed_grams[formula], precision),
                    "consumed_moles": _round(consumed, precision),
                    "remaining_moles": _round(remaining, precision),
                    "remaining_grams": _round(remaining * masses[formula], precision),
                    "limiting": formula in limiting,
                }
            )

        return {
            "equation": balanced["equation"],
            "coefficients": balanced["coefficients"],
            "reactant_coefficients": r_coeffs,
            "product_coefficients": p_coeffs,
            "precision": precision,
            "extent": _round(extent, precision),
            "limiting_reactants": limiting,
            "products": product_rows,
            "reactants": reactant_rows,
        }
