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
    from .schemas import FeedItem, MaterialSpec


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


def compute_batch_plan(
    reactants: list[str],
    products: list[str],
    target_product: str,
    target_mass: Decimal,
    yield_fraction: Decimal,
    materials: list[MaterialSpec],
    precision: int,
) -> dict:
    """批次投料规划：目标产量与收率 → 理论反应进度 → 逐原料称量。

    流程：配平 → 目标物与参数校验 → 原料配置校验 → 由目标纯品质量与收率
    反推理论反应进度 → 按配平系数求各反应物的化学计量摩尔数 → 叠加过量
    百分比 → 按纯度折算实际称量克数。
    FormulaError / BalanceError 继续向上抛，由路由层映射为 4xx。
    """
    balanced = balance_equation(reactants, products)
    r_coeffs = balanced["reactant_coefficients"]
    p_coeffs = balanced["product_coefficients"]

    # ---- 目标物与全局参数校验 ----
    if target_product not in products:
        raise StoichiometryError(
            f"目标生成物 {target_product!r} 不在生成物列表中"
        )
    if not target_mass.is_finite() or target_mass <= 0:
        raise StoichiometryError(
            f"目标纯品质量必须为正数: target_mass 为 {target_mass}"
        )
    if not yield_fraction.is_finite() or yield_fraction <= 0 or yield_fraction > 1:
        raise StoichiometryError(
            f"预计收率必须在 (0,1] 区间内: yield_fraction 为 {yield_fraction}"
        )

    # ---- 原料配置校验：重复 / 非反应物 / 纯度与过量越界 / 缺失 ----
    specs: dict[str, MaterialSpec] = {}
    for spec in materials:
        formula = spec.formula.strip()
        if formula in specs:
            raise StoichiometryError(f"同一原料重复配置: {formula!r}")
        if formula not in reactants:
            raise StoichiometryError(f"原料配置 {formula!r} 不在反应物列表中")
        if not spec.purity.is_finite() or spec.purity <= 0 or spec.purity > 1:
            raise StoichiometryError(
                f"纯度必须在 (0,1] 区间内: {formula!r} 的 purity 为 {spec.purity}"
            )
        if not spec.excess_percent.is_finite() or spec.excess_percent < 0:
            raise StoichiometryError(
                f"过量百分比必须为非负数: "
                f"{formula!r} 的 excess_percent 为 {spec.excess_percent}"
            )
        specs[formula] = spec
    missing = [r for r in reactants if r not in specs]
    if missing:
        raise StoichiometryError(f"缺少原料配置: {', '.join(missing)}")

    # ---- 摩尔质量：任一物质含表外元素即报错 ----
    masses = {f: molar_mass(f) for f in list(reactants) + list(products)}
    target_molar_mass = masses[target_product]
    if target_molar_mass == 0:
        raise StoichiometryError(
            f"目标生成物 {target_product!r} 摩尔质量为 0（如电子 e），不能按质量规划"
        )

    with localcontext() as ctx:
        # 内部精度远高于输出精度（2～8 位），避免中间舍入误差
        ctx.prec = 50

        # 由目标纯品质量与收率反推理论产量，再除以配平系数得反应进度
        target_theoretical_grams = target_mass / yield_fraction
        target_theoretical_moles = target_theoretical_grams / target_molar_mass
        extent = target_theoretical_moles / p_coeffs[products.index(target_product)]

        material_rows = []
        for i, formula in enumerate(reactants):
            spec = specs[formula]
            stoich_moles = r_coeffs[i] * extent
            excess_moles = stoich_moles * (1 + spec.excess_percent / 100)
            weigh = excess_moles * masses[formula] / spec.purity
            material_rows.append(
                {
                    "formula": formula,
                    "coefficient": r_coeffs[i],
                    "purity": str(spec.purity),
                    "excess_percent": str(spec.excess_percent),
                    "stoichiometric_moles": _round(stoich_moles, precision),
                    "moles_with_excess": _round(excess_moles, precision),
                    "weigh_grams": _round(weigh, precision),
                }
            )

        return {
            "equation": balanced["equation"],
            "coefficients": balanced["coefficients"],
            "reactant_coefficients": r_coeffs,
            "product_coefficients": p_coeffs,
            "precision": precision,
            "target_product": target_product,
            "target_mass": str(target_mass),
            "yield_fraction": str(yield_fraction),
            "target_theoretical_moles": _round(target_theoretical_moles, precision),
            "target_theoretical_grams": _round(target_theoretical_grams, precision),
            "extent": _round(extent, precision),
            "materials": material_rows,
        }
