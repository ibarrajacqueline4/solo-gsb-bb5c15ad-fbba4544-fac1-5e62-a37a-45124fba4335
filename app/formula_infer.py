"""经验式 / 分子式推断：由元素质量百分比反推最简式。

流程：
  1. 校验百分比（均为正、总和与 100 的偏差不超过 tolerance）；
  2. 复用 stoichiometry 的内置原子量表，把质量百分比换算为摩尔比
     （各元素百分比 / 原子量，再除以最小值归一化）；
  3. 在 1～12 倍尺度内把摩尔比四舍五入（ROUND_HALF_UP）为正整数，
     约去最大公约数得到互质下标候选（去重）；
  4. 用候选式反算各元素理论质量百分比，保留最大绝对偏差
     不超过 tolerance 的候选；偏差最小者胜出，并列时依次取
     总原子数较少、Hill 顺序字符串较小者；
  5. 提供实测摩尔质量时，取最接近实测值的正整数倍作为分子式，
     其相对偏差（%）须不超过 tolerance，否则判为无解。

全程 Decimal 计算（内部精度 50），输出按请求 precision 四舍五入
为定长小数字符串。无合格候选或摩尔质量不匹配时抛 NoSolutionError
（对应 HTTP 422）；输入数据本身不合法时抛 FormulaInferError
（对应 HTTP 400）。
"""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP, localcontext
from functools import reduce
from math import gcd

from .stoichiometry import ATOMIC_WEIGHTS, _round, _SUPPORTED_ELEMENTS

# 摩尔比放大倍率上限（含）
MAX_SCALE = 12


class FormulaInferError(ValueError):
    """百分比 / 摩尔质量等输入数据不合法（对应 HTTP 400）。"""


class NoSolutionError(ValueError):
    """容差内无合格候选式，或摩尔质量无匹配的正整数倍（对应 HTTP 422）。"""

    def __init__(self, message: str, **context: object):
        super().__init__(message)
        self.message = message
        self.context = context


def hill_order(elements: list[str]) -> list[str]:
    """Hill 顺序：含 C 时 C 优先、H 次之，其余元素按字母序。"""
    rest = sorted(elements)
    ordered = []
    if "C" in rest:
        ordered.append("C")
        rest.remove("C")
        if "H" in rest:
            ordered.append("H")
            rest.remove("H")
    return ordered + rest


def hill_string(subscripts: dict[str, int]) -> str:
    """按 Hill 顺序拼化学式字符串（下标为 1 时省略）。"""
    parts = []
    for el in hill_order(list(subscripts)):
        n = subscripts[el]
        parts.append(el if n == 1 else f"{el}{n}")
    return "".join(parts)


def _candidate_subscripts(
    ratios: dict[str, Decimal], elements: list[str]
) -> list[dict[str, int]]:
    """1～MAX_SCALE 倍尺度内生成互质正整数下标候选（按尺度升序去重）。

    摩尔比已按最小值归一化（均 ≥ 1），故任意尺度下四舍五入结果
    必为正整数；仍保留防御性检查。
    """
    seen: set[tuple[int, ...]] = set()
    candidates: list[dict[str, int]] = []
    for scale in range(1, MAX_SCALE + 1):
        subs = []
        for el in elements:
            n = int((ratios[el] * scale).to_integral_value(rounding=ROUND_HALF_UP))
            if n < 1:
                break
            subs.append(n)
        else:
            g = reduce(gcd, subs)
            key = tuple(s // g for s in subs)
            if key not in seen:
                seen.add(key)
                candidates.append(dict(zip(elements, key)))
    return candidates


def _validate_percentages(percentages: dict[str, Decimal], tolerance: Decimal) -> None:
    for el, pct in percentages.items():
        if el not in ATOMIC_WEIGHTS:
            raise FormulaInferError(
                f"元素 {el} 不在内置原子量表（支持：{_SUPPORTED_ELEMENTS}）"
            )
        if not pct.is_finite() or pct <= 0:
            raise FormulaInferError(f"元素 {el} 的质量百分比必须为正数: {pct}")
    total = sum(percentages.values(), Decimal(0))
    if abs(total - 100) > tolerance:
        raise FormulaInferError(
            f"质量百分比总和 {total} 与 100 的偏差超过 tolerance={tolerance}"
        )


def infer_formula(
    percentages: dict[str, Decimal],
    tolerance: Decimal,
    molar_mass: Decimal | None = None,
    precision: int = 4,
) -> dict:
    """主流程：百分比校验 → 摩尔比 → 候选生成 → 偏差筛选 → 分子式推断。"""
    _validate_percentages(percentages, tolerance)
    if molar_mass is not None and (not molar_mass.is_finite() or molar_mass <= 0):
        raise FormulaInferError(f"实测摩尔质量必须为正数: {molar_mass}")

    with localcontext() as ctx:
        # 内部精度远高于输出精度（2～8 位），避免中间舍入误差
        ctx.prec = 50

        elements = sorted(percentages)
        moles = {el: percentages[el] / ATOMIC_WEIGHTS[el] for el in elements}
        min_moles = min(moles.values())
        ratios = {el: moles[el] / min_moles for el in elements}

        best = None  # (排序键, 下标, 理论百分比, 偏差, 经验式摩尔质量, Hill 字符串)
        nearest = None  # 偏差最小的候选（用于 422 报错信息）
        for subs in _candidate_subscripts(ratios, elements):
            emp_mass = sum(ATOMIC_WEIGHTS[el] * n for el, n in subs.items())
            theo = {
                el: ATOMIC_WEIGHTS[el] * subs[el] / emp_mass * 100 for el in elements
            }
            devs = {el: abs(theo[el] - percentages[el]) for el in elements}
            max_dev = max(devs.values())
            hill = hill_string(subs)
            if nearest is None or max_dev < nearest[0]:
                nearest = (max_dev, hill)
            if max_dev > tolerance:
                continue
            # 偏差最小优先；并列时依次取总原子数较少、Hill 字符串较小者
            key = (max_dev, sum(subs.values()), hill)
            if best is None or key < best[0]:
                best = (key, subs, theo, devs, emp_mass, hill)

        if best is None:
            context: dict[str, object] = {"tolerance": str(tolerance)}
            if nearest is not None:
                context["nearest_candidate"] = nearest[1]
                context["nearest_max_deviation"] = _round(nearest[0], 4)
            raise NoSolutionError(
                f"1～{MAX_SCALE} 倍尺度内无候选式满足最大绝对偏差 ≤ tolerance={tolerance}",
                **context,
            )

        _, subs, theo, devs, emp_mass, emp_hill = best

        molecular_formula = None
        multiplier = None
        molecular_mass = None
        relative_deviation = None
        if molar_mass is not None:
            # 最接近实测值的正整数倍
            k = int(
                (molar_mass / emp_mass).to_integral_value(rounding=ROUND_HALF_UP)
            )
            k = max(k, 1)
            molecular_mass = emp_mass * k
            relative_deviation = abs(molecular_mass - molar_mass) / molar_mass * 100
            if relative_deviation > tolerance:
                raise NoSolutionError(
                    f"经验式 {emp_hill} 的 {k} 倍摩尔质量 {molecular_mass} "
                    f"与实测值 {molar_mass} 的相对偏差超过 tolerance={tolerance}%",
                    empirical_formula=emp_hill,
                    multiplier=k,
                    theoretical_molar_mass=_round(molecular_mass, 4),
                    relative_deviation_percent=_round(relative_deviation, 4),
                    tolerance=str(tolerance),
                )
            multiplier = k
            molecular_formula = hill_string({el: n * k for el, n in subs.items()})

        ordered = hill_order(elements)
        return {
            "empirical_formula": emp_hill,
            "subscripts": {el: subs[el] for el in ordered},
            "theoretical_percentages": {
                el: _round(theo[el], precision) for el in ordered
            },
            "deviations": {el: _round(devs[el], precision) for el in ordered},
            "max_deviation": _round(max(devs.values()), precision),
            "empirical_molar_mass": _round(emp_mass, precision),
            "molecular_formula": molecular_formula,
            "multiplier": multiplier,
            "molecular_molar_mass": (
                _round(molecular_mass, precision) if molecular_mass is not None else None
            ),
            "molar_mass_relative_deviation": (
                _round(relative_deviation, precision)
                if relative_deviation is not None
                else None
            ),
        }
