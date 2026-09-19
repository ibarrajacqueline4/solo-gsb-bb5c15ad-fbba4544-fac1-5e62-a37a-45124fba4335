"""化学式解析器。

支持语法：
  - 元素符号与下标：H2O、Fe2O3
  - 嵌套括号（()、[]、{} 可混用嵌套）：Ca(OH)2、K4[Fe(CN)6]
  - 点号结晶水（.、·、* 均可）：CuSO4.5H2O、Na2CO3·10H2O
  - 末尾离子电荷：
      * 插入符形式（推荐，无歧义）：SO4^2-、Fe^3+、NH4^+
      * 纯符号尾（每个 +/- 计 1 单位电荷）：Na+、Cl-、Ca++、SO4--
      * 注意：不带 ^ 时，符号前的数字属于前一个元素/基团，
        即 NH4+ 解析为 NH4 带 +1 电荷；要表示 Ca2+ 请写 Ca^2+。
  - 电子（用于半反应）：e、e-
"""
from __future__ import annotations

import re
from collections import Counter

# 118 种元素符号，用于校验拼写
ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}


class FormulaError(ValueError):
    """化学式无法解析。"""


_OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}"}
_CLOSES = set(_OPEN_TO_CLOSE.values())
_DOTS = {".", "·", "*"}

_TOKEN_RE = re.compile(r"[A-Z][a-z]?|\d+|[()\[\]{}]|[.·*]")
_CARET_CHARGE_RE = re.compile(r"\^(\d*)([+-])$")
_TRAILING_SIGN_RE = re.compile(r"([+-]+)$")


def _extract_charge(text: str) -> tuple[int, str]:
    """剥离末尾电荷，返回 (电荷数, 剩余化学式主体)。"""
    m = _CARET_CHARGE_RE.search(text)
    if m:
        magnitude = int(m.group(1)) if m.group(1) else 1
        sign = 1 if m.group(2) == "+" else -1
        return sign * magnitude, text[: m.start()]
    m = _TRAILING_SIGN_RE.search(text)
    if m:
        signs = m.group(1)
        if set(signs) == {"+"}:
            return len(signs), text[: m.start()]
        if set(signs) == {"-"}:
            return -len(signs), text[: m.start()]
        raise FormulaError(f"电荷符号混乱: {signs!r}")
    return 0, text


def _tokenize(text: str) -> list[str]:
    tokens = []
    pos = 0
    for m in _TOKEN_RE.finditer(text):
        if m.start() != pos:
            raise FormulaError(f"无法识别的字符 {text[pos]!r}（位置 {pos}）")
        tokens.append(m.group(0))
        pos = m.end()
    if pos != len(text):
        raise FormulaError(f"无法识别的字符 {text[pos]!r}（位置 {pos}）")
    return tokens


def _read_number(tokens: list[str], i: int) -> tuple[int, int]:
    if i < len(tokens) and tokens[i][0].isdigit():
        n = int(tokens[i])
        if n == 0:
            raise FormulaError("下标/倍数不能为 0")
        return n, i + 1
    return 1, i


def _parse_group(tokens: list[str], i: int) -> tuple[Counter, int]:
    """解析一个基团，遇到未匹配的闭合括号或 token 耗尽时返回。"""
    counts: Counter = Counter()
    while i < len(tokens):
        tok = tokens[i]
        if tok in _CLOSES:
            return counts, i
        if tok in _OPEN_TO_CLOSE:
            sub, j = _parse_group(tokens, i + 1)
            if j >= len(tokens) or tokens[j] != _OPEN_TO_CLOSE[tok]:
                raise FormulaError(f"括号 {tok!r} 未正确闭合")
            n, i = _read_number(tokens, j + 1)
            for el, c in sub.items():
                counts[el] += c * n
        elif tok in _DOTS:
            # 结晶水：点号后的可选倍数作用于紧随其后的整个基团
            n, j = _read_number(tokens, i + 1)
            sub, i = _parse_group(tokens, j)
            for el, c in sub.items():
                counts[el] += c * n
            return counts, i
        elif tok[0].isdigit():
            raise FormulaError(f"数字 {tok!r} 前面没有元素或括号")
        else:
            if tok not in ELEMENTS:
                raise FormulaError(f"未知元素符号: {tok!r}")
            n, i = _read_number(tokens, i + 1)
            counts[tok] += n
    return counts, i


def parse_formula(text: str) -> tuple[Counter, int]:
    """解析化学式，返回 (元素计数 Counter, 电荷数)。"""
    if not isinstance(text, str) or not text.strip():
        raise FormulaError("化学式不能为空")
    text = text.strip()
    charge, body = _extract_charge(text)
    if body == "e":  # 电子（半反应用）
        return Counter(), charge
    tokens = _tokenize(body)
    if not tokens:
        raise FormulaError("化学式不能为空")
    counts, i = _parse_group(tokens, 0)
    if i != len(tokens):
        raise FormulaError(f"多余的闭合括号 {tokens[i]!r}")
    if not counts:
        raise FormulaError("化学式中不含任何元素")
    return counts, charge
