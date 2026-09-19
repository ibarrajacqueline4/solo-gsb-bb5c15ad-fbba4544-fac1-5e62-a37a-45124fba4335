"""请求 / 响应模型。"""
from __future__ import annotations

from decimal import Decimal
from typing import Literal

from pydantic import BaseModel, Field


class BalanceRequest(BaseModel):
    reactants: list[str] = Field(
        min_length=1,
        description="反应物分子式列表",
        examples=[["KMnO4", "HCl"]],
    )
    products: list[str] = Field(
        min_length=1,
        description="生成物分子式列表",
        examples=[["KCl", "MnCl2", "H2O", "Cl2"]],
    )


class ElementAccount(BaseModel):
    element: str
    reactants: int = Field(description="反应物侧该元素原子总数")
    products: int = Field(description="生成物侧该元素原子总数")
    balanced: bool


class ChargeAccount(BaseModel):
    reactants: int = Field(description="反应物侧总电荷")
    products: int = Field(description="生成物侧总电荷")
    balanced: bool


class BalanceResponse(BaseModel):
    equation: str
    coefficients: list[int] = Field(description="全部系数（先反应物后生成物）")
    reactant_coefficients: list[int]
    product_coefficients: list[int]
    elements: list[ElementAccount] = Field(description="逐元素核算")
    charge: ChargeAccount = Field(description="电荷核算")
    balanced: bool


class VerifyRequest(BaseModel):
    reactants: list[str] = Field(min_length=1, examples=[["H2", "O2"]])
    products: list[str] = Field(min_length=1, examples=[["H2O"]])
    reactant_coefficients: list[int] = Field(min_length=1, examples=[[2, 1]])
    product_coefficients: list[int] = Field(min_length=1, examples=[[2]])


class VerifyResponse(BaseModel):
    balanced: bool
    equation: str
    elements: list[ElementAccount]
    charge: ChargeAccount
    gcd: int = Field(description="系数的最大公约数")
    is_reduced: bool = Field(description="系数是否已约简（gcd == 1）")
    message: str


class FeedItem(BaseModel):
    formula: str = Field(
        description="投料物质分子式（必须是反应物之一）", examples=["H2"]
    )
    amount: Decimal = Field(description="投料量，必须为正数", examples=["2.5"])
    unit: Literal["mol", "g"] = Field(description="投料单位，仅支持 mol 或 g")


class StoichiometryRequest(BaseModel):
    reactants: list[str] = Field(
        min_length=1,
        description="反应物分子式列表",
        examples=[["H2", "O2"]],
    )
    products: list[str] = Field(
        min_length=1,
        description="生成物分子式列表",
        examples=[["H2O"]],
    )
    feed: list[FeedItem] = Field(
        description="各反应物的投料（每种反应物恰好一条，单位 mol 或 g）"
    )
    precision: int = Field(
        default=4,
        ge=2,
        le=8,
        description="输出小数位数（四舍五入），2～8，默认 4",
    )


class ProductYield(BaseModel):
    formula: str
    coefficient: int = Field(description="配平系数")
    theoretical_moles: str = Field(description="理论产量（mol，定精度小数字符串）")
    theoretical_grams: str = Field(description="理论产量（g，定精度小数字符串）")


class ReactantRemainder(BaseModel):
    formula: str
    coefficient: int = Field(description="配平系数")
    feed_moles: str = Field(description="投料量（mol）")
    feed_grams: str = Field(description="投料量（g）")
    consumed_moles: str = Field(description="消耗量（mol）")
    remaining_moles: str = Field(description="剩余量（mol）")
    remaining_grams: str = Field(description="剩余量（g）")
    limiting: bool = Field(description="是否为限制反应物")


class StoichiometryResponse(BaseModel):
    equation: str
    coefficients: list[int] = Field(description="全部系数（先反应物后生成物）")
    reactant_coefficients: list[int]
    product_coefficients: list[int]
    precision: int = Field(description="实际采用的小数位数")
    extent: str = Field(description="理论反应进度 ξ（mol，定精度小数字符串）")
    limiting_reactants: list[str] = Field(
        description="限制反应物（反应进度最小，可并列多个）"
    )
    products: list[ProductYield] = Field(description="各生成物理论产量")
    reactants: list[ReactantRemainder] = Field(description="各反应物投料/消耗/剩余")


class FormulaInferRequest(BaseModel):
    percentages: dict[str, Decimal] = Field(
        min_length=1,
        description=(
            "元素质量百分比（元素符号 -> 百分比）；须均为正数，"
            "且总和与 100 的偏差不超过 tolerance"
        ),
        examples=[{"C": "40.0", "H": "6.7", "O": "53.3"}],
    )
    tolerance: Decimal = Field(
        default=Decimal("0.5"),
        ge=0,
        description=(
            "容差：百分比总和与 100 的允许偏差、候选式理论百分比的最大绝对偏差上限"
            "（均为百分点），以及摩尔质量相对偏差上限（%）；默认 0.5"
        ),
    )
    molar_mass: Decimal | None = Field(
        default=None,
        description="实测摩尔质量（g/mol）；提供时按最接近的正整数倍给出分子式",
        examples=["180.0"],
    )
    precision: int = Field(
        default=4,
        ge=2,
        le=8,
        description="输出小数位数（四舍五入），2～8，默认 4",
    )


class FormulaInferResponse(BaseModel):
    empirical_formula: str = Field(description="经验式（Hill 顺序）")
    subscripts: dict[str, int] = Field(description="经验式下标（Hill 顺序）")
    theoretical_percentages: dict[str, str] = Field(
        description="候选式反算的各元素质量百分比（定精度小数字符串）"
    )
    deviations: dict[str, str] = Field(
        description="逐元素绝对偏差（百分点，定精度小数字符串）"
    )
    max_deviation: str = Field(description="最大绝对偏差（百分点）")
    empirical_molar_mass: str = Field(description="经验式摩尔质量（g/mol）")
    molecular_formula: str | None = Field(
        description="分子式（仅在提供 molar_mass 且偏差合格时给出）"
    )
    multiplier: int | None = Field(description="分子式相对经验式的正整数倍")
    molecular_molar_mass: str | None = Field(description="分子式理论摩尔质量（g/mol）")
    molar_mass_relative_deviation: str | None = Field(
        description="分子式理论摩尔质量与实测值的相对偏差（%）"
    )
