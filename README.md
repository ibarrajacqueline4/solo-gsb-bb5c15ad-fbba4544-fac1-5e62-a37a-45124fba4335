# 化学方程式精确配平 API

基于 FastAPI 的化学方程式配平服务。服务端为反应物/生成物构建**元素守恒矩阵（含电荷行）**，用**有理数高斯消元**求零空间；仅当零空间为一维且基向量可化为全正时，返回最大公约数约简后的**最小正整数系数**及逐元素、电荷核算，否则返回 `422` 并说明原因（零空间维数或无法全正）。

另提供**化学计量接口**：在配平结果基础上，按内置原子量表（H、C、N、O、Na、Mg、Al、S、Cl、K、Ca、Fe、Cu、Mn）换算投料物质的量，以反应进度最小值判定（可并列的）限制反应物，返回理论反应进度、各生成物理论产量与各反应物剩余量。全程 `Decimal` 计算，输出按请求 `precision` 四舍五入。

以及**经验式推断接口**：由元素质量百分比反推经验式——先按内置原子量表换算摩尔比，再在 1～12 倍尺度内生成互质正整数下标候选，用候选式反算理论百分比并保留最大绝对偏差不超过 `tolerance` 者（并列时依次取总原子数较少、Hill 顺序字符串较小者）；提供实测摩尔质量时，按最接近的正整数倍给出分子式。

## 快速启动（单条命令）

```bash
docker compose up --build
```

依赖随容器启动自动安装（`pip install -r requirements.txt` 后再拉起 uvicorn）。

- 服务地址：http://localhost:8000
- 交互式文档（Swagger UI）：http://localhost:8000/docs
- 健康检查：http://localhost:8000/health

## 配置

| 环境变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PORT` | `8000` | 容器内服务监听端口 |

宿主机端口映射在 `docker-compose.yml` 的 `ports` 中调整（默认 `8000:8000`）。

## 分子式语法

| 语法 | 示例 |
| --- | --- |
| 元素与下标 | `H2O`、`Fe2O3` |
| 嵌套括号（`()` `[]` `{}` 可混用） | `Ca(OH)2`、`K4[Fe(CN)6]` |
| 点号结晶水（`.` `·` `*`） | `CuSO4.5H2O`、`Na2CO3·10H2O` |
| 插入符电荷（推荐，无歧义） | `SO4^2-`、`Fe^3+`、`NH4^+` |
| 纯符号尾电荷（每个 `+`/`-` 计 1） | `Na+`、`Cl-`、`Ca++`、`SO4--` |
| 电子（半反应） | `e-` |

> 注意：不带 `^` 时，符号前的数字属于前一个元素/基团——`NH4+` 解析为 NH₄⁺（电荷 +1）；要表示 Ca²⁺ 请写 `Ca^2+`。

## API

### `POST /api/v1/balance` — 配平

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["KMnO4", "HCl"], "products": ["KCl", "MnCl2", "H2O", "Cl2"]}'
```

响应（`200`）：

```json
{
  "equation": "2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 8H2O + 5Cl2",
  "coefficients": [2, 16, 2, 2, 8, 5],
  "reactant_coefficients": [2, 16],
  "product_coefficients": [2, 2, 8, 5],
  "elements": [
    {"element": "Cl", "reactants": 16, "products": 16, "balanced": true},
    {"element": "H",  "reactants": 16, "products": 16, "balanced": true},
    {"element": "K",  "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "Mn", "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "O",  "reactants": 8,  "products": 8,  "balanced": true}
  ],
  "charge": {"reactants": 0, "products": 0, "balanced": true},
  "balanced": true
}
```

离子方程式（含电荷行守恒）：

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4-", "H+", "e-"], "products": ["Mn^2+", "H2O"]}'
# -> MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O
```

无法配平时返回 `422`，并说明零空间维数或无法全正：

```json
{
  "detail": {
    "message": "零空间维数为 2（必须恰好为 1 才能唯一配平）",
    "null_space_dim": 2
  }
}
```

化学式语法错误返回 `400`。

### `POST /api/v1/verify` — 校验系数

```bash
curl -s http://localhost:8000/api/v1/verify \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "reactant_coefficients": [2, 1],
        "product_coefficients": [2]
      }'
```

响应包含逐元素与电荷核算、`balanced` 结论，以及系数的 `gcd` 与是否已约简（`is_reduced`）。

### `POST /api/v1/stoichiometry` — 化学计量

在配平结果基础上，根据各反应物投料（单位仅支持 `mol` 或 `g`）计算限制反应物、理论反应进度、各生成物理论产量与各反应物剩余量。

- 内置原子量表（g/mol）：H 1.008、C 12.011、N 14.007、O 15.999、Na 22.98977、Mg 24.305、Al 26.982、S 32.06、Cl 35.45、K 39.0983、Ca 40.078、Fe 55.845、Cu 63.546、Mn 54.938（离子电荷的电子质量忽略不计）。
- 每种反应物必须恰好投料一次；以反应进度 ξ = n/ν 的最小值判定限制反应物（可并列多个）。
- 全程 `Decimal` 计算，所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/stoichiometry \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "feed": [
          {"formula": "H2", "amount": 10,  "unit": "mol"},
          {"formula": "O2", "amount": 100, "unit": "g"}
        ],
        "precision": 4
      }'
```

响应（`200`）：

```json
{
  "equation": "2H2 + O2 -> 2H2O",
  "coefficients": [2, 1, 2],
  "reactant_coefficients": [2, 1],
  "product_coefficients": [2],
  "precision": 4,
  "extent": "3.1252",
  "limiting_reactants": ["O2"],
  "products": [
    {"formula": "H2O", "coefficient": 2, "theoretical_moles": "6.2504", "theoretical_grams": "112.6008"}
  ],
  "reactants": [
    {"formula": "H2", "coefficient": 2, "feed_moles": "10.0000", "feed_grams": "20.1600",
     "consumed_moles": "6.2504", "remaining_moles": "3.7496", "remaining_grams": "7.5592", "limiting": false},
    {"formula": "O2", "coefficient": 1, "feed_moles": "3.1252", "feed_grams": "100.0000",
     "consumed_moles": "3.1252", "remaining_moles": "0.0000", "remaining_grams": "0.0000", "limiting": true}
  ]
}
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析；缺失投料；同一物质重复投料；投料量非正数；投料物质不是反应物；元素无内置原子量 |
| `422` | 方程无法配平（零空间维数不为 1 或无法全正）；单位不是 `mol`/`g`；`precision` 超出 2～8 |

### `POST /api/v1/formula/infer` — 经验式/分子式推断

由元素质量百分比反推经验式：按内置原子量表把百分比换算为摩尔比，在 **1～12 倍尺度**内四舍五入并约简为**互质正整数下标**候选；用候选式反算各元素理论百分比，保留**最大绝对偏差 ≤ `tolerance`** 的候选，偏差最小者胜出（并列时依次取总原子数较少、Hill 顺序字符串较小者）。提供 `molar_mass` 时，取**最接近实测值的正整数倍**作为分子式，其相对偏差（%）同样须 ≤ `tolerance`。

- 百分比须均为正数，且总和与 100 的偏差不超过 `tolerance`（默认 `0.5`，单位：百分点）。
- 元素须在内置原子量表内（同化学计量接口）。
- 所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{
        "percentages": {"C": 40.0, "H": 6.7, "O": 53.3},
        "molar_mass": 180.0
      }'
```

响应（`200`）：

```json
{
  "empirical_formula": "CH2O",
  "subscripts": {"C": 1, "H": 2, "O": 1},
  "theoretical_percentages": {"C": "40.0020", "H": "6.7142", "O": "53.2838"},
  "deviations": {"C": "0.0020", "H": "0.0142", "O": "0.0162"},
  "max_deviation": "0.0162",
  "empirical_molar_mass": "30.0260",
  "molecular_formula": "C6H12O6",
  "multiplier": 6,
  "molecular_molar_mass": "180.1560",
  "molar_mass_relative_deviation": "0.0867"
}
```

不提供 `molar_mass` 时，分子式相关字段为 `null`：

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{"percentages": {"N": 36.86, "O": 63.14}, "tolerance": 0.5}'
# -> empirical_formula: N2O3（尺度 2 的候选比尺度 1 的 NO2 偏差更小）
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 百分比非正数；百分比总和与 100 的偏差超过 `tolerance`；元素无内置原子量；`molar_mass` 非正数 |
| `422` | 1～12 倍尺度内无候选式满足容差；摩尔质量无满足相对偏差容差的正整数倍；`tolerance` 为负；`precision` 超出 2～8 |

`422` 响应的 `detail` 附带诊断信息（最接近的候选式及其偏差，或摩尔质量倍数与相对偏差），例如：

```json
{
  "detail": {
    "message": "1～12 倍尺度内无候选式满足最大绝对偏差 ≤ tolerance=0.5",
    "tolerance": "0.5",
    "nearest_candidate": "CO",
    "nearest_max_deviation": "0.8811"
  }
}
```

## 本地开发

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload --port 8000
pytest                 # 运行测试
```

## 项目结构

```
app/
  parser.py         # 分子式解析：嵌套括号、结晶水、离子电荷
  balancer.py       # 守恒矩阵（含电荷行）、有理数 RREF、零空间、整数化
  stoichiometry.py  # 内置原子量表、投料换算、限制反应物、产量与剩余量（Decimal）
  formula_infer.py  # 质量百分比 → 摩尔比 → 互质下标候选 → 经验式/分子式推断（Decimal）
  schemas.py        # 请求/响应模型
  main.py           # FastAPI 路由
tests/test_api.py / tests/test_stoichiometry.py / tests/test_formula_infer.py
Dockerfile / docker-compose.yml
```
