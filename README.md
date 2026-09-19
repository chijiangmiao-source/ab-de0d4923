# 家系相位分析服务（Pedigree Phasing Service）

纯后端分析服务：对一次提交的核心家系（父母 + 1–6 名子女）未定相双等位基因型，
计算**孟德尔一致且相邻标记间传递链切换总数（交叉数）全局最小**的家系相位，
并给出规范解与逐项歧义掩码。无前端、无任何在线服务调用，全部计算在本地完成。

## 快速开始

```bash
docker compose up --build          # 默认宿主机端口 8000
API_PORT=9000 docker compose up    # 宿主机端口由环境变量 API_PORT 配置
```

健康检查：

```bash
curl http://localhost:${API_PORT:-8000}/health     # {"status":"ok"}
```

Dockerfile 内置 `HEALTHCHECK`，Compose 服务亦配置了健康检查，均探测 `GET /health`。

交互式 API 文档（OpenAPI/Swagger）：`http://localhost:${API_PORT:-8000}/docs`。

## API

### `POST /api/v1/phase`

请求体：

```json
{
  "markers": ["rs01", "rs02", "rs03"],
  "father":  ["0/1", "0/1", "0/0"],
  "mother":  ["0/0", "0/1", "0/0"],
  "children": [
    {"id": "proband", "genotypes": ["0/0", "0/1", "0/0"]}
  ]
}
```

| 字段 | 约束 |
| --- | --- |
| `markers` | 2–1500 个非空标记 id，**按顺序排列** |
| `father` / `mother` | 基因型数组，长度必须等于标记数 |
| `children` | 1–6 名子女；`id` 可缺省（自动命名为 `child1..child6`），`genotypes` 长度必须等于标记数 |

基因型写法（双等位，未定向）：`0/0`、`0/1`、`1/1`、`./.`（完全缺失）、
`0/.`、`./1`（半缺失）；`|` 可作为 `/` 的别名（仍按未定向处理）。

成功响应（`200`，字段顺序固定、重复调用逐字节一致）：

```json
{
  "status": "ok",
  "marker_count": 3,
  "child_count": 1,
  "min_crossovers": 0,
  "all_fixed": false,
  "solution": [
    {
      "index": 0, "id": "rs01",
      "father": {"alleles": [0, 1], "fixed": [true, true]},
      "mother": {"alleles": [0, 0], "fixed": [true, true]},
      "children": [
        {"index": 0, "id": "proband",
         "paternal": 0, "maternal": 0,
         "paternal_fixed": true, "maternal_fixed": false}
      ]
    }
  ]
}
```

语义：

- `father.alleles` / `mother.alleles`：该标记处父/母两条单倍型（0 号、1 号）携带的等位基因。
- `paternal` / `maternal`：该子女在该标记处分别继承父/母的哪一条单倍型（0 或 1）。
- `min_crossovers`：所有子女、双亲两侧、相邻标记间传递链切换的总数的最小值。
- `fixed` / `*_fixed`：**歧义掩码**——该项在**所有**最优解（满足方向固定规则、
  交叉数最小的全部解）中是否取同一值；`all_fixed` 为全表汇总。
- 规范解：全部最优解中，按 **标记 → 父亲 → 母亲 → 子女（按提交顺序，每人先父源后母源）**
  的相位与传递位序列字典序最小者。

### 错误返回

- `422` 非法结构（可定位）：`detail` 为数组，每项含 `loc`（出错位置路径，如
  `["body","father",3]`）、`msg`、`type`。涵盖：标记数越界、基因型格式非法、
  各体基因型数与标记数不一致、子女数越界、多余字段等。
- `409` 结构合法但无孟德尔一致解释：`detail` 指出**最早无解标记**及其基因型：

```json
{
  "detail": {
    "error": "mendelian_inconsistency",
    "message": "no Mendelian-consistent explanation exists; earliest unsolvable marker: index 1 (id 'm2')",
    "marker_index": 1,
    "marker_id": "m2",
    "genotypes": {"father": "0/0", "mother": "0/0", "children": ["1/1"]}
  }
}
```

## 判定规则精确定义

1. **可行性（孟德尔一致）**：每个标记独立判定——存在父母各自的有序等位基因对
   （与父/母该位点基因型一致）及每名子女的传递位 `(p, m)`，使子女实际继承的
   有序等位基因对与其基因型一致。父母相位与子女传递逐标记自由选择，因此
   **整体可行 ⟺ 每个标记各自可行**；最早不可行标记即最早无解标记。
2. **目标函数**：`Σ_子女 Σ_相邻标记 (父源传递位变化 + 母源传递位变化)` 最小化。
3. **方向固定**：父（母）两条单倍型以其**首个确定杂合位点**（双等位均被检出且不同，
   即 `0/1`/`1/0`；半缺失不算）固定 0 号方向——该位点 0 号单倍型携带等位基因 0。
   无杂合位点的亲本不施加该约束。
4. **规范解裁决**：在满足 1–3 的全部最优解中，取位序列
   `标记1: [父0,父1,母0,母1,子女1父源,子女1母源,...], 标记2: ...`
   的字典序最小者（0 < 1）。
5. **歧义掩码**：对每项（父/母各单倍型等位、每名子女两个传递位），若全部最优解
   在该项取值相同则 `fixed=true`，否则 `false`。
6. **确定性**：所有裁决均为全序域上的显式取最小，无随机、无时间戳、无哈希序依赖，
   相同输入重复调用逐字节一致。

## 算法与复杂度（不枚举解集合）

- 状态：每名子女 2 个传递位，联合状态至多 `4^6 = 4096`。
- 前向/后向 DP：相邻状态间代价为汉明距离，min-plus 转移用超立方体距离变换，
  每标记 `O(2K·4^K)`；整体 `O(M·2K·4^K)`，1500×6 规模实测 < 1 秒。
- 规范解：沿标记贪心，每步用前向/后向表判定“已选前缀是否仍可扩展为最优解”。
- 歧义掩码：由前向/后向表推出“位于某条最优路径上的状态集合”及其兼容的父母相位，
  逐项检查取值是否唯一。**全程不枚举最优解集合**（其数量可指数级大）。

## 本地开发与测试

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/uvicorn app.main:app --reload --port 8000
.venv/bin/pytest            # 54 项测试
```

测试包含：手工推算用例、与**枚举全部最优解**的独立 oracle 在小规模随机/模拟家系上的
全量交叉验证（最优值、规范解、掩码逐项一致）、与朴素 O(16^K) 转移 DP 的中规模对照、
1500×6 满负荷冒烟、422/409 错误形状、重复调用字节一致性。

## 目录结构

```
app/
  main.py      # FastAPI 应用：路由、健康检查、异常处理（422/409）
  schemas.py   # 请求/响应模型（Pydantic 校验）
  solver.py    # 核心求解器：DP、规范解、歧义掩码
  errors.py    # 领域异常
tests/         # 单元测试 + 独立 oracle 交叉验证 + API 测试
examples/      # 示例请求体
Dockerfile / docker-compose.yml
```
