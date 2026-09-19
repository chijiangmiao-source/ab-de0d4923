# 家系相位分析服务（Family Phasing Service）

纯后端 HTTP 服务：对含缺失位点的双亲—多名子女家系，求**全局最少交叉**的孟德尔
一致定相，并给出唯一可复算的规范解与逐字段歧义掩码。无前端、不调用任何在线服务。

## 快速开始

```bash
# 宿主机端口由环境变量 HOST_PORT 配置（默认 8080）
HOST_PORT=9000 docker compose up --build
```

容器内固定监听 8080；compose 通过 `${HOST_PORT:-8080}:8080` 发布到宿主机。
健康检查同时定义在 Dockerfile（`HEALTHCHECK`）与 compose（`healthcheck`）中。

```bash
curl http://localhost:9000/health
# {"status":"ok"}
```

## 请求格式 `POST /phase`

等位基因为 `A/C/G/T`，缺失为 `.`（也接受 `-` 与两元素数组形式如 `["A","."]`）。
标记按数组顺序排列，需 **2–1500** 个；子女 **1–6** 名，各基因型列表长度必须与
标记数一致。

```json
{
  "markers": ["m1", "m2"],
  "father": {"genotypes": ["AC", "TT"]},
  "mother": {"genotypes": ["GG", "AC"]},
  "children": [
    {"name": "k", "genotypes": ["AG", "TC"]}
  ]
}
```

成功（200）：

```json
{
  "min_crossovers": 0,
  "child_crossovers": {"k": 0},
  "solution": [
    {
      "marker": "m1",
      "father": {"haplotype": ["A", "C"], "ambiguous": [false, false]},
      "mother": {"haplotype": ["G", "G"], "ambiguous": [false, false]},
      "children": [
        {"name": "k", "transmitted_haplotype": [0, 1], "ambiguous": [false, false]}
      ]
    }
  ]
}
```

- `min_crossovers`：相邻标记间所有子女父源/母源链切换总数的全局最小值。
- 每个字段的 `ambiguous`（歧义掩码）表示该值是否在**所有最优解中不固定**
  （`true` = 至少存在两个同为最优的解在此不同）。
- `transmitted_haplotype`：`[父源链位, 母源链位]`，0/1 分别指向该标记父亲/母亲
  的第 0/1 条单倍型。
- 响应以固定分隔符、排序键序列化，**相同输入重复调用逐字节一致**。

### 方向锚定与规范解

全局标签翻转（交换某亲本两条同源体编号并同步翻转所有子女的对应传递位）不改变
代价。规范约定：父亲、母亲分别以其**首个真杂合位点**（`A.` 这类半缺失不算）
固定 0 号方向——该位点 0 号单倍型取字典序较小的等位基因；此后字典序按
「标记 → 父亲两条 → 母亲两条 → 各子女（父源位, 母源位）」逐字段裁决。全程
**不枚举完整解集合**：前向/后向代价层做最小和 DP，用超立方体汉明距离变换在
`O(M·K·2^K)`（K ≤ 12）内完成；规范解由 fwd/bwd 证书下逐字段贪心得到；固定性
由能参与某条全局最优路径的局部取值集合判定。若某亲本在全部标记都无真杂合
位点，其物理标签不可识别，对应传递位/方向会在掩码中如实报为歧义，但仍输出
唯一的字典序代表。

### 错误（422）

- 结构非法：`{"error":{"type":"invalid_structure","message":...,"path":"$.father.genotypes[1]"}}`
  ，`path` 可定位到具体字段/下标。
- 结构合法但无孟德尔一致解释：`type=mendel_inconsistent`，并给出**最早无解
  标记**的下标与名称，如 `"marker_index":1,"marker":"m2"`。

## 本地开发与测试

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
pip install pytest httpx
pytest -q
uvicorn app.main:app --host 0.0.0.0 --port 8080
```

`tests/test_phasing.py` 在数百个随机小家系上用**全路径暴力枚举**独立校验：
最小交叉数、字典序规范解、每个字段的固定性掩码、锚定方向，以及孟德尔无解时
的最早标记；`tests/test_api.py` 覆盖 HTTP 200/422、错误定位与响应逐字节确定性。

## 项目结构

```
app/model.py     输入解析与可定位的结构校验（422）
app/phasing.py   位点可行性、汉明距离最小和 DP、锚定、规范解与歧义掩码
app/main.py      FastAPI：/health 与 /phase，确定性 JSON 序列化
Dockerfile       python:3.11-slim，含 HEALTHCHECK
docker-compose.yml  HOST_PORT 环境变量配置宿主机端口
```

