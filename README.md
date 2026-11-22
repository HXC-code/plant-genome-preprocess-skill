# plant-wgdi-preprocess（Agent Skill）

植物比较基因组 WGDI 流程预处理的 Agent 调用入口，按实验室《整体流程》覆盖
**第 1–4 步**：数据收集整理 → 预处理四件套 → DIAMOND 同源比对 → WGDI 共线性。
Ks 分析为第 4 步内可选分支。第 5–8 步（系统发育树、文献核对、多重共线性、圈图）
不在范围内。

| 项 | 值 |
| --- | --- |
| Skill 名 | `plant-wgdi-preprocess` |
| 依赖程序仓 | [plant-genome-preprocess](https://github.com/HXC-code/plant-genome-preprocess)（当前 v3.7.0） |
| 运行副本 | `pipeline/`（与程序仓镜像同步，含源码、测试、0–6 脚本、资源表） |
| 范围文档 | `pipeline/code/docs/WORKFLOW.md` |
| 发布门禁 | `pipeline/VALIDATION.md` |

## 安装

把**本仓库根目录**（包含 `SKILL.md` 和 `pipeline/`）复制到 Agent 的 skills 目录：

```powershell
# 按所用 Agent 的约定选择目标目录
xcopy /E /I . "$env:USERPROFILE\.codex\skills\plant-wgdi-preprocess"
# 或
xcopy /E /I . "$env:USERPROFILE\.workbuddy\skills\plant-wgdi-preprocess"
```

`pipeline/` 是可运行副本，不需要再从程序仓单独下载。

## 使用

```text
使用 plant-wgdi-preprocess 把 D:\raw 下的这批植物基因组整理成 WGDI 四件套
```

Agent 会按 `SKILL.md` 执行：注册/校验本地数据 → `run` → `validate` → `check-names`
→（可选）DIAMOND 与 WGDI 共线性。分析开始前必须先声明第 5–8 步不支持。

手动运行（不经 Agent）也可以直接使用 `pipeline/`：

```powershell
cd pipeline\code
python -m pip install -e ".[gbff]"
python -m pgci_preprocess init-project --project-root D:\plant_project
python -m pgci_preprocess run      --project-root D:\plant_project
python -m pgci_preprocess validate --project-root D:\plant_project --refresh-links
```

详细命令见 `pipeline/QUICKSTART.md` 与 `pipeline/README.md`。

## 强制执行规则

Skill 依赖程序引擎落实、而非口头约定：

1. 代表转录本：CDS 最长优先，等长按稳定转录本 ID 字典序；不以蛋白长度或出现顺序为依据。
2. 默认只保留染色体挂载序列；无染色体证据 fail-closed，碎片化组装用
   `--sequence-mode adaptive`。
3. 剔除必须留痕（`removed_list.tsv`）。
4. Ks 默认关闭，仅明确要求时 `WGDI_KS=true`。
5. 审计链保留（输入 SHA-256、manifest、绑定 QC 的裁决记录）。

## 维护

`pipeline/` 由程序仓的 `code/tools/sync_skill_package.py` 生成，**不要手工修改**；
改动先在程序仓落地，再同步并在本仓库提交。CI 会拉取程序仓比对镜像，不一致即失败。

```powershell
cd pipeline\code   # 在程序仓中执行
python tools/sync_skill_package.py --pipeline-target ../../plant-genome-preprocess-skill/pipeline --check
```

## 许可

MIT License（见 `LICENSE`）。第三方工具与脚本按其上游许可使用，
详见 `pipeline/THIRD_PARTY_NOTICES.md`。
