# IQ2MC — IQ-TREE 3 + MCMCTree 贝叶斯系统发育定年 pipeline

结合 **IQ-TREE 3** 的 `--dating mcmctree` 模式与 **MCMCTree**，从多序列比对和有根校准树出发，自动化推断带时间信息的系统发育树。

## 目录

```
9-IQ2MC/
├── pipe/run_pipeline.sh      # 主控脚本
├── python/                    # Python 工具（双模块）
│   ├── validate_inputs.py     #   输入校验
│   ├── build_iqtree_cmd.py    #   组装 IQ-TREE 命令
│   ├── prepare_mcmctree_ctl.py#   生成 MCMCTree 控制文件
│   └── collect_summary.py     #   运行汇总
├── script/                    # Bash 辅助（环境检查、配置加载）
├── conf/Config.yaml           # 配置文件
├── bin/iqtree3                # IQ-TREE 3 可执行文件
├── Linux/                     # PAML/MCMCTree 软件包
└── data/                      # 示例输入数据
```

## 三步流程

| Step | 说明 |
|------|------|
| 1 | IQ-TREE 3 最大似然树推断（可选，`run_step1`） |
| 2 | IQ-TREE 3 的 `--dating mcmctree` 模式，生成 Hessian 矩阵、有根树、dummy 比对和 `.ctl` 控制文件 |
| 3 | MCMCTree 贝叶斯 MCMC 定年，支持断点续跑 |

## 用法

1. 编辑 `conf/Config.yaml`，填入输入文件路径和参数
2. 运行：`bash pipe/run_pipeline.sh conf/Config.yaml`
3. 结果在 `output/` 目录下，`logs/summary.md` 为运行摘要
