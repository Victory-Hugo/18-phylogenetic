# `baseml_llt` 多线程使用说明

`baseml_llt` 是当前仓库里为 `baseml` 维护的 OpenMP 多线程版本。
它保留原来的 `baseml` 用法，但多了可选的线程并行路径。

## 1. 编译

在源码目录执行：

```bash
cd /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/src
make baseml_llt
```

编译完成后会生成：

- `src/baseml_llt`

仓库里的可执行文件建议放在：

- `bin/baseml_llt`

## 2. 运行

最基本的运行方式：

```bash
/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

## 3. 如何开启多线程

`baseml_llt` 的线程数由 OpenMP 环境变量控制。

例如使用 4 线程：

```bash
OMP_NUM_THREADS=4 /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

## 4. 如何开启第三阶段并行内核

当前更快的并行路径默认是关闭的。
如果要启用第三阶段内核，请加上 `BASEML_PAR3=1`。

推荐用法：

```bash
BASEML_PAR3=1 OMP_NUM_THREADS=4 /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

可选参数：

- `BASEML_PAR3_BLOCK=512`
  - 控制并行 block 大小
  - 一般保持默认即可

例如：

```bash
BASEML_PAR3=1 BASEML_PAR3_BLOCK=512 OMP_NUM_THREADS=4 /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

## 5. 推荐设置

在当前机器和现有实现下，通常建议先试：

```bash
BASEML_PAR3=1 OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

说明：

- `OMP_NUM_THREADS=4`
  - 目前通常比 `8` 线程更稳
- `OMP_PROC_BIND=close`
  - 减少线程漂移
- `OMP_PLACES=cores`
  - 尽量按物理核心绑定

## 6. 关闭并行内核

如果你只想保留 OpenMP 线程框架，但不启用第三阶段内核：

```bash
OMP_NUM_THREADS=4 /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

如果你想显式关闭第三阶段：

```bash
BASEML_PAR3=0 OMP_NUM_THREADS=4 /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

## 7. 注意事项

- `baseml_llt` 是 `baseml` 的多线程变体，不修改 `baseml.ctl` 格式
- 线程数不是在 ctl 文件里设置，而是通过环境变量设置
- 小数据集上，多线程不一定更快
- 如果多个任务同时跑，请给每个任务单独工作目录，避免输出文件互相覆盖
