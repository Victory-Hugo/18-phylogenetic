# 环境

请使用名为`IQ2MC`的conda环境进行分析。


软件安装方式：

若需计算群体遗传相关指标，**优先使用成熟的生物信息学工具包或者软件**，如`scikit-allel`、`pysam`等；若不存在相关工具，可以自己编写代码实现相关功能；但是**自己编写代码不作为首选方案**。

1. 优先使用`/home/luolintao/miniconda3/condabin/conda`安装。新创建一个conda环境，命名为`IQ2MC`。
2. 若没有，则通过`apt`安装。

我的sudo密码是`luolintao`。

# 编写代码的要求

1. 你需要编写的是清晰的生物信息学**pipeline代码**，注释要清晰，变量命名要有意义（注释要中文）。
2. 所有的python代码均需使用**双模块**。
3. 所有的文件均需放在`/mnt/c/Users/Administrator/Desktop/9-IQ2MC/`目录下，不要去删除、修改其他目录下的文件。
4. 严禁在`//mnt/c/Users/Administrator/Desktop/9-IQ2MC/`目录之外生成任何中间文件、测试文件、缓存文件、日志文件或结果文件；如需临时产物，也必须写入该项目目录内部并在任务结束后清理。

