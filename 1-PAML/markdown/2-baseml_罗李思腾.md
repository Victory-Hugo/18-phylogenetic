# 说明
一个心血来潮的研究生觉得`baseml`的速度太慢，所以添加了`OpenMP`来进行并行化处理。
# 编译
进入工作目录：
```sh
cd /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/src
```
开始编译：
```sh
# 注意，在/src/ 内需要treesub_LLST.c 文件，在编译时默认需要
gcc -O3 -std=c11 -fopenmp -o baseml_LLST baseml_LLST.c tools.c -lm
mv baseml_LLST ../bin/
```

# 运行
运行的代码与之前相同，只需确保使用新编译的可执行文件，例如`pipe/2_运行.sh`。

>[!question] 性能问题：
>罗李思腾的程序的确能够使用超多的线程，但是我却觉得速度没有明显提升，不知道为什么……