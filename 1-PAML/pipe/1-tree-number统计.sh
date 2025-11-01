#!/usr/bin/env bash
set -euo pipefail

# 1) 指定树文件路径
tree_file='/mnt/c/Users/Administrator/Desktop/BIG.treefile'

if [[ ! -f "$tree_file" ]]; then
  echo "找不到文件：$tree_file" >&2
  exit 1
fi

# 定义“有效树行”的判定：同一行内同时出现 "("、")"、";"
tree_line_regex='\(.*\).*;.*'

# 2) 从文件中裁剪出：自第一棵有效树行起直到文件结束的内容（保证树在第二行）
trimmed_content="$(
  awk -v RS='\n' '
    BEGIN{flag=0}
    {
      if (!flag) {
        if ($0 ~ /\(.*\).*;.*$/) { flag=1; print $0 }
      } else {
        print $0
      }
    }
  ' "$tree_file"
)"

# 若整文件没有找到有效树行
if [[ -z "${trimmed_content}" ]]; then
  num_trees=0
  num_samples=0
  tmp_file="$(mktemp)"
  printf " %d %d;\n" "$num_samples" "$num_trees" > "$tmp_file"
  mv "$tmp_file" "$tree_file"
  echo "已更新 ${tree_file}：未发现有效树行。写入头部 '0 0;'"
  exit 0
fi

# 3) 基于裁剪后的内容统计树的数量
num_trees="$(printf "%s\n" "$trimmed_content" | grep -E -c "$tree_line_regex" || true)"

# 4) 从第一棵树行估算样本数
first_line="$(printf "%s\n" "$trimmed_content" | grep -m1 -E "$tree_line_regex" || true)"
if [[ -n "$first_line" ]]; then
  commas="$(printf "%s" "$first_line" | tr -cd ',' | wc -c | awk '{print $1}')"
  num_samples=$((commas + 1))
else
  num_samples=0
fi

# 5) 写回（第一行数字带分号，第二行就是第一棵树）
tmp_file="$(mktemp)"
{
  printf " %d %d;\n" "$num_samples" "$num_trees"
  printf "%s\n" "$trimmed_content"
} > "$tmp_file"

mv "$tmp_file" "$tree_file"

echo "已成功更新 ${tree_file}：找到 ${num_trees} 棵树，样本数为 ${num_samples}（第一行数字已带分号）"
