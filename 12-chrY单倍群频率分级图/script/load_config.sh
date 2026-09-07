#!/usr/bin/env bash
# 将两层结构的 YAML 配置解析为大写 shell 变量赋值语句，供 pipe 脚本 eval 使用
parse_yaml() {
    awk '
        /^[[:space:]]*#/ { next }
        /^[[:space:]]+[A-Za-z_]+:[[:space:]]*[^[:space:]]/ {
            line = $0
            sub(/[[:space:]]*#.*$/, "", line)
            idx = index(line, ":")
            key = substr(line, 1, idx - 1)
            val = substr(line, idx + 1)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", val)
            gsub(/"/, "", val)
            print toupper(key) "=\"" val "\""
        }
    ' "$1"
}
