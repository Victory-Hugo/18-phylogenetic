#!/usr/bin/env bash
# 读取两层 YAML 配置中的标量值。

set -euo pipefail

yaml_get() {
    local config_path="$1"
    local section="$2"
    local key="$3"

    awk -v wanted_section="$section" -v wanted_key="$key" '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        /^[^[:space:]][^:]*:[[:space:]]*$/ {
            current_section = $0
            sub(/:[[:space:]]*$/, "", current_section)
            next
        }
        current_section == wanted_section {
            line = $0
            sub(/^[[:space:]]+/, "", line)
            split(line, parts, ":")
            if (parts[1] == wanted_key) {
                sub(/^[^:]+:[[:space:]]*/, "", line)
                gsub(/^"|"$/, "", line)
                print line
                exit
            }
        }
    ' "$config_path"
}
