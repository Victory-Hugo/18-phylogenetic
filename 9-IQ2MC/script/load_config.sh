#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "[ERROR] Source this script instead: source script/load_config.sh <config_path>" >&2
    exit 1
fi

_config_path="${1:-}"
if [[ -z "$_config_path" ]]; then
    echo "[ERROR] Missing config path. Usage: source script/load_config.sh <config_path>" >&2
    return 1
fi

if [[ ! -f "$_config_path" ]]; then
    echo "[ERROR] Config file not found: $_config_path" >&2
    return 1
fi

mapfile -t _config_exports < <(
    awk '
    function trim(value) {
        sub(/^[[:space:]]+/, "", value)
        sub(/[[:space:]]+$/, "", value)
        return value
    }

    function unquote(value, first_char, last_char) {
        first_char = substr(value, 1, 1)
        last_char = substr(value, length(value), 1)
        if ((first_char == "\"" && last_char == "\"") || (first_char == "\047" && last_char == "\047")) {
            return substr(value, 2, length(value) - 2)
        }
        return value
    }

    function shq(value) {
        gsub(/\047/, "\047\"\047\"\047", value)
        return "\047" value "\047"
    }

    /^[[:space:]]*($|#)/ {
        next
    }

    match($0, /^([A-Za-z0-9_-]+):[[:space:]]*$/, groups) {
        section = groups[1]
        next
    }

    match($0, /^  ([A-Za-z0-9_-]+):[[:space:]]*(.*)[[:space:]]*$/, groups) {
        if (section == "") {
            printf "[ERROR] Invalid YAML layout near line %d: %s\n", NR, $0 > "/dev/stderr"
            exit 1
        }

        key = groups[1]
        value = trim(groups[2])

        if (value ~ /^[^"\047]*[[:space:]]+#/) {
            sub(/[[:space:]]+#.*/, "", value)
            value = trim(value)
        }

        value = unquote(value)
        var_name = toupper(section "_" key)
        gsub(/-/, "_", var_name)
        printf "export CFG_%s=%s\n", var_name, shq(value)
        next
    }

    {
        printf "[ERROR] Unsupported YAML syntax in %s at line %d: %s\n", FILENAME, NR, $0 > "/dev/stderr"
        exit 1
    }
    ' "$_config_path"
)

for _config_export in "${_config_exports[@]}"; do
    eval "$_config_export"
done

unset _config_path
unset _config_export
unset _config_exports
