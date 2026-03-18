#!/usr/bin/env bash

if [[ -n "${_PHYLO_CONSOLE_UI_SH_LOADED:-}" ]]; then
    return 0
fi
_PHYLO_CONSOLE_UI_SH_LOADED=1

ui_init() {
    if [[ "${_PHYLO_UI_READY:-0}" -eq 1 ]]; then
        return 0
    fi
    _PHYLO_UI_READY=1

    local stderr_tty=0
    local color_count=0
    local locale_hint="${LC_ALL:-${LC_CTYPE:-${LANG:-}}}"

    if [[ -t 2 ]]; then
        stderr_tty=1
    fi

    UI_USE_COLOR=0
    if [[ -z "${NO_COLOR:-}" && "${TERM:-}" != "dumb" && "$stderr_tty" -eq 1 ]]; then
        if command -v tput >/dev/null 2>&1; then
            color_count=$(tput colors 2>/dev/null || printf '0')
        else
            color_count=8
        fi
        if [[ "$color_count" =~ ^[0-9]+$ ]] && (( color_count >= 8 )); then
            UI_USE_COLOR=1
        fi
    fi

    UI_USE_UNICODE=0
    case "$locale_hint" in
        *UTF-8*|*utf8*|*utf-8*)
            UI_USE_UNICODE=1
            ;;
    esac
    if [[ -n "${PHYLO_UI_ASCII_ONLY:-}" ]]; then
        UI_USE_UNICODE=0
    fi

    if [[ "$UI_USE_COLOR" -eq 1 ]]; then
        UI_RESET=$'\033[0m'
        UI_BOLD=$'\033[1m'
        UI_DIM=$'\033[2m'
        UI_BLUE=$'\033[38;5;39m'
        UI_CYAN=$'\033[38;5;45m'
        UI_GREEN=$'\033[38;5;42m'
        UI_YELLOW=$'\033[38;5;220m'
        UI_RED=$'\033[38;5;196m'
        UI_MAGENTA=$'\033[38;5;141m'
        UI_WHITE=$'\033[38;5;255m'
    else
        UI_RESET=""
        UI_BOLD=""
        UI_DIM=""
        UI_BLUE=""
        UI_CYAN=""
        UI_GREEN=""
        UI_YELLOW=""
        UI_RED=""
        UI_MAGENTA=""
        UI_WHITE=""
    fi

    if [[ "$UI_USE_UNICODE" -eq 1 ]]; then
        UI_ICON_INFO="•"
        UI_ICON_OK="✔"
        UI_ICON_WARN="⚠"
        UI_ICON_ERROR="✖"
        UI_ICON_STAGE="▶"
        UI_RULE_CHAR="─"
        UI_SECTION_EDGE="│"
        UI_SECTION_JOIN="┄"
    else
        UI_ICON_INFO="*"
        UI_ICON_OK="+"
        UI_ICON_WARN="!"
        UI_ICON_ERROR="x"
        UI_ICON_STAGE=">"
        UI_RULE_CHAR="-"
        UI_SECTION_EDGE="|"
        UI_SECTION_JOIN="-"
    fi
}

ui_repeat() {
    local char="$1"
    local count="$2"
    local out=""
    local i=0
    while (( i < count )); do
        out+="$char"
        ((i += 1))
    done
    printf '%s' "$out"
}

ui_print_line() {
    local color="$1"
    local label="$2"
    local message="$3"
    ui_init
    if [[ -n "$color" ]]; then
        printf '%b[%s]%b %s\n' "$color" "$label" "$UI_RESET" "$message" >&2
    else
        printf '[%s] %s\n' "$label" "$message" >&2
    fi
}

ui_info() {
    ui_print_line "$UI_CYAN" "INFO" "$1"
}

ui_ok() {
    ui_print_line "$UI_GREEN" " OK " "$1"
}

ui_warn() {
    ui_print_line "$UI_YELLOW" "WARN" "$1"
}

ui_error() {
    ui_print_line "$UI_RED" "ERR " "$1"
}

ui_rule() {
    ui_init
    local rule
    rule=$(ui_repeat "$UI_RULE_CHAR" 72)
    printf '%b%s%b\n' "$UI_DIM" "$rule" "$UI_RESET" >&2
}

ui_section() {
    ui_init
    local title="$1"
    local subtitle="${2:-}"
    ui_rule
    printf '%b%s %s%b\n' "$UI_BOLD$UI_BLUE" "$UI_SECTION_EDGE" "$title" "$UI_RESET" >&2
    if [[ -n "$subtitle" ]]; then
        printf '%b%s%b %s\n' "$UI_DIM" "$UI_SECTION_EDGE" "$UI_RESET" "$subtitle" >&2
    fi
    ui_rule
}

ui_logo() {
    ui_init
    if [[ "$UI_USE_UNICODE" -eq 1 ]]; then
        printf '%b%s%b\n' "$UI_CYAN$UI_BOLD"    "             ┌─●" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_CYAN$UI_BOLD"    "         ┌───┤" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_BLUE$UI_BOLD"    "     ┌───┤   └─●" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_BLUE$UI_BOLD"    " ┌───┤" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_MAGENTA$UI_BOLD" "─┤   │       ┌─●" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_MAGENTA$UI_BOLD" " │   └───────┤" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_BLUE$UI_BOLD"    "─┤           └─●" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_CYAN$UI_BOLD"    " │" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_CYAN$UI_BOLD"    " └─────────────┬─●" "$UI_RESET" >&2
        printf '%b%s%b\n' "$UI_MAGENTA$UI_BOLD" "               └─●" "$UI_RESET" >&2
    else
        printf '%s\n' "                o" >&2
        printf '%s\n' "            ,---|" >&2
        printf '%s\n' "        ,---|   '-o" >&2
        printf '%s\n' "    ,---|" >&2
        printf '%s\n' "---|   |       ,-o" >&2
        printf '%s\n' "   |   '-------|" >&2
        printf '%s\n' "---|           '-o" >&2
        printf '%s\n' "   |" >&2
        printf '%s\n' "   '-------------,-o" >&2
        printf '%s\n' "                 '-o" >&2
    fi
}

ui_stage_start() {
    local step="$1"
    local title="$2"
    ui_section "$UI_ICON_STAGE $step" "$title"
}

ui_stage_end() {
    local step="$1"
    local title="$2"
    ui_ok "Completed | $step | $title"
}

ui_kv() {
    ui_init
    local key="$1"
    local value="$2"
    printf '%b%s%b %b%-16s%b %s\n' "$UI_DIM" "$UI_SECTION_JOIN" "$UI_RESET" "$UI_BOLD" "$key" "$UI_RESET" "$value" >&2
}

ui_summary_line() {
    ui_init
    local label="$1"
    local value="$2"
    printf '%b%s%b %s %b%s%b\n' "$UI_DIM" "$UI_ICON_INFO" "$UI_RESET" "$label" "$UI_BOLD" "$value" "$UI_RESET" >&2
}
