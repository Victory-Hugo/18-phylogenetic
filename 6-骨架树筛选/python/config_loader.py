#!/usr/bin/env python3

import argparse
import ast
import json
from pathlib import Path


def parse_scalar(token):
    if token in {"null", "Null", "NULL", "~"}:
        return None
    if token in {"true", "True", "TRUE"}:
        return True
    if token in {"false", "False", "FALSE"}:
        return False
    if token.startswith('"') and token.endswith('"'):
        return json.loads(token)
    if token.startswith("'") and token.endswith("'"):
        return token[1:-1].replace("''", "'")
    try:
        return ast.literal_eval(token)
    except (SyntaxError, ValueError):
        return token


def preprocess_yaml(text):
    items = []
    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        if "\t" in raw_line[:indent]:
            raise ValueError(f"Tabs are not supported in YAML indentation (line {lineno}).")
        items.append((indent, raw_line[indent:], lineno))
    return items


def split_mapping_entry(content, lineno):
    if ":" not in content:
        raise ValueError(f"Invalid YAML mapping entry on line {lineno}: {content}")
    key, value = content.split(":", 1)
    key = key.strip()
    if not key:
        raise ValueError(f"Missing YAML key on line {lineno}.")
    return key, value.strip()


def parse_yaml_block(items, index, indent):
    if index >= len(items):
        return None, index
    _, content, _ = items[index]
    if content.startswith("- "):
        return parse_yaml_sequence(items, index, indent)
    return parse_yaml_mapping(items, index, indent)


def parse_yaml_mapping(items, index, indent):
    data = {}
    while index < len(items):
        current_indent, content, lineno = items[index]
        if current_indent < indent:
            break
        if current_indent != indent:
            raise ValueError(f"Unexpected indentation on line {lineno}.")
        if content.startswith("- "):
            raise ValueError(f"Unexpected list item inside mapping on line {lineno}.")
        key, value_part = split_mapping_entry(content, lineno)
        index += 1
        if value_part:
            data[key] = parse_scalar(value_part)
            continue
        if index < len(items) and items[index][0] > current_indent:
            value, index = parse_yaml_block(items, index, items[index][0])
            data[key] = value
        else:
            data[key] = None
    return data, index


def parse_yaml_sequence(items, index, indent):
    sequence = []
    while index < len(items):
        current_indent, content, lineno = items[index]
        if current_indent < indent:
            break
        if current_indent != indent:
            raise ValueError(f"Unexpected indentation on line {lineno}.")
        if not content.startswith("- "):
            break
        value_part = content[2:].strip()
        index += 1
        if value_part:
            sequence.append(parse_scalar(value_part))
            continue
        if index < len(items) and items[index][0] > current_indent:
            value, index = parse_yaml_block(items, index, items[index][0])
            sequence.append(value)
        else:
            sequence.append(None)
    return sequence, index


def load_yaml(text):
    items = preprocess_yaml(text)
    if not items:
        return {}
    data, index = parse_yaml_block(items, 0, items[0][0])
    if index != len(items):
        _, _, lineno = items[index]
        raise ValueError(f"Failed to parse YAML near line {lineno}.")
    return data


def load_config(path):
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        if config_path.suffix.lower() == ".json":
            return json.load(handle)
        return load_yaml(handle.read())


def lookup_key(data, dotted_key):
    current = data
    for part in dotted_key.split("."):
        if not isinstance(current, dict) or part not in current:
            raise KeyError(f"Missing config key: {dotted_key}")
        current = current[part]
    return current


def run(config, key, default=None, use_default_on_missing=False):
    data = load_config(config)
    try:
        value = lookup_key(data, key)
    except KeyError:
        if not use_default_on_missing:
            raise
        value = default
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    if value is None:
        return "null"
    return str(value)


def build_parser():
    parser = argparse.ArgumentParser(description="Read a value from a JSON or YAML config file.")
    parser.add_argument("--config", required=True, help="Path to a config file (.json/.yaml/.yml)")
    parser.add_argument("--key", required=True, help="Dotted key such as paths.meta_tsv")
    parser.add_argument("--default", default=None, help="Fallback value when the key is missing.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    print(run(args.config, args.key, default=args.default, use_default_on_missing=args.default is not None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
