#!/usr/bin/env python3

import argparse
import json


def lookup_key(data, dotted_key):
    current = data
    for part in dotted_key.split("."):
        if not isinstance(current, dict) or part not in current:
            raise KeyError(f"Missing config key: {dotted_key}")
        current = current[part]
    return current


def run(config, key):
    with open(config, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    value = lookup_key(data, key)
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    if value is None:
        return "null"
    return str(value)


def build_parser():
    parser = argparse.ArgumentParser(description="Read a value from Config.json.")
    parser.add_argument("--config", required=True, help="Path to Config.json")
    parser.add_argument("--key", required=True, help="Dotted key such as paths.input_tree")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    print(run(args.config, args.key))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
