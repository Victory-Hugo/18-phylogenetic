#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

python3 "$SCRIPT_DIR/merge_using_fixed_backbone_scale.py" \
  --split-output-dir "$PROJECT_ROOT/output/split" \
  --backbone-tree "$PROJECT_ROOT/output/paml/backbone_analysis/backbone_ultrametric_tree.nwk" \
  --analysis-tree-dir "$PROJECT_ROOT/output/paml/analysis_trees" \
  --output-dir "$SCRIPT_DIR/output_fixed_backbone_merge"
