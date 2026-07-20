"""Architectural test: oncocartograph.scoring must have zero cross-imports.

This package is meant to be extractable to a standalone, independently
publishable PyPI package (see src/oncocartograph/scoring/README.md and
this package's __init__.py docstring). That promise only holds if it
never imports oncocartograph.data_ingestion, .preprocessing,
.integration, or .drug_targets -- this test enforces that mechanically
rather than relying on developers remembering the rule.
"""

from __future__ import annotations

import ast
from pathlib import Path

_FORBIDDEN_MODULE_PREFIXES = (
    "oncocartograph.data_ingestion",
    "oncocartograph.preprocessing",
    "oncocartograph.integration",
    "oncocartograph.drug_targets",
    "oncocartograph.reporting",
    "oncocartograph.validation",
)

_SCORING_SRC_DIR = Path(__file__).resolve().parents[2] / "src" / "oncocartograph" / "scoring"


def _imported_module_names(source_path: Path) -> set[str]:
    tree = ast.parse(source_path.read_text(), filename=str(source_path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def test_scoring_package_has_no_cross_imports_from_other_submodules() -> None:
    """No .py file under src/oncocartograph/scoring may import a forbidden submodule."""
    violations: dict[str, set[str]] = {}
    for path in _SCORING_SRC_DIR.rglob("*.py"):
        imported = _imported_module_names(path)
        forbidden = {
            name
            for name in imported
            if any(
                name == prefix or name.startswith(prefix + ".")
                for prefix in _FORBIDDEN_MODULE_PREFIXES
            )
        }
        if forbidden:
            violations[str(path)] = forbidden

    assert not violations, f"Forbidden cross-imports found: {violations}"
