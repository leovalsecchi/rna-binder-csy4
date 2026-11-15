"""ck.site — load site.yaml and expand ${key} placeholders in loaded configs.

This is the portability layer: cluster/*.yaml and models/*.yaml use ${rb_root},
${shared}, ${user}; a new user sets those once in site.yaml (from site.yaml.example)
and the whole kit resolves to their paths. No absolute paths are baked into any
committed file.
"""
from __future__ import annotations
import os
import re

try:
    import yaml
except ImportError:                     # pragma: no cover
    yaml = None

_VAR = re.compile(r"\$\{(\w+)\}")


def load_site(kit_root: str) -> dict:
    """Read site.yaml (env CK_SITE overrides). Clear error if it's missing."""
    path = os.environ.get("CK_SITE") or os.path.join(kit_root, "site.yaml")
    if not os.path.exists(path):
        raise SystemExit(
            f"no site.yaml at {path}\n"
            f"  -> copy site.yaml.example to site.yaml and edit it (user/rb_root/shared)")
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def resolve(obj, site: dict):
    """Recursively expand ${key} in every string, using site (supports nesting)."""
    if isinstance(obj, str):
        def sub(m):
            k = m.group(1)
            if k not in site:
                raise KeyError(f"site.yaml has no '{k}' (referenced in a config)")
            return str(site[k])
        prev = None
        while prev != obj:               # allow ${a} whose value contains ${b}
            prev, obj = obj, _VAR.sub(sub, obj)
        return obj
    if isinstance(obj, list):
        return [resolve(x, site) for x in obj]
    if isinstance(obj, dict):
        return {k: resolve(v, site) for k, v in obj.items()}
    return obj
