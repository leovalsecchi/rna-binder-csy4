#!/usr/bin/env python
"""
Walk the campaigns/ tree, read MANIFEST.json files, and rebuild _index/.

Updates:
    _index/designs.db        SQLite database (source of truth for queries)
    _index/campaigns.csv     Human-readable summary
    _index/designs.parquet   For pandas/duckdb analysis
    _index/last_updated.txt

Not yet implemented.
"""

if __name__ == "__main__":
    raise NotImplementedError("TODO: implement reindex.py")
