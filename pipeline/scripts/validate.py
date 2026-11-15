#!/usr/bin/env python
"""
Walk the campaigns/ tree and check MANIFEST integrity.

Reports:
    - campaigns with missing MANIFEST.json
    - campaigns with stages marked complete but missing files
    - archive files that don't match their .index
    - orphaned files not referenced by any manifest

Not yet implemented.
"""

if __name__ == "__main__":
    raise NotImplementedError("TODO: implement validate.py")
