import gzip
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


def fixture_path(name: str) -> Path:
    p = FIXTURES / f"{name}.cif.gz"
    if not p.exists():
        pytest.skip(f"fixture {name} missing; run scripts/fetch_fixtures.py")
    return p


def fixture_text(name: str) -> str:
    return gzip.decompress(fixture_path(name).read_bytes()).decode()


@pytest.fixture
def rfd3_design() -> Path:
    return fixture_path("rfd3_design")


@pytest.fixture
def rfd3_design_seq7() -> Path:
    return fixture_path("rfd3_design_seq7")


@pytest.fixture
def bgbif_design() -> Path:
    return fixture_path("bgbif_design")


@pytest.fixture
def bd1_design() -> Path:
    return fixture_path("bd1_design")


@pytest.fixture
def af3_refold() -> Path:
    return fixture_path("af3_refold")


@pytest.fixture
def esm_refold() -> Path:
    return fixture_path("esm_refold")
