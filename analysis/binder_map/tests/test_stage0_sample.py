import numpy as np
import pandas as pd
import pytest

from stage0_sample import draw_sample


def _nodes(n: int) -> pd.DataFrame:
    return pd.DataFrame({"bb": [f"bb{i:06d}" for i in range(n)]})


def test_draw_is_exact_size_and_without_replacement():
    mask = draw_sample(_nodes(10_000), n=1_000, seed=42)
    assert mask.sum() == 1_000
    assert mask.dtype == bool
    assert len(mask) == 10_000


def test_draw_is_reproducible_for_a_seed():
    a = draw_sample(_nodes(10_000), n=1_000, seed=42)
    b = draw_sample(_nodes(10_000), n=1_000, seed=42)
    assert np.array_equal(a, b)


def test_a_different_seed_gives_a_different_draw():
    a = draw_sample(_nodes(10_000), n=1_000, seed=42)
    b = draw_sample(_nodes(10_000), n=1_000, seed=43)
    assert not np.array_equal(a, b)


def test_draw_larger_than_universe_raises():
    with pytest.raises(ValueError):
        draw_sample(_nodes(10), n=100, seed=42)


def test_draw_of_the_whole_universe_takes_everything():
    mask = draw_sample(_nodes(500), n=500, seed=42)
    assert mask.all()


def test_draw_is_uniform_across_an_arbitrary_label():
    nodes = _nodes(100_000)
    grp = np.arange(len(nodes)) % 4                 # four equal groups
    mask = draw_sample(nodes, n=20_000, seed=42)
    shares = pd.Series(grp[mask]).value_counts() / 20_000
    assert shares.between(0.24, 0.26).all()
