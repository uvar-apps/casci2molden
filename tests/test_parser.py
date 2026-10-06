"""Fast unit tests for the CASCI.log parser (no PySCF run required)."""
from __future__ import annotations

import textwrap

import numpy as np
import pytest

from casci2molden.parser import parse_casci_log


def _write(tmp_path, text):
    p = tmp_path / "CASCI.log"
    p.write_text(textwrap.dedent(text))
    return p


def test_parses_single_block_two_columns(tmp_path):
    log = _write(
        tmp_path,
        """\
        Basis: 6-31g*
          CAS (1e+1e, 2o), ncore = 3, nvir = 9
          Natural occ [1.6 0.4]
          Natural orbital (expansion on meta-Lowdin AOs) in CAS space

        #1          #2
        0 H 1s     0.70711    0.70711
        1 H 1s     0.70711   -0.70711

        Mulliken spin population:
        """,
    )
    d = parse_casci_log(log)
    assert d["basis"] == "6-31g*"
    assert d["nelec"] == (1, 1)
    assert d["ncas"] == 2
    assert d["ncore"] == 3
    assert d["nvir"] == 9
    np.testing.assert_allclose(d["occ"], [1.6, 0.4])
    assert d["labels"] == [(0, "H", "1s"), (1, "H", "1s")]
    assert d["c_ml"].shape == (2, 2)
    np.testing.assert_allclose(d["c_ml"][0], [0.70711, 0.70711])
    np.testing.assert_allclose(d["c_ml"][1], [0.70711, -0.70711])


def test_stitches_horizontal_column_blocks(tmp_path):
    # Three blocks of widths 2 + 1 = ncas 3, same row labels.
    log = _write(
        tmp_path,
        """\
        Basis: sto-3g
          CAS (2e+2e, 3o), ncore = 1, nvir = 1
          Natural occ [1.9 1.0 0.1]
          Natural orbital (expansion on meta-Lowdin AOs) in CAS space

        #1          #2
        0 C 2s     0.10000    0.20000
        1 C 2p     0.30000    0.40000

        #3
        0 C 2s     0.50000
        1 C 2p     0.60000
        """,
    )
    d = parse_casci_log(log)
    assert d["c_ml"].shape == (2, 3)
    np.testing.assert_allclose(d["c_ml"][:, 0], [0.1, 0.3])
    np.testing.assert_allclose(d["c_ml"][:, 1], [0.2, 0.4])
    np.testing.assert_allclose(d["c_ml"][:, 2], [0.5, 0.6])


def test_missing_basis_defaults(tmp_path):
    log = _write(
        tmp_path,
        """\
          CAS (1e+1e, 2o), ncore = 3, nvir = 9
          Natural occ [1.5 0.5]
          Natural orbital (expansion on meta-Lowdin AOs) in CAS space

        0 H 1s     0.70711    0.70711
        1 H 1s     0.70711   -0.70711
        """,
    )
    assert parse_casci_log(log)["basis"] == "6-31g*"


def test_errors(tmp_path):
    with pytest.raises(ValueError, match="CAS line not found"):
        parse_casci_log(_write(tmp_path, "Natural occ [1.0]\n"))
    with pytest.raises(ValueError, match="occupations but ncas"):
        parse_casci_log(
            _write(
                tmp_path,
                "CAS (1e+1e, 2o), ncore = 1, nvir = 1\nNatural occ [1.0]\n",
            )
        )
    with pytest.raises(ValueError, match="meta-Lowdin natural orbital table not found"):
        parse_casci_log(
            _write(
                tmp_path,
                "CAS (1e+1e, 2o), ncore = 1, nvir = 1\nNatural occ [1.0 1.0]\n",
            )
        )
