"""Output-naming logic (pure, no PySCF)."""
from __future__ import annotations

from pathlib import Path

from casci2molden.converter import _output_name


def test_dirs_joined_with_underscore():
    root = Path("/data/Electronic_H")
    log = root / "F43" / "Ni_2" / "singlet" / "8e_5o" / "CASCI.log"
    assert _output_name(log, root, ()) == "F43_Ni_2_singlet_8e_5o.molden"


def test_strip_removes_components():
    root = Path("/data")
    log = root / "Electronic_H" / "F43" / "8e_5o" / "CASCI.log"
    assert _output_name(log, root, ("Electronic_H",)) == "F43_8e_5o.molden"


def test_log_directly_in_root_uses_root_name():
    root = Path("/data/mycalc")
    log = root / "CASCI.log"
    assert _output_name(log, root, ()) == "mycalc.molden"


def test_log_outside_root_uses_its_own_path():
    root = Path("/data/root")
    log = Path("/elsewhere/a/b/CASCI.log")
    # Not under root -> the log's own (absolute) directory path is used.
    assert _output_name(log, root, ()) == "elsewhere_a_b.molden"
