"""End-to-end round-trip: generate a real PySCF CASCI.log, convert it, verify.

This is the strongest test of the meta-Lowdin inversion: it produces a genuine
CASCI.log with PySCF, converts it, and checks that the recovered AO coefficients
are orthonormal under the AO overlap matrix (which they must be, since the natural
orbitals are orthonormal).
"""
from __future__ import annotations

import io
from pathlib import Path

import pytest

pytest.importorskip("pyscf")


def _generate(tmp_path):
    from pyscf import gto, mcscf, scf

    mol = gto.M(
        atom="O 0 0 0; H 0 0.96 0; H 0 -0.48 0.83",
        basis="6-31g*",
        spin=0,
        unit="Angstrom",
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.kernel()
    mc = mcscf.CASCI(mf, 4, 2)  # CAS(2e, 4o)
    mc.verbose = 5
    buf = io.StringIO()
    # PySCF logs to mc.stdout (bound at construction), not sys.stdout.
    mc.stdout = buf
    mc.kernel()
    mc.analyze()
    log = buf.getvalue()

    calc = tmp_path / "calc"
    calc.mkdir()
    logp = calc / "CASCI.log"
    logp.write_text(log)
    # Write the geometry in Angstrom (mol.tofile has no unit control).
    coords = mol.atom_coords(unit="Angstrom")
    symbols = [mol.atom_symbol(i) for i in range(mol.natm)]
    xyz = f"{mol.natm}\nH2O test geometry\n" + "".join(
        f"{s} {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}\n" for s, c in zip(symbols, coords)
    )
    (calc / "h2o.xyz").write_text(xyz)
    (calc / "RHF.log").write_text(
        "[INPUT] charge = 0\n[INPUT] spin (= nelec alpha-beta = 2S) = 0\n"
    )
    return logp, calc, mol.nao_nr()


def test_roundtrip_orthonormality(tmp_path):
    from casci2molden import convert_one

    logp, calc, nao = _generate(tmp_path)
    text = logp.read_text()
    assert "meta-Lowdin" in text, "PySCF did not print the meta-Lowdin table"

    out = tmp_path / "out"
    rep = convert_one(logp, out, root=calc)

    assert rep["nao"] == nao
    assert rep["ncas"] == 4
    # Recovered AO coefficients must be orthonormal under S (natural orbitals are).
    assert rep["max_s_orth_err_ao"] < 1e-3
    # NOONs sum to the active-electron count.
    assert rep["occ_sum_err"] < 1e-6

    molden = out / "calc.molden"  # log sits directly in root -> root name
    assert molden.exists()
    body = molden.read_text()
    assert body.startswith("[Molden Format]")
    assert "Spin=Alpha" in body
    assert "# converted from" in body


def test_convert_tree_and_report(tmp_path):
    from casci2molden import convert_tree

    logp, calc, _ = _generate(tmp_path)
    out = tmp_path / "out"
    reports, failures = convert_tree(calc, out)
    assert failures == []
    assert len(reports) == 1
    assert (out / "calc.molden").exists()
