"""Parsing of PySCF CASCI.log files and their sibling context files.

A CASCI.log produced by PySCF (tested with 2.14.0) contains, per active-space
calculation:

* ``Basis: 6-31g*``
* ``CAS (n_a e+n_b e, ncas o), ncore = ..., nvir = ...``
* ``Natural occ [...]``            -- natural orbital occupation numbers (NOONs)
* ``Natural orbital (expansion on meta-Lowdin AOs) in CAS space``
  a table of CAS-space natural-orbital coefficients printed in the
  *meta-Lowdin orthonormalised AO basis* (rows labelled ``<atom> <elem> <ao>``),
  laid out in horizontal column blocks of five.

PySCF prints the natural orbitals as ``C_ml = X^T S C_ao`` where ``S`` is the AO
overlap and ``X = lo.orth_ao(mol, 'meta_lowdin', s=S)`` satisfies ``X^T S X = I``.
Because ``X`` is square, ``X^T S = X^{-1}`` and the raw contracted-Gaussian AO
coefficients are recovered exactly by ``C_ao = X @ C_ml`` (see :mod:`.converter`).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

RE_BASIS = re.compile(r"^\s*Basis:\s*(\S+)", re.MULTILINE)
RE_CAS = re.compile(
    r"^\s*CAS \((\d+)e\+(\d+)e,\s*(\d+)o\),\s*ncore\s*=\s*(\d+),\s*nvir\s*=\s*(\d+)",
    re.MULTILINE,
)
RE_NAT_OCC = re.compile(r"Natural occ \[([^\]]+)\]", re.DOTALL)
RE_NO_HEADER = re.compile(r"Natural orbital \(expansion on meta-Lowdin AOs\) in CAS space")
RE_ROW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+(\S+)((?:\s+[-+0-9.eE-]+)*)$")
RE_INPUT = re.compile(r"^\[INPUT\] (charge|spin).*=\s*(-?\d+)\s*$", re.MULTILINE)

DEFAULT_BASIS = "6-31g*"


def parse_casci_log(path: Path) -> dict:
    """Parse a CASCI.log into basis, CAS size, occupations and meta-Lowdin coefficients.

    Returns a dict with keys: ``basis``, ``nelec`` (a, b), ``ncas``, ``ncore``,
    ``nvir``, ``occ`` (np.ndarray of NOONs), ``labels`` (list of (idx, elem, ao))
    and ``c_ml`` (nao x ncas array of meta-Lowdin coefficients).
    """
    text = Path(path).read_text()

    m = RE_BASIS.search(text)
    basis = m.group(1) if m else DEFAULT_BASIS

    m = RE_CAS.search(text)
    if not m:
        raise ValueError(f"{path}: CAS line not found")
    nelec_a, nelec_b, ncas, ncore, nvir = (int(g) for g in m.groups())

    m = RE_NAT_OCC.search(text)
    if not m:
        raise ValueError(f"{path}: Natural occ list not found")
    occ = np.array([float(t) for t in m.group(1).split()])
    if occ.size != ncas:
        raise ValueError(f"{path}: {occ.size} occupations but ncas={ncas}")

    # Locate the meta-Lowdin natural-orbital table and read its column blocks.
    start = RE_NO_HEADER.search(text)
    if not start:
        raise ValueError(f"{path}: meta-Lowdin natural orbital table not found")
    tail = text[start.end():]
    # The table ends at the next section header (e.g. Mulliken spin population).
    end = re.search(r"^\S.*:$", tail, re.MULTILINE)
    body = tail[: end.start()] if end else tail

    blocks: list[list[tuple[tuple[int, str, str], list[float]]]] = []
    cur_rows: list[tuple[tuple[int, str, str], list[float]]] = []
    for line in body.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            if cur_rows:
                blocks.append(cur_rows)
                cur_rows = []
            continue
        m = RE_ROW.match(line)
        if not m:
            continue
        lab = (int(m.group(1)), m.group(2), m.group(3))
        vals = [float(t) for t in m.group(4).split()]
        cur_rows.append((lab, vals))
    if cur_rows:
        blocks.append(cur_rows)
    if not blocks:
        raise ValueError(f"{path}: no rows parsed from natural orbital table")

    labels = [lab for lab, _ in blocks[0]]
    c_ml = np.zeros((len(labels), ncas))
    col = 0
    for block in blocks:
        if col >= ncas:
            break
        blk_labels = [lab for lab, _ in block]
        if blk_labels != labels:
            raise ValueError(f"{path}: inconsistent row labels between column blocks")
        widths = {len(vals) for _, vals in block}
        if len(widths) != 1:
            raise ValueError(f"{path}: ragged column block")
        w = widths.pop()
        if col + w > ncas:
            raise ValueError(f"{path}: table wider than ncas={ncas}")
        for i, (_, vals) in enumerate(block):
            c_ml[i, col : col + w] = vals
        col += w
    if col != ncas:
        raise ValueError(f"{path}: parsed {col} columns, expected {ncas}")

    return {
        "basis": basis,
        "nelec": (nelec_a, nelec_b),
        "ncas": ncas,
        "ncore": ncore,
        "nvir": nvir,
        "occ": occ,
        "labels": labels,
        "c_ml": c_ml,
    }


def find_context(casci_log: Path, max_up: int = 2) -> tuple[Path, int, int]:
    """Return ``(xyz_path, charge, spin)``.

    Searches the CASCI.log's parent directory, then its grandparent (up to
    ``max_up`` levels), for a ``*.xyz`` plus a ``RHF.log``/``ROHF.log`` that
    carries ``[INPUT] charge`` and ``[INPUT] spin`` lines.
    """
    casci_log = Path(casci_log)
    for k in range(max_up):
        parent = casci_log.parent if k == 0 else casci_log.parents[k]
        xyzs = sorted(parent.glob("*.xyz"))
        mf_logs = sorted(parent.glob("RHF.log")) + sorted(parent.glob("ROHF.log"))
        if xyzs and mf_logs:
            vals = dict(RE_INPUT.findall(mf_logs[0].read_text()))
            if "charge" in vals and "spin" in vals:
                return xyzs[0], int(vals["charge"]), int(vals["spin"])
    raise FileNotFoundError(
        f"{casci_log}: could not find .xyz plus RHF/ROHF.log with charge/spin "
        f"within {max_up} level(s)"
    )
