"""Rebuild the molecule in PySCF and write CAS-space natural orbitals as Molden.

The log stores natural orbitals in the meta-Lowdin AO basis (``C_ml``). With the
geometry (``.xyz``), the basis-set name and the charge/spin (from ``RHF.log`` /
``ROHF.log``) the molecule is rebuilt in PySCF, the meta-Lowdin transform ``X`` is
recomputed, and the raw AO coefficients are recovered as ``C_ao = X @ C_ml``.
These are written as a standard Molden ``[MO]`` section (spherical harmonics,
``Spin= Alpha`` with fractional ``Occup=`` values = NOONs, which EigenVista treats
as spatial orbitals).

Precision note: the log prints coefficients with 5 decimals, so reconstructed
orbitals carry ~1e-5 rounding noise (column norms deviate from 1 by ~1e-4). This
is irrelevant for visualisation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .parser import find_context, parse_casci_log


def _output_name(log_path: Path, root: Path, strip: tuple[str, ...]) -> str:
    """Derive the flat output filename from the log path relative to ``root``.

    Directory components (the ``CASCI.log`` filename itself excluded) are joined
    with ``_``; any component listed in ``strip`` is dropped. Falls back to the
    root directory name when the log sits directly in ``root``.
    """
    try:
        rel = log_path.relative_to(root)
    except ValueError:
        rel = log_path
    dirs = [p for p in rel.parts[:-1] if p not in ("/", ".", "..")]
    if not dirs:
        dirs = [root.name]
    parts = [p for p in dirs if p not in strip]
    if not parts:
        parts = [root.name]
    return "_".join(parts) + ".molden"


def convert_one(
    casci_log: Path,
    out_dir: Path,
    root: Path | None = None,
    strip: tuple[str, ...] = (),
    dry: bool = False,
) -> dict:
    """Convert a single ``CASCI.log`` to a Molden file.

    ``root`` anchors the output naming (defaults to the log's parent directory).
    ``strip`` removes path components from the generated name. Returns a report
    dict with the conversion metrics.
    """
    from pyscf import gto, lo
    from pyscf.tools import molden

    casci_log = Path(casci_log)
    out_dir = Path(out_dir)
    root = Path(root) if root is not None else casci_log.parent

    data = parse_casci_log(casci_log)
    xyz, charge, spin = find_context(casci_log)

    mol = gto.M(
        atom=str(xyz),
        basis=data["basis"],
        charge=charge,
        spin=spin,
        unit="Angstrom",
        cart=False,
        verbose=0,
    )
    nao = mol.nao_nr()
    if nao != len(data["labels"]):
        raise ValueError(
            f"{casci_log}: rebuilt mol has nao={nao} but log table has "
            f"{len(data['labels'])} rows (basis mismatch?)"
        )

    c_ml = data["c_ml"]
    occ = data["occ"]

    # Sanity: meta-Lowdin functions are orthonormal, so printed columns are ~unit norm.
    gram = c_ml.T @ c_ml
    orth_err = float(np.abs(gram - np.eye(data["ncas"])).max())
    trace_err = float(abs(occ.sum() - sum(data["nelec"])))

    # Invert the meta-Lowdin representation: C_ao = X @ C_ml
    s_ao = mol.intor("int1e_ovlp")
    x_ml = lo.orth_ao(mol, "meta_lowdin", s=s_ao)
    c_ao = x_ml @ c_ml
    s_orth = float(np.abs(c_ao.T @ s_ao @ c_ao - np.eye(data["ncas"])).max())

    name = _output_name(casci_log, root, tuple(strip))
    out_path = out_dir / name
    if dry:
        print(f"[dry] {casci_log} -> {out_path}")
        return {"source": str(casci_log), "out": str(out_path), "nao": nao}

    out_dir.mkdir(parents=True, exist_ok=True)
    molden.from_mo(
        mol,
        str(out_path),
        c_ao,
        spin="Alpha",
        ene=np.zeros(data["ncas"]),
        occ=occ,
        ignore_h=False,
    )
    # Prepend provenance comments (EigenVista skips '#' lines).
    lines = out_path.read_text().splitlines()
    prov = [
        f"# converted from {casci_log} by casci2molden",
        f"# CAS({data['nelec'][0]}e+{data['nelec'][1]}e,{data['ncas']}o) "
        f"ncore={data['ncore']} nvir={data['nvir']} basis={data['basis']}",
        "# orbitals: CAS-space natural orbitals (meta-Lowdin coefficients inverted to AO basis)",
        "# occupations: NOONs; Spin=Alpha set treated as spatial orbitals by EigenVista",
    ]
    out_path.write_text("\n".join(lines[:2] + prov + lines[2:]) + "\n")

    return {
        "source": str(casci_log),
        "out": str(out_path),
        "basis": data["basis"],
        "nao": nao,
        "ncas": data["ncas"],
        "nelec": list(data["nelec"]),
        "max_col_orth_err_ml": orth_err,
        "occ_sum_err": trace_err,
        "max_s_orth_err_ao": s_orth,
    }


def convert_tree(
    root: Path,
    out_dir: Path,
    pattern: str = "CASCI.log",
    strip: tuple[str, ...] = (),
    dry: bool = False,
) -> tuple[list[dict], list[tuple[str, str]]]:
    """Recursively convert every ``pattern`` file under ``root``.

    Returns ``(reports, failures)`` where ``failures`` is a list of
    ``(source, error)`` pairs. Individual failures never abort the run.
    """
    root = Path(root)
    logs = sorted(root.rglob(pattern))
    reports: list[dict] = []
    failures: list[tuple[str, str]] = []
    for log in logs:
        try:
            reports.append(convert_one(log, out_dir, root=root, strip=strip, dry=dry))
        except Exception as exc:  # noqa: BLE001 - report and continue
            failures.append((str(log), f"{type(exc).__name__}: {exc}"))
    return reports, failures
