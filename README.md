# casci2molden

Convert [PySCF](https://pyscf.org) `CASCI.log` files into **Molden** files.

It was mainly intended to convert file to be loaded from the 
[EigenVista](https://github.com/uvar-apps/eigenvista) program but 
any visualiszation tool support the molden format should work.

The active-space (CAS/CASSCF) logs printed by PySCF store their natural orbitals in
the *meta-Löwdin orthonormalised AO basis*. `casci2molden` rebuilds the molecule
from the sibling `.xyz` and basis name, inverts the meta-Löwdin representation back
to raw contracted-Gaussian AO coefficients, and writes a standard Molden `[MO]`
section with the natural-orbital occupation numbers (NOONs) as fractional
occupations.

## Install

```sh
# from a checkout
pip install .

# or straight from the repo
pip install git+https://github.com/uvar-apps/casci2molden.git
```

This pulls in `numpy` and `pyscf`. Python ≥ 3.10. The logs are version-tolerant but
were produced with PySCF 2.14.0.

## Usage

Point it at any directory tree containing `CASCI.log` files:

```sh
casci2molden /path/to/data                 # writes ./molden/*.molden
casci2molden /path/to/data -o /tmp/out    # custom output dir
casci2molden /path/to/data --dry-run      # list what would be written
```

Or as a module: `python -m casci2molden ...`.

### What each `CASCI.log` needs nearby

For every `CASCI.log`, the tool looks in the **same directory, then the parent** for:

| File | Provides |
|---|---|
| `*.xyz` | geometry (Ångström) |
| `RHF.log` or `ROHF.log` | `[INPUT] charge` and `[INPUT] spin` |

The basis name is read from the log (`Basis: 6-31g*`); if absent it falls back to
`6-31g*`.

### Output naming

The output filename is the log's **directory path relative to the scan root**, with
components joined by `_`:

```
<data>/F43/Ni_2/singlet/8e_5o/CASCI.log   ->  F43_Ni_2_singlet_8e_5o.molden
```

Use `--strip COMPONENT` (repeatable) to drop noisy path components, e.g.

```sh
casci2molden inputs/chemistry --strip Electronic_H
```

A JSON report (`<out>/conversion_report.json`, override with `--report`) records
per-file metrics: `nao`, `ncas`, `nelec`, and the orthonormality / occupation-sum
residuals.

### CLI reference

```
casci2molden [ROOT] [-o OUT] [--pattern GLOB] [--strip COMPONENT]
                     [--report PATH] [--dry-run] [--version]
```

## How it works (the math)

PySCF prints the CAS-space natural orbitals as

```
C_ml = Xᵀ S C_ao      with   X = lo.orth_ao(mol, 'meta_lowdin', s=S),  Xᵀ S X = I
```

Because `X` is square, `Xᵀ S = X⁻¹`, so the raw AO coefficients are recovered
exactly by

```
C_ao = X @ C_ml
```

`C_ao` is then written through `pyscf.tools.molden.from_mo` as a single
`Spin=Alpha` set with fractional `Occup=` values (the NOONs). EigenVista treats a
single Alpha set as **spatial orbitals** using those occupations — exactly the
natural-orbital picture.

**Precision:** the log prints coefficients to 5 decimals, so reconstructed
orbitals carry ~1e-5 rounding noise (column norms deviate from 1 by ~1e-4). This
is irrelevant for visualisation.

**Scope:** only CAS-space natural orbitals are present in the logs (core RHF
coefficients are not printed), and `*_FCIDUMP.txt` files are *not* convertible to
Molden (they carry no AO basis). They remain useful as sidecar data.

## Python API

```python
from casci2molden import convert_one, convert_tree, parse_casci_log, find_context

report = convert_one("CASCI.log", "out", root=".", strip=("Electronic_H",))
reports, failures = convert_tree("inputs/chemistry", "molden")
```

## Development

```sh
python -m venv .venv && source .venv/bin/activate
pip install -e . pytest
pytest
```

The test suite includes a real PySCF round-trip (H₂O CAS(2e,4o)) that verifies the
recovered AO coefficients are orthonormal under the overlap matrix.

## License

GPL-2.0-or-later. See [LICENSE](LICENSE). Third-party components (PySCF, NumPy)
retain their own licenses.
