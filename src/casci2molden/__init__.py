"""casci2molden — convert PySCF CASCI.log files into Molden files for EigenVista.

Public API:
    parse_casci_log(path)   -> dict   (basis, CAS size, NOONs, meta-Lowdin coefficients)
    find_context(log_path)  -> (xyz_path, charge, spin)
    convert_one(log_path, out_dir, root=..., strip=()) -> report dict
    convert_tree(root, out_dir, ...) -> (reports, failures)
"""
from __future__ import annotations

from .converter import convert_one, convert_tree
from .parser import find_context, parse_casci_log

__version__ = "0.1.0"

__all__ = [
    "parse_casci_log",
    "find_context",
    "convert_one",
    "convert_tree",
    "__version__",
]
