"""Check the paper's page budget, references, and TeX layout diagnostics.

Uses the installed Poppler tools; visual inspection of rendered pages is
still required. Run from paper/ with make check.
"""
from pathlib import Path
import re
import subprocess
import sys

pdf, log = map(Path, sys.argv[1:])
info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
pages_match = re.search(r"^Pages:\s+(\d+)", info, re.MULTILINE)
assert pages_match, "Could not determine the PDF page count"
pages = int(pages_match[1])
assert 1 <= pages <= 15, f"Paper has {pages} pages; maximum is 15"
# TeX can emit raw font-encoding bytes when reporting problematic glyphs.
diagnostics = log.read_text(errors="replace")
for bad in ("Overfull", "undefined references", "Citation " + chr(96), "Missing character:"):
    assert bad not in diagnostics, f"TeX diagnostic needs review: {bad}"
text = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
assert "??" not in text, "Unresolved reference in rendered text"
assert "HQP project" in text and "MSc thesis" in text
print(f"Paper checks passed: {pages} pages; references and layout diagnostics clean.")
