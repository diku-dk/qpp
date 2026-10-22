"""Check a paper's TeX layout diagnostics, references, page budget and required text.

Usage: check_paper.py PDF LOG [--max-pages N] [--require TEXT]...

Every paper gets the diagnostics and reference checks; the page budget and
the required text are per-paper options, which the Makefile passes. Uses the
installed Poppler tools; visual inspection of rendered pages is still
required. Run from paper/ with make check.
"""
import argparse
from pathlib import Path
import re
import subprocess

parser = argparse.ArgumentParser(description="Check a built paper.")
parser.add_argument("pdf", type=Path)
parser.add_argument("log", type=Path, help="the latexmk log of the same build")
parser.add_argument("--max-pages", type=int, metavar="N", help="fail above N pages")
parser.add_argument("--require", action="append", default=[], metavar="TEXT",
                    help="text that must appear in the rendered paper; repeatable")
args = parser.parse_args()

info = subprocess.check_output(["pdfinfo", str(args.pdf)], text=True)
pages_match = re.search(r"^Pages:\s+(\d+)", info, re.MULTILINE)
assert pages_match, "Could not determine the PDF page count"
pages = int(pages_match[1])
assert pages >= 1, "Paper has no pages"
if args.max_pages is not None:
    assert pages <= args.max_pages, f"Paper has {pages} pages; maximum is {args.max_pages}"
# TeX can emit raw font-encoding bytes when reporting problematic glyphs.
diagnostics = args.log.read_text(errors="replace")
for bad in ("Overfull", "undefined references", "Citation " + chr(96), "Missing character:"):
    assert bad not in diagnostics, f"TeX diagnostic needs review: {bad}"
text = subprocess.check_output(["pdftotext", "-layout", str(args.pdf), "-"], text=True)
assert "??" not in text, "Unresolved reference in rendered text"
for required in args.require:
    assert required in text, f"Required text missing: {required}"
budget = f" of at most {args.max_pages}" if args.max_pages is not None else ""
print(f"{args.pdf.name}: {pages} pages{budget}; references and layout diagnostics clean.")
