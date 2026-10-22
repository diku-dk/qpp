SHELL := /bin/sh
ROOT := $(abspath ..)
BUILD := $(ROOT)/tmp/pdfs
OUTPUT := $(ROOT)/output/pdf
SOURCES := $(wildcard *.tex)
STEMS := $(SOURCES:.tex=)
TARGETS := $(addsuffix .pdf,$(STEMS))
PDFS := $(addprefix $(OUTPUT)/,$(TARGETS))
CHECKS := $(addprefix check-,$(STEMS))
RENDERS := $(addprefix render-,$(STEMS))

# Per-paper options for check_paper.py; every paper gets the TeX diagnostics
# and reference checks. exact-qop was budgeted at 15 pages and now has 16.
CHECK_ARGS_exact-qop := --max-pages 16 --require "HQP project" --require "MSc thesis"

# Check the build of paper $* in its build directory.
CHECK = python3 check_paper.py "$(BUILD)/$*/$*.pdf" "$(BUILD)/$*/$*.log" $(CHECK_ARGS_$*)

.PHONY: all check render $(TARGETS) $(CHECKS) $(RENDERS)
all: $(PDFS)

# A PDF reaches output/pdf only after it passes its checks.
$(OUTPUT)/%.pdf: %.tex Makefile check_paper.py
	mkdir -p "$(BUILD)/$*" "$(OUTPUT)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD)/$*" "$<"
	$(CHECK)
	cp "$(BUILD)/$*/$*.pdf" "$@"

# Allow e.g. `make uop-algebra.pdf`, keeping PDFs in output/pdf.
$(TARGETS): %.pdf: $(OUTPUT)/%.pdf

# Re-run the checks on every built paper, and the paper's worked examples.
check: $(CHECKS)
	python3 check_examples.py

$(CHECKS): check-%: $(OUTPUT)/%.pdf
	$(CHECK)

# Render pages for visual inspection, e.g. `make render-exact-qop`.
render: $(RENDERS)

$(RENDERS): render-%: $(OUTPUT)/%.pdf
	mkdir -p "$(BUILD)/$*/pages"
	pdftoppm -r 110 -png "$<" "$(BUILD)/$*/pages/page"
