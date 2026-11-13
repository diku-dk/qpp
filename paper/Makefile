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

# Every paper gets the TeX diagnostics and reference checks of check_paper.py;
# options for one paper go in CHECK_ARGS_<paper>.

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

# Re-run the checks on every built paper.
check: $(CHECKS)

$(CHECKS): check-%: $(OUTPUT)/%.pdf
	$(CHECK)

# Render pages for visual inspection, e.g. `make render-uop-algebra`.
render: $(RENDERS)

$(RENDERS): render-%: $(OUTPUT)/%.pdf
	mkdir -p "$(BUILD)/$*/pages"
	pdftoppm -r 110 -png "$<" "$(BUILD)/$*/pages/page"
