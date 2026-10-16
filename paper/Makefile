SHELL := /bin/sh
ROOT := $(abspath ..)
BUILD := $(ROOT)/tmp/pdfs/exact-qop
OUTPUT := $(ROOT)/output/pdf/exact-qop.pdf
SOURCE := $(CURDIR)/exact-qop.tex

.PHONY: all check render
all: $(OUTPUT)

$(BUILD)/exact-qop.pdf: $(SOURCE) Makefile check_paper.py
	mkdir -p "$(BUILD)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD)" "$(SOURCE)"
	python3 check_paper.py "$(BUILD)/exact-qop.pdf" "$(BUILD)/exact-qop.log"

$(OUTPUT): $(BUILD)/exact-qop.pdf
	mkdir -p "$(dir $(OUTPUT))"
	cp "$(BUILD)/exact-qop.pdf" "$(OUTPUT)"

check: $(OUTPUT)
	python3 check_examples.py
	python3 check_paper.py "$(OUTPUT)" "$(BUILD)/exact-qop.log"

render: $(OUTPUT)
	mkdir -p "$(BUILD)/pages"
	pdftoppm -r 110 -png "$(OUTPUT)" "$(BUILD)/pages/page"
