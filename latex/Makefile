SHELL := /bin/sh
ROOT := $(abspath ..)
OUTPUT := $(ROOT)/output
BUILD := $(ROOT)/tmp/pdfs/tikz
SOURCES := qpp-tikz.sty operator-diagrams.tex ../exercises/TikZOperators.hs ../src/QPP/PrettyPrint/TikZ.hs
SOURCES += ../src/Programs/QFT.hs ../src/QPP/Syntax.hs $(wildcard ../src/QPP/UOp/*.hs) ../qpp.cabal
SOURCES += ../src/QPP/PrettyPrint/QCirc.hs $(wildcard ../src/QPP/QCirc/*.hs)

.PHONY: all
all: $(OUTPUT)/pdf/operator-diagrams.pdf

$(OUTPUT)/pdf/operator-diagrams.pdf: $(SOURCES) Makefile
	cd "$(ROOT)" && cabal run TikZOperators
	mkdir -p "$(BUILD)" "$(OUTPUT)/pdf"
	cd "$(OUTPUT)/tikz" && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD)" operator-diagrams.tex
	cp "$(BUILD)/operator-diagrams.pdf" "$@"
