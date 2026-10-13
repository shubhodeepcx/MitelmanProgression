# MitelmanProgression (Shubhodeep Chatterjee)
PY ?= python3

.PHONY: all data tables models validate figures report pdf test lint clean

all:            ## full pipeline and tests
	bash run_all.sh

tables:         ## steps 1-5 (needs data/raw)
	$(PY) -m karyoprog run --from 1 --to 5

models:         ## baseline, CBN, MHN, clone timing
	$(PY) -m karyoprog run --from 6 --to 9

validate:
	$(PY) -m karyoprog run --only 10

figures:
	$(PY) -m karyoprog run --only 11

report:         ## LaTeX and DOCX from the results
	$(PY) report/build_report.py

pdf: report     ## compile the IEEE PDF (needs pdflatex)
	cd report && pdflatex -interaction=nonstopmode main.tex >/dev/null && pdflatex -interaction=nonstopmode main.tex >/dev/null \
	  && mv main.pdf MitelmanProgression_report.pdf && rm -f main.aux main.log main.out

test:
	$(PY) -m pytest -q tests

lint:
	ruff check karyoprog scripts tests report

clean:
	rm -rf .pytest_cache **/__pycache__ report/main.aux report/main.log report/main.out
