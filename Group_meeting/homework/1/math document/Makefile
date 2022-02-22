#include LaTeX.mk

# commands
LATEX   := xelatex
BIBTEX  := bibtex
RM      := rm -rf

default: math_document.pdf

%.pdf : %.tex bib png
	$(LATEX) $<
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<


clean:
	$(RM) *.aux *.bbl *.blg *.dvi *.log *.spl *~ */*~

purge:
	$(MAKE) clean
	$(RM) *.pdf

%.tar.gz : %.tex bib png
	mkdir $* $*/bib $*/png
	cp $*.tex Makefile $*/
	cp bib/$*.bib $*/bib/
	cp png/*.png $*/png/
	tar cvf $*.tar $*
	gzip $*.tar
	rm -rf $* $*.tar
