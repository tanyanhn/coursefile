#include LaTeX.mk

BIBTEX  := bibtex
DVIPS   := dvips
LATEX   := pdflatex
PSPDF   := ps2pdf
RM      := rm -rf

default : COVID.pdf

%.pdf : %.tex bib sec tab
#	$(MAKE) --directory=pst
	$(LATEX) $<
	$(BIBTEX) $*
	$(LATEX) $<
	$(LATEX) $<

#%.ps  : %.dvi
#	dvips $< -t letter -o $@

#%.pdf : %.ps
#	$(PSPDF) $< $@

clean:
#	$(MAKE) --directory=pst clean
	$(RM) *.aux *.bbl *.blg *.log *.dvi ./eps/*.pdf

realclean:
#	$(MAKE) --directory=pst realclean
	$(RM) *.aux *.bbl *.blg *.log *.dvi *.pdf *~ ./eps/*.pdf

%.zip : %.tex bib pdf pst sec tab
	$(MAKE)
	mkdir -p $* $*/bib $*/eps $*/pdf $*/pst $*/sec $*/tab
	cp $*.tex Makefile $*/
	cp bib/*.bib $*/bib/
	cp -r pdf/*.pdf $*/pdf/
	cp -r eps/*.eps $*/eps/
	cp -r pst/Makefile pst/*.pdf $*/pst/
	cp sec/*.tex $*/sec/
	cp tab/*.tex $*/tab/
	zip -r $*.zip $*
	$(RM) $*

