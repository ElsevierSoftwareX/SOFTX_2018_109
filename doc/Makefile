srcs = niget

.PHONY: all clean cleanall pdfdoc htmldoc

all: pdfdoc htmldoc

clean:
	-$(RM) *.pdf *.xml
	-$(RM) -r html4
	-$(RM) -r html5

cleanall: clean
	-$(RM) *.out *.aux *.toc *.log *.bbl *.blg
	-$(RM) *.css *.cache

pdfdoc:
	pdflatex $(srcs)
	bibtex $(srcs)
	pdflatex $(srcs)
	pdflatex $(srcs)

# predelat, aby dokazal zpracovat vic vstupnich souboru
htmldoc:
	latexml --destination=niget.xml $(srcs).tex
	latexmlpost --destination=html4/index.html --format=html4 --split niget.xml
	latexmlpost --destination=html5/index.html --format=html5 --split niget.xml
