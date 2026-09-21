PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVER  := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: prep build check clean

prep:
	R -q -e "Rcpp::compileAttributes(); devtools::document()"
	# roxygen2 7.3.2 emits a help-file name here instead of a topic alias.
	R -q -e "p <- 'man/reexports.Rd'; l <- readLines(p); l <- gsub('[BiocGenerics:row_colnames]', '[BiocGenerics:rownames]', l, fixed = TRUE); writeLines(l, p)"

build:
	cd ..;\
	R CMD build $(PKGSRC)

check:
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVER).tar.gz

bioc-check:
	cd ..;\
	R -q -e "BiocCheck::BiocCheck('$(PKGNAME)_$(PKGVER).tar.gz', 'quit-with-status'=TRUE)"

coverage:
	[[ -f "coverage.R" ]] && ./coverage.R || echo "missing coverage.R file"

install:
	cd ..;\
	R CMD install $(PKGNAME)_$(PKGVER).tar.gz

clean:
	mkdir -p builds;\
	mkdir -p Rcheck;\
	rm -f builds/$(PKGNAME)_$(PKGVER).tar.gz;\
	rm -rf Rcheck/$(PKGNAME).Rcheck_$(PKGVER);\
	mv -f ../$(PKGNAME).Rcheck Rcheck/$(PKGNAME).Rcheck_$(PKGVER);\
	mv -f ../$(PKGNAME)_$(PKGVER).tar.gz builds/$(PKGNAME)_$(PKGVER).tar.gz

.PHONY: vignettes benchmarks

vignettes:
	cd vignettes;\
	R -q -e "rmarkdown::render('MotifComparisonAndPvalues.Rmd')" &&\
	R -q -e "rmarkdown::render('Introduction.Rmd')" &&\
	R -q -e "rmarkdown::render('IntroductionToSequenceMotifs.Rmd')" &&\
	R -q -e "rmarkdown::render('MotifManipulation.Rmd')" &&\
	R -q -e "rmarkdown::render('SequenceSearches.Rmd')" &&\
	rm -rf *.log &&\
	rm -rf *_files

benchmarks:
	mkdir -p benchmarks/results;\
	cd benchmarks;\
	./benchmarks.R
