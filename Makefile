PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVER  := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: prep build check clean

prep:
	R -q -e "Rcpp::compileAttributes(); devtools::document()";\
	R -q -e "l <- readLines(con <- file('man/reexports.Rd')); close(con); l <- sub('[BiocGenerics]{ncol}', '[BiocGenerics:nrow]{ncol}', l, fixed = TRUE); l <- sub('[BiocGenerics]{rownames}', '[BiocGenerics:row_colnames]{rownames}', l, fixed = TRUE); l <- sub('[BiocGenerics]{colnames}', '[BiocGenerics:row_colnames]{colnames}', l, fixed = TRUE); l <- sub('[BiocGenerics]{rowSums}', '[BiocGenerics:colSums]{rowSums}', l, fixed = TRUE); l <- sub('[BiocGenerics]{colSums}', '[BiocGenerics:colSums]{colSums}', l, fixed = TRUE); l <- sub('[BiocGenerics]{rowMeans}', '[BiocGenerics:colSums]{rowMeans}', l, fixed = TRUE); l <- sub('[BiocGenerics]{colMeans}', '[BiocGenerics:colSums]{colMeans}', l, fixed = TRUE); writeLines(l, 'man/reexports.Rd')";\

build:
	cd ..;\
	R CMD BUILD $(PKGSRC)

check:
	cd ..;\
	R CMD CHECK $(PKGNAME)_$(PKGVER).tar.gz

bioc-check:
	cd ..;\
	R -q -e "BiocCheck::BiocCheck('$(PKGNAME)_$(PKGVER).tar.gz', 'quit-with-status'=TRUE)"

coverage:
	[[ -f "coverage.R" ]] && ./coverage.R || echo "missing coverage.R file"

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVER).tar.gz

clean:
	mkdir -p builds;\
	mkdir -p Rcheck;\
	cd ..;\
	mv -f $(PKGNAME).Rcheck $(PKGSRC)/Rcheck/$(PKGNAME).Rcheck_$(PKGVER); \
	mv -f $(PKGNAME)_$(PKGVER).tar.gz -t $(PKGNAME)/builds

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

