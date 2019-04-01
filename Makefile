PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVER  := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: prep build check clean

prep:
	R -q -e "Rcpp::compileAttributes(); devtools::document()"

build:
	cd ..;\
	R CMD BUILD $(PKGSRC)

check:
	cd ..;\
	R CMD CHECK $(PKGNAME)_$(PKGVER).tar.gz

bioc-check:
	cd ..;\
	R -q -e "BiocCheck::BiocCheck('$(PKGNAME)_$(PKGVER).tar.gz', 'quit-with-status'=TRUE)"

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVER).tar.gz

clean:
	mkdir -p builds;\
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck;\
	mv -f $(PKGNAME)_$(PKGVER).tar.gz -t $(PKGNAME)/builds

.PHONY: vignettes benchmarks

vignettes:
	cd vignettes;\
	R -q -e "rmarkdown::render('AdvancedUsage.Rmd')" &&\
	R -q -e "rmarkdown::render('Introduction.Rmd')" &&\
	R -q -e "rmarkdown::render('IntroductionToSequenceMotifs.Rmd')" &&\
	R -q -e "rmarkdown::render('MotifManipulation.Rmd')" &&\
	R -q -e "rmarkdown::render('SequenceSearches.Rmd')" &&\
	$(RM) -rf *.log &&\
	$(RM) -rf *_files

benchmarks:
	mkdir -p benchmarks/results;\
	cd benchmarks;\
	./benchmarks.R
