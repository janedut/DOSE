PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
BIOCVER := RELEASE_3_19

all: rd check clean

alldocs: rd readme

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd", encoding="UTF-8")'

build:
	#cd ..;\
	#R CMD build $(PKGSRC)
	Rscript -e 'devtools::build()'
	
build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: 
	##cd ..;\
	## Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'
	Rscript -e 'devtools::check()'
	
check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

gitmaintain:
	git gc --auto;\
	git prune -v;\
	git fsck --full

update:
	git fetch --all;\
	git checkout devel;\
	git merge upstream/devel;\
	git merge origin/devel


push:
	git push upstream devel;\
	git push origin devel

rmrelease:
	git branch -D $(BIOCVER)

release:
	git checkout $(BIOCVER);\
	git fetch --all

biocinit:
	git remote add upstream git@git.bioconductor.org:packages/$(PKGNAME).git;\
	git fetch --all
