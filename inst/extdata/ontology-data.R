library(rappdirs)



dbfile <- "HDO.sqlite"
getwd()
file.exists(dbfile)
dbconn <- dbFileConnect(dbfile)
dbListTables(dbconn)      
dbReadTable(dbconn, 'offspring') -> xx
head(xx)
head(xx)
dbconn



library(AnnotationDbi)
.keys <- getFromNamespace(".keys", "AnnotationDbi")
.cols <- getFromNamespace(".cols", "AnnotationDbi")
smartKeys <- getFromNamespace("smartKeys", "AnnotationDbi")

.queryForKeys <- getFromNamespace(".queryForKeys", "AnnotationDbi")
dbQuery <- getFromNamespace("dbQuery", "AnnotationDbi")



head(keys(x, 'id'))
columns(x)

select(x, keys(x, 'id')[1:6], keytype='id', columns='parent') 


get_ont_info <- function(ontology) {
    ## selected columns of `genemap`
    cols <- c(2, 1)
    if (ontology == "HDO" || ontology == "DO") {
        check_pkg("HDO.db")
        # not used in current version
    } else if (ontology == "HPO") {
        check_pkg("HPO.db")
        genemap <- get_fun_from_pkg("HPO.db", "HPOGENE")
        ancmap <- get_fun_from_pkg("HPO.db", "HPOANCESTOR")
        termmap <- get_fun_from_pkg("HPO.db", "HPOTERM")
    } else if (ontology == "MPO") {
        check_pkg("MPO.db")
        genemap <- get_fun_from_pkg("MPO.db", "MPOMPMGI")
        ancmap <- get_fun_from_pkg("MPO.db", "MPOANCESTOR")
        termmap <- get_fun_from_pkg("MPO.db", "MPOTERM")
    } else if (ontology == "MDO") {
        check_pkg("MPO.db")
        genemap <- get_fun_from_pkg("MPO.db", "MPOMGIDO")
        cols <- c(1, 2)

        check_pkg("HDO.db")
        ancmap <- get_fun_from_pkg("HDO.db", "HDOANCESTOR")
        termmap <- get_fun_from_pkg("HDO.db", "HDOTERM")
    }
    # toTable(genemap)[, cols]
    res <- list(genemap = genemap,
            cols = cols,
            ancmap = ancmap,
            termmap = termmap
        )
}





head(xx)



columns(x)

xx <- function() {
x <- dbReadTable(con=dbconn, 'do_offspring')
#head(x)
xx <- split(x$offspring, x$doid)
}

head(xx)
xx['DOID:0040041']

yy <- function() {
Offsprings <- AnnotationDbi::as.list(HDO.db::HDOOFFSPRING)
}

Offsprings['DOID:0040041']


microbenchmark::microbenchmark(xx(), yy())



### create sqlite

library(obolite)
date <- '20240628'
name <- "Disease Ontology"
url <- "https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"

create_sqlite("Downloads/HumanDO.obo", "HDO.sqlite", name, date, url)


