library(rappdirs)
dir <- rappdirs::user_data_dir()

dir <- file.path(dir, 'ontology')
if (!dir.exists(dir)) dir.create(dir)

library(RSQLite)
library(AnnotationDbi)

dbfile <- file.path(dir, 'HDO.sqlite')
dbfile <- "HDO.sqlite"
getwd()
file.exists(dbfile)
dbconn <- dbFileConnect(dbfile)
dbListTables(dbconn)      


library(AnnotationDbi)

x = loadDb(dbfile)
x
columns(x)




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





date <- '20240628'
name <- "Disease Ontology"
url <- "https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"

create_sqlite("Downloads/HumanDO.obo", "HDO.sqlite", name, date, url)


