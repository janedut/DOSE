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



