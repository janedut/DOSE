##' @importClassesFrom AnnotationDbi AnnotationDb
setRefClass("OntDb", contains="AnnotationDb")

#' @importMethodsFrom AnnotationDbi keys
setMethod("keys", "OntDb",
    function(x, keytype, ...){
        if(missing(keytype)) keytype <- "id"
        term <- toTable(x)
        term[, keytype]
    }
)

#' @importMethodsFrom AnnotationDbi keytypes
setMethod("keytypes", "OntDb",
    function(x) {
        c("id", "term")
    }

)


#' @importMethodsFrom BiocGenerics toTable
setMethod("toTable", "OntDb",
    function(x) {
        dbReadTable(dbconn(x), 'term') |>
        setNames(c("id", "term"))
    }
)



#' @importMethodsFrom AnnotationDbi select
setMethod("select", "OntDb",
    function(x, keys, columns, keytype, ...){
        if (missing(keytype)) keytype <- "id"
        keytype <- match.arg(keytype, c("id","term"))
        strKeys <- paste0("\"", keys, "\"", collapse = ",")
        if (keytype == "term") {
            sql_key <- paste("SELECT doid FROM do_term WHERE term in (",
                strKeys, ")")
            doids <- AnnotationDbi:::dbQuery(dbconn(x), sql_key)[, 1]
            strKeys <- paste0("\"", doids, "\"", collapse = ",")
        }
        columns <- unique(c("id", columns))

        sqls <- paste("SELECT ", paste(columns, collapse = ","),
            " FROM term")
        columns2 <- setdiff(columns, c("id", "term"))
        for (col in columns2) {
            leftJoin <- paste0("LEFT JOIN  ", col, " USING (id)")
            sqls <- c(sqls, leftJoin)
        }
        sqls <- c(sqls, paste0("WHERE term.id in (", strKeys, ")"))
        sqls <- paste(sqls, collapse = " ")
        res <- dbQuery(dbconn(x), sqls)
        res
    }
)

#' @importMethodsFrom AnnotationDbi columns
setMethod("columns", "OntDb",
    function(x) {
        c("id","term", "alias", "synonym", "parent", "children",
            "ancestor", "offspring")
    }
)

get_do_offspring <- function(output = "list") {
    x <- load_onto('HDO')
    get_onto_data(x, output, 'offspring')
}

get_do_parent <- function(output = "list") {
    x <- load_onto('HDO')
    get_onto_data(x, output, 'parent')
}


get_do_ancestor <- function(output = "list") {
    x <- load_onto('HDO')
    get_onto_data(x, output, 'ancestor')
}

get_onto_data <- function(x, output='list', table="offspring") {
    output <- match.arg(output, c("data.frame", "list"))
    res <- dbReadTable(dbconn(x), table)
    if (output == 'data.frame') return(res)

    # column 1 is ID, column 2 is the related term
    split(res[,2], res[,1]) 
}

load_onto <- function(onto = "HDO") {
    .env <- get_dose_env()
    .onto <- sprintf(".onto_%s", onto)
    if (exists(.onto, envir=.env)) {
        db <- get(.onto, envir=.env)
        return(db)
    }

    dir <- rappdirs::user_data_dir()

    dir <- file.path(dir, 'ontology')
    if (!dir.exists(dir)) dir.create(dir)

    dbfile <- file.path(dir, sprintf("%s.sqlite", onto))

    if (!file.exists(dbfile)) {
        # download the file
    }

    setRefClass("OntDb", contains="AnnotationDb")
    db <- loadDb(dbfile)
    assign(.onto, db, envir = .env)
    return(db)
}

