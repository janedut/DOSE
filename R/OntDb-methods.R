##' @importClassesFrom AnnotationDbi AnnotationDb
setRefClass("OntDb", contains="AnnotationDb")

#' @importMethodsFrom AnnotationDbi keys
setMethod("keys", "OntDb",
    function(x, keytype, ...){
        if(missing(keytype)) keytype <- "id"
        term[, keytype]
    }
)

#' @importMethodsFrom AnnotationDbi keytypes
setMethod("keytypes", "OntDb",
    function(x) {
        c("id", "term")
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

