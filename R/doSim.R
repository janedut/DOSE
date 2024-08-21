##' measuring similarities between two DO term vectors.
##'
##' provide two term vectors, this function will calculate their similarities.
##' @title doseSim
##' @rdname doseSim
##' @param DOID1 DO term, MPO term or HPO term vector
##' @param DOID2 DO term, MPO term or HPO term vector
##' @param ont one of "HDO", "HPO" and "MPO"
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", "Lin", and "TCSS".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{https://yulab-smu.top}
doseSim <- function(DOID1,
                  DOID2,
                  measure="Wang",
                  ont = "HDO") {
    ont <- match.arg(ont, c("DO", "HDO", "MPO", "HPO"))                

    if (ont == "DO") ont <- 'HDO'

    processTCSS <- FALSE
    if (measure == "TCSS") {
        processTCSS <- TRUE
    } 

    scores <- GOSemSim::termSim(
        DOID1,
        DOID2, 
        semdata2(processTCSS = processTCSS, ont = ont), 
        measure
    )    

    if(length(scores) == 1)
        scores <- as.numeric(scores)

    return(scores)
}

##' @rdname doseSim
##' @export
doSim <- doseSim
