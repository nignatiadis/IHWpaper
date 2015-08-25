# generics related to multiple testing

# given the output of a method, this function returns a boolen vector of rejected hypotheses
rejected_hypotheses <- function(object,...) UseMethod("rejected_hypotheses")
