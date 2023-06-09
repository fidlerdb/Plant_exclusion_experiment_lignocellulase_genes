# Create a function to capitalise every letter after a " " or a "-"
capStr <- function(text){
  st <- strsplit(text, "(?=[ -])", perl = TRUE)
  st <- lapply(st, FUN = function(x)paste(toupper(substring(x, 1,1)), substring(x, 2),
                                          sep = "", collapse = ""))
  st <- unlist(st)
  return(st)
}

# Test
#text <- c("glucose"
#          , "fucose", "xylose", "3,6-anhydro-D-galactose"
#          , "vanillic acid", "4-hydroxybenzoic acid", "benzoic acid")
#capStr(text)