my.get_eolid <- function (sciname, ask = TRUE, verbose = TRUE, key = NULL) {
    fun <- function(sciname, ask, verbose) {
      taxize:::mssg(verbose, "\nRetrieving data for taxon '", sciname, 
                    "'\n")
      tmp <- eol_search(terms = sciname, key = key)
      if (all(is.na(tmp))) {
        taxize:::mssg(verbose, "Not found. Consider checking the spelling or alternate classification")
        id <- NA
      } else {
        pageids <- tmp[grep(sciname, tmp$name), "pageid"]
        dfs <- compact(lapply(pageids, function(x) eol_pages(x, key = key)$scinames))
        dfs <- ldply(dfs[!sapply(dfs, nrow) == 0])
        df <- try(dfs[, c("identifier", "scientificname", "nameaccordingto")])
        if(class(df) %in% 'try-error') return(NA)
        names(df) <- c("eolid", "name", "source")
        df <- taxize:::getsourceshortnames(df)
        if (nrow(df) == 0) {
          taxize:::mssg(verbose, "Not found. Consider checking the spelling or alternate classification")
          id <- NA
        } else {
          id <- df$eolid
        }
      }
      if (length(id) == 0) {
        message("Not found. Consider checking the spelling or alternate classification")
        id <- NA
      }
      if (length(id) > 1) {
        id <- df$eolid[1]
#        if (ask) {
#          rownames(df) <- 1:nrow(df)
#          message("\n\n")
#          message("\nMore than one eolid found for taxon '", 
#                  sciname, "'!\n\n Enter rownumber of taxon (other inputs will return 'NA'):\n")
#          print(df)
#          take <- scan(n = 1, quiet = TRUE, what = "raw")
#          if (length(take) == 0) 
#            take <- "notake"
#          if (take %in% seq_len(nrow(df))) {
#            take <- as.numeric(take)
#            message("Input accepted, took eolid '", as.character(df$eolid[take]), "'.\n")
#            id <- as.character(df$eolid[take])
#          } else {
#            id <- NA
#            taxize:::mssg(verbose, "\nReturned 'NA'!\n\n")
#          }
#        } else {
#          id <- NA
#        }
      }
      return(id)
    }
  out <- laply(sciname, fun, ask, verbose)
  class(out) <- "eolid"
  return(out)
}
