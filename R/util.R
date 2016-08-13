# Function to locate the full path to the helper application
# Returns:
#   "" - if the application is not found in the path.
#   <full_path> - that includes the name of the application
findApplication <- function(appname="") {
   if (appname=="") {
      return (FALSE)
   } else {
      #system("which plink", intern=T) # capture the full path of program into a variable
      result = tryCatch({
         system(paste0("which ", appname), intern=T)
      }, warning = function(w) {
         #cat(paste0("Cannot find ", appname, " in $PATH.\n"))
         return ("")
      }, error = function(e) {
         #cat(paste0("Cannot find ", appname, " in $PATH.\n"))
         return ("")
      }, finally = {
         #cleanup-code
      })
   }
} # findApplication()
