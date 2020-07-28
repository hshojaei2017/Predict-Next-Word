# initialize
setwd("~/Documents/data/CapstoneProject/final/en_US")

# set options
options(stringsAsFactors = FALSE)

# get list of files
flist <- list.files()

# define a function to extract important file information
file_stat <- function(f) {
  fsize <- file.info(f)[1]/1024/1024
  lines <- readLines(f)
  nwords <- sum(sapply(strsplit(lines, "\\s+"), length))
  return(c(f, round(fsize, 2), length(lines), nwords))
}

# apply the function to all the files 
stats_list <- lapply(flist, file_stat)

# put the results in a data frame
df <- data.frame(matrix(unlist(stats_list), 
                        nrow=length(stats_list), byrow=T))
colnames(df) <- c("File", "Size(MB)", "Num_of_Lines", "Num_of_Words")
df

