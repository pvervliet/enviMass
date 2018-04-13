LXBUW5535V <- function(infile, folderout, msconvert_path, notintern = FALSE, use_format = "mzXML", sim_as_scan = FALSE){
if(nchar(Sys.which("msconvert")[[1]])==0){
cat("msconvert not in system path - ok if msconvert_path correct")
}
if(
sum(substr(infile, nchar(infile) - 3, nchar(infile)) != ".RAW", substr(infile, nchar(infile) - 3, nchar(infile)) != ".raw") == 1
){cat("running .RAW file conversion.")}
if(!sim_as_scan){
there2 <- paste(" -o ",shQuote(folderout), sep = "")
filtered0 <- paste(shQuote("--"), use_format, sep = "")
filtered1 <- paste(shQuote("--32"), sep = "")
filtered2 <- paste(shQuote("--zlib"), sep = "")
filtered3 <- paste(" --filter ", shQuote("peakPicking true 1-2"), sep = "")
system(
paste(
shQuote(msconvert_path),
shQuote(infile),
filtered1,
filtered2,
filtered0,
filtered3,
there2
)
,intern = notintern)
}
if(sim_as_scan){
there2 <- paste(" -o ",shQuote(folderout), sep = "")
filtered0 <- paste(shQuote("--"),use_format, sep = "")
filtered1 <- paste(shQuote("--32"), sep = "")
filtered2 <- paste(shQuote("--zlib"), sep = "")
filtered3 <- paste(" --filter ", shQuote("peakPicking true 1-2"), sep = "")
filtered5 <- paste(shQuote("--simAsSpectra"), sep = "")
system(
paste(
shQuote(msconvert_path),
shQuote(infile),
filtered1,
filtered2,
filtered0,
filtered3,
filtered5,
there2
)
,intern = notintern)
}
}
