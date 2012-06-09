## Adapted from Rcpp

## Make sure system.file returns an absolute path
Rplinkseq.system.file <- function(...){
	tools:::file_path_as_absolute( base:::system.file( ..., package = "Rplinkseq" ) )
}

## Use R's internal knowledge of path settings to find the lib/ directory
## plus optinally an arch-specific directory on system building multi-arch
RplinkseqLdPath <- function() {
	if (nzchar(.Platform$r_arch)) {	## eg amd64, ia64, mips
		path <- Rplinkseq.system.file("lib",.Platform$r_arch)
	} else {
		path <- Rplinkseq.system.file("lib")
	}
	path
}

# Provide linker flags -- i.e. -L/path/to/libRplinkseq -- as well as an
# optional rpath call needed to tell the dynamic linker about the
# location.  
RplinkseqLdFlags <- function() {
	rplinkseqdir <- RplinkseqLdPath()
	
	flags <- paste("-L", rplinkseqdir, " -lRplinkseq", sep="")  # Baseline setting
	if ((.Platform$OS.type == "unix") && (length(grep("^linux",R.version$os)))) {
			flags <- paste(flags, " -Wl,-rpath,", rplinkseqdir, sep="")
	}
	
	invisible(flags)
}

LdFlags <- function() {
	cat(RplinkseqLdFlags())
}

