.onLoad <- function(libname, pkgname) {
	invisible()
}

.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is eicm ",
                          utils::packageDescription("eicm", fields="Version"))
                          
	packageStartupMessage("Be sure to check the vignette:\nvignette(\"eicm\")")
}

