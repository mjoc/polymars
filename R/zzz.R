.First.lib <- function(lib, pkg) {
  library.dynam("polymars", pkg, lib)
}

if (version$minor < "62.0")
  library.dynam("polymars")

