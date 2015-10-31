# PolyMARS

An Splus package, writen in S and C that I wrote while studying for an M.S. at the University of Washington.

It was written under the supervision of Charles Kooperberg, now at The Fred Hutchison Cancer Research Center.

It is not as old as the dates on the github repository would suggest (I have forked someone else's repo version), but is was about the first implementations of a version of Jerome Friedman's MARS procedure, free in wild.
It is not exactly MARS, as it uses different basis functions for the splines, and can do Polychotomous Regression (multiple Y values to create a classifier) for Classification. It was submitted to the Statlib Splus open-source package repository.

I didn't have a home computer for a number of years after leaving college (it was back in the day) so I had essentially orphaned it; happily others have taken good care of it.

The package was ported to R by Guido Masarotto and was eventually combined into another R package (polspline) by Charles Kooperberg, so is still available although newer better, fancier versions are now available in R packages.
