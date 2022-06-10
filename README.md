# gonad-development
R code for analysing gonad length relative to larval stage or pharynx-vulva length in C. elegans
We have measured gonad length from tip to tip using a lag-2::gfp (qIs56) marker for all larval stages (L1 to L4) in C. elegans, to be able to create a model where we 
can predict gonad length for a given larval stage or pharynx-vulva length. The R code can be used to analyse both types of data. First it will determine if the data are
normally distributed and then compare gonad length at each stage vs. every other. Then it will use the "segmented" package to determine how many breakpoints best fit 
the gonad length relative to pharynx-vulva length data, generating both the number of breakpoints, their location and the slope of the lines. This package was developed
by Dr. Robert Page at Texas A&M University-San Antonio and the model will (hopefully) be published at micropublication.org. Included is a plain text file of the data 
used in this analysis, which you are welcome to use to test the code. Thanks!
