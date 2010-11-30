This Application uses openMP so you need g++ 4.2 version.
If not, just remove the -fopenmp option from Makefile

mstatx takes a multiple alignment in input with the option -ma.
You can specifie a statistic measure by the -s flag. The default statistic is wentropy which means a Shannon entropy + sequence weighting defined by Henikoff & Henikoff.