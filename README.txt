The aim of MstatX project is to provide an easy-to-use and easy-to-extend application to derive statistics from a multiple alignment file.

MstatX only takes a multiple alignment in input with the option -ma.

The statistic measure used is provided by the -s flag. 

The default statistic is a weighted-entropy which means a Shannon entropy + sequence weighting defined by Henikoff & Henikoff.

This application does not validate a multiple alignment. The given multiple alignment is supposed to be exact (that is obviously not true)