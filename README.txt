The aim of MstatX project is to provide an easy-to-use and easy-to-extend 
application to derive statistics from a multiple alignment file.

MstatX only takes a multiple alignment in input with the option -ma.
Then you can specify the statistic with the flag -s.

The default statistic is a weighted-entropy which means a Shannon 
entropy based on probabilities computed with the sequence weighting 
scheme defined by Henikoff & Henikoff.

Statistics currently available are :
 - wentropy (1)
 - trident (1)
 - kabat (1)
 - jensen (1)
 - mvector (1)
 - mlc (2)
 - mutual (2)
 - sca (2)

These statistics are divided in two. in one hand, statistics on one 
column (1), in the other hand, statistics on pairs of columns (2).
The description of each statistic measure is provided in the code.

This application is not designed to validate a multiple alignment. 
The given multiple alignment is supposed to be exact (obviously,
this assumption is not true).

Mstatx is distributed under the term of the MIT licence. For any bug 
report or information, contact me at : gcollet [AT] ouvaton.org

