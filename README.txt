The aim of the MstatX project is to provide an easy-to-use and easy-to-extend application to derive statistics from a multiple alignment file.

MstatX only takes a multiple alignment in input with the option flag -i or --input. The default statistic is the famous sum-of-pairs based on a blosum62 matrix. The description of the available statistics is in the manual and also in the code.

The output of MstatX consists in two files:
 - The output file with the statistics asked by the user.
 - The information file (.info) which gives basic information about the alignment.

Statistics currently available are :
 - sumofpairs (1)

This application is not designed to validate a multiple alignment but only to calculate a statistical score. In consequence, the multiple alignment, given in input, is supposed to be exact (obviously, this assumption is not true). MstatX was meant to compute statistics for each columns but with the flag -g, you can also output a global score of a multiple alignment (the mean of the column scores).

Mstatx is distributed under the term of the MIT licence. For any bug 
report or information, contact me at gcollet[AT]ouvaton.org or through my github account : https://github.com/gcollet

License:

Copyright (c) 2010 Guillaume Collet

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
