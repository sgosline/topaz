=====
TOPAZ: Identifying transcriptional regulators from changing microRNA and  mRNA levels.
=====
Contact: Sara JC Gosline sgosline@mit.edu

Copyright (c) 2014-2015 Sara JC Gosline

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in th
e Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so, subj
ect to the following conditions:

The above copyright notice and this permission notice shall be included in all c
opies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLI
ED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYR
IGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WIT
H THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=====
To use:
________

1. Download SAMNet [from here](http://www.github.com/sgosline/SAMNet) and Garnet
from [OmicsIntegrator package](http://www.github.com/sgosline/OmicsIntegrator)
2. Use GARNet to create transcriptional regulatory networks from BED-formatted
files of histone or other chromatin accessibility data
3. Collect the required files for topaz (type `topaz.py --h` for file
descriptions)
4. Make sure you have optimization code (ampl) required for SAMNet
5. Run!
