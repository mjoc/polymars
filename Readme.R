Changes to the original package:
1. Re-organisation of the files.
2. Conversion of the docs.
3. New Makefile. No problem with the C code.
4. Minor modifications to the interpreter code:
   a) persp -> image
   b) removed xtitle and ytitle to contour. They don't work
      on my machine. 
   c) converted multiline strings 
        "xxxxxxxxxx'
	 yyyyyyyyyy"
      to
       paste("xxxxxxxxxx",
	     "yyyyyyyyyy")
	     
guido masarotto
Dept. Statistical Sciences
University of Padua
Italy
http://sirio.stat.unipd.it
guido@sirio.stat.unipd.it
	 
