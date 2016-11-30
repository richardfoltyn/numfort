                       README FILE FOR CDFLIB90

CDFLIB90 is a  set of Fortran 90 routines  that compute the cumulative
distribution  function,  inverse,  and  more  for  common  statistical
distributions.  The complete description  of functionality is found in
the documentation files.

============================== FILES ==============================

The following files should be present when you unpack this code.

INSTALLATION: Instructions on compiling  and installing the code.  The
instructions are  necessarily vague  because almost any  computer with
any operating system could contain these programs.

-----

LEGALITIES: Describes the legalities of  use of the code.  Briefly: We
place  our  efforts  in  the  public  domain, but  much  of  the  code
originates from ACM publications.   ACM retains copyright to this code
but  allows  the  use  in  any computer  for  any  purpose.   However,
distribution for commercial gain  requires permission from the ACM and
perhaps payment.

-----

README: This file.

============================== DIRECTORIES ==============================

The following directories should be present when you unpack this code.

-----

doc: Documentation for the  package.  The source for the documentation
is the LaTeX  document 'cdflib_doc.tex'.  Postscript ('cdflib_doc.ps')
was obtained from dvips.   The Acrobat pdf file, 'cdflib_doc.pdf', was
obtained from pdflatex.  The text file was obtained using HeVea.

Free Acrobat pdf readers for most computer systems can be obtained from
                         http://www.adobe.com

HeVea can be obtained (free) from
        http://pauillac.inria.fr/~maranget/hevea/index.html

-----

SOURCE:   Source  code   for   the  library.    Also  'Makefile'   and
'compile.cdflib90' that may be useful in installing this code.  (See
INSTALL.)

          ++++++++++++++++++++++++++++++++++++++++++++++++++

COMPARISON TO  DCDFLIB: The capabilities are generally  the same.  The
packaging  is better.   This is  one of  the major  advantages  of f95
MODULEs.  The  accuracy may be  a bit better:  we now always  invert a
function on  the minimum of CUM  and CCUM.  The speed  may be slightly
improved,  since we  use  a later  and  better zero  finder than  does
dcdflib.

          ++++++++++++++++++++++++++++++++++++++++++++++++++

C/C++ VERSION.   Alas, there is  none.  We have  a f77 to  C converter
which allowed  us to produce  DCDFLIB in  C, but we  have no f95  to C
converter.  (Please don't  tell us to use Sierra-Pacific's  f90 to f77
converter  and then  the f77  to C  converter.  The  odds  of anything
resembling  readable  C  code  following the  double  conversion  seem
miniscule.)   If anyone  knows of  a f95  to C  converter or  wants to
convert this code manually, we will  be glad to post the result on our
web site if the quality is adequate.
