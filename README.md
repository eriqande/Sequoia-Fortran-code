# Sequoia-Fortran-code

_No suffix:_ 
To be compiled and run under linux or windows + cygwin

This Fortran script is functionally identical to the script in the R package, but has partly optional extra output, useful for debugging. The code has been developped on a Windows 8 machine with Cygwin64 and gfortran. 

_suffix R:_
Requires Rtools. After renaming file to sequoia.f90, on the command line type: R CMD SHLIB sequoia.f90

## Compiling With a Mac Running El Capitan

I (Eric Anderson) struggled with this for a while, but I think to have gotten it, now.  Compiling 
Fortran on a Mac is always a bit of a drag.  The fortran compiler `gfortran-4.8` that is distributed 
from CRAN did not work due to a `kern.osversion` error. So, I installed a newer version, 
`gfortran 6.1`, of the compiler that I downloaded from [https://gcc.gnu.org/wiki/GFortranBinaries#Mac](https://gcc.gnu.org/wiki/GFortranBinaries#Mac).
Note that the installer for that compiler puts everything into `/usr/local/gfortran`.  That will become
important later.

After installing that I was able to compile the non-R version like so:
```sh 
gfortran -std=f95 -Wall -pedantic -fbounds-check Sequoia.f90 -o Seq
```
with no error messages or anything.  So, I knew I was on the right track.

But then the question was whether I could compile Sequoia_R.f90 using
`R CMD SHLIB`.  Doing so required some more finagling.  The Makeconf that
R comes with wants to use `gfortan-4.8`.  So, I modified my file,
`~/.R/Makevars` by adding these lines:
```make
FC = /usr/local/gfortran/bin/gfortran
F77 = /usr/local/gfortran/bin/gfortran
FLIBS = -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm
```