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
After which I did like so:
```sh
R CMD SHLIB Sequoia_R.f90 --output=sequoia.so
```
which worked, but gave a few errors:
```
/usr/local/gfortran/bin/gfortran  -fPIC -Wall -g -O2  -c  Sequoia_R.f90 -o Sequoia_R.o
Sequoia_R.f90:3093:65:

 double precision :: LLA(2,7,7), dLL, TopLL, LLcp(3,2), LLINOUT(2), LLU(3), LLtmp(3), ALR(3)
                                                                 1
Warning: Unused variable ‘llinout’ declared at (1) [-Wunused-variable]
Sequoia_R.f90:3623:71:

 double precision :: LLtmp(2), ALR(7), LLx(6), LLz(2), LRHS, LLM(6), LUX, LLHA(3), &
                                                                       1
Warning: Unused variable ‘lux’ declared at (1) [-Wunused-variable]
Sequoia_R.f90:6293:0:

 integer :: x, y, z, gcat
 
Warning: ‘gcat’ may be used uninitialized in this function [-Wmaybe-uninitialized]
Sequoia_R.f90:2472:0:

                 do j=1,nB
 
Warning: ‘nb’ may be used uninitialized in this function [-Wmaybe-uninitialized]
Sequoia_R.f90:1740:0:

 integer :: l, x, y, curGP(2), m, z, AncB(2, mxA), v
 
Warning: ‘m’ may be used uninitialized in this function [-Wmaybe-uninitialized]
/usr/local/gfortran/bin/gfortran -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o sequoia.so Sequoia_R.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
```
So, now I am going to see if I can swap that in for the lib in the package...
