## LAPACK
Instruction for using LAPACK in c on linux

### Installation
* Download latest version from http://www.netlib.org/lapack/.
* Uncompress tar file.
* `$cd lapack-*/`
* `$cmake ./`
### Test
* Get test program from http://www.netlib.org/lapack/lapacke.html.
* Compile with `gcc test_lapack.c -L lib/ -llapacke -llapack -lrefblas -lgfortran -lm -I include/ -o test_lapack` 
      
### Documentation
* http://www.netlib.org/lapack/explore-html/index.html.
