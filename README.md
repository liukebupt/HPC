Instruction for using LAPACK in c on linux

1.  Installation
      Download latest version from http://www.netlib.org/lapack/.
      Uncompress tar file.
      $cd lapack-*/
      $cmake ./
      
2.  Test
      Get test program from http://www.netlib.org/lapack/lapacke.html.
      Compile with 
        $gcc test_lapack.c -L lib/ -llapacke -llapack -lrefblas -lgfortran -lm -I include/ -o test_lapack 
      
3.  Documentation
      http://www.netlib.org/lapack/explore-html/index.html.
