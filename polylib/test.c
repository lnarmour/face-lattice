#include <stdio.h>
#include <polylib/polylib64.h>


int main() {
  Matrix *a;
  Polyhedron *A;

  a = Matrix_Read(); /* 1 */

  A = Constraints2Polyhedron(a, 500); /* 3 */ 

  Matrix_Free(a); /* 5 */ 

  a = Polyhedron2Rays(A); /* 7 */ 

  Matrix_Print(stdout,P_VALUE_FMT,a); /* 9 */ 

  Matrix_Free(a); /* 11 */ 

  return 0;
}

