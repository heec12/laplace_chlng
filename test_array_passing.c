#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void func(int** array, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++)
  {
    for (j=0; j<cols; j++)
    {
      array[i][j] = i*j;
    }
  }
}

int main()
{
  int rows, cols, i;
  int **x;

  /* obtain values for rows & cols */

  /* allocate the array */
  x = malloc(rows * sizeof *x);
  for (i=0; i<rows; i++)
  {
    x[i] = malloc(cols * sizeof *x[i]);
  }

  /* use the array */
  func(x, rows, cols);

  /* deallocate the array */
  for (i=0; i<rows; i++)
  {
    free(x[i]);
  }
}
