#include <stdio.h>
#include "sort.h"

void sort_lt(int n, int list[], double f[])
{
  int i, j, l, ir, ira;

  l = n/2+1;
  ir = n;
  while (n > 1)
  {
    if (l > 1)
    {
      l = l-1;
      ira = list[l-1];
    } else
    {
      ira = list[ir-1];
      list[ir-1] = list[0];
      ir = ir-1;
      if (ir == 1)
      {
        list[0] = ira;
        break;
      }
    }
    i = l;
    j = l*2;
    while (j <= ir)
    {
      if (j < ir)
        if ( f[list[j-1]] < f[list[j]] )
          j = j+1;
      if ( f[ira] < f[list[j-1]] )
      {
        list[i-1] = list[j-1];
        i = j;
        j = j+i;
      } else
        j = ir+1;
    }
    list[i-1] = ira;
  }
 
}

void sort_gt(int n, int list[], double f[])
{
  int i, j, l, ir, ira;

  l = n/2+1;
  ir = n;
  while (n > 1)
  {
    if (l > 1)
    {
      l = l-1;
      ira = list[l-1];
    } else
    {
      ira = list[ir-1];
      list[ir-1] = list[0];
      ir = ir-1;
      if (ir == 1)
      {
        list[0] = ira;
        break;
      }
    }
    i = l;
    j = l*2;
    while (j <= ir)
    {
      if (j < ir)
        if ( f[list[j-1]] > f[list[j]] )
          j = j+1;
      if ( f[ira] > f[list[j-1]] )
      {
        list[i-1] = list[j-1];
        i = j;
        j = j+i;
      } else
        j = ir+1;
    }
    list[i-1] = ira;
  }
 
}
