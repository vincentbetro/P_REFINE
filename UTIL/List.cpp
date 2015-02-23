#include <stdio.h>
#include "List.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

int List::Redimension(int size)
{
  if (size <= 0)
  {
    if (list != 0)
      free(list);
    list = 0;
    dim = max = 0;
    return(1);
  } else if (size != dim) {
    if (list != 0)
      list=(int*)realloc((void*)list,size*sizeof(int));
    else
      list=(int*)malloc(size*sizeof(int));
    if (list == 0)
    {
      fprintf(stderr,"Could not allocate space for list.");
      return(0);
    }
    //for (int i=max; i < size; i++)
    //  list[i] = -1;
    dim = size;
    if (dim < max)
      max = dim;
  }
  return(1);
}

int List::Is_In_List(int n)
{
  int i, y;

  for (y=i=0; i < max && !y; i++)
    if (list[i] == n)
      y=1;

  return(y);
}

int List::Times_In_List(int n)
{
  int i, y;

  for (y=i=0; i < max; i++)
    if (list[i] == n)
      y++;

  return(y);
}

int List::Index(int n)
{
  int i, j;

  j = -1;
  for (i=0; i < max && j < 0; i++)
    if (list[i] == n)
      j = i;

  return(j);
}

int List::Check_List(int n)
{
  const int INC = 5;
  int new_dim;
  //new_dim = MIN(dim+100000,MAX(INC,dim*10));
  new_dim = dim+INC;
  if (!Is_In_List(n))
  {
    if (max >= dim)
    {
      if (!Redimension(new_dim))
      {
        fprintf(stderr,"Could not add to list.");
        return(0);
      }
    }
    list[max++] = n;
  }
  return(1);
}

int List::Add_To_List(int n)
{
  const int INC = 5;
  int new_dim;
  //new_dim = MIN(dim+100000,MAX(INC,dim*10));
  new_dim = dim+INC;
  if (max < dim)
    list[max++] = n;
  else if (Redimension(new_dim))
  {
    list[max++] = n;
    return(1);
  } else {
    fprintf(stderr,"Could not add to list.");
    return(0);
  }
  return(1);
}

int List::Delete_From_List(int n)
{
  int i, j, flag = 0;

  for (i=0; i < max; i++)
    if (list[i] == n)
    {
      for (j=i; j < max-1; j++)
        list[j] = list[j+1];
      max--;
      flag = 1;
      break;
    }

  return(flag);
}

int List::Replace(int m, int n)
{
  int i, flag = 0;

  for (i=0; i < max; i++)
    if (list[i] == m)
    {
      list[i] = n;
      flag = 1;
    }

  return(flag);
}
