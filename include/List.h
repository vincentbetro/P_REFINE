#include <stdio.h>
#include <stdlib.h>

#ifndef List_h
#define List_h

class List
{
  public:
  List()
  { dim = max = 0; list = 0; }
  ~List()
  {
    dim = max = 0;
    if (list != 0)
    {
      free(list);
      list = 0;
    }
  }
  void construct()
  { dim = max = 0; list = 0; }
  void destruct()
  {
    dim = max = 0;
    if (list != 0)
    {
      free(list);
      list = 0;
    }
  }
  List(List* s) // copy constructor
  {
    int i;
    Redimension(s->dim);
    max = 0;
    for (i=0; i < s->max; i++)
      Add_To_List(s->list[i]);
  }
  List(List& s) // copy constructor
  {
    int i;
    Redimension(s.dim);
    max = 0;
    for (i=0; i < s.max; i++)
      Add_To_List(s.list[i]);
  }
  List& operator=(List& s) // operator = for use when LHS already exists
  {
    int i;
    Redimension(s.dim);
    max = 0;
    for (i=0; i < s.max; i++)
      Add_To_List(s.list[i]);
    return *this;
  }
  void print_dim(FILE *outf)
  {
    fprintf(outf,"\nInteger list dimension =%d",dim);
  }
  void print(FILE *outf)
  { 
    print_dim(outf);
    fprintf(outf,"\nInteger list maximum index =%d",max);
    for (int i=0; i < max; i++)
      fprintf(outf,"\nlist(%d)= %d",i,list[i]);
  }
  int Redimension(int size);
  int Is_In_List(int n);
  int Check_List(int n);
  int Add_To_List(int n);
  int Delete_From_List(int n);
  int Replace(int m, int n);
  int Times_In_List(int n);
  int Index(int n);
  
  int *list;
  int max;
  private:
  int dim;
};

#endif
