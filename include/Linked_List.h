#include <stdio.h>
#include <stdlib.h>

#ifndef Linked_List_h
#define Linked_List_h

class Linked_Node
{
  public:
  int data;
  Linked_Node *next;
  Linked_Node(int i, Linked_Node *n)
  { data = i, next = n; }
};

class Linked_List
{
  public:
  Linked_Node *head;
  Linked_List()
  { head = 0; }
  ~Linked_List()
  {
    Linked_Node *current = head;
    while (current)
    {
      Linked_Node *nxt = current->next;
      delete current;
      current = nxt;
    }
    head = 0;
  }
  void Insert(int i)
  { head = new Linked_Node(i, head); }

  int Length()
  {
    int l = 0;
    Linked_Node *current = head;
    while (current)
    {
      l++;
      current = current->next;
    }
    return(l);
  }

  void Remove(int i)
  {
    Linked_Node *current = head;
    Linked_Node *prev = 0;
    while (current && current->data != i)
    {
      prev = current;
      current = current->next;
    }
    if (current && current->data == i)
    {
      Linked_Node *nxt = current->next;
      if (prev == 0)
        head = nxt;
      else
        prev->next = nxt;
      delete current;
    }
  }

  int Replace(int i, int j)
  {
    int flag = 0;
    Linked_Node *current = head;
    while (current && current->data != i)
    {
      current = current->next;
    }
    if (current && current->data == i)
    {
      current->data = j;
      flag = 1;
    }

    return (flag);
  }

  int In_list(int i)
  {
    int flag = 0;
    Linked_Node *current = head;
    while (current && current->data != i)
    {
      current = current->next;
    }
    if (current && current->data == i)
      flag = 1;

    return (flag);
  }

};

#endif
