#include <stdio.h>
#include "List.h"

#ifndef Smesh_obj_h
#define Smesh_obj_h

class mesh_obj
{
 public:
  mesh_obj()
  { nn=ntet=npyr=npri=nhex=nb=0;
    node=0;
    t_n=0;
    q_n=0;
    tet_n=0;
    pyr_n=0;
    pri_n=0;
    hex_n=0;
    nt=0;
    nq=0;
    bname=0;
    pmap = 0;
    tri_map = 0;
    quad_map = 0;
    tet_map = 0;
    pyr_map = 0;
    pri_map = 0;
    hex_map = 0;
    tri_emap = 0;
    quad_emap = 0;
    tet_emap = 0;
    pyr_emap = 0;
    pri_emap = 0;
    hex_emap = 0;
    
  }
  ~mesh_obj()
  {
    int b,i;
    for (b=0; b < nb; b++)
    {
      if (nt[b] > 0 && t_n != 0)
      {
        for (i=0; i < nt[b]; i++)
          free(t_n[b][i]);
        free(t_n[b]);
      }
      if (nq[b] > 0 && q_n != 0)
      {
        for (i=0; i < nq[b]; i++)
          free(q_n[b][i]);
        free(q_n[b]);
      }
      if (nt[b] > 0 && tri_map != 0)
      {
        free(tri_map[b]);
      }
      if (nq[b] > 0 && quad_map != 0)
      {
        free(quad_map[b]);
      }
      if (nt[b] > 0 && tri_emap != 0)
      {
        for (i=0; i < nt[b]; i++)
          free(tri_emap[b][i]);
        free(tri_emap[b]);
      }
      if (nq[b] > 0 && quad_emap != 0)
      {
        for (i=0; i < nq[b]; i++)
          free(quad_emap[b][i]);
        free(quad_emap[b]);
      }
      free(bname[b]);
    }
    if (ntet > 0)
    {
      for (i=0; i < ntet; i++)
        free(tet_n[i]);
      free(tet_n);
    }
    if (npyr > 0)
    {
      for (i=0; i < npyr; i++)
        free(pyr_n[i]);
      free(pyr_n);
    }
    if (npri > 0)
    {
      for (i=0; i < npri; i++)
        free(pri_n[i]);
      free(pri_n);
    }
    if (nhex > 0)
    {
      for (i=0; i < nhex; i++)
        free(hex_n[i]);
      free(hex_n);
    }
    free(node);
    if (t_n != 0)
      free(t_n);
    if (q_n != 0)
      free(q_n);
    free(bname);
    if (tri_map != 0)
      free(tri_map);
    if (quad_map != 0)
      free(quad_map);
    if (tet_map != 0)
      free(tet_map);
    if (pyr_map != 0)
      free(pyr_map);
    if (pri_map != 0)
      free(pri_map);
    if (hex_map != 0)
      free(hex_map);
    if (tri_emap != 0)
      free(tri_emap);
    if (quad_emap != 0)
      free(quad_emap);
    if (nt != 0)
      free(nt);
    if (nq != 0)
      free(nq);
    if (pmap != 0)
    {
      for (i = 0; i < nn; i++)
        free(pmap[i]);
      free(pmap);
    }
    if (tet_emap != 0)
    {
      for (i = 0; i < ntet; i++)
        free(tet_emap[i]);
      free(tet_emap);
    }
    if (pyr_emap != 0)
    {
      for (i = 0; i < npyr; i++)
        free(pyr_emap[i]);
      free(pyr_emap);
    }
    if (pri_emap != 0)
    {
      for (i = 0; i < npri; i++)
        free(pri_emap[i]);
      free(pri_emap);
    }
    if (hex_emap != 0)
    {
      for (i = 0; i < nhex; i++)
        free(hex_emap[i]);
      free(hex_emap);
    }
    nn=ntet=npyr=npri=nhex=nb=0;
  }
  int mesh_stats();
  int refine(int target, int tri_flag, int pri_stack, double threshold, char aname[], int mltype, int bnd_flag);
  int refine_io(int mode, char sname[]);
  int make_nabors( List **tethash, List **pyrhash, List **prihash, List **hexhash, int timing);
  void element_pi();
  void exchange_emap(int **emap, int nelem, List **hash, int type, int *elem_map, int bndy);
  void create_emap(int nelem, int type, int **elem_node, List **hash, int **emap);

  int nn, ntet, npyr, npri, nhex, nb;
  double **node;
  int ***t_n;
  int ***q_n;
  int **tet_n;
  int **pyr_n;
  int **pri_n;
  int **hex_n;
  int *nt;
  int *nq;
  char **bname;
  int **pmap;
  int **tri_map, **quad_map, *tet_map, *pyr_map, *pri_map, *hex_map;
  int ***tri_emap, ***quad_emap, **tet_emap, **pyr_emap, **pri_emap, **hex_emap;
};
#endif

