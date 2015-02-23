#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/param.h>
#include "CGNS.h"
#include "List.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))

#include "cgnslib.h"

double dist(double *p1, double *p2)
{
  double ds = 0.0;
  for (int i=0; i < 3; i++)
    ds += (p2[i]-p1[i])*(p2[i]-p1[i]);
  return(sqrt(ds));
}

void CGNS_Errorcheck(int ier)
{
  char *error_message;
  if (ier != 0)
  {
    fprintf(stderr,"\nCGNS Error detected. Error code %d",ier);
    error_message = (char *)cg_get_error();
    fprintf(stderr,"\n      Error message = %s",error_message);
    fflush(stderr);
    //cg_error_exit();
  }
}

//
// These routines assume the CGNS files connectivities that start at node 1
// Internally the indexing is convert to "C" style, starting at node 0
//

// Serial wrapper for parallel CGNS read function
// Pointers to pointers are passed so the allocated memory can be passed out of the routine
// THIS ROUTINE ASSUMES THE POINTERS ARE NOT ALREADY ALLLOCATED!!!!!!
int CGNS_read(char *fname, int &nn, double ***p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn)
{
  int **nmap, **tri_map, **quad_map;
  int *tet_map, *pyr_map, *pri_map, *hex_map;
  int parallel = 0;
  nmap = tri_map = quad_map = 0;
  tet_map = pyr_map = pri_map = hex_map = 0;

  int flag = P_CGNS_read(parallel,fname,nn,p,nb,b_name,nt,tri_conn,nq,quad_conn,
                         ntet,tet_conn,npyr,pyr_conn,npri,pri_conn,nhex,hex_conn,
                         &nmap,&tri_map,&quad_map,&tet_map,&pyr_map,&pri_map,&hex_map);

  return(flag);
}

// Serial wrapper for parallel CGNS write function
int CGNS_write(char *fname, int nn, double **p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn)
{
  int **nmap, **tri_map, **quad_map;
  int *tet_map, *pyr_map, *pri_map, *hex_map;
  int parallel = 0;
  nmap = tri_map = quad_map = 0;
  tet_map = pyr_map = pri_map = hex_map = 0;

  int flag = P_CGNS_write(parallel,fname,nn,p,nb,b_name,nt,tri_conn,nq,quad_conn,
               ntet,tet_conn,npyr,pyr_conn,npri,pri_conn,nhex,hex_conn,
               nmap,tri_map,quad_map,tet_map,pyr_map,pri_map,hex_map);

  return(flag);
}

// Parallel CGNS read function
// Pointers to pointers are passed so the allocated memory can be passed out of the routine
// THIS ROUTINE ASSUMES THE POINTERS ARE NOT ALREADY ALLLOCATED!!!!!!
int P_CGNS_read(int parallel, char *fname, int &nn, double ***p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn,
              int ***nmap, int ***tri_map, int ***quad_map,
              int **tet_map, int **pyr_map, int **pri_map, int **hex_map)
{
  int b, c, i, j, m, n, s, ns, nstart, nend, nbndry, pflag;
  int icelldim, iphysdim, flag, nzones, nn_zone, ntet_zone, npyr_zone, npri_zone, nhex_zone;
  int index_file, index_base, index_zone, narrays, nb_zone;
  int isize[3], ndim, idim;
  int ier;
  DataType_t itype;
  ZoneType_t zonetype;
  BCType_t ibocotype;
  int npts,normalindex[3],normallistflag,ndataset,normallist;
  PointSetType_t iptset;
  DataType_t normaldatatype;
  int *conn;
  char *error_message;
  char basename[32];
  char zonename[32];
  char sectionname[74];
  char array_name[32];
  ElementType_t etype;
  double *x;
  List *zone_nodes = 0;

  pflag = 0;
  flag = 0;

  conn = 0;

  // initialize counters to 0
  nn = nb = ntet = npyr = npri = nhex = 0;
  // initialize pointer to 0
  (*b_name) = 0;
  (*nt) = 0;
  (*tri_conn) = 0;
  (*nq) = 0;
  (*quad_conn) = 0;
  (*tet_conn) = 0;
  (*pyr_conn) = 0;
  (*pri_conn) = 0;
  (*hex_conn) = 0;

  ier = cg_open(fname,CG_MODE_READ,&index_file);
  CGNS_Errorcheck(ier);

  ier = cg_nbases(index_file,&index_base);
  CGNS_Errorcheck(ier);
  if (index_base != 1)
    fprintf(stderr,"\nP_CGNS_read: Incorrect number of bases!");
  index_base=1;

  ier = cg_base_read(index_file,index_base,basename,&icelldim,&iphysdim);
  CGNS_Errorcheck(ier);

  ier = cg_nzones(index_file,index_base,&nzones);
  CGNS_Errorcheck(ier);

  if (nzones > 1)
    zone_nodes = new List();

  for (index_zone=1; index_zone <= nzones; index_zone++)
  {
    // initialize counters to 0
    nn_zone = nb_zone = ntet_zone = npyr_zone = npri_zone = nhex_zone = 0;
    ier=cg_zone_read(index_file,index_base,index_zone,zonename,isize);
    CGNS_Errorcheck(ier);
    ier=cg_zone_type(index_file,index_base,index_zone,&zonetype);
    CGNS_Errorcheck(ier);
    if (zonetype != Unstructured)
    {
      fprintf(stderr,"\nP_CGNS_read: Incorrect zone type!");
      error_message = (char *)cg_get_error();
      cg_error_exit();
    }

    if (isize[1] == 0)  // Are there elements in the zone?
      continue;

    nn_zone = isize[0];
    if (nn_zone > 0)
    {
      (*p) = (double **)realloc((void*)(*p),(nn+nn_zone)*sizeof(double *));
      if ((*p) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for point array!\n");
        fflush(stderr);
        exit(0);
      }
      for (n=nn; n < nn+nn_zone; n++)
        (*p)[n] = (double *)malloc(3*sizeof(double));

      x = (double*)malloc(nn_zone*sizeof(double));
      if (x == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to allocate memory for X array!\n");
        fflush(stderr);
        exit(0);
      }
      i=1;
      j=nn_zone;
      ier=cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
                        RealDouble,&i,&j,x);
      CGNS_Errorcheck(ier);
      for (n=0; n < nn_zone; n++)
        (*p)[n+nn][0] = x[n];
      ier=cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
                        RealDouble,&i,&j,x);
      CGNS_Errorcheck(ier);
      for (n=0; n < nn_zone; n++)
        (*p)[n+nn][1] = x[n];
      ier=cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
                        RealDouble,&i,&j,x);
      CGNS_Errorcheck(ier);
      for (n=0; n < nn_zone; n++)
        (*p)[n+nn][2] = x[n];

      free(x);
    }

    int (*bc_range)[2];

    cg_goto(index_file,index_base,"Zone_t",1,"ZoneBC_t",1,"end");
    ier=cg_nbocos(index_file,index_base,index_zone,&nb_zone);
    CGNS_Errorcheck(ier);
    if (nb_zone > 0)
    {
      bc_range = new int[nb_zone][2];
      if (bc_range == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to allocate memory for bc_range array!\n");
        fflush(stderr);
        exit(0);
      }
      (*b_name) = (char**)realloc((void*)(*b_name),(nb+nb_zone)*sizeof(char*));
      if ((*b_name) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to allocate memory for *b_name array!\n");
        fflush(stderr);
        exit(0);
      }
      for (b=1; b <= nb_zone; b++)
      {
        (*b_name)[nb+b-1] = (char*)malloc(33*sizeof(char));
        if ((*b_name)[nb+b-1] == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to allocate memory for (*b_name)[%d] array!\n",nb+b-1);
          fflush(stderr);
          exit(0);
        }

        ier=cg_boco_info(index_file,index_base,index_zone,b,sectionname,&ibocotype,
                         &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
        CGNS_Errorcheck(ier);
        strcpy((*b_name)[nb+b-1],sectionname);
        if (iptset != ElementRange)
        {
          fprintf(stderr,"\nError.  For this program, BCs must be set up as ElementRange \
                 type %s\n",PointSetTypeName[iptset]);
          exit(0);
        }
        
        ier=cg_boco_read(index_file,index_base,index_zone,b,bc_range[b-1],&normallist);
        CGNS_Errorcheck(ier);
      }

      (*nt) = (int*)realloc((void*)(*nt),(nb+nb_zone)*sizeof(int));
      if ((*nt) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *nt array!\n");
        fflush(stderr);
        exit(0);
      }
      (*tri_conn) = (int***)realloc((void*)(*tri_conn),(nb+nb_zone)*sizeof(int**));
      if ((*tri_conn) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *tri_conn array!\n");
        fflush(stderr);
        exit(0);
      }
      (*nq) = (int*)realloc((void*)(*nq),(nb+nb_zone)*sizeof(int));
      if ((*nq) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *nq array!\n");
        fflush(stderr);
        exit(0);
      }
      (*quad_conn) = (int***)realloc((void*)(*quad_conn),(nb+nb_zone)*sizeof(int**));
      if ((*quad_conn) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *quad_conn array!\n");
        fflush(stderr);
        exit(0);
      }
      for (n=0; n < nb_zone; n++)
      {
        (*nt)[nb+n] = 0;
        (*tri_conn)[nb+n] = 0;
        (*nq)[nb+n] = 0;
        (*quad_conn)[nb+n] = 0;
      }
    }

    cg_goto(index_file,index_base,"Zone_t",1,"end");
    ier = cg_nsections(index_file, index_base, index_zone, &ns);
    CGNS_Errorcheck(ier);

    m = 0;
    for (s=1; s <= ns; s++)
    {
      ier = cg_section_read(index_file,index_base,index_zone,s,sectionname,
                             &etype,&nstart,&nend,&nbndry,&pflag);
      CGNS_Errorcheck(ier);

      // if boundary section determine boundary number
      b=-1;
      for (j=0; j < nb_zone && b < 0; j++)
        if (nstart == bc_range[j][0] && nend == bc_range[j][1])
          b = j+nb;

      if (etype == MIXED)
      {
        conn = (int*)realloc((void*)conn,9*(nend-nstart+1)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);

        // now read each element and process accordingly
        j=0;
        for(c=0;c < (nend-nstart+1); c++)
        {
          switch(conn[j++])
          {
            case TRI_3:
              if (b < 0)
              {
                //fprintf(stderr,"\nTri elements and no boundary identified!");
                //fprintf(stderr,"\nAssuming these are zone connectivity elements!");
                zone_nodes->Check_List(conn[j++]-1+nn);
                zone_nodes->Check_List(conn[j++]-1+nn);
                zone_nodes->Check_List(conn[j++]-1+nn);
              } else
              {
                (*tri_conn)[b] = (int**)realloc((void*)(*tri_conn)[b],((*nt)[b]+1)*sizeof(int*));
                if ((*tri_conn)[b] == NULL)
                {
                  fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tri_conn)[%d] array!\n",b);
                  fflush(stderr);
                  exit(0);
                }
                for (i=(*nt)[b]; i < (*nt)[b]+1; i++)
                {
                  (*tri_conn)[b][i] = (int*)malloc(3*sizeof(int));
                  (*tri_conn)[b][i][0] = conn[j++]-1 + nn;
                  (*tri_conn)[b][i][1] = conn[j++]-1 + nn;
                  (*tri_conn)[b][i][2] = conn[j++]-1 + nn;
                }
                (*nt)[b]++;
              }
              break;
            case QUAD_4:
              if (b < 0)
              {
                //fprintf(stderr,"\nQuad elements and no boundary identified!");
                //fprintf(stderr,"\nAssuming these are zone connectivity elements!");
                zone_nodes->Check_List(conn[j++]-1+nn);
                zone_nodes->Check_List(conn[j++]-1+nn);
                zone_nodes->Check_List(conn[j++]-1+nn);
                zone_nodes->Check_List(conn[j++]-1+nn);
              } else
              {
                (*quad_conn)[b] = (int**)realloc((void*)(*quad_conn)[b],((*nq)[b]+1)*sizeof(int*));
                if ((*quad_conn)[b] == NULL)
                {
                  fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*quad_conn)[%d] array!\n",b);
                  fflush(stderr);
                  exit(0);
                }
                for (i=(*nq)[b]; i < (*nq)[b]+1; i++)
                {
                  (*quad_conn)[b][i] = (int*)malloc(4*sizeof(int));
                  (*quad_conn)[b][i][0] = conn[j++]-1 + nn;
                  (*quad_conn)[b][i][1] = conn[j++]-1 + nn;
                  (*quad_conn)[b][i][2] = conn[j++]-1 + nn;
                  (*quad_conn)[b][i][3] = conn[j++]-1 + nn;
                }
                (*nq)[b]++;
              }
              break;
            case TETRA_4:
              (*tet_conn) = (int**)realloc((void*)(*tet_conn),(ntet+ntet_zone+1)*sizeof(int*));
              if ((*tet_conn) == NULL)
              {
                fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tet_conn) array!\n");
                fflush(stderr);
                exit(0);
              }
              (*tet_conn)[ntet+ntet_zone] = (int*)malloc(4*sizeof(int));
              (*tet_conn)[ntet+ntet_zone][0] = conn[j++]-1 + nn;
              (*tet_conn)[ntet+ntet_zone][1] = conn[j++]-1 + nn;
              (*tet_conn)[ntet+ntet_zone][2] = conn[j++]-1 + nn;
              (*tet_conn)[ntet+ntet_zone][3] = conn[j++]-1 + nn;
              ntet_zone++;
              break;
            case PYRA_5:
              (*pyr_conn) = (int**)realloc((void*)(*pyr_conn),(npyr+npyr_zone+1)*sizeof(int*));
              if ((*pyr_conn) == NULL)
              {
                fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pyr_conn) array!\n");
                fflush(stderr);
                exit(0);
              }
              (*pyr_conn)[npyr+npyr_zone] = (int*)malloc(5*sizeof(int));
              (*pyr_conn)[npyr+npyr_zone][0] = conn[j++]-1 + nn;
              (*pyr_conn)[npyr+npyr_zone][1] = conn[j++]-1 + nn;
              (*pyr_conn)[npyr+npyr_zone][2] = conn[j++]-1 + nn;
              (*pyr_conn)[npyr+npyr_zone][3] = conn[j++]-1 + nn;
              (*pyr_conn)[npyr+npyr_zone][4] = conn[j++]-1 + nn;
              npyr_zone++;
              break;
            case PENTA_6:
              (*pri_conn) = (int**)realloc((void*)(*pri_conn),(npri+npri_zone+1)*sizeof(int*));
              if ((*pri_conn) == NULL)
              {
                fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pri_conn) array!\n");
                fflush(stderr);
                exit(0);
              }
              (*pri_conn)[npri+npri_zone] = (int*)malloc(6*sizeof(int));
              (*pri_conn)[npri+npri_zone][0] = conn[j++]-1 + nn;
              (*pri_conn)[npri+npri_zone][1] = conn[j++]-1 + nn;
              (*pri_conn)[npri+npri_zone][2] = conn[j++]-1 + nn;
              (*pri_conn)[npri+npri_zone][3] = conn[j++]-1 + nn;
              (*pri_conn)[npri+npri_zone][4] = conn[j++]-1 + nn;
              (*pri_conn)[npri+npri_zone][5] = conn[j++]-1 + nn;
              npri_zone++;
              break;
            case HEXA_8:
              (*hex_conn) = (int**)realloc((void*)(*hex_conn),(nhex+nhex_zone+1)*sizeof(int*));
              if ((*hex_conn) == NULL)
              {
                fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*hex_conn) array!\n");
                fflush(stderr);
                exit(0);
              }
              (*hex_conn)[nhex+nhex_zone] = (int*)malloc(8*sizeof(int));
              (*hex_conn)[nhex+nhex_zone][0] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][1] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][2] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][3] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][4] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][5] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][6] = conn[j++]-1 + nn;
              (*hex_conn)[nhex+nhex_zone][7] = conn[j++]-1 + nn;
              nhex_zone++;
              break;
            default:
              fprintf(stderr,"\nElement type %d not supported!",conn[j]);
              break;
          }
        }
      }

      if (etype == TRI_3)
      {
        conn = (int*)realloc((void*)conn,3*(nend-nstart+1)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);

        if (b < 0)
        {
          fprintf(stderr,"\nTri elements and no boundary identified!");
          fprintf(stderr,"\nAssuming these are zone connectivity elements!");
          for (j=0,i=0; i < nend-nstart+1; i++)
          {
            zone_nodes->Check_List(conn[j++]-1+nn);
            zone_nodes->Check_List(conn[j++]-1+nn);
            zone_nodes->Check_List(conn[j++]-1+nn);
          }
        } else
        {
          (*tri_conn)[b] = (int**)realloc((void*)(*tri_conn)[b],((*nt)[b]+nend-nstart+1)*sizeof(int*));
          if ((*tri_conn)[b] == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tri_conn)[%d] array!\n",b);
            fflush(stderr);
            exit(0);
          }
          for (j=0,i=(*nt)[b]; i < (*nt)[b]+nend-nstart+1; i++)
          {
            (*tri_conn)[b][i] = (int*)malloc(3*sizeof(int));
            (*tri_conn)[b][i][0] = conn[j++]-1 + nn;
            (*tri_conn)[b][i][1] = conn[j++]-1 + nn;
            (*tri_conn)[b][i][2] = conn[j++]-1 + nn;
          }
          (*nt)[b] += nend-nstart+1;
        }
      }

      if (etype == QUAD_4)
      {
        conn = (int*)realloc((void*)conn,4*(nend-nstart+1)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);

        if (b < 0)
        {
          fprintf(stderr,"\nQuad elements and no boundary identified!");
          fprintf(stderr,"\nAssuming these are zone connectivity elements!");
          for (j=0,i=0; i < nend-nstart+1; i++)
          {
            zone_nodes->Check_List(conn[j++]-1+nn);
            zone_nodes->Check_List(conn[j++]-1+nn);
            zone_nodes->Check_List(conn[j++]-1+nn);
            zone_nodes->Check_List(conn[j++]-1+nn);
          }
        } else
        {
          (*quad_conn)[b] = (int**)realloc((void*)(*quad_conn)[b],((*nq)[b]+nend-nstart+1)*sizeof(int*));
          if ((*quad_conn)[b] == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*quad_conn)[%d] array!\n",b);
            fflush(stderr);
            exit(0);
          }
          for (j=0,i=(*nq)[b]; i < (*nq)[b]+nend-nstart+1; i++)
          {
            (*quad_conn)[b][i] = (int*)malloc(4*sizeof(int));
            (*quad_conn)[b][i][0] = conn[j++]-1 + nn;
            (*quad_conn)[b][i][1] = conn[j++]-1 + nn;
            (*quad_conn)[b][i][2] = conn[j++]-1 + nn;
            (*quad_conn)[b][i][3] = conn[j++]-1 + nn;
          }
          (*nq)[b] += nend-nstart+1;
        }
      }

      if (etype == TETRA_4)
      {
        (*tet_conn) = (int**)realloc((void*)(*tet_conn),(ntet+ntet_zone+nend-nstart+1)*sizeof(int*));
        if ((*tet_conn) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tet_conn) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,((nend-nstart+1)*4)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);
        for (j=0,c=ntet_zone; c < ntet_zone+nend-nstart+1; c++)
        {
          (*tet_conn)[c+ntet] = (int*)malloc(4*sizeof(int));
          (*tet_conn)[c+ntet][0] = conn[j++]-1 + nn;
          (*tet_conn)[c+ntet][1] = conn[j++]-1 + nn;
          (*tet_conn)[c+ntet][2] = conn[j++]-1 + nn;
          (*tet_conn)[c+ntet][3] = conn[j++]-1 + nn;
        }
        ntet_zone += (nend-nstart+1);
      }
    
      if (etype == PYRA_5)
      {
        (*pyr_conn) = (int**)realloc((void*)(*pyr_conn),(npyr+npyr_zone+nend-nstart+1)*sizeof(int*));
        if ((*pyr_conn) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pyr_conn) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,((nend-nstart+1)*5)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);
        for (j=0,c=npyr_zone; c < npyr_zone+nend-nstart+1; c++)
        {
          (*pyr_conn)[c+npyr] = (int*)malloc(5*sizeof(int));
          (*pyr_conn)[c+npyr][0] = conn[j++]-1 + nn;
          (*pyr_conn)[c+npyr][1] = conn[j++]-1 + nn;
          (*pyr_conn)[c+npyr][2] = conn[j++]-1 + nn;
          (*pyr_conn)[c+npyr][3] = conn[j++]-1 + nn;
          (*pyr_conn)[c+npyr][4] = conn[j++]-1 + nn;
        }
        npyr_zone += (nend-nstart+1);
      }
      
      if (etype == PENTA_6)
      {
        (*pri_conn) = (int**)realloc((void*)(*pri_conn),(npri+npri_zone+nend-nstart+1)*sizeof(int*));
        if ((*pri_conn) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pri_conn) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,((nend-nstart+1)*6)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);
        for (j=0,c=npri_zone; c < npri_zone+nend-nstart+1; c++)
        {
          (*pri_conn)[c+npri] = (int*)malloc(6*sizeof(int));
          (*pri_conn)[c+npri][0] = conn[j++]-1 + nn;
          (*pri_conn)[c+npri][1] = conn[j++]-1 + nn;
          (*pri_conn)[c+npri][2] = conn[j++]-1 + nn;
          (*pri_conn)[c+npri][3] = conn[j++]-1 + nn;
          (*pri_conn)[c+npri][4] = conn[j++]-1 + nn;
          (*pri_conn)[c+npri][5] = conn[j++]-1 + nn;
        }
        npri_zone += (nend-nstart+1);
      }
    
      if (etype == HEXA_8)
      {
        (*hex_conn) = (int**)realloc((void*)(*hex_conn),(nhex+nhex_zone+nend-nstart+1)*sizeof(int*));
        if ((*hex_conn) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*hex_conn) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,((nend-nstart+1)*8)*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        ier = cg_elements_read(index_file,index_base,index_zone,s,conn,&pflag);
        CGNS_Errorcheck(ier);
        for (j=0,c=nhex_zone; c < nhex_zone+nend-nstart+1; c++)
        {
          (*hex_conn)[c+nhex] = (int*)malloc(8*sizeof(int));
          (*hex_conn)[c+nhex][0] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][1] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][2] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][3] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][4] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][5] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][6] = conn[j++]-1 + nn;
          (*hex_conn)[c+nhex][7] = conn[j++]-1 + nn;
        }
        nhex_zone += (nend-nstart+1);
      }
    }
    nn += nn_zone;
    ntet += ntet_zone;
    npyr += npyr_zone;
    npri += npri_zone;
    nhex += nhex_zone;
    nb += nb_zone;
    delete[] bc_range;
  }

  if ((ntet + npyr + npri + nhex) > 0)
  {
    // if multiple zones were read then eliminate duplicate nodes and grid-connectivity elements
    if (nzones > 1)
    {
      int *map = new int[nn];
      for (n=0; n < nn; n++)
        map[n] = n;

      // compute a tolerance
      double tol = 1e20;

      for (c=0; c < ntet; c++)
        for (i=0; i < 3; i++)
        {
          n = (*tet_conn)[c][i];
          for (j=i+1; j < 4; j++)
          {
            m = (*tet_conn)[c][j];
            tol = MIN(tol,dist((*p)[(*tet_conn)[c][i]],(*p)[(*tet_conn)[c][j]]));
          }
        }
      for (c=0; c < npyr; c++)
        for (i=0; i < 4; i++)
          for (j=i+1; j < 5; j++)
            tol = MIN(tol,dist((*p)[(*pyr_conn)[c][i]],(*p)[(*pyr_conn)[c][j]]));
      for (c=0; c < npri; c++)
        for (i=0; i < 5; i++)
          for (j=i+1; j < 6; j++)
            tol = MIN(tol,dist((*p)[(*pri_conn)[c][i]],(*p)[(*pri_conn)[c][j]]));
      for (c=0; c < nhex; c++)
        for (i=0; i < 7; i++)
          for (j=i+1; j < 8; j++)
            tol = MIN(tol,dist((*p)[(*hex_conn)[c][i]],(*p)[(*hex_conn)[c][j]]));
      tol *= 0.001;
      if (tol < 1e-15)
      {
        fprintf(stderr,"Computed tolerance for node check is too small, tol = %lg",tol);
        exit(0);
      }

      double p1[3], p2[3];
      for (i=0; i < zone_nodes->max-1; i++)
      {
        n = zone_nodes->list[i];
        p1[0] = (*p)[n][0];
        p1[1] = (*p)[n][1];
        p1[2] = (*p)[n][2];
        if (map[n] != n) continue;
        for (j=i+1; j < zone_nodes->max; j++)
        {
          m = zone_nodes->list[j];
          if (map[m] != m) continue;
          p2[0] = (*p)[m][0];
          p2[1] = (*p)[m][1];
          p2[2] = (*p)[m][2];
          if (dist(p1,p2) < tol)
            map[m] = n;
        }
      }

      // apply current node map
      for (c=0; c < ntet; c++)
        for (i=0; i < 4; i++)
          (*tet_conn)[c][i] = map[(*tet_conn)[c][i]];
      for (c=0; c < npyr; c++)
        for (i=0; i < 5; i++)
          (*pyr_conn)[c][i] = map[(*pyr_conn)[c][i]];
      for (c=0; c < npri; c++)
        for (i=0; i < 6; i++)
          (*pri_conn)[c][i] = map[(*pri_conn)[c][i]];
      for (c=0; c < nhex; c++)
        for (i=0; i < 8; i++)
          (*hex_conn)[c][i] = map[(*hex_conn)[c][i]];
      for (b=0; b < nb; b++)
      {
        for (c=0; c < (*nt)[b]; c++)
          for (i=0; i < 3; i++)
            (*tri_conn)[b][c][i] = map[(*tri_conn)[b][c][i]];
        for (c=0; c < (*nq)[b]; c++)
          for (i=0; i < 4; i++)
            (*quad_conn)[b][c][i] = map[(*quad_conn)[b][c][i]];
      }

      // now shift nodes down to eliminate duplicates
      m = 0;
      for (n=0; n < nn; n++)
        if (map[n] == n)
        {
          (*p)[m] = (*p)[n];
          map[n] = m++;
        } else
          map[n] = -1;

      // apply current node map
      for (c=0; c < ntet; c++)
        for (i=0; i < 4; i++)
          (*tet_conn)[c][i] = map[(*tet_conn)[c][i]];
      for (c=0; c < npyr; c++)
        for (i=0; i < 5; i++)
          (*pyr_conn)[c][i] = map[(*pyr_conn)[c][i]];
      for (c=0; c < npri; c++)
        for (i=0; i < 6; i++)
          (*pri_conn)[c][i] = map[(*pri_conn)[c][i]];
      for (c=0; c < nhex; c++)
        for (i=0; i < 8; i++)
          (*hex_conn)[c][i] = map[(*hex_conn)[c][i]];
      for (b=0; b < nb; b++)
      {
        for (c=0; c < (*nt)[b]; c++)
          for (i=0; i < 3; i++)
            (*tri_conn)[b][c][i] = map[(*tri_conn)[b][c][i]];
        for (c=0; c < (*nq)[b]; c++)
          for (i=0; i < 4; i++)
            (*quad_conn)[b][c][i] = map[(*quad_conn)[b][c][i]];
      }
      for (n=m; n < nn; n++)
        free((*p)[n]);
      nn = m;
      (*p) = (double **)realloc((void*)(*p),nn*sizeof(double *));

      delete[] map;
      delete zone_nodes;
    }
    
    // check for duplicate boundary names
    // if found transfer elements to lower boundary and delete higher boundary
    for (b=0; b < nb-1; b++)
    {
      c=b+1;
      while (c < nb)
      {
        if (strcmp((*b_name)[b],(*b_name)[c]) == 0)
        {
          fprintf(stderr,"\nDuplicate boundary names found for boundaries %d & %d",b+1,c+1);
          fprintf(stderr,"\nBoundary %d is %s",b+1,(*b_name)[b]);
          fprintf(stderr,"\nBoundary %d is %s",c+1,(*b_name)[c]);
          fprintf(stderr,"\nTransferring elements from boundary %d to boundary %d",c+1,b+1);
          fflush(stderr);

          // transfer boundary c to b
          if ((*nt)[c] > 0)
          {
            (*tri_conn)[b] = (int**)realloc((void*)(*tri_conn)[b],((*nt)[b]+(*nt)[c])*sizeof(int*));
            for (j=0,i=(*nt)[b]; i < (*nt)[b]+(*nt)[c]; i++, j++)
            {
              (*tri_conn)[b][i] = (int*)malloc(3*sizeof(int));
              (*tri_conn)[b][i][0] = (*tri_conn)[c][j][0];
              (*tri_conn)[b][i][1] = (*tri_conn)[c][j][1];
              (*tri_conn)[b][i][2] = (*tri_conn)[c][j][2];
            }
            (*nt)[b] += (*nt)[c];
          }
          if ((*nq)[c] > 0)
          {
            (*quad_conn)[b] = (int**)realloc((void*)(*quad_conn)[b],((*nq)[b]+(*nq)[c])*sizeof(int*));
            for (j=0,i=(*nq)[b]; i < (*nq)[b]+(*nq)[c]; i++, j++)
            {
              (*quad_conn)[b][i] = (int*)malloc(4*sizeof(int));
              (*quad_conn)[b][i][0] = (*quad_conn)[c][j][0];
              (*quad_conn)[b][i][1] = (*quad_conn)[c][j][1];
              (*quad_conn)[b][i][2] = (*quad_conn)[c][j][2];
              (*quad_conn)[b][i][3] = (*quad_conn)[c][j][3];
            }
            (*nq)[b] += (*nq)[c];
          }
          
          // shift remaining boundaries up
          for (m=c; m < nb-1; m++)
          {
            fprintf(stderr,"\nShifting boundary %d to boundary %d",m+2,m+1);
            fflush(stderr);
            (*b_name)[m][0] = '\0';
            strcpy((*b_name)[m],(*b_name)[m+1]);
            for (i=(*nt)[m]; i < (*nt)[m]; i++)
              free((*tri_conn)[m][i]);
            for (i=(*nq)[m]; i < (*nq)[m]; i++)
              free((*quad_conn)[m][i]);
            if ((*nt)[m+1] > 0)
            {
              (*tri_conn)[m] = (int**)realloc((void*)(*tri_conn)[m],((*nt)[m+1])*sizeof(int*));
              for (i=0; i < (*nt)[m+1]; i++)
              {
                (*tri_conn)[m][i] = (int*)malloc(3*sizeof(int));
                (*tri_conn)[m][i][0] = (*tri_conn)[m+1][i][0];
                (*tri_conn)[m][i][1] = (*tri_conn)[m+1][i][1];
                (*tri_conn)[m][i][2] = (*tri_conn)[m+1][i][2];
              }
              (*nt)[m] = (*nt)[m+1];
            } else if ((*nt)[m] > 0)
            {
              free((*tri_conn)[m]);
              (*tri_conn)[m] = 0;
            }
            if ((*nq)[m+1] > 0)
            {
              (*quad_conn)[m] = (int**)realloc((void*)(*quad_conn)[m],((*nq)[m+1])*sizeof(int*));
              for (i=0; i < (*nq)[m+1]; i++)
              {
                (*quad_conn)[m][i] = (int*)malloc(4*sizeof(int));
                (*quad_conn)[m][i][0] = (*quad_conn)[m+1][i][0];
                (*quad_conn)[m][i][1] = (*quad_conn)[m+1][i][1];
                (*quad_conn)[m][i][2] = (*quad_conn)[m+1][i][2];
                (*quad_conn)[m][i][3] = (*quad_conn)[m+1][i][3];
              }
              (*nq)[m] = (*nq)[m+1];
            } else if ((*nq)[m] > 0)
            {
              free((*quad_conn)[m]);
              (*quad_conn)[m] = 0;
            }
          }
          // delete last boundary
          m = nb-1;
          for (i=(*nt)[m]; i < (*nt)[m]; i++)
            free((*tri_conn)[m][i]);
          for (i=(*nq)[m]; i < (*nq)[m]; i++)
            free((*quad_conn)[m][i]);
          if ((*nt)[m] > 0)
             free((*tri_conn)[m]);
          if ((*nq)[m] > 0)
             free((*quad_conn)[m]);
          free((*b_name)[m]);
          nb--;
          (*b_name) = (char**)realloc((void*)(*b_name),nb*sizeof(char*));
          (*nt) = (int*)realloc((void*)(*nt),nb*sizeof(int));
          (*tri_conn) = (int***)realloc((void*)(*tri_conn),nb*sizeof(int**));
          (*nq) = (int*)realloc((void*)(*nq),nb*sizeof(int));
          (*quad_conn) = (int***)realloc((void*)(*quad_conn),nb*sizeof(int**));
        } else
          c++;
      }
    }

    if (parallel)
    {
      (*nmap) = (*tri_map) = (*quad_map) = 0;
      (*tet_map) = (*pyr_map) = (*pri_map) = (*hex_map) = 0;

      int nud = 0;
      // Read node map triplet
      nud++;
      ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
      CGNS_Errorcheck(ier);
      cg_narrays(&narrays);
      if (narrays != 3)
      {
        fprintf(stderr,"\nNumber of parallel node map arrays = %d",narrays);
        fprintf(stderr,"\nExpecting three!");
        fflush(stderr);
        cg_error_exit();
      }
      (*nmap) = (int**)malloc(nn*sizeof(int*));
      if ((*nmap) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *nmap array!\n");
        fflush(stderr);
        exit(0);
      }
      for (n=0; n < nn; n++)
        (*nmap)[n] = (int*)malloc(3*sizeof(int));
      conn = (int*)realloc((void*)conn,nn*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      for (i=1; i <= narrays; i++)
      {
        cg_array_info(i,array_name,&itype,&ndim,&idim);
        if (idim != nn)
        {
          fprintf(stderr,"\nP_CGNS_read: nmap dimension incorrect!");
          fprintf(stderr,"\n           expected dimension = %d",nn);
          fprintf(stderr,"\n               file dimension = %d",idim);
          fflush(stderr);
          cg_error_exit();
        }
        switch(i)
        {
          case 1:
            if (strncmp(array_name,"Global_Node_Number",18) != 0)
            {
              fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
              fflush(stderr);
              cg_error_exit();
            }
            break;
          case 2:
            if (strncmp(array_name,"Node_Processor_Number",21) != 0)
            {
              fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
              fflush(stderr);
              cg_error_exit();
            }
            break;
          case 3:
            if (strncmp(array_name,"Node_Index_Number",17) != 0)
            {
              fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
              fflush(stderr);
              cg_error_exit();
            }
            break;
        }

        ier = cg_array_read(i,conn);
        CGNS_Errorcheck(ier);
        for (n=0; n < nn; n++)
          (*nmap)[n][i-1] = conn[n]-1; // decrement for C indexing
      }

      // Read boundary element global maps
      nud++;
      ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
      CGNS_Errorcheck(ier);
      cg_narrays(&narrays);
      (*tri_map) = (int**)malloc(nb*sizeof(int*));
      if ((*tri_map) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *tri_map array!\n");
        fflush(stderr);
        exit(0);
      }
      (*quad_map) = (int**)malloc(nb*sizeof(int*));
      if ((*quad_map) == NULL)
      {
        fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for *quad_map array!\n");
        fflush(stderr);
        exit(0);
      }
      for (j=0; j < nb; j++)
      {
        (*tri_map)[j] = 0;
        (*quad_map)[j] = 0;
      }
      for (j=1; j <= narrays; j++)
      {
        cg_array_info(j,array_name,&itype,&ndim,&idim);
        if (strstr(array_name,"Global_Tri_Number") != NULL)
        {
          sscanf(array_name,"Boundary_%d_Global_Tri_Number",&b);
          b--;
          if (b < 0 || b >= nb)
          {
            fprintf(stderr,"\nP_CGNS_read: boundary number for tri map out of bounds!");
            fprintf(stderr,"\n     input value = %d. Maximum boundary = %d",b+1,nb);
            fflush(stderr);
            cg_error_exit();
          }
          if (idim != (*nt)[b])
          {
            fprintf(stderr,"\nP_CGNS_read: tri_map dimension incorrect!");
            fprintf(stderr,"\n           expected dimension = %d",(*nt)[b]);
            fprintf(stderr,"\n               file dimension = %d",idim);
            fflush(stderr);
            cg_error_exit();
          }
          (*tri_map)[b] = (int*)malloc((*nt)[b]*sizeof(int));
          if ((*tri_map)[b] == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tri_map)[%d] array!\n",b);
            fflush(stderr);
            exit(0);
          }
          conn = (int*)realloc((void*)conn,(*nt)[b]*sizeof(int));
          if (conn == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
            fflush(stderr);
            exit(0);
          }
          cg_array_read(j,conn);
          for (c=0; c < (*nt)[b]; c++)
            (*tri_map)[b][c] = conn[c]-1; // decrement for C indexing
        }
        if (strstr(array_name,"Global_Quad_Number") != NULL)
        {
          sscanf(array_name,"Boundary_%d_Global_Quad_Number",&b);
          b--;
          if (b < 0 || b >= nb)
          {
            fprintf(stderr,"\nP_CGNS_read: boundary number for quad map out of bounds!");
            fprintf(stderr,"\n     input value = %d. Maximum boundary = %d",b+1,nb);
            fflush(stderr);
            cg_error_exit();
          }
          if (idim != (*nq)[b])
          {
            fprintf(stderr,"\nP_CGNS_read: quad_map dimension incorrect!");
            fprintf(stderr,"\n           expected dimension = %d",(*nq)[b]);
            fprintf(stderr,"\n               file dimension = %d",idim);
            fflush(stderr);
            cg_error_exit();
          }
          (*quad_map)[b] = (int*)malloc((*nq)[b]*sizeof(int));
          if ((*quad_map)[b] == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*quad_map)[%d] array!\n",b);
            fflush(stderr);
            exit(0);
          }
          conn = (int*)realloc((void*)conn,(*nq)[b]*sizeof(int));
          if (conn == NULL)
          {
            fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
            fflush(stderr);
            exit(0);
          }
          cg_array_read(j,conn);
          for (c=0; c < (*nq)[b]; c++)
            (*quad_map)[b][c] = conn[c]-1; // decrement for C indexing
        }
      }
        
      if (ntet > 0)
      {
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        cg_narrays(&narrays);
        if (narrays != 1)
        {
          fprintf(stderr,"\nNumber of parallel tet map arrays = %d",narrays);
          fprintf(stderr,"\nExpecting 1!");
          fflush(stderr);
          cg_error_exit();
        }
        cg_array_info(narrays,array_name,&itype,&ndim,&idim);
        if (idim != ntet)
        {
          fprintf(stderr,"\nP_CGNS_read: tet_map dimension incorrect!");
          fprintf(stderr,"\n           expected dimension = %d",ntet);
          fprintf(stderr,"\n               file dimension = %d",idim);
          fflush(stderr);
          cg_error_exit();
        }
        if (strcmp(array_name,"Global_Tet_Number") != 0)
        {
          fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
          fflush(stderr);
          cg_error_exit();
        }
        (*tet_map) = (int*)malloc(ntet*sizeof(int));
        if ((*tet_map) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*tet_map) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,ntet*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_array_read(1,conn);
        for (c=0; c < ntet; c++)
          (*tet_map)[c] = conn[c]-1; // decrement for C indexing
      }

      if (npyr > 0)
      {
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        cg_narrays(&narrays);
        if (narrays != 1)
        {
          fprintf(stderr,"\nNumber of parallel Pyramid map arrays = %d",narrays);
          fprintf(stderr,"\nExpecting 1!");
          fflush(stderr);
          cg_error_exit();
        }
        cg_array_info(narrays,array_name,&itype,&ndim,&idim);
        if (idim != npyr)
        {
          fprintf(stderr,"\nP_CGNS_read: pyr_map dimension incorrect!");
          fprintf(stderr,"\n           expected dimension = %d",npyr);
          fprintf(stderr,"\n               file dimension = %d",idim);
          fflush(stderr);
          cg_error_exit();
        }
        if (strcmp(array_name,"Global_Pyramid_Number") != 0)
        {
          fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
          fflush(stderr);
          cg_error_exit();
        }
        (*pyr_map) = (int*)malloc(npyr*sizeof(int));
        if ((*pyr_map) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pyr_map) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,npyr*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_array_read(1,conn);
        for (c=0; c < npyr; c++)
          (*pyr_map)[c] = conn[c]-1; // decrement for C indexing
      }

      if (npri > 0)
      {
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        cg_narrays(&narrays);
        if (narrays != 1)
        {
          fprintf(stderr,"\nNumber of parallel Prism map arrays = %d",narrays);
          fprintf(stderr,"\nExpecting 1!");
          fflush(stderr);
          cg_error_exit();
        }
        cg_array_info(narrays,array_name,&itype,&ndim,&idim);
        if (idim != npri)
        {
          fprintf(stderr,"\nP_CGNS_read: pri_map dimension incorrect!");
          fprintf(stderr,"\n           expected dimension = %d",npri);
          fprintf(stderr,"\n               file dimension = %d",idim);
          fflush(stderr);
          cg_error_exit();
        }
        if (strcmp(array_name,"Global_Prism_Number") != 0)
        {
          fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
          fflush(stderr);
          cg_error_exit();
        }
        (*pri_map) = (int*)malloc(npri*sizeof(int));
        if ((*pri_map) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*pri_map) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,npri*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_array_read(1,conn);
        for (c=0; c < npri; c++)
          (*pri_map)[c] = conn[c]-1; // decrement for C indexing
      }

      if (nhex > 0)
      {
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        cg_narrays(&narrays);
        if (narrays != 1)
        {
          fprintf(stderr,"\nNumber of parallel Hex map arrays = %d",narrays);
          fprintf(stderr,"\nExpecting 1!");
          fflush(stderr);
          cg_error_exit();
        }
        cg_array_info(narrays,array_name,&itype,&ndim,&idim);
        if (idim != nhex)
        {
          fprintf(stderr,"\nP_CGNS_read: hex_map dimension incorrect!");
          fprintf(stderr,"\n           expected dimension = %d",nhex);
          fprintf(stderr,"\n               file dimension = %d",idim);
          fflush(stderr);
          cg_error_exit();
        }
        if (strcmp(array_name,"Global_Hex_Number") != 0)
        {
          fprintf(stderr,"\nP_CGNS_read: Incorrect array name!");
          fflush(stderr);
          cg_error_exit();
        }
        (*hex_map) = (int*)malloc(nhex*sizeof(int));
        if ((*hex_map) == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for (*hex_map) array!\n");
          fflush(stderr);
          exit(0);
        }
        conn = (int*)realloc((void*)conn,nhex*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_read: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_array_read(1,conn);
        for (c=0; c < nhex; c++)
          (*hex_map)[c] = conn[c]-1; // decrement for C indexing
      }
    }
  }

  ier = cg_close(index_file);
  CGNS_Errorcheck(ier);

  if (conn)
    free(conn);

  return(flag);

}

// Parallel CGNS write function
int P_CGNS_write(int parallel, char *fname, int nn, double **p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn,
               int **nmap, int **tri_map, int **quad_map,
               int *tet_map, int *pyr_map, int *pri_map, int *hex_map)
{
  int b, c, i, j, n, nstart, nend;
  int nc, nf, ier;
  int index_file, index_base, index_zone, index_coord, index_section;
  int icelldim, iphysdim;
  int isize[3], *bc_section[2];
  int *conn;
  double *x;
  char sectionname[74];
  char descriptorname[74];

  conn = 0;

  ier = cg_open(fname,CG_MODE_WRITE,&index_file);
  CGNS_Errorcheck(ier);

  if (nb > 0)
  {
    bc_section[0] = (int*)malloc((nb+1)*sizeof(int));
    if (bc_section[0] == NULL)
    {
      fprintf(stderr,"\nP_CGNS_write: Unable to allocate memory for bc_section[0] array!\n");
      fflush(stderr);
      exit(0);
    }
    bc_section[1] = (int*)malloc((nb+1)*sizeof(int));
    if (bc_section[1] == NULL)
    {
      fprintf(stderr,"\nP_CGNS_write: Unable to allocate memory for bc_section[1] array!\n");
      fflush(stderr);
      exit(0);
    }
  }

  char basename[]="Base";
  icelldim=3;
  iphysdim=3;
  ier = cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
  CGNS_Errorcheck(ier);
  char zonename[] = "Zone 1";

  time_t tm;
  char *t_char;

  time(&tm);
  t_char = ctime(&tm);
 
  cg_goto(index_file,index_base,"end");
  sprintf(descriptorname,"Exported from P_CGNS\n%s",t_char);
  cg_descriptor_write("Information",descriptorname);

  nc = ntet + npyr + npri + nhex;
  for (nf=n=0; n < nb; n++)
    nf += nt[n] + nq[n];

  isize[0]=nn;
  isize[1]=nc;
  isize[2]=0;
  ier=cg_zone_write(index_file,index_base,zonename,isize,Unstructured,&index_zone);
  CGNS_Errorcheck(ier);

  if (nc > 0)
  {
    x = (double*)malloc(nn*sizeof(double));
    if (x == NULL)
    {
      fprintf(stderr,"\nP_CGNS_write: Unable to allocate memory for x array!\n");
      fflush(stderr);
      exit(0);
    }
    for (n=0; n < nn; n++)
      x[n] = p[n][0];
    index_coord = 1;
    ier=cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateX",x,
                       &index_coord);
    CGNS_Errorcheck(ier);
    for (n=0; n < nn; n++)
      x[n] = p[n][1];
    index_coord = 2;
    ier=cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateY",x,
                       &index_coord);
    CGNS_Errorcheck(ier);
    for (n=0; n < nn; n++)
      x[n] = p[n][2];
    index_coord = 3;
    ier=cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateZ",x,
                       &index_coord);
    CGNS_Errorcheck(ier);
    free(x);

    int (*bc_range)[2];
    if (nb > 0)
      bc_range = new int[nb][2];

    nc = 0;
    if (ntet > 0)
    {
      i = 4*ntet;
      conn = (int*)realloc((void*)conn,i*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      nstart = nc+1;
      for (j=c=0; c < ntet; c++)
      {
        nc++;
        conn[j++] = tet_conn[c][0]+1;
        conn[j++] = tet_conn[c][1]+1;
        conn[j++] = tet_conn[c][2]+1;
        conn[j++] = tet_conn[c][3]+1;
      }
      nend = nc;
      ier = cg_section_write(index_file,index_base,index_zone,"Tetra_Elements",
                             TETRA_4,nstart,nend,0,conn,&index_section);
      CGNS_Errorcheck(ier);
    }
    if (npyr > 0)
    {
      i = 5*npyr;
      conn = (int*)realloc((void*)conn,i*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      nstart = nc+1;
      for (j=c=0; c < npyr; c++)
      {
        nc++;
        conn[j++] = pyr_conn[c][0]+1;
        conn[j++] = pyr_conn[c][1]+1;
        conn[j++] = pyr_conn[c][2]+1;
        conn[j++] = pyr_conn[c][3]+1;
        conn[j++] = pyr_conn[c][4]+1;
      }
      nend = nc;
      ier = cg_section_write(index_file,index_base,index_zone,"Pyramid_Elements",
                             PYRA_5,nstart,nend,0,conn,&index_section);
      CGNS_Errorcheck(ier);
    }
    if (npri > 0)
    {
      i = 6*npri;
      conn = (int*)realloc((void*)conn,i*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      nstart = nc+1;
      for (j=c=0; c < npri; c++)
      {
        nc++;
        conn[j++] = pri_conn[c][0]+1;
        conn[j++] = pri_conn[c][1]+1;
        conn[j++] = pri_conn[c][2]+1;
        conn[j++] = pri_conn[c][3]+1;
        conn[j++] = pri_conn[c][4]+1;
        conn[j++] = pri_conn[c][5]+1;
      }
      nend = nc;
      ier = cg_section_write(index_file,index_base,index_zone,"Prism_Elements",
                             PENTA_6,nstart,nend,0,conn,&index_section);
      CGNS_Errorcheck(ier);
    }
    if (nhex > 0)
    {
      i = 8*nhex;
      conn = (int*)realloc((void*)conn,i*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      nstart = nc+1;
      for (j=c=0; c < nhex; c++)
      {
        nc++;
        conn[j++] = hex_conn[c][0]+1;
        conn[j++] = hex_conn[c][1]+1;
        conn[j++] = hex_conn[c][2]+1;
        conn[j++] = hex_conn[c][3]+1;
        conn[j++] = hex_conn[c][4]+1;
        conn[j++] = hex_conn[c][5]+1;
        conn[j++] = hex_conn[c][6]+1;
        conn[j++] = hex_conn[c][7]+1;
      }
      nend = nc;
      ier = cg_section_write(index_file,index_base,index_zone,"Hex_Elements",
                             HEXA_8,nstart,nend,0,conn,&index_section);
      CGNS_Errorcheck(ier);
    }
    
    // Boundary faces
    nf = nc;
    for (n=0; n < nb; n++)
    {
      bc_range[n][0] = nf;
      bc_range[n][1] = nf-1;
      sprintf(sectionname,"%s",b_name[n]);
      if (nt[n] > 0 && nq[n] > 0)
      {
        i = 4*nt[n]+5*nq[n];
        conn = (int*)realloc((void*)conn,i*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        nstart = nf+1;
        i=0;
        for (c=0; c < nt[n]; c++)
        {
          nf++;
          conn[i++]=TRI_3;
          conn[i++]=tri_conn[n][c][0]+1;
          conn[i++]=tri_conn[n][c][1]+1;
          conn[i++]=tri_conn[n][c][2]+1;
        }
        for (c=0; c < nq[n]; c++)
        {
          nf++;
          conn[i++]=QUAD_4;
          conn[i++]=quad_conn[n][c][0]+1;
          conn[i++]=quad_conn[n][c][1]+1;
          conn[i++]=quad_conn[n][c][2]+1;
          conn[i++]=quad_conn[n][c][3]+1;
        }
        nend = nf;
        ier = cg_section_write(index_file,index_base,index_zone,sectionname,
                         MIXED,nstart,nend,0,conn,&bc_section[0][n]);
        CGNS_Errorcheck(ier);
        bc_range[n][0] = nstart;
        bc_range[n][1] = nend;
      } else if (nt[n] > 0)
      {
        i = 3*nt[n];
        conn = (int*)realloc((void*)conn,i*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        nstart = nf+1;
        for (i=c=0; c < nt[n]; c++)
        {
          nf++;
          conn[i++]=tri_conn[n][c][0]+1;
          conn[i++]=tri_conn[n][c][1]+1;
          conn[i++]=tri_conn[n][c][2]+1;
        }
        nend = nf;
        ier = cg_section_write(index_file,index_base,index_zone,sectionname,
                         TRI_3,nstart,nend,0,conn,&bc_section[0][n]);
        CGNS_Errorcheck(ier);
        bc_range[n][0] = nstart;
        bc_range[n][1] = nend;
      } else if (nq[n] > 0)
      {
        i = 4*nq[n];
        conn = (int*)realloc((void*)conn,i*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        nstart = nf+1;
        for (i=c=0; c < nq[n]; c++)
        {
          nf++;
          conn[i++]=quad_conn[n][c][0]+1;
          conn[i++]=quad_conn[n][c][1]+1;
          conn[i++]=quad_conn[n][c][2]+1;
          conn[i++]=quad_conn[n][c][3]+1;
        }
        nend = nf;
        ier = cg_section_write(index_file,index_base,index_zone,sectionname,
                         QUAD_4,nstart,nend,0,conn,&bc_section[1][n]);
        CGNS_Errorcheck(ier);
        bc_range[n][0] = nstart;
        bc_range[n][1] = nend;
      }
    }

    for (n=0; n < nb; n++)
    {
      sprintf(sectionname,"%s",b_name[n]);
      ier = cg_boco_write(index_file,index_base,index_zone,sectionname,(BCType_t)CG_Null,ElementRange,2,bc_range[n],&bc_section[0][n]);
      CGNS_Errorcheck(ier);
    }

    if (nb > 0)
    {
      delete[] bc_range;

      free(bc_section[0]);
      free(bc_section[1]);
    }

    if (parallel) // Write node and element maps to global numbers
    {
      int nud = 0;
      conn = (int*)realloc((void*)conn,nn*sizeof(int));
      if (conn == NULL)
      {
        fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
        fflush(stderr);
        exit(0);
      }
      ier = cg_goto(index_file,index_base,"Zone_t",1,"end");
      CGNS_Errorcheck(ier);
      ier = cg_user_data_write("Parallel_Node_Map");
      CGNS_Errorcheck(ier);
      nud++;
      ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
      CGNS_Errorcheck(ier);
      for (c=0; c < nn; c++)
        conn[c] = nmap[c][0]+1;
      ier = cg_array_write("Global_Node_Number",Integer,1,&nn,conn);
      CGNS_Errorcheck(ier);
      for (c=0; c < nn; c++)
        conn[c] = nmap[c][1]+1;
      ier = cg_array_write("Node_Processor_Number",Integer,1,&nn,conn);
      CGNS_Errorcheck(ier);
      for (c=0; c < nn; c++)
        conn[c] = nmap[c][2]+1;
      ier = cg_array_write("Node_Index_Number",Integer,1,&nn,conn);
      CGNS_Errorcheck(ier);

      cg_goto(index_file,index_base,"Zone_t",1,"end");
      cg_user_data_write("Parallel_Boundary_Maps");
      nud++;
      ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
      CGNS_Errorcheck(ier);
      for (b=0; b < nb; b++)
      {
        if (nt[b] > 0)
        {
          conn = (int*)realloc((void*)conn,nt[b]*sizeof(int));
          if (conn == NULL)
          {
            fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
            fflush(stderr);
            exit(0);
          }
          //sprintf(sectionname,"Boundary_%d_Parallel_Tri_Map",b+1);
          //cg_user_data_write(sectionname);
          for (c=0; c < nt[b]; c++)
            conn[c] = tri_map[b][c]+1;
          sprintf(sectionname,"Boundary_%d_Global_Tri_Number",b+1);
          ier = cg_array_write(sectionname,Integer,1,&nt[b],conn);
          CGNS_Errorcheck(ier);
        }
        if (nq[b] > 0)
        {
          conn = (int*)realloc((void*)conn,nq[b]*sizeof(int));
          if (conn == NULL)
          {
            fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
            fflush(stderr);
            exit(0);
          }
          //sprintf(sectionname,"Boundary_%d_Parallel_Quad_Map",b+1);
          //cg_user_data_write(sectionname);
          for (c=0; c < nq[b]; c++)
            conn[c] = quad_map[b][c]+1;
          sprintf(sectionname,"Boundary_%d_Global_Quad_Number",b+1);
          ier = cg_array_write(sectionname,Integer,1,&nq[b],conn);
          CGNS_Errorcheck(ier);
        }
      }

      if (ntet > 0)
      {
        conn = (int*)realloc((void*)conn,ntet*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_goto(index_file,index_base,"Zone_t",1,"end");
        cg_user_data_write("Parallel_Tet_Map");
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        for (c=0; c < ntet; c++)
          conn[c] = tet_map[c]+1;
        ier = cg_array_write("Global_Tet_Number",Integer,1,&ntet,conn);
        CGNS_Errorcheck(ier);
      }
      if (npyr > 0)
      {
        conn = (int*)realloc((void*)conn,npyr*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_goto(index_file,index_base,"Zone_t",1,"end");
        cg_user_data_write("Parallel_Pyramid_Map");
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        for (c=0; c < npyr; c++)
          conn[c] = pyr_map[c]+1;
        ier = cg_array_write("Global_Pyramid_Number",Integer,1,&npyr,conn);
        CGNS_Errorcheck(ier);
      }
      if (npri > 0)
      {
        conn = (int*)realloc((void*)conn,npri*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_goto(index_file,index_base,"Zone_t",1,"end");
        cg_user_data_write("Parallel_Prism_Map");
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        for (c=0; c < npri; c++)
          conn[c] = pri_map[c]+1;
        ier = cg_array_write("Global_Prism_Number",Integer,1,&npri,conn);
        CGNS_Errorcheck(ier);
      }
      if (nhex > 0)
      {
        conn = (int*)realloc((void*)conn,nhex*sizeof(int));
        if (conn == NULL)
        {
          fprintf(stderr,"\nP_CGNS_write: Unable to reallocate memory for conn array!\n");
          fflush(stderr);
          exit(0);
        }
        cg_goto(index_file,index_base,"Zone_t",1,"end");
        cg_user_data_write("Parallel_Hex_Map");
        nud++;
        ier = cg_goto(index_file,index_base,"Zone_t",1,"UserDefinedData_t",nud,"end");
        CGNS_Errorcheck(ier);
        for (c=0; c < nhex; c++)
          conn[c] = hex_map[c]+1;
        ier = cg_array_write("Global_Hex_Number",Integer,1,&nhex,conn);
        CGNS_Errorcheck(ier);
      }
    }
  }

  ier = cg_close(index_file);
  CGNS_Errorcheck(ier);

  if (conn)
    free(conn);

  return(0);
}
