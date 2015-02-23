#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "mesh.h"
#include "refine_elements.h"
#include "Spacing_Field.h"
#include "List.h"
#include "CGNS.h"
#include "sort.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

void mesh_obj::create_emap(int nelem, int type, int **elem_node, List **hash, int **emap)
{
  int i, j, k, n, flag;

  int temp[num_procs];

  //set to zero
  for (i = 0; i < num_procs; i++)
    temp[i] = 0;

  //look through elements, if local, then number, else, give -1, also create hash   
  for (n = 0; n < nelem; n++)
  {
    for (i = 0; i < type; i++)
    {
      j = elem_node[n][i];
      hash[j]->Add_To_List(n);
      temp[pmap[j][1]]++;
    }
    //loop thru temp to determine element ownership, if any node in fight, lowest wins
    flag = 0;
    for (i = 0; i < num_procs && !flag; i++)
    {
      if (temp[i] > 0)
      {
        k = i;
        flag = 1;
      }
    }
    /*//set final proc in temp to be value to beat
    k = num_procs - 1;
    j = temp[num_procs - 1];  
    //use temp to determine element ownership
    for (i = num_procs - 2; i > -1; i--)
    {
      //start at highest proc, since in a tie, lowest wins
      if (temp[i] >= j)
      {
        k = i;
        j = temp[i];
      }
    }*/
    //set map
    emap[n][0] = k;
    if (k == my_rank)
      emap[n][1] = n;
    else
      emap[n][1] = -1;

    //reset temp
    for (i = 0; i < num_procs; i++)
      temp[i] = 0;
  }
  
  return;
}

void mesh_obj::exchange_emap(int **emap, int nelem, List **hash, int type, int *elem_map, int bndy)
{
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses;
  int *bdim, *bdim2, *bdim21;
  char **sbuff, **rbuff;
  int n, p, i, j, k, flag, q;

  sendcnt = new int [num_procs];
  recvcnt = new int [num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sendcnt[n] = 0;
    recvcnt[n] = 0;
  }
  srequest = new MPI_Request[num_procs];
  rrequest = new MPI_Request[num_procs];
  statuses = new MPI_Status[num_procs];
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nAllocated MPI.\n");
  //fflush(debug_f);
  #endif

  //cycle thru elements	
  for (n = 0; n < nelem; n++)
    if (emap[n][0] != my_rank)
      sendcnt[emap[n][0]]++; 

  #ifdef _DEBUG
  //fprintf(debug_f,"\nDetermined sendcounts.\n");
  //fflush(debug_f);
  #endif

  //now, send and recv count for all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (p == my_rank)
      continue;

    MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
    nreq_s++;

    MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
    nreq_r++;
  }
 
  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished first comm.\n");
  //fflush(debug_f);
  #endif

  //allocate buffers to send/recv element global and local and node to look in hash table for
  bdim = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim[n] = MAX(sendcnt[n],recvcnt[n])*3*sizeof(int); //indices
  sbuff = new char*[num_procs];
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
  {
    sbuff[n] = new char[bdim[n]];
    rbuff[n] = new char[bdim[n]];
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nAllocated buffers.\n");
  //fflush(debug_f);
  #endif

  //now, pack indices
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] == 0 || p == my_rank)
      continue;

    //set position to 0
    sposition=0;
 
    //pack node indices needed, as would be seen by other proc and the corresponding index on the recv proc
    for (n = 0; n < nelem; n++)
    {
      flag = 0;
      if (emap[n][0] == p)
      {
        for (i = 0; i < type && !flag; i++)
          {
          if (type == 4 && bndy == -1)
            j = tet_n[n][i];
          if (type == 5 && bndy == -1)
            j = pyr_n[n][i];
          if (type == 6 && bndy == -1)
            j = pri_n[n][i];
          if (type == 8 && bndy == -1)
            j = hex_n[n][i];
          if (type == 3 && bndy > -1)
            j = t_n[bndy][n][i];
          if (type == 4 && bndy > -1)
            j = q_n[bndy][n][i];
          if (pmap[j][1] == p)
            {
            flag = 1;
            MPI_Pack(&(pmap[j][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            }
          }
        MPI_Pack(&(elem_map[n]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
      }
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nPacked info.\n");
  //fflush(debug_f);
  #endif

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished second comm.\n");
  //fflush(debug_f);
  #endif

  //now, delete sbuff and resize
  for (n = 0; n < num_procs; n++)
    if (sbuff[n] != 0) delete [] sbuff[n];
  delete [] sbuff;
 
  //resize
  bdim21 = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim21[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //marking and indices
  sbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
    sbuff[n] = new char[bdim21[n]];

  //now, unpack node, look for triangle in hash, find local for global sent, send back owning local with non-owning local
  for (p = 0; p < num_procs; p++)
  {
    #ifdef _DEBUG
    //fprintf(debug_f,"\nrecvcnt[%d] = %d.\n",p,recvcnt[p]);
    //fflush(debug_f);
    #endif
    sposition = rposition = 0; //reset pos
    for (n = 0; n < recvcnt[p]; n++)
    {
      //unpack node index initially
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
      
      //unpack global index
      MPI_Unpack(rbuff[p],bdim[p],&rposition,&k,1,MPI_INT,MPI_COMM_WORLD);
      
      //look through hash
      flag = 0;
      for (j = 0; j < hash[i]->max && !flag; j++)
      {
        if (elem_map[hash[i]->list[j]] == k)
        {
          flag = 1;
          //set aside this owning local index
          k = hash[i]->list[j];
          //unpack other procs local index to send back into q
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&q,1,MPI_INT,MPI_COMM_WORLD);
          //pack back up other procs local and then owning local
          MPI_Pack(&q,1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
          MPI_Pack(&k,1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
        }
      }
      if (flag == 0)
      {
        fprintf(stderr,"\nUnable to find matching element for global %d.\n",k);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }

  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished unpack and pack.\n");
  //fflush(debug_f);
  #endif

  //now, delete rbuff and resize
  for (n = 0; n < num_procs; n++)
    if (rbuff[n] != 0) delete [] rbuff[n];
  delete [] rbuff;
 
  //resize
  bdim2 = new int[num_procs];
  for (n = 0; n < num_procs; n++)
    bdim2[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //marking and indices
  rbuff = new char*[num_procs];
  for (n = 0; n < num_procs; n++)
    rbuff[n] = new char[bdim2[n]];

  //now, send and recv packets from all procs
  nreq_s = nreq_r = 0;
  for (p = 0; p < num_procs; p++)
  {
    if (recvcnt[p] > 0 && p != my_rank) 
    {
      MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
    }

    if (sendcnt[p] > 0 && p != my_rank) 
    {
      MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
  }

  //now, wait for finish before unpacking
  MPI_Waitall(nreq_s,srequest,statuses);
  MPI_Waitall(nreq_r,rrequest,statuses);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished third comm.\n");
  //fflush(debug_f);
  #endif

  //finally, unpack nodes in proper position
  for (p = 0; p < num_procs; p++)
  {
    rposition = 0; //reset pos
    for (n = 0; n < sendcnt[p]; n++)
    {
      //unpack local and owning local
      MPI_Unpack(rbuff[p],bdim2[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
      MPI_Unpack(rbuff[p],bdim2[p],&rposition,&(emap[i][1]),1,MPI_INT,MPI_COMM_WORLD);
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished final unpack.\n");
  //fflush(debug_f);
  #endif

  //finally, free MPI mem
  delete[] sendcnt;
  delete[] recvcnt;
  delete[] srequest;
  delete[] rrequest;
  delete[] statuses;
  for (n = 0; n < num_procs; n++)
  {
    if (sbuff[n] != 0) delete [] sbuff[n];
    if (rbuff[n] != 0) delete [] rbuff[n];
  }
  delete[] sbuff;
  delete[] rbuff;
  delete[] bdim;
  delete[] bdim2;
  delete[] bdim21;
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nFreed mem.\n");
  //fflush(debug_f);
  #endif

  return;
}

void mesh_obj::element_pi()
{
  int i, j, n, k;

  #ifdef _DEBUG
  //fprintf(debug_f,"\nAllocating emap mem.\n");
  //fflush(debug_f);
  #endif

  //allocate for all emaps
  tri_emap = (int***)calloc(nb,sizeof(int**));
  quad_emap = (int***)calloc(nb,sizeof(int**));
  for (i = 0; i < nb; i++)
  {
    tri_emap[i] = (int**)calloc(nt[i],sizeof(int*));
    quad_emap[i] = (int**)calloc(nq[i],sizeof(int*));
    for (j = 0; j < nt[i]; j++)
    {
      tri_emap[i][j] = (int*)calloc(2,sizeof(int));
      for (k = 0; k < 2; k++)
        tri_emap[i][j][k] = -1;
    }
    for (j = 0; j < nq[i]; j++)
    {
      quad_emap[i][j] = (int*)calloc(2,sizeof(int));
      for (k = 0; k < 2; k++)
        quad_emap[i][j][k] = -1;
    }
  }
  
  tet_emap = (int**)calloc(ntet,sizeof(int*));
  for (i = 0; i < ntet; i++)
  {
    tet_emap[i] = (int*)calloc(2,sizeof(int));
    for (k = 0; k < 2; k++)
      tet_emap[i][k] = -1;
  }
  pri_emap = (int**)calloc(npri,sizeof(int*));
  for (i = 0; i < npri; i++)
  {
    pri_emap[i] = (int*)calloc(2,sizeof(int));
    for (k = 0; k < 2; k++)
      pri_emap[i][k] = -1;
  }
  pyr_emap = (int**)calloc(npyr,sizeof(int*));
  for (i = 0; i < npyr; i++)
  {
    pyr_emap[i] = (int*)calloc(2,sizeof(int));
    for (k = 0; k < 2; k++)
      pyr_emap[i][k] = -1;
  }
  hex_emap = (int**)calloc(nhex,sizeof(int*));
  for (i = 0; i < nhex; i++)
  {
    hex_emap[i] = (int*)calloc(2,sizeof(int));
    for (k = 0; k < 2; k++)
      hex_emap[i][k] = -1;
  }

  //create hash that can be reused for each element type
  List **hash;
  hash = new List*[nn];
  for (n=0; n < nn; n++)
    hash[n] = new List();

  //DO 4 BASIC TYPES

  #ifdef _DEBUG
  //fprintf(debug_f,"\nFinished allocating emap mem.\n");
  //fflush(debug_f);
  #endif

  //create hash table and element map for local elements
  create_emap(ntet, 4, tet_n, hash, tet_emap);

  #ifdef _DEBUG
  //fprintf(debug_f,"\nCreated tet emap.\n");
  //fflush(debug_f);
  #endif
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  //now, comm with other procs to determine local element number of unowned elements by searching 
  //nodes on that proc hash, finding all possible local and searching for that with matching global
  exchange_emap(tet_emap, ntet, hash, 4, tet_map, -1);

  //redim hash    
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);

  #ifdef _DEBUG
  //fprintf(debug_f,"\nExchanged tet emap, redimmed hash.\n");
  //fflush(debug_f);
  #endif
    
  //create hash table and element map for local elements
  create_emap(npyr, 5, pyr_n, hash, pyr_emap);

  #ifdef _DEBUG
  //fprintf(debug_f,"\nCreated pyr emap.\n");
  //fflush(debug_f);
  #endif
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  //now, comm with other procs to determine local element number of unowned elements by searching 
  //nodes on that proc hash, finding all possible local and searching for that with matching global
  exchange_emap(pyr_emap, npyr, hash, 5, pyr_map, -1);

  //redim hash    
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nExchanged pyr emap, redimmed hash.\n");
  //fflush(debug_f);
  #endif
  
  //create hash table and element map for local elements
  create_emap(npri, 6, pri_n, hash, pri_emap);
 
  #ifdef _DEBUG
  //fprintf(debug_f,"\nCreated pri emap.\n");
  //fflush(debug_f);
  #endif
  
  MPI_Barrier(MPI_COMM_WORLD);
    
  //now, comm with other procs to determine local element number of unowned elements by searching 
  //nodes on that proc hash, finding all possible local and searching for that with matching global
  exchange_emap(pri_emap, npri, hash, 6, pri_map, -1);

  //redim hash    
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nExchanged pri emap, redimmed hash.\n");
  //fflush(debug_f);
  #endif
  
  //create hash table and element map for local elements
  create_emap(nhex, 8, hex_n, hash, hex_emap);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nCreated hex emap.\n");
  //fflush(debug_f);
  #endif
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  //now, comm with other procs to determine local element number of unowned elements by searching 
  //nodes on that proc hash, finding all possible local and searching for that with matching global
  exchange_emap(hex_emap, nhex, hash, 8, hex_map, -1);

  //redim hash    
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nExchanged hex emap, redimmed hash.\n");
  //fflush(debug_f);
  #endif
 
  //DO BOUNDARY FACETS
  for (i = 0 ; i < nb; i++)
  {
    //create hash table and element map for local elements
    create_emap(nt[i], 3, t_n[i], hash, tri_emap[i]);
  
    #ifdef _DEBUG
    //fprintf(debug_f,"\nCreated bndry %d tri emap.\n",i);
    //fflush(debug_f);
    #endif

    MPI_Barrier(MPI_COMM_WORLD);
  
    //now, comm with other procs to determine local element number of unowned elements by searching 
    //nodes on that proc hash, finding all possible local and searching for that with matching global
    exchange_emap(tri_emap[i], nt[i], hash, 3, tri_map[i], i);

    //redim hash    
    for (n=0; n < nn; n++)
      hash[n]->Redimension(0);

    #ifdef _DEBUG
    //fprintf(debug_f,"\nExchanged tri %d emap, redimmed hash.\n",i);
    //fflush(debug_f);
    #endif
    
    //create hash table and element map for local elements
    create_emap(nq[i], 4, q_n[i], hash, quad_emap[i]);
  
    #ifdef _DEBUG
    //fprintf(debug_f,"\nCreated bndry %d quad emap.\n",i);
    //fflush(debug_f);
    #endif

    MPI_Barrier(MPI_COMM_WORLD);
  
    //now, comm with other procs to determine local element number of unowned elements by searching 
    //nodes on that proc hash, finding all possible local and searching for that with matching global
    exchange_emap(quad_emap[i], nq[i], hash, 4, quad_map[i], i);

    //redim hash    
    for (n=0; n < nn; n++)
      hash[n]->Redimension(0);
    
    #ifdef _DEBUG
    //fprintf(debug_f,"\nExchanged quad %d emap, redimmed hash.\n",i);
    //fflush(debug_f);
    #endif
  }

  #ifdef _DEBUG
  //(THIS IS A DEBUGGING TOOL):  run check to assure all filled in...will print off any proc so when it occurs, we can see at least one reason why abourt occurred.
  for (i = 0; i < ntet; i++)
  {
    if (tet_emap[i][1] < 0)
    {
      fprintf(debug_f,"\nEmap for tet %d = %d.  Mapping not correct.  Exiting....\n",i,tet_emap[i][1]);
      fflush(debug_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (i = 0; i < npyr; i++)
  {
    if (pyr_emap[i][1] < 0)
    {
      fprintf(debug_f,"\nEmap for pyr %d = %d.  Mapping not correct.  Exiting....\n",i,pyr_emap[i][1]);
      fflush(debug_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (i = 0; i < npri; i++)
  {
    if (pri_emap[i][1] < 0)
    {
      fprintf(debug_f,"\nEmap for pri %d = %d.  Mapping not correct.  Exiting....\n",i,pri_emap[i][1]);
      fflush(debug_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (i = 0; i < nhex; i++)
  {
    if (hex_emap[i][1] < 0)
    {
      fprintf(debug_f,"\nEmap for hex %d = %d.  Mapping not correct.  Exiting....\n",i,hex_emap[i][1]);
      fflush(debug_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }

  for (i = 0; i < nb; i++)
  {
    for (j = 0; j < nt[i]; j++)
    {
      if (tri_emap[i][j][1] < 0)
      {
        fprintf(debug_f,"\nEmap for tri %d on bd %d = %d.  Mapping not correct.  Exiting....\n",j,i,tri_emap[i][j][1]);
        fflush(debug_f);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
    for (j = 0; j < nq[i]; j++)
    {
      if (quad_emap[i][j][1] < 0)
      {
        fprintf(debug_f,"\nEmap for quad %d on bd %d = %d.  Mapping not correct.  Exiting....\n",j,i,quad_emap[i][j][1]);
        fflush(debug_f);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }
  
  #endif
  
  //free mem
  for (n=0; n < nn; n++)
    delete hash[n];
  delete hash;

  return;
}

//will print off any proc so when it occurs, we can see at least one reason why abourt occurred.
int retrieve_face_nodes(int n0, int n1, int n2, int n3, int tet, int pyr, int pri, int hex,
              int *conn, int *fnodes)
{
  int i, j, m;
  int corner[4];

  m = 0;
  if (tet >= 0) m++;
  if (pyr >= 0) m++;
  if (pri >= 0) m++;
  if (hex >= 0) m++;
  if (m != 1)
  {
    fprintf(stderr,"\nRETRIEVE_FACE_NODES: Zero or multiple elements passed in!");
    fprintf(stderr,"\n          Only one allowed!");
    fprintf(stderr,"\n          tet = %d",tet);
    fprintf(stderr,"\n          pyr = %d",pyr);
    fprintf(stderr,"\n          pri = %d",pri);
    fprintf(stderr,"\n          hex = %d\n",hex);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }

  corner[0] = n0;
  corner[1] = n1;
  corner[2] = n2;
  corner[3] = n3;
  for (i=0; i < 9; i++)
    fnodes[i] = -1;

  if (tet >= 0 && corner[3] < 0)
  {
    // face 1
    j=0;
    for (i=0; i < 3 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 3 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 3 && j == 2; i++)
      if (corner[i] == conn[2]) j++;
    if (j == 3)
    {
      fnodes[0] = conn[0];
      fnodes[1] = conn[1];
      fnodes[2] = conn[2];
      fnodes[3] = conn[4];
      fnodes[4] = conn[5];
      fnodes[5] = conn[6];
      return(0);
    }
    // face 2
    j=0;
    for (i=0; i < 3 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 3 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 3 && j == 2; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 3)
    {
      fnodes[0] = conn[0];
      fnodes[1] = conn[3];
      fnodes[2] = conn[1];
      fnodes[3] = conn[7];
      fnodes[4] = conn[8];
      fnodes[5] = conn[4];
      return(1);
    }
    // face 3
    j=0;
    for (i=0; i < 3 && j == 0; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 3 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 3 && j == 2; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 3)
    {
      fnodes[0] = conn[1];
      fnodes[1] = conn[3];
      fnodes[2] = conn[2];
      fnodes[3] = conn[8];
      fnodes[4] = conn[9];
      fnodes[5] = conn[5];
      return(2);
    }
    // face 4
    j=0;
    for (i=0; i < 3 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 3 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 3 && j == 2; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 3)
    {
      fnodes[0] = conn[2];
      fnodes[1] = conn[3];
      fnodes[2] = conn[0];
      fnodes[3] = conn[9];
      fnodes[4] = conn[7];
      fnodes[5] = conn[6];
      return(3);
    }
  }

  if (pyr >= 0)
  {
    // face 1
    if (corner[3] >= 0)
    {
      j=0;
      for (i=0; i < 4 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 4 && j == 1; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 4 && j == 2; i++)
        if (corner[i] == conn[2]) j++;
      for (i=0; i < 4 && j == 3; i++)
        if (corner[i] == conn[3]) j++;
      if (j == 4)
      {
        fnodes[0] = conn[0];
        fnodes[1] = conn[1];
        fnodes[2] = conn[2];
        fnodes[3] = conn[3];
        fnodes[4] = conn[5];
        fnodes[5] = conn[6];
        fnodes[6] = conn[7];
        fnodes[7] = conn[8];
        fnodes[8] = conn[13];
        return(0);
      }
    }
    // face 2
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[4]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[0];
        fnodes[1] = conn[4];
        fnodes[2] = conn[1];
        fnodes[3] = conn[9];
        fnodes[4] = conn[10];
        fnodes[5] = conn[5];
        return(1);
      }
    }
    // face 3
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[2]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[4]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[1];
        fnodes[1] = conn[4];
        fnodes[2] = conn[2];
        fnodes[3] = conn[10];
        fnodes[4] = conn[11];
        fnodes[5] = conn[6];
        return(2);
      }
    }
    // face 4
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[2]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[3]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[4]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[2];
        fnodes[1] = conn[4];
        fnodes[2] = conn[3];
        fnodes[3] = conn[11];
        fnodes[4] = conn[12];
        fnodes[5] = conn[7];
        return(3);
      }
    }
    // face 5
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[3]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[4]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[3];
        fnodes[1] = conn[4];
        fnodes[2] = conn[0];
        fnodes[3] = conn[12];
        fnodes[4] = conn[9];
        fnodes[5] = conn[8];
        return(4);
      }
    }
  }

  if (pri >= 0)
  {
    // face 1
    if (corner[3] >= 0)
    {
      j=0;
      for (i=0; i < 4 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 4 && j == 1; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 4 && j == 2; i++)
        if (corner[i] == conn[3]) j++;
      for (i=0; i < 4 && j == 3; i++)
        if (corner[i] == conn[4]) j++;
      if (j == 4)
      {
        fnodes[0] = conn[0];
        fnodes[1] = conn[3];
        fnodes[2] = conn[4];
        fnodes[3] = conn[1];
        fnodes[4] = conn[9];
        fnodes[5] = conn[12];
        fnodes[6] = conn[10];
        fnodes[7] = conn[6];
        fnodes[8] = conn[15];
        return(0);
      }
    }
    // face 2
    if (corner[3] >= 0)
    {
      j=0;
      for (i=0; i < 4 && j == 0; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 4 && j == 1; i++)
        if (corner[i] == conn[2]) j++;
      for (i=0; i < 4 && j == 2; i++)
        if (corner[i] == conn[4]) j++;
      for (i=0; i < 4 && j == 3; i++)
        if (corner[i] == conn[5]) j++;
      if (j == 4)
      {
        fnodes[0] = conn[1];
        fnodes[1] = conn[4];
        fnodes[2] = conn[5];
        fnodes[3] = conn[2];
        fnodes[4] = conn[10];
        fnodes[5] = conn[13];
        fnodes[6] = conn[11];
        fnodes[7] = conn[7];
        fnodes[8] = conn[16];
        return(1);
      }
    }
    // face 3
    if (corner[3] >= 0)
    {
      j=0;
      for (i=0; i < 4 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 4 && j == 1; i++)
        if (corner[i] == conn[2]) j++;
      for (i=0; i < 4 && j == 2; i++)
        if (corner[i] == conn[3]) j++;
      for (i=0; i < 4 && j == 3; i++)
        if (corner[i] == conn[5]) j++;
      if (j == 4)
      {
        fnodes[0] = conn[2];
        fnodes[1] = conn[5];
        fnodes[2] = conn[3];
        fnodes[3] = conn[0];
        fnodes[4] = conn[11];
        fnodes[5] = conn[14];
        fnodes[6] = conn[9];
        fnodes[7] = conn[8];
        fnodes[8] = conn[17];
        return(2);
      }
    }
    // face 4
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[0]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[1]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[2]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[0];
        fnodes[1] = conn[1];
        fnodes[2] = conn[2];
        fnodes[3] = conn[6];
        fnodes[4] = conn[7];
        fnodes[5] = conn[8];
        return(3);
      }
    }
    // face 5
    if (corner[3] < 0)
    {
      j=0;
      for (i=0; i < 3 && j == 0; i++)
        if (corner[i] == conn[3]) j++;
      for (i=0; i < 3 && j == 1; i++)
        if (corner[i] == conn[4]) j++;
      for (i=0; i < 3 && j == 2; i++)
        if (corner[i] == conn[5]) j++;
      if (j == 3)
      {
        fnodes[0] = conn[3];
        fnodes[1] = conn[5];
        fnodes[2] = conn[4];
        fnodes[3] = conn[14];
        fnodes[4] = conn[13];
        fnodes[5] = conn[12];
        return(4);
      }
    }
  }

  if (hex >= 0 && corner[3] >= 0)
  {
    // face 1
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[0];
      fnodes[1] = conn[1];
      fnodes[2] = conn[2];
      fnodes[3] = conn[3];
      fnodes[4] = conn[8];
      fnodes[5] = conn[9];
      fnodes[6] = conn[10];
      fnodes[7] = conn[11];
      fnodes[8] = conn[20];
      return(0);
    }
    // face 2
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[5]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[0];
      fnodes[1] = conn[4];
      fnodes[2] = conn[5];
      fnodes[3] = conn[1];
      fnodes[4] = conn[12];
      fnodes[5] = conn[16];
      fnodes[6] = conn[13];
      fnodes[7] = conn[8];
      fnodes[8] = conn[21];
      return(1);
    }
    // face 3
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[5]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[6]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[1];
      fnodes[1] = conn[5];
      fnodes[2] = conn[6];
      fnodes[3] = conn[2];
      fnodes[4] = conn[13];
      fnodes[5] = conn[17];
      fnodes[6] = conn[14];
      fnodes[7] = conn[9];
      fnodes[8] = conn[22];
      return(2);
    }
    // face 4
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[6]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[2];
      fnodes[1] = conn[6];
      fnodes[2] = conn[7];
      fnodes[3] = conn[3];
      fnodes[4] = conn[14];
      fnodes[5] = conn[18];
      fnodes[6] = conn[15];
      fnodes[7] = conn[10];
      fnodes[8] = conn[23];
      return(3);
    }
    // face 5
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[3];
      fnodes[1] = conn[7];
      fnodes[2] = conn[4];
      fnodes[3] = conn[0];
      fnodes[4] = conn[15];
      fnodes[5] = conn[19];
      fnodes[6] = conn[12];
      fnodes[7] = conn[11];
      fnodes[8] = conn[24];
      return(4);
    }
    // face 6
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[5]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[6]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      fnodes[0] = conn[7];
      fnodes[1] = conn[6];
      fnodes[2] = conn[5];
      fnodes[3] = conn[4];
      fnodes[4] = conn[18];
      fnodes[5] = conn[17];
      fnodes[6] = conn[16];
      fnodes[7] = conn[19];
      fnodes[8] = conn[25];
      return(5);
    }
  }

  return(-1);
}

void store_face_node(int n0, int n1, int n2, int n3, int n, int pyr, int pri, int hex,
              int *conn)
{
  int i, j, m;
  int corner[4];

  m = 0;
  if (pyr >= 0) m++;
  if (pri >= 0) m++;
  if (hex >= 0) m++;
  if (m != 1)
  {
    fprintf(stderr,"\nSTORE_FACE_NODE: Zero or multiple elements passed in!");
    fprintf(stderr,"\n          Only one allowed!");
    fprintf(stderr,"\n          pyr = %d",pyr);
    fprintf(stderr,"\n          pri = %d",pri);
    fprintf(stderr,"\n          hex = %d\n",hex);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }

  corner[0] = n0;
  corner[1] = n1;
  corner[2] = n2;
  corner[3] = n3;

  if (hex >= 0)
  {
    // face 1
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 4)
    {
      if (conn[20] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 1 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[20] = n;
      return;
    }
    // face 2
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[5]) j++;
    if (j == 4)
    {
      if (conn[21] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 2 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[21] = n;
      return;
    }
    // face 3
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[5]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[6]) j++;
    if (j == 4)
    {
      if (conn[22] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 3 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[22] = n;
      return;
    }
    // face 4
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[6]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      if (conn[23] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 4 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[23] = n;
      return;
    }
    // face 5
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      if (conn[24] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 5 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[24] = n;
      return;
    }
    // face 6
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[5]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[6]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[7]) j++;
    if (j == 4)
    {
      if (conn[25] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Hexahedron face 6 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[25] = n;
      return;
    }
  }

  if (pri >= 0)
  {
    // face 1
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[4]) j++;
    if (j == 4)
    {
      if (conn[15] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Prism face 1 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[15] = n;
      return;
    }
    // face 2
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[4]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[5]) j++;
    if (j == 4)
    {
      if (conn[16] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Prism face 2 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[16] = n;
      return;
    }
    // face 3
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[3]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[5]) j++;
    if (j == 4)
    {
      if (conn[17] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Prism face 3 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[17] = n;
      return;
    }
  }

  if (pyr >= 0)
  {
    // face 1
    j=0;
    for (i=0; i < 4 && j == 0; i++)
      if (corner[i] == conn[0]) j++;
    for (i=0; i < 4 && j == 1; i++)
      if (corner[i] == conn[1]) j++;
    for (i=0; i < 4 && j == 2; i++)
      if (corner[i] == conn[2]) j++;
    for (i=0; i < 4 && j == 3; i++)
      if (corner[i] == conn[3]) j++;
    if (j == 4)
    {
      if (conn[13] >= 0)
      {
        fprintf(stderr,"FACE_NODE: Pyramid face 1 mid-node defined!");
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      conn[13] = n;
      return;
    }
  }

  fprintf(stderr,"FACE_NODE: unable to identify neighbor face");
  fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD,0);
  exit(0);

}

//parallel allowances made in make_nabor!!!
//will print off any proc so when it occurs, we can see at least one reason why abort occurred.
void neighbor(int n1, int n2, int n3, int n4,
             int &tet, int &pyr, int &pri, int &hex,
             List **tethash,
             List **pyrhash,
             List **prihash,
             List **hexhash)
{
  int i, tetin, pyrin, priin, hexin, m;

  m = 0;
  if (tet >= 0) m++;
  if (pyr >= 0) m++;
  if (pri >= 0) m++;
  if (hex >= 0) m++;
  if (m > 1)
  {
    fprintf(stderr,"\nNEIGHBOR: Multiple elements passed in!");
    fprintf(stderr,"\n          Maximum of one allowed!");
    fprintf(stderr,"\n          tet = %d",tet);
    fprintf(stderr,"\n          pyr = %d",pyr);
    fprintf(stderr,"\n          pri = %d",pri);
    fprintf(stderr,"\n          hex = %d\n",hex);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }

  tetin = tet;
  pyrin = pyr;
  priin = pri;
  hexin = hex;
  tet = pyr = pri = hex = -1;

  // check for neighbor hexahedron
  if (n4 >= 0)
  {
    for (i=0; i < hexhash[n1]->max && hex < 0; i++)
    {
      m = hexhash[n1]->list[i];
      if (hexin < 0 || m != hexin)
      {
        if (hexhash[n2]->Is_In_List(m))
          if (hexhash[n3]->Is_In_List(m))
            if (hexhash[n4]->Is_In_List(m))
              hex = m;
      }
    }
  }

  if (hex >= 0)
    return;

  // check for neighbor prism
  for (i=0; i < prihash[n1]->max && pri < 0; i++)
  {
    m = prihash[n1]->list[i];
    if (priin < 0 || m != priin)
    {
      if (prihash[n2]->Is_In_List(m))
        if (prihash[n3]->Is_In_List(m))
        {
          if (n4 < 0)
            pri = m;
          else if (prihash[n4]->Is_In_List(m))
            pri = m;
        }
    }
  }

  if (pri >= 0)
    return;

  // check for neighbor pyramid
  for (i=0; i < pyrhash[n1]->max && pyr < 0; i++)
  {
    m = pyrhash[n1]->list[i];
    if (pyrin < 0 || m != pyrin)
    {
      if (pyrhash[n2]->Is_In_List(m))
        if (pyrhash[n3]->Is_In_List(m))
        {
          if (n4 < 0)
            pyr = m;
          else if (pyrhash[n4]->Is_In_List(m))
            pyr = m;
        }
    }
  }

  if (pyr >= 0)
    return;

  // check for neighbor tetrahedra
  if (n4 < 0)
  {
    for (i=0; i < tethash[n1]->max && tet < 0; i++)
    {
      m = tethash[n1]->list[i];
      if (tetin < 0 || m != tetin)
      {
        if (tethash[n2]->Is_In_List(m))
          if (tethash[n3]->Is_In_List(m))
            tet = m;
      }
    }
  }

  return;

}

int tet_side(int conn[4], int n0, int n1, int n2)
{
  // return side of tet with proper winding
  if ((conn[0] == n0 && conn[2] == n1 && conn[1] == n2) ||
      (conn[0] == n2 && conn[2] == n0 && conn[1] == n1) ||
      (conn[0] == n1 && conn[2] == n2 && conn[1] == n0))
    return(0);
  if ((conn[0] == n0 && conn[1] == n1 && conn[3] == n2) ||
      (conn[0] == n2 && conn[1] == n0 && conn[3] == n1) ||
      (conn[0] == n1 && conn[1] == n2 && conn[3] == n0))
    return(1);
  if ((conn[1] == n0 && conn[2] == n1 && conn[3] == n2) ||
      (conn[1] == n2 && conn[2] == n0 && conn[3] == n1) ||
      (conn[1] == n1 && conn[2] == n2 && conn[3] == n0))
    return(2);
  if ((conn[2] == n0 && conn[0] == n1 && conn[3] == n2) ||
      (conn[2] == n2 && conn[0] == n0 && conn[3] == n1) ||
      (conn[2] == n1 && conn[0] == n2 && conn[3] == n0))
    return(3);

  return(-1);
}

int pyramid_side(int conn[5], int n0, int n1, int n2, int n3)
{
  // return side of pyramid with proper winding
  if (n3 >= 0)
  {
    if ((conn[0] == n0 && conn[3] == n1 && conn[2] == n2 && conn[1] == n3) ||
        (conn[0] == n3 && conn[3] == n0 && conn[2] == n1 && conn[1] == n2) ||
        (conn[0] == n2 && conn[3] == n3 && conn[2] == n0 && conn[1] == n1) ||
        (conn[0] == n1 && conn[3] == n2 && conn[2] == n3 && conn[1] == n0))
      return(0);
  } else
  {
    if ((conn[0] == n0 && conn[1] == n1 && conn[4] == n2) ||
        (conn[0] == n2 && conn[1] == n0 && conn[4] == n1) ||
        (conn[0] == n1 && conn[1] == n2 && conn[4] == n0))
      return(1);
    if ((conn[1] == n0 && conn[2] == n1 && conn[4] == n2) ||
        (conn[1] == n2 && conn[2] == n0 && conn[4] == n1) ||
        (conn[1] == n1 && conn[2] == n2 && conn[4] == n0))
      return(2);
    if ((conn[2] == n0 && conn[3] == n1 && conn[4] == n2) ||
        (conn[2] == n2 && conn[3] == n0 && conn[4] == n1) ||
        (conn[2] == n1 && conn[3] == n2 && conn[4] == n0))
      return(3);
    if ((conn[3] == n0 && conn[0] == n1 && conn[4] == n2) ||
        (conn[3] == n2 && conn[0] == n0 && conn[4] == n1) ||
        (conn[3] == n1 && conn[0] == n2 && conn[4] == n0))
      return(4);
  }

  return(-1);
}

int prism_side(int conn[6], int n0, int n1, int n2, int n3)
{
  // return side of prism with proper winding
  if (n3 >= 0)
  {
    if ((conn[0] == n0 && conn[1] == n1 && conn[4] == n2 && conn[3] == n3) ||
        (conn[0] == n3 && conn[1] == n0 && conn[4] == n1 && conn[3] == n2) ||
        (conn[0] == n2 && conn[1] == n3 && conn[4] == n0 && conn[3] == n1) ||
        (conn[0] == n1 && conn[1] == n2 && conn[4] == n3 && conn[3] == n0))
      return(0);
    if ((conn[1] == n0 && conn[2] == n1 && conn[5] == n2 && conn[4] == n3) ||
        (conn[1] == n3 && conn[2] == n0 && conn[5] == n1 && conn[4] == n2) ||
        (conn[1] == n2 && conn[2] == n3 && conn[5] == n0 && conn[4] == n1) ||
        (conn[1] == n1 && conn[2] == n2 && conn[5] == n3 && conn[4] == n0))
      return(1);
    if ((conn[2] == n0 && conn[0] == n1 && conn[3] == n2 && conn[5] == n3) ||
        (conn[2] == n3 && conn[0] == n0 && conn[3] == n1 && conn[5] == n2) ||
        (conn[2] == n2 && conn[0] == n3 && conn[3] == n0 && conn[5] == n1) ||
        (conn[2] == n1 && conn[0] == n2 && conn[3] == n3 && conn[5] == n0))
      return(2);
  } else
  {
    if ((conn[0] == n0 && conn[2] == n1 && conn[1] == n2) ||
        (conn[0] == n2 && conn[2] == n0 && conn[1] == n1) ||
        (conn[0] == n1 && conn[2] == n2 && conn[1] == n0))
      return(3);
    if ((conn[3] == n0 && conn[4] == n1 && conn[5] == n2) ||
        (conn[3] == n2 && conn[4] == n0 && conn[5] == n1) ||
        (conn[3] == n1 && conn[4] == n2 && conn[5] == n0))
      return(4);
  }

  return(-1);
}

int hex_side(int conn[8], int n0, int n1, int n2, int n3)
{
  // return side of hex with proper winding
  if ((conn[0] == n0 && conn[3] == n1 && conn[2] == n2 && conn[1] == n3) ||
      (conn[0] == n3 && conn[3] == n0 && conn[2] == n1 && conn[1] == n2) ||
      (conn[0] == n2 && conn[3] == n3 && conn[2] == n0 && conn[1] == n1) ||
      (conn[0] == n1 && conn[3] == n2 && conn[2] == n3 && conn[1] == n0))
    return(0);
  if ((conn[0] == n0 && conn[1] == n1 && conn[5] == n2 && conn[4] == n3) ||
      (conn[0] == n3 && conn[1] == n0 && conn[5] == n1 && conn[4] == n2) ||
      (conn[0] == n2 && conn[1] == n3 && conn[5] == n0 && conn[4] == n1) ||
      (conn[0] == n1 && conn[1] == n2 && conn[5] == n3 && conn[4] == n0))
    return(1);
  if ((conn[1] == n0 && conn[2] == n1 && conn[6] == n2 && conn[5] == n3) ||
      (conn[1] == n3 && conn[2] == n0 && conn[6] == n1 && conn[5] == n2) ||
      (conn[1] == n2 && conn[2] == n3 && conn[6] == n0 && conn[5] == n1) ||
      (conn[1] == n1 && conn[2] == n2 && conn[6] == n3 && conn[5] == n0))
    return(2);
  if ((conn[2] == n0 && conn[3] == n1 && conn[7] == n2 && conn[6] == n3) ||
      (conn[2] == n3 && conn[3] == n0 && conn[7] == n1 && conn[6] == n2) ||
      (conn[2] == n2 && conn[3] == n3 && conn[7] == n0 && conn[6] == n1) ||
      (conn[2] == n1 && conn[3] == n2 && conn[7] == n3 && conn[6] == n0))
    return(3);
  if ((conn[0] == n0 && conn[4] == n1 && conn[7] == n2 && conn[3] == n3) ||
      (conn[0] == n3 && conn[4] == n0 && conn[7] == n1 && conn[3] == n2) ||
      (conn[0] == n2 && conn[4] == n3 && conn[7] == n0 && conn[3] == n1) ||
      (conn[0] == n1 && conn[4] == n2 && conn[7] == n3 && conn[3] == n0))
    return(4);
  if ((conn[4] == n0 && conn[5] == n1 && conn[6] == n2 && conn[7] == n3) ||
      (conn[4] == n3 && conn[5] == n0 && conn[6] == n1 && conn[7] == n2) ||
      (conn[4] == n2 && conn[5] == n3 && conn[6] == n0 && conn[7] == n1) ||
      (conn[4] == n1 && conn[5] == n2 && conn[6] == n3 && conn[7] == n0))
    return(5);

  return(-1);
}

int mesh_obj::refine_io(int mode, char sname[])
{
  int flag = 0; 
  int parallel = 1;
  int gnn, gntet, gnpyr, gnpri, gnhex, gnt, gnq, n, b;
  
  switch (mode)
  {
    case -1:
      if (num_procs > 1)
      {
        flag = P_CGNS_read(parallel,sname,nn,&node,nb,&bname,&nt,&t_n,&nq,&q_n,ntet,&tet_n,
                           npyr,&pyr_n,npri,&pri_n,nhex,&hex_n,&pmap,&tri_map,&quad_map,
                           &tet_map,&pyr_map,&pri_map,&hex_map);
      } else
      {
        //no need for maps since no need to reassemble
        flag = CGNS_read(sname,nn,&node,nb,&bname,&nt,&t_n,&nq,&q_n,ntet,&tet_n,npyr,&pyr_n,
                         npri,&pri_n,nhex,&hex_n);
      }
      break;
    case 1:

      if (num_procs > 1)
      {
        flag=P_CGNS_write(parallel,sname,nn,node,nb,bname,nt,t_n,nq,q_n,ntet,tet_n,npyr,pyr_n,
                          npri,pri_n,nhex,hex_n,pmap,tri_map,quad_map,tet_map,pyr_map,pri_map,hex_map);
      } else
      {
        //no need for maps since no need to reassemble
        flag=CGNS_write(sname,nn,node,nb,bname,nt,t_n,nq,q_n,ntet,tet_n,npyr,pyr_n,
                        npri,pri_n,nhex,hex_n); 
      }
      break;
  }

  gnn=nn;
  gntet=ntet;
  gnpyr=npyr;
  gnpri=npri;
  gnhex=nhex;

  int local;

  if (num_procs > 1)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    gnn=local=0;
    for (n=0; n < nn; n++)
      local = MAX(local,pmap[n][0]+1);
    MPI_Allreduce(&local,&gnn,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  
    gntet=local=0;
    for (n=0; n < ntet; n++)
      local = MAX(local,tet_map[n]+1);
    MPI_Allreduce(&local,&gntet,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  
    gnpyr=local=0;
    for (n=0; n < npyr; n++)
      local = MAX(local,pyr_map[n]+1);
    MPI_Allreduce(&local,&gnpyr,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  
    gnpri=local=0;
    for (n=0; n < npri; n++)
      local = MAX(local,pri_map[n]+1);
    MPI_Allreduce(&local,&gnpri,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  
    gnhex=local=0;
    for (n=0; n < nhex; n++)
      local = MAX(local,hex_map[n]+1);
    MPI_Allreduce(&local,&gnhex,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  }

  if (my_rank == 0)
  {
    fprintf(out_f,"\nTotal number of nodes = %d",gnn);
    fprintf(out_f,"\nTotal number of tetrahedra = %d",gntet);
    fprintf(out_f,"\nTotal number of pyramids = %d",gnpyr);
    fprintf(out_f,"\nTotal number of prisms = %d",gnpri);
    fprintf(out_f,"\nTotal number of hexahedra = %d",gnhex);
    fprintf(out_f,"\nTotal number of boundaries = %d",nb);
    fflush(out_f);
  }
  for (b=0; b < nb; b++)
  {
    gnt=nt[b];
    gnq=nq[b];
    if (num_procs > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      gnt=local=0;
      for (n=0; n < nt[b]; n++)
        local = MAX(local,tri_map[b][n]+1);
      MPI_Allreduce(&local,&gnt,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      gnq=local=0;
      for (n=0; n < nq[b]; n++)
        local = MAX(local,quad_map[b][n]+1);
      MPI_Allreduce(&local,&gnq,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    }
    if (my_rank == 0)
    {
      fprintf(out_f,"\nBoundary %d, name = %s",b+1,bname[b]);
      fprintf(out_f,"\nTotal number of triangles      = %d",gnt);
      fprintf(out_f,"\nTotal number of quadrilaterals = %d",gnq);
      fflush(out_f);
    }
  }

  return(flag);

}

int mesh_obj::refine(int target, int tri_flag, int pri_stack, double threshold, char aname[], int mltype, int bnd_flag)
{
  int b, c, e, i, j, k, m, n, nnhash, t, q, p, r, z, flag, changes, total_changes, pass;
  long int tester;
  int n0, n1, n2, n3, n4, n5, n6, n7;
  List **nhash, **tethash, **pyrhash, **prihash, **hexhash;
  int const cdim = 24;
  int tet[cdim][4], pyr[cdim][5], pri[cdim][6], hex[cdim][8], tri[cdim][3], quad[cdim][4];
  int new_tet, new_pyr, new_pri, new_hex, new_t, new_q;
  int ne, lnn, ltet, lpyr, lpri, lhex;
  int (*edge)[3], (*tet_edges)[6], (*pyr_edges)[8], (*pri_edges)[9], (*hex_edges)[12];
  int conn[9];
  double *ml;
  int gcent, gface, gtet, ghex, gpyr, gpri;
  
  int nreq_s, nreq_r, sposition, rposition, temp_node, temp_node2, temp_node3, temp_node6, globalmax = 0, fgnn = 0;
  long int temp_node4, temp_node5;
  int *sendcnt, *recvcnt;
  MPI_Request *srequest, *rrequest;
  MPI_Status *statuses, single_status;
  int *bdim, *bdim2, *bdim21, sdim;
  char **sbuff, **rbuff, *sendbuff, *recvbuff;
  int *oldnt, *oldnq;
  //alloc vars to hold elem type
  int tt = -10;
  int prt = -20;
  int pyt = -30;
  int ht = -40;
  
  oldnt = (int*)calloc(nb,sizeof(int));
  oldnq = (int*)calloc(nb,sizeof(int));

  MPI_Barrier(MPI_COMM_WORLD);

  double lo[3], hi[3];
  //lo[0]=lo[1]=lo[2]=1.0e20;
  //hi[0]=hi[1]=hi[2]=-1.0e20;
  //for (int i=0; i < 3; i++)
  //  for (int n=0; n < nn; n++)
  //  {
  //    lo[i] = MIN(lo[i],node[n][i]);
  //    hi[i] = MAX(hi[i],node[n][i]);
  //  }
  lo[0]=lo[1]=lo[2]=-1.0e20;
  hi[0]=hi[1]=hi[2]=1.0e20;

  if (my_rank == 0)
  {
    fprintf(out_f,"\n\n");
    fflush(out_f);
  }

  SF_Initialize(aname,my_rank+1,num_procs,lo,hi);

  int nsp = SF_Number_of_Entries();
    
  MPI_Barrier(MPI_COMM_WORLD);

  if (num_procs > 1)
  {
    int nspmn, nspmx;
    MPI_Allreduce(&nsp,&nspmn,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&nsp,&nspmx,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (my_rank == 0)
    {
      fprintf(out_f,"\nMinimum number of entries from adaptive file stored = %d",nspmn);
      fprintf(out_f,"\nMaximum number of entries from adaptive file stored = %d\n",nspmx);
      fflush(out_f);
    }
  } else
    if (my_rank == 0) fprintf(out_f,"\n\nNumber of entries in adaptive file = %d\n",nsp);

  nnhash = nn;
  nhash = new List*[nnhash];
  for (n=0; n < nnhash; n++)
    nhash[n] = new List();

  for (c=0; c < ntet; c++)
  {
    for (i=0; i < 4; i++)
    {
      n = tet_n[c][i];
      for (j=i+1; j < 4; j++)
      {
        m = tet_n[c][j];
        nhash[n]->Check_List(m);
        nhash[m]->Check_List(n);
      }
    }
  }
  for (c=0; c < npyr; c++)
  {
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n4);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n4);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n2);
    nhash[n4]->Check_List(n3);
  }
  for (c=0; c < npri; c++)
  {
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n2);
    nhash[n0]->Check_List(n3);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n4);
    nhash[n2]->Check_List(n0);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n5);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n4);
    nhash[n3]->Check_List(n5);
    nhash[n4]->Check_List(n1);
    nhash[n4]->Check_List(n3);
    nhash[n4]->Check_List(n5);
    nhash[n5]->Check_List(n2);
    nhash[n5]->Check_List(n3);
    nhash[n5]->Check_List(n4);
  }
  for (c=0; c < nhex; c++)
  {
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];
    nhash[n0]->Check_List(n1);
    nhash[n0]->Check_List(n3);
    nhash[n0]->Check_List(n4);
    nhash[n1]->Check_List(n0);
    nhash[n1]->Check_List(n2);
    nhash[n1]->Check_List(n5);
    nhash[n2]->Check_List(n1);
    nhash[n2]->Check_List(n3);
    nhash[n2]->Check_List(n6);
    nhash[n3]->Check_List(n0);
    nhash[n3]->Check_List(n2);
    nhash[n3]->Check_List(n7);
    nhash[n4]->Check_List(n0);
    nhash[n4]->Check_List(n5);
    nhash[n4]->Check_List(n7);
    nhash[n5]->Check_List(n1);
    nhash[n5]->Check_List(n4);
    nhash[n5]->Check_List(n6);
    nhash[n6]->Check_List(n2);
    nhash[n6]->Check_List(n5);
    nhash[n6]->Check_List(n7);
    nhash[n7]->Check_List(n3);
    nhash[n7]->Check_List(n4);
    nhash[n7]->Check_List(n6);
  }

  // create unique edges
  ne=0;
  for (n=0; n < nn; n++)
    for (j=0; j < nhash[n]->max; j++)
      if (nhash[n]->list[j] > n)
        ne++;
  edge = new int[ne][3];
  ml = new double[ne];
  ne=0;
  for (n=0; n < nn; n++)
  {
    for (j=0; j < nhash[n]->max; j++)
    {
      if ((m=nhash[n]->list[j]) > n)
      {
        edge[ne][0] = n;
        edge[ne][1] = m;
        edge[ne][2] = -1;
        ne++;
      }
    }
  }
  for (n=0; n < nn; n++)
    nhash[n]->Redimension(0);

  int tenth, hundredth;

  // mark edges for refinement using Spacing Field

  int mark = 0;
  int newpts = 0;
  int nep;

  List *elist;
  elist = new List();
  int *tag;
  tag = new int[nn];
  for (n=0; n < nn; n++)
    tag[n] = 0;
  if (bnd_flag == 1)
  {
    for (b=0; b < nb; b++)
    {
      for (i=0; i < nt[b]; i++)
        for (j=0; j < 3; j++)
          tag[t_n[b][i][j]] = 1;
      for (i=0; i < nq[b]; i++)
        for (j=0; j < 4; j++)
          tag[q_n[b][i][j]] = 1;
    }
  }


  double mlocal, mglobal, mlmin, mlmax, mlavg;
  double p1[3], p2[3];
  mlmin = 1.0e20;
  mlmax = -1.0e20;
  mlavg = 0.0;
  int count = 0;
  for (p=0; p < num_procs; p++)
  {
    if (my_rank == p)
    {
      for (e=0; e < ne; e++)
      {
        ml[e] = 0.0;
        if (bnd_flag == 0)
          elist->Add_To_List(e);
        else if (tag[edge[e][0]] == 1 && tag[edge[e][1]] == 1)
          elist->Add_To_List(e);
      }
      nep = elist->max;
    }

    if (num_procs > 1)
    {
      MPI_Bcast(&nep, 1, MPI_INT, p, MPI_COMM_WORLD);

      if (nep > 0)
      {
        sdim = nep*6*sizeof(double);
        sendbuff= new char[sdim];
        recvbuff= new char[sdim];

        if (p == my_rank)
        {
          sposition=0;
          for (i=0; i < elist->max; i++)
          {
            e = elist->list[i];
            n = edge[e][0];
            m = edge[e][1];
            MPI_Pack(&(node[n][0]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(node[n][1]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(node[n][2]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(node[m][0]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(node[m][1]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(node[m][2]),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //now, send and recv packets from all procs
        if (p == my_rank)
          for (n = 0; n < num_procs; n++)
          {
            if (n == p) continue;

            MPI_Send(sendbuff,sdim,MPI_CHAR,n,p,MPI_COMM_WORLD);
          }
        if (p != my_rank)
          MPI_Recv(recvbuff,sdim,MPI_CHAR,p,p,MPI_COMM_WORLD,&single_status);
      }
    }

    tenth = (int)((double)nep/10.0+1);
    hundredth = MAX(1,(int)((double)nep/100.0));
    if (my_rank == 0)
    {
      fprintf(out_f,"\nProcessor %d, number of edges = %d",p,nep);
      fprintf(out_f,"\nComputing edge metric lengths. (dot every %d edges)\n",hundredth);
      fflush(out_f);
    }

    sposition=0;
    rposition=0;
    for (i=0; i < nep; i++)
    {
      //MPI_Barrier(MPI_COMM_WORLD);  // this is ensure the progress bar is realistic
      if (p == my_rank)
      {
        e = elist->list[i];
        n = edge[e][0];
        m = edge[e][1];
        p1[0] = node[n][0];
        p1[1] = node[n][1];
        p1[2] = node[n][2];
        p2[0] = node[m][0];
        p2[1] = node[m][1];
        p2[2] = node[m][2];
      } else
      {
        MPI_Unpack(recvbuff,sdim,&rposition,&(p1[0]),1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(recvbuff,sdim,&rposition,&(p1[1]),1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(recvbuff,sdim,&rposition,&(p1[2]),1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(recvbuff,sdim,&rposition,&(p2[0]),1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(recvbuff,sdim,&rposition,&(p2[1]),1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(recvbuff,sdim,&rposition,&(p2[2]),1,MPI_DOUBLE,MPI_COMM_WORLD);
      }

      mlocal=SF_edge_metric_length(p1,p2,mltype);

      if (p == my_rank)
        ml[e] = mlocal;
      else
        MPI_Pack(&(mlocal),1,MPI_DOUBLE,sendbuff,sdim,&sposition,MPI_COMM_WORLD);

      if (my_rank == 0 && (i+1) % hundredth == 0)
      {
        fprintf(out_f,".");
        fflush(out_f);
      }
      if (my_rank == 0 && ((i+1) % tenth == 0 || i == nep-1))
      {
        fprintf(out_f," %d\n",i+1);
        fflush(out_f);
      }
    }

    if (num_procs > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      for (n=0; n < num_procs; n++)
      {
        if (p == n)
          continue;

        MPI_Barrier(MPI_COMM_WORLD);

        if (my_rank == 0)
        {
          fprintf(out_f,"\nTransferring metric lengths from process %d to process %d",n,p);
          fflush(out_f);
        }

        if (n == my_rank) 
          MPI_Send(sendbuff,sdim,MPI_CHAR,p,n,MPI_COMM_WORLD);
        if (p == my_rank)
        {
          MPI_Recv(recvbuff,sdim,MPI_CHAR,n,n,MPI_COMM_WORLD, &single_status);

          rposition=0;
          for (i=0; i < nep; i++)
          {
            e = elist->list[i];
            MPI_Unpack(recvbuff,sdim,&rposition,&mlocal,1,MPI_DOUBLE,MPI_COMM_WORLD);
            if (mltype == 0)
              ml[e] = MAX(ml[e],mlocal);
            else
              ml[e] += mlocal;
          }
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (my_rank == 0)
      {
        fprintf(out_f,"\nTransfer of metric lengths to process %d completed.\n",p);
        fflush(out_f);
      }
    }

    if (p == my_rank)
    {
      for (e=0; e < ne; e++)
      {
        if (ml[e] > threshold)
          mark++;
        mlmin = MIN(mlmin,ml[e]);
        mlmax = MAX(mlmax,ml[e]);
        mlavg += ml[e];
        count++;
      }
    }

    if (num_procs > 1)
    {
      delete[] sendbuff;
      delete[] recvbuff;
      sendbuff = 0;
      recvbuff = 0;
    }
  }

  delete elist;
  delete[] tag;
  
  SF_Finalize();

  int gmark = mark;
  if (num_procs > 1)
  {
    MPI_Allreduce(&count,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    count = gmark;
    gmark=0;
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    double gdouble;
    MPI_Allreduce(&mlmin,&gdouble,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    mlmin = gdouble;
    MPI_Allreduce(&mlmax,&gdouble,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    mlmax = gdouble;
    MPI_Allreduce(&mlavg,&gdouble,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    mlavg = gdouble;
  }
  mlavg /= count;

  if (my_rank == 0)
  {
    fprintf(out_f,"\nMinimum metric length = %lg",mlmin);
    fprintf(out_f,"\nAverage metric length = %lg",mlavg);
    fprintf(out_f,"\nMaximum metric length = %lg",mlmax);
    fflush(out_f);
  }
  
  if (gmark > 0 && my_rank == 0)
  {
    fprintf(out_f,"\nOriginal GLOBAL number of edges marked (includes dups)= %d",gmark);
    fflush(out_f);
  }

  if (target < gmark)
  {
    if (num_procs == 1)
    {
      // sort edges and unmark lower valued metric length edges
      int *ed = new int[ne];
      for (e=0; e < ne; e++)
        ed[e] = e;

      sort_gt(ne,ed,ml);

      for (i=target; i < ne; i++)
        ml[ed[i]] = 0.0;

      delete[] ed;
    } else
    {
      // create buckets to count number of marked edges per bucket
      const int nbuckets = 1000;
      int bucket[nbuckets];
      double range[nbuckets+1];
      double div = (mlmax-mlmin)/nbuckets;
      i=nbuckets;
      range[i] = mlmax+1e-10;
      while (i > 0)
      {
        range[i-1] = range[i] - div;
        i--;
      }
      for (i=0; i < nbuckets; i++)
        bucket[i] = 0;
      for (e=0; e < ne; e++)
      {
        for (i=0; i < nbuckets; i++)
        {
          if (ml[e] >= range[i] && ml[e] <= range[i+1])
          {
            bucket[i]++;
            break;
          }
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);

      int count = 0;
      for (i=nbuckets-1; i >= 0 && count < target; i--)
      {
        int bcount;
        MPI_Allreduce(&bucket[i],&bcount,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        count += bcount;
        if (my_rank == 0)
        {
          fprintf(out_f,"\nBucket [%lg : %lg] contains %d entries, total = %d.",range[i],range[i+1],bcount,count);
          fflush(out_f);
        }
      }
      i = MAX(0,i);
      if (my_rank == 0)
      {
        fprintf(out_f,"\nNew threshold from buckets = %lg",range[i]);
        fflush(out_f);
      }
      for (e=0; e < ne; e++)
        if (ml[e] < range[i])
          ml[e] = 0.0;
    }
  }

  // now mark edges
  mark = newpts = 0;
  for (e=0; e < ne; e++)
  {
    n = edge[e][0];
    m = edge[e][1];
    if (ml[e] > threshold)
    {
      mark++;
      edge[e][2] = nn+newpts++;
    }
  }

  // create node to edge hash table
  for (e=0; e < ne; e++)
  {
    nhash[edge[e][0]]->Add_To_List(e);
    nhash[edge[e][1]]->Add_To_List(e);
  }

  // create element to edge arrays
  if (ntet > 0)
    tet_edges = new int[ntet][6];
  if (npyr > 0)
    pyr_edges = new int[npyr][8];
  if (npri > 0)
    pri_edges = new int[npri][9];
  if (nhex > 0)
    hex_edges = new int[nhex][12];

  //
  // create hash tables for elements
  //
  tethash = new List*[nnhash];
  pyrhash = new List*[nnhash];
  prihash = new List*[nnhash];
  hexhash = new List*[nnhash];
  for (n=0; n < nnhash; n++)
  {
    tethash[n] = new List();
    pyrhash[n] = new List();
    prihash[n] = new List();
    hexhash[n] = new List();
  }

  //
  // create node-to-tetrahedra hash list
  //
  for (c=0; c < ntet; c++)
    for (i=0; i < 4; i++)
      tethash[tet_n[c][i]]->Add_To_List(c);

  //
  // create node-to-pyramid hash list
  //
  for (c=0; c < npyr; c++)
    for (i=0; i < 5; i++)
      pyrhash[pyr_n[c][i]]->Add_To_List(c);

  //
  // create node-to-prism hash list
  //
  for (c=0; c < npri; c++)
    for (i=0; i < 6; i++)
      prihash[pri_n[c][i]]->Add_To_List(c);

  //
  // create node-to-hexahedra hash list
  //
  for (c=0; c < nhex; c++)
    for (i=0; i < 8; i++)
      hexhash[hex_n[c][i]]->Add_To_List(c);

  if (my_rank == 0)
  {
    fprintf(out_f,"\nNode-to-node and node-to-element hash tables created.");
    fflush(out_f);
  }

  // check current connectivity of mesh (THIS IS A DEBUGGING TOOL...the version at the end only works in serial)
  make_nabors(tethash,pyrhash,prihash,hexhash,1);


  // create tet-edge connectivity
  for (c=0; c < ntet; c++)
  {
    for (i=0; i < 6; i++)
    {
      switch (i)
      {
        case 0: n0 = tet_n[c][0]; n1 = tet_n[c][1]; break;
        case 1: n0 = tet_n[c][1]; n1 = tet_n[c][2]; break;
        case 2: n0 = tet_n[c][2]; n1 = tet_n[c][0]; break;
        case 3: n0 = tet_n[c][0]; n1 = tet_n[c][3]; break;
        case 4: n0 = tet_n[c][1]; n1 = tet_n[c][3]; break;
        case 5: n0 = tet_n[c][2]; n1 = tet_n[c][3]; break;
      }

      for (j=0; j < nhash[n0]->max; j++)
      {
        e = nhash[n0]->list[j];
        if (nhash[n1]->Is_In_List(e))
        {
          tet_edges[c][i] = e;
          break;
        }
      }
    }
  }

  // create pyr-edge connectivity
  for (c=0; c < npyr; c++)
  {
    for (i=0; i < 8; i++)
    {
      switch (i)
      {
        case 0: n0 = pyr_n[c][0]; n1 = pyr_n[c][1]; break;
        case 1: n0 = pyr_n[c][1]; n1 = pyr_n[c][2]; break;
        case 2: n0 = pyr_n[c][2]; n1 = pyr_n[c][3]; break;
        case 3: n0 = pyr_n[c][3]; n1 = pyr_n[c][0]; break;
        case 4: n0 = pyr_n[c][0]; n1 = pyr_n[c][4]; break;
        case 5: n0 = pyr_n[c][1]; n1 = pyr_n[c][4]; break;
        case 6: n0 = pyr_n[c][2]; n1 = pyr_n[c][4]; break;
        case 7: n0 = pyr_n[c][3]; n1 = pyr_n[c][4]; break;
      }

      for (j=0; j < nhash[n0]->max; j++)
      {
        e = nhash[n0]->list[j];
        if (nhash[n1]->Is_In_List(e))
        {
          pyr_edges[c][i] = e;
          break;
        }
      }
    }
  }

  // create pri-edge connectivity
  for (c=0; c < npri; c++)
  {
    for (i=0; i < 9; i++)
    {
      switch (i)
      {
        case 0: n0 = pri_n[c][0]; n1 = pri_n[c][1]; break;
        case 1: n0 = pri_n[c][1]; n1 = pri_n[c][2]; break;
        case 2: n0 = pri_n[c][2]; n1 = pri_n[c][0]; break;
        case 3: n0 = pri_n[c][0]; n1 = pri_n[c][3]; break;
        case 4: n0 = pri_n[c][1]; n1 = pri_n[c][4]; break;
        case 5: n0 = pri_n[c][2]; n1 = pri_n[c][5]; break;
        case 6: n0 = pri_n[c][3]; n1 = pri_n[c][4]; break;
        case 7: n0 = pri_n[c][4]; n1 = pri_n[c][5]; break;
        case 8: n0 = pri_n[c][5]; n1 = pri_n[c][3]; break;
      }

      for (j=0; j < nhash[n0]->max; j++)
      {
        e = nhash[n0]->list[j];
        if (nhash[n1]->Is_In_List(e))
        {
          pri_edges[c][i] = e;
          break;
        }
      }
    }
  }

  // create hex-edge connectivity
  for (c=0; c < nhex; c++)
  {
    for (i=0; i < 12; i++)
    {
      switch (i)
      {
        case  0: n0 = hex_n[c][0]; n1 = hex_n[c][1]; break;
        case  1: n0 = hex_n[c][1]; n1 = hex_n[c][2]; break;
        case  2: n0 = hex_n[c][2]; n1 = hex_n[c][3]; break;
        case  3: n0 = hex_n[c][3]; n1 = hex_n[c][0]; break;
        case  4: n0 = hex_n[c][0]; n1 = hex_n[c][4]; break;
        case  5: n0 = hex_n[c][1]; n1 = hex_n[c][5]; break;
        case  6: n0 = hex_n[c][2]; n1 = hex_n[c][6]; break;
        case  7: n0 = hex_n[c][3]; n1 = hex_n[c][7]; break;
        case  8: n0 = hex_n[c][4]; n1 = hex_n[c][5]; break;
        case  9: n0 = hex_n[c][5]; n1 = hex_n[c][6]; break;
        case 10: n0 = hex_n[c][6]; n1 = hex_n[c][7]; break;
        case 11: n0 = hex_n[c][7]; n1 = hex_n[c][4]; break;
      }

      for (j=0; j < nhash[n0]->max; j++)
      {
        e = nhash[n0]->list[j];
        if (nhash[n1]->Is_In_List(e))
        {
          hex_edges[c][i] = e;
          break;
        }
      }
    }
  }


  //
  // expand element connectivities to allow for mid-edge nodes
  //
  for (c=0; c < ntet; c++)
  {
    tet_n[c] = (int*)realloc((void*)tet_n[c],10*sizeof(int));
    for (i=4; i < 10; i++)
      tet_n[c][i] = -1;
  }
  for (c=0; c < npyr; c++)
  {
    pyr_n[c] = (int*)realloc((void*)pyr_n[c],15*sizeof(int));
    for (i=5; i < 15; i++)
      pyr_n[c][i] = -1;
  }
  for (c=0; c < npri; c++)
  {
    pri_n[c] = (int*)realloc((void*)pri_n[c],19*sizeof(int));
    for (i=6; i < 19; i++)
      pri_n[c][i] = -1;
  }
  for (c=0; c < nhex; c++)
  {
    hex_n[c] = (int*)realloc((void*)hex_n[c],27*sizeof(int));
    for (i=8; i < 27; i++)
      hex_n[c][i] = -1;
  }

  if (my_rank == 0)
  {
    fprintf(out_f,"\nElement connectivities expanded.");
    fflush(out_f);
  }

  // mark all edges in element connectivities
  for (c=0; c < ntet; c++)
    for (i=0; i < 6; i++)
      tet_n[c][i+4] = edge[tet_edges[c][i]][2];
  for (c=0; c < npyr; c++)
    for (i=0; i < 8; i++)
      pyr_n[c][i+5] = edge[pyr_edges[c][i]][2];
  for (c=0; c < npri; c++)
    for (i=0; i < 9; i++)
      pri_n[c][i+6] = edge[pri_edges[c][i]][2];
  for (c=0; c < nhex; c++)
    for (i=0; i < 12; i++)
      hex_n[c][i+8] = edge[hex_edges[c][i]][2];

  gmark=mark;
  int gpts=newpts;
  if (num_procs > 1)
  {
    gmark=0;
    gpts=0;
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&newpts,&gpts,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (gmark > 0 && my_rank == 0)
  {
    fprintf(out_f,"\nNumber of GLOBAL edges marked after reduction (includes dups)= %d",gmark);
    fprintf(out_f,"\nNumber of GLOBAL new points (includes dups)= %d",gpts);
    fflush(out_f);
  }

  delete[] ml;

  #ifdef _DEBUG
  fprintf(debug_f,"\nAbout to begin edge marking loop.\n");
  fflush(debug_f);
  #endif

  // loop to check for refinement marking constraints until no changes
  int e0, e1, e2, i0, i1, i2, m0, m1, m2;
  int tet_count, pyr_count, pri_count;
  int md[4], ed[4];
  int sweep = 0;
  total_changes = 1;
  changes = 1; //to be sure goes through once
  do
  {
    sweep++;
    changes = 0;
    total_changes = 0;

    // check triangle faces for 2 marked nodes and force refinement of 3rd
    if (tri_flag == 1)
    {
      for (c=0; c < ntet; c++)
      {
        for (k=0; k < 4; k++)
        {
          switch (k)
          {
            case 0: md[0] = 4; md[1] = 5; md[2] = 6; ed[0] = 0; ed[1] = 1; ed[2] = 2; break;
            case 1: md[0] = 4; md[1] = 7; md[2] = 8; ed[0] = 0; ed[1] = 3; ed[2] = 4; break;
            case 2: md[0] = 5; md[1] = 8; md[2] = 9; ed[0] = 1; ed[1] = 4; ed[2] = 5; break;
            case 3: md[0] = 6; md[1] = 7; md[2] = 9; ed[0] = 2; ed[1] = 3; ed[2] = 5; break;
          }
          j = 0;
          if (tet_n[c][md[0]] >= 0) j++;
          if (tet_n[c][md[1]] >= 0) j++;
          if (tet_n[c][md[2]] >= 0) j++;
          if (j == 2)
          {
            for (i=0; i < 3; i++)
            {
              if (tet_n[c][md[i]] < 0)
              {
                e = tet_edges[c][ed[i]];
                if (edge[e][2] < 0)
                {
                  edge[e][2] = nn+newpts++;
                  changes++;
                }
                tet_n[c][md[i]] = edge[e][2];
              }
            }
          }
        }
      }
      for (c=0; c < npyr; c++)
      {
        for (k=1; k < 5; k++)
        {
          switch (k)
          {
            case 1: md[0] = 5; md[1] =  9; md[2] = 10; ed[0] = 0; ed[1] = 4; ed[2] = 5; break;
            case 2: md[0] = 6; md[1] = 10; md[2] = 11; ed[0] = 1; ed[1] = 5; ed[2] = 6; break;
            case 3: md[0] = 7; md[1] = 11; md[2] = 12; ed[0] = 2; ed[1] = 6; ed[2] = 7; break;
            case 4: md[0] = 8; md[1] =  9; md[2] = 12; ed[0] = 3; ed[1] = 4; ed[2] = 7; break;
          }
          j = 0;
          if (pyr_n[c][md[0]] >= 0) j++;
          if (pyr_n[c][md[1]] >= 0) j++;
          if (pyr_n[c][md[2]] >= 0) j++;
          if (j == 2)
          {
            for (i=0; i < 3; i++)
            {
              if (pyr_n[c][md[i]] < 0)
              {
                e = pyr_edges[c][ed[i]];
                if (edge[e][2] < 0)
                {
                  edge[e][2] = nn+newpts++;
                  changes++;
                }
                pyr_n[c][md[i]] = edge[e][2];
              }
            }
          }
        }
      }
      for (c=0; c < npri; c++)
      {
        for (k=3; k < 5; k++)
        {
          switch (k)
          {
            case 3: md[0] =  6; md[1] =  7; md[2] =  8; ed[0] = 0; ed[1] = 1; ed[2] = 2; break;
            case 4: md[0] = 12; md[1] = 13; md[2] = 14; ed[0] = 6; ed[1] = 7; ed[2] = 8; break;
          }
          j = 0;
          if (pri_n[c][md[0]] >= 0) j++;
          if (pri_n[c][md[1]] >= 0) j++;
          if (pri_n[c][md[2]] >= 0) j++;
          if (j == 2)
          {
            for (i=0; i < 3; i++)
            {
              if (pri_n[c][md[i]] < 0)
              {
                e = pri_edges[c][ed[i]];
                if (edge[e][2] < 0)
                {
                  edge[e][2] = nn+newpts++;
                  changes++;
                }
                pri_n[c][md[i]] = edge[e][2];
              }
            }
          }
        }
      }
    }

    if (pri_stack == 1)
    {
      for (c=0; c < npri; c++)
      {
        for (i=6; i < 9; i++)
        {
          if ((pri_n[c][i] >= 0 && pri_n[c][i+6] < 0) || (pri_n[c][i] < 0 && pri_n[c][i+6] >= 0))
          {
            e = pri_edges[c][i-6];
            if (edge[e][2] < 0)
            {
              edge[e][2] = nn+newpts++;
              changes++;
            }
            pri_n[c][i] = edge[e][2];
            e = pri_edges[c][i];
            if (edge[e][2] < 0)
            {
              edge[e][2] = nn+newpts++;
              changes++;
            }
            pri_n[c][i+6] = edge[e][2];
          }
        }
      }
    }

    // mark all edges in element connectivities
    for (c=0; c < ntet; c++)
    {
      for (i=0; i < 6; i++)
      {
        if (tet_n[c][i+4] != edge[tet_edges[c][i]][2])
        {
          if (tet_n[c][i+4] >= 0 || edge[tet_edges[c][i]][2] < 0)
          {
            fprintf(stderr,"\nProcess %d: Edge and tet markings inconsistent.",my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
          tet_n[c][i+4] = edge[tet_edges[c][i]][2];
          changes++;
        }
      }
    }
    for (c=0; c < npyr; c++)
    {
      for (i=0; i < 8; i++)
      {
        if (pyr_n[c][i+5] != edge[pyr_edges[c][i]][2])
        {
          if (pyr_n[c][i+5] >= 0 || edge[pyr_edges[c][i]][2] < 0)
          {
            fprintf(stderr,"\nProcess %d: Edge and pyramid markings inconsistent.",my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
          pyr_n[c][i+5] = edge[pyr_edges[c][i]][2];
          changes++;
        }
      }
    }
    for (c=0; c < npri; c++)
    {
      for (i=0; i < 9; i++)
      {
        if (pri_n[c][i+6] != edge[pri_edges[c][i]][2])
        {
          if (pri_n[c][i+6] >= 0 || edge[pri_edges[c][i]][2] < 0)
          {
            fprintf(stderr,"\nProcess %d: Edge and prism markings inconsistent.",my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
          pri_n[c][i+6] = edge[pri_edges[c][i]][2];
          changes++;
        }
      }
    }
    for (c=0; c < nhex; c++)
    {
      for (i=0; i < 12; i++)
      {
        if (hex_n[c][i+8] != edge[hex_edges[c][i]][2])
        {
          if (hex_n[c][i+8] >= 0 || edge[hex_edges[c][i]][2] < 0)
          {
            fprintf(stderr,"\nProcess %d: Edge and hex markings inconsistent.",my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
          hex_n[c][i+8] = edge[hex_edges[c][i]][2];
          changes++;
        }
      }
    }

    // mark mid face nodes where appropriate
    for (c=0; c < npyr; c++)
    {
      n0 = pyr_n[c][0];
      n1 = pyr_n[c][3];
      n2 = pyr_n[c][2];
      n3 = pyr_n[c][1];
      if ((pyr_n[c][5] < 0 || pyr_n[c][6] < 0 || pyr_n[c][7] < 0 || pyr_n[c][8] < 0) && pyr_n[c][13] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE CENTROID FOR PYRAMID %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pyr_n[c][5] >= 0 && pyr_n[c][6] >= 0 && pyr_n[c][7] >= 0 && pyr_n[c][8] >= 0 && pyr_n[c][13] < 0)
      {
        pyr_n[c][13] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lpyr=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
    }
    for (c=0; c < npri; c++)
    {
      n0 = pri_n[c][0];
      n1 = pri_n[c][1];
      n2 = pri_n[c][4];
      n3 = pri_n[c][3];
      if ((pri_n[c][6] < 0 || pri_n[c][9] < 0 || pri_n[c][12] < 0 || pri_n[c][10] < 0) && pri_n[c][15] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 0 CENTROID FOR PRISM %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pri_n[c][6] >= 0 && pri_n[c][9] >= 0 && pri_n[c][12] >= 0 && pri_n[c][10] >= 0 && pri_n[c][15] < 0)
      {
        pri_n[c][15] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lpri=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = pri_n[c][1];
      n1 = pri_n[c][2];
      n2 = pri_n[c][5];
      n3 = pri_n[c][4];
      if ((pri_n[c][7] < 0 || pri_n[c][10] < 0 || pri_n[c][13] < 0 || pri_n[c][11] < 0) && pri_n[c][16] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 1 CENTROID FOR PRISM %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pri_n[c][7] >= 0 && pri_n[c][10] >= 0 && pri_n[c][13] >= 0 && pri_n[c][11] >= 0 && pri_n[c][16] < 0)
      {
        pri_n[c][16] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lpri=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = pri_n[c][2];
      n1 = pri_n[c][0];
      n2 = pri_n[c][3];
      n3 = pri_n[c][5];
      if ((pri_n[c][8] < 0 || pri_n[c][11] < 0 || pri_n[c][14] < 0 || pri_n[c][9] < 0) && pri_n[c][17] >= 0)
      {
         fprintf(stderr,"\nProcess %d: ILLEGAL FACE 2 CENTROID FOR PRISM %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pri_n[c][8] >= 0 && pri_n[c][11] >= 0 && pri_n[c][14] >= 0 && pri_n[c][9] >= 0 && pri_n[c][17] < 0)
      {
        pri_n[c][17] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lpri=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
    }
    for (c=0; c < nhex; c++)
    {
      n0 = hex_n[c][0];
      n1 = hex_n[c][3];
      n2 = hex_n[c][2];
      n3 = hex_n[c][1];
      if ((hex_n[c][8] < 0 || hex_n[c][9] < 0 || hex_n[c][10] < 0 || hex_n[c][11] < 0) && hex_n[c][20] >= 0)
      {
       fprintf(stderr,"\nProcess %d: ILLEGAL FACE 0 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (hex_n[c][8] >= 0 && hex_n[c][9] >= 0 && hex_n[c][10] >= 0 && hex_n[c][11] >= 0 && hex_n[c][20] < 0)
      {
        hex_n[c][20] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = hex_n[c][0];
      n1 = hex_n[c][1];
      n2 = hex_n[c][5];
      n3 = hex_n[c][4];
      if ((hex_n[c][8] < 0 || hex_n[c][12] < 0 || hex_n[c][16] < 0 || hex_n[c][13] < 0) && hex_n[c][21] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 1 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);

      }
      if (hex_n[c][8] >= 0 && hex_n[c][12] >= 0 && hex_n[c][16] >= 0 && hex_n[c][13] >= 0 && hex_n[c][21] < 0)
      {
        hex_n[c][21] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = hex_n[c][1];
      n1 = hex_n[c][2];
      n2 = hex_n[c][6];
      n3 = hex_n[c][5];
      if ((hex_n[c][9] < 0 || hex_n[c][13] < 0 || hex_n[c][17] < 0 || hex_n[c][14] < 0) && hex_n[c][22] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 2 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (hex_n[c][9] >= 0 && hex_n[c][13] >= 0 && hex_n[c][17] >= 0 && hex_n[c][14] >= 0 && hex_n[c][22] < 0)
      {
        hex_n[c][22] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = hex_n[c][2];
      n1 = hex_n[c][3];
      n2 = hex_n[c][7];
      n3 = hex_n[c][6];
      if ((hex_n[c][10] < 0 || hex_n[c][14] < 0 || hex_n[c][18] < 0 || hex_n[c][15] < 0) && hex_n[c][23] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 3 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (hex_n[c][10] >= 0 && hex_n[c][14] >= 0 && hex_n[c][18] >= 0 && hex_n[c][15] >= 0 && hex_n[c][23] < 0)
      {
        hex_n[c][23] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = hex_n[c][0];
      n1 = hex_n[c][4];
      n2 = hex_n[c][7];
      n3 = hex_n[c][3];
      if ((hex_n[c][11] < 0 || hex_n[c][15] < 0 || hex_n[c][19] < 0 || hex_n[c][12] < 0) && hex_n[c][24] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 4 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (hex_n[c][11] >= 0 && hex_n[c][15] >= 0 && hex_n[c][19] >= 0 && hex_n[c][12] >= 0 && hex_n[c][24] < 0)
      {
        hex_n[c][24] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
      n0 = hex_n[c][4];
      n1 = hex_n[c][5];
      n2 = hex_n[c][6];
      n3 = hex_n[c][7];
      if ((hex_n[c][16] < 0 || hex_n[c][17] < 0 || hex_n[c][18] < 0 || hex_n[c][19] < 0) && hex_n[c][25] >= 0)
      {
        fprintf(stderr,"\nProcess %d: ILLEGAL FACE 5 CENTROID FOR HEX %d",my_rank,c);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (hex_n[c][16] >= 0 && hex_n[c][17] >= 0 && hex_n[c][18] >= 0 && hex_n[c][19] >= 0 && hex_n[c][25] < 0)
      {
        hex_n[c][25] = nn+newpts;

        // notify neighbor element
        ltet=lpyr=lpri=lhex=-1;
        lhex=c;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,lpyr,-1,-1,pyr_n[lpyr]);
        if (lpri >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,lpri,-1,pri_n[lpri]);
        if (lhex >= 0) store_face_node(n0,n1,n2,n3,nn+newpts,-1,-1,lhex,hex_n[lhex]);

        newpts++;
        changes++;
      }
    }
    for (c=0; c < nhex; c++)
    {
      if (hex_n[c][20] >= 0 && hex_n[c][21] >= 0 && hex_n[c][22] >= 0 &&
          hex_n[c][23] >= 0 && hex_n[c][24] >= 0 && hex_n[c][25] >= 0 && hex_n[c][26] < 0)
      {
        hex_n[c][26] = nn+newpts++;
        changes++;
      }
    }
    gmark=newpts;
    if (num_procs > 1)
    {
      gmark=0;
      MPI_Allreduce(&newpts,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }
    if (gmark > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nFace sweep %d, global # of new points = %d",sweep,gmark);
      fflush(out_f);
    }
   
    //check after all marking done that it is consistent between procs (re sweep if not, to be sure all uniform)
    
    total_changes = changes;

    if (num_procs > 1)
    {
      sendcnt = new int [num_procs];
      recvcnt = new int [num_procs];
      for (n = 0; n < num_procs; n++)
      {
        sendcnt[n] = 0;
        recvcnt[n] = 0;
      }
      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];
      bdim = new int[num_procs];
      sbuff = new char*[num_procs];
      rbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
      {  
        bdim[n] = 0;
        sbuff[n] = 0;
        rbuff[n] = 0;
      }
       
      // exchange edge markings between processors
      // send all edges with nodes on other processors to compare refinement
      for (n=0; n < nn; n++)
      {
        if ((p=pmap[n][1]) == my_rank)
          continue;
        sendcnt[p] += nhash[n]->max;
      }

      //now, send and recv count for all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;

        MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
 
      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //allocate buffers to send/recv node indices
      for (n = 0; n < num_procs; n++)
        bdim[n] = MAX(sendcnt[n],recvcnt[n])*9*sizeof(int);
      for (n = 0; n < num_procs; n++)
      {
        sbuff[n] = new char[bdim[n]];
        rbuff[n] = new char[bdim[n]];
      }

      //now, pack node indices
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;

        //set position to 0
        sposition=0;
 
        //pack node indices needed
        for (n = 0; n < nn; n++)
        {
          if (pmap[n][1] == p)
          {
            for (i=0; i < nhash[n]->max; i++)
            {
              e=nhash[n]->list[i];
              n0 = edge[e][0];
              n1 = edge[e][1];
              n2 = edge[e][2];

              MPI_Pack(&n0,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n0][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n0][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n0][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&n1,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n1][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n1][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&(pmap[n1][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              MPI_Pack(&n2,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            }
          }
        }
      }

      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }

        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      int m0, m0g, m0p, m0i, m1, m1g, m1p, m1i, m2;

      //now, unpack compare and repack for return
      for (p = 0; p < num_procs; p++)
      {
        sposition = rposition = 0; //reset pos
        for (n = 0; n < recvcnt[p]; n++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m0,1,MPI_INT,MPI_COMM_WORLD); 
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m0g,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m0p,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m0i,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m1,1,MPI_INT,MPI_COMM_WORLD); 
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m1g,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m1p,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m1i,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&m2,1,MPI_INT,MPI_COMM_WORLD); 

          // send back end nodes & middle marking
          MPI_Pack(&m0,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          MPI_Pack(&m1,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);

          if (m0p == my_rank)
          {
            m=m0i;
          } else if (m1p == my_rank)
          {
            m=m1i;
          } else
          {
            fprintf(stderr,"\nEdge nodes passed to proc %d not owned by proc!",my_rank);
            fprintf(stderr,"\nnode 0: %d, %d, %d, %d",m0,m0g,m0p,m0i,my_rank);
            fprintf(stderr,"\nnode 1: %d, %d, %d, %d",m1,m1g,m1p,m1i,my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
          for (k=0,i=0; i < nhash[m]->max && !k; i++)
          {
            e = nhash[m]->list[i];
            n0 = edge[e][0];
            n1 = edge[e][1];
            n2 = edge[e][2];
            if ((pmap[n0][0] == m0g && pmap[n1][0] == m1g) ||
                (pmap[n0][0] == m1g && pmap[n1][0] == m0g))
            {
              if (n2 < 0 && m2 >= 0)
              {
                edge[e][2] = nn+newpts++;
                changes++;
              } else if (m2 < 0 && n2 >= 0)
              {
                m2 = n2;
                changes++;
              }
              MPI_Pack(&m2,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              k=1;
            }
          }
        }
      }

      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }

        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //now, unpack return messages
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;

        //set position to 0
        rposition=0;
 
        for (n = 0; n < nn; n++)
        {
          if (pmap[n][1] == p)
          {
            for (i=0; i < nhash[n]->max; i++)
            {
              e=nhash[n]->list[i];
              n0 = edge[e][0];
              n1 = edge[e][1];
              n2 = edge[e][2];
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&m0,1,MPI_INT,MPI_COMM_WORLD); 
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&m1,1,MPI_INT,MPI_COMM_WORLD); 
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&m2,1,MPI_INT,MPI_COMM_WORLD); 
              if (n0 != m0 || n1 != m1)
              {
                fprintf(stderr,"\nReturn edge nodes from to proc %d not correct!",p);
                fprintf(stderr,"\nnode 0 = %d, %d",n0,m0);
                fprintf(stderr,"\nnode 1 = %d, %d",n1,m1);
                fprintf(stderr,"\nnode 2 = %d, %d",n2,m2);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              if (n2 < 0 && m2 >= 0)
              {
                changes++;
                edge[e][2] = nn+newpts++;
              }
            }
          }
        }
      }
     
      delete[] sendcnt;
      delete[] recvcnt;
      delete[] srequest;
      delete[] rrequest;
      delete[] statuses;
      for (n = 0; n < num_procs; n++)
      {
        if (sbuff[n] != 0) delete [] sbuff[n];
        if (rbuff[n] != 0) delete [] rbuff[n];
      }
      delete[] sbuff;
      delete[] rbuff;
      delete[] bdim;
    }
    
    //cull together changes...in serial, use only changes flag, else in parallel, use total_changes to be sure if any proc behind, willl catch up
    total_changes = 0;
    MPI_Allreduce(&changes,&total_changes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  } while (total_changes > 0);

  #ifdef _DEBUG
  fprintf(debug_f,"\nFinished edge marking loop.\n");
  fflush(debug_f);
  #endif

  // create centroid node where needed...will follow checked edge refinement to be same on all procs 
  // while fp calc might differ, the owning proc will write out final node coords and settle
  //more importantly, pmapping is done after pmap is expanded, later
  sweep=0;
  do
  {
    sweep++;
    changes = 0;
    // mark centroid nodes for pyramid
    for (c=0; c < npyr; c++)
    {
      // check for centroid needed for some pyramid refine cases
      tester = -1;
      if (!refine_pyr(pyr_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester))
      {
        pyr_n[c][14] = nn+newpts++;
        changes++;
      }
    }
    // mark centroid nodes for prisms
    for (c=0; c < npri; c++)
    {
      // check for centroid needed for some prism refine cases
      tester = -1;
      if (!refine_pri(pri_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester))
      {
        pri_n[c][18] = nn+newpts++;
        changes++;
      }
    }
    // mark centroid nodes for hexes
    for (c=0; c < nhex; c++)
    {
      // check for centroid needed for some hex refine cases
      tester = -1;
      if (!refine_hex(hex_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester))
      {
        hex_n[c][26] = nn+newpts++;
        changes++;
      }
    }
    gmark=changes;
    if (num_procs > 1)
    {
      gmark=0;
      MPI_Allreduce(&changes,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }
    if (gmark > 0 && my_rank == 0)
    {
      fprintf(out_f,"\nCentroid sweep %d, number of centroids marked = %d",sweep,gmark);
      fflush(out_f);
    }
  } while (gmark > 0);
  
  gmark=newpts;
  if (num_procs > 1)
  {
    gmark=0;
    MPI_Allreduce(&newpts,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  //since we need element to code correspondence, we alloc for all elem, init to -1
  long int *optet = 0, *oppyr = 0, *oppri = 0, *ophex = 0;
  long int **optri = 0, **opquad = 0;

  optet = (long int*)calloc(ntet,sizeof(long int));
  oppyr = (long int*)calloc(npyr,sizeof(long int));
  oppri = (long int*)calloc(npri,sizeof(long int));
  ophex = (long int*)calloc(nhex,sizeof(long int));

  optri = (long int**)calloc(nb,sizeof(long int*));
  opquad = (long int**)calloc(nb,sizeof(long int*));
  for (b = 0; b < nb; b++)
  {
    optri[b] = (long int*)calloc(nt[b],sizeof(long int));
    opquad[b] = (long int*)calloc(nq[b],sizeof(long int));
  }
  
  //init
  for (t = 0; t < ntet; t++)
    optet[t] = -1;
  for (t = 0; t < npyr; t++)
    oppyr[t] = -1;
  for (t = 0; t < npri; t++)
    oppri[t] = -1;
  for (t = 0; t < nhex; t++)
    ophex[t] = -1;

  for (b = 0; b < nb; b++)
  {
    for (t = 0; t < nt[b]; t++)
    {
      optri[b][t] = -1;
    }
    for (t = 0; t < nq[b]; t++)
    {
      opquad[b][t] = -1;
    }
  }

  // update boundaries
  if (gmark > 0)
  {

    #ifdef _DEBUG
    fprintf(debug_f,"\nAbout to begin boundary loop.\n");
    fflush(debug_f);
    #endif
  
    int ntnew, nqnew;
    // count new triangles and quads per boundary
    for (b=0; b < nb; b++)
    {
      total_changes = 1;
      pass = 0;
      do
      {
        changes = 0;
        total_changes = 0;
        pass++;
        ntnew = nqnew = 0;
        for (t=0; t < nt[b]; t++)
        {
          n0 = t_n[b][t][0];
          n1 = t_n[b][t][1];
          n2 = t_n[b][t][2];
          n3 = -1;
          ltet=lpyr=lpri=lhex=-1;
          neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
          if (num_procs > 1)
          {
            if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
            {
              double cg[3];
              for (i=0; i < 3; i++)
                cg[i] = (node[n0][i]+node[n1][i]+node[n2][i])/3.0;
              //VCB: debug only to see all triangles with issues          
              //fprintf(debug_f,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR TRIANGLE %d (%d-%d)!",my_rank,b,t,tri_emap[b][t][0],tri_emap[b][t][1]);
              //fprintf(debug_f,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
              //for (i=0; i < 3; i++)
                //fprintf(debug_f,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",t_n[b][t][i],pmap[t_n[b][t][i]][0],pmap[t_n[b][t][i]][1],pmap[t_n[b][t][i]][2],node[t_n[b][t][i]][0],node[t_n[b][t][i]][1],node[t_n[b][t][i]][2]);
              //fflush(debug_f);
              fprintf(stderr,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR TRIANGLE %d (%d-%d)!",my_rank,b,t,tri_emap[b][t][0],tri_emap[b][t][1]);
              fprintf(stderr,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
              for (i=0; i < 3; i++)
                fprintf(stderr,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",t_n[b][t][i],pmap[t_n[b][t][i]][0],pmap[t_n[b][t][i]][1],pmap[t_n[b][t][i]][2],node[t_n[b][t][i]][0],node[t_n[b][t][i]][1],node[t_n[b][t][i]][2]);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
          }
          if (num_procs == 1)
          {
            if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
            {
              fprintf(stderr,"\nProcess %d: NO MATCHING VOLUME ELEMENT FOUND FOR TRIANGLE!",my_rank);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
          }
          if (ltet >= 0) retrieve_face_nodes(n0,n1,n2,n3,ltet,-1,-1,-1,tet_n[ltet],conn);
          if (lpyr >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
          if (lpri >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
          tester = optri[b][t];
          refine_tri(conn,cdim,new_t,tri,new_q,quad,tester);
          //create list of refinements on other procs
          //do for all for comm
          optri[b][t] = tester;
          ntnew += MAX(0,new_t-1);
          nqnew += new_q;
        }
        for (q=0; q < nq[b]; q++)
        {
          n0 = q_n[b][q][0];
          n1 = q_n[b][q][1];
          n2 = q_n[b][q][2];
          n3 = q_n[b][q][3];
          ltet=lpyr=lpri=lhex=-1;
          neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
          if (num_procs > 1)
          {
            if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
            {
              double cg[3];
              for (i=0; i < 3; i++)
                cg[i] = (node[n0][i]+node[n1][i]+node[n2][i]+node[n3][i])/4.0;
              
              //VCB: debug only to see all quads with issues          
              //fprintf(debug_f,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR QUADRILATERAL %d (%d-%d)!",my_rank,b,q,quad_emap[b][q][0],quad_emap[b][q][1]);
              //fprintf(debug_f,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
              //for (i=0; i < 4; i++)
                //fprintf(debug_f,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",q_n[b][q][i],pmap[q_n[b][q][i]][0],pmap[q_n[b][q][i]][1],pmap[q_n[b][q][i]][2],node[q_n[b][q][i]][0],node[q_n[b][q][i]][1],node[q_n[b][q][i]][2]);
              //fflush(debug_f);
              fprintf(stderr,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR QUADRILATERAL %d (%d-%d)!",my_rank,b,q,quad_emap[b][q][0],quad_emap[b][q][1]);
              fprintf(stderr,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
              for (i=0; i < 4; i++)
                fprintf(stderr,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",q_n[b][q][i],pmap[q_n[b][q][i]][0],pmap[q_n[b][q][i]][1],pmap[q_n[b][q][i]][2],node[q_n[b][q][i]][0],node[q_n[b][q][i]][1],node[q_n[b][q][i]][2]);
              fflush(stderr);
              //look to see if two triangles by nodes
              List trifind;
              trifind.Redimension(0);
              trifind.Add_To_List(n0);
              trifind.Add_To_List(n1);
              trifind.Add_To_List(n2);
              trifind.Add_To_List(n3);
              m = 0;
              for (i = 0; i < nt[b] && m < 2; i++)
              {
                k = 0; //cntr
                for (j = 0; j < 3; j++)
                {
                  if (trifind.Is_In_List(t_n[b][i][j]))
                    k++;
                }
                if (k == 3)
                {
                  m++;
                }
              }
              if (m == 2)
              {
                fprintf(stderr,"\nNODALLY FOUND TWO TRIANGLES COVERING QUADRILATERAL %d (%d-%d)!",q,quad_emap[b][q][0],quad_emap[b][q][1]);
                fflush(stderr);
              }
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
          }
          if (num_procs == 1)
          {
            if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
            {
              fprintf(stderr,"\nProcess %d: NO MATCHING VOLUME ELEMENT FOUND FOR QUADRILATERAL!",my_rank);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
          }
          if (lpyr >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
          if (lpri >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
          if (lhex >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,-1,lhex,hex_n[lhex],conn);
          tester = opquad[b][q];
          refine_quad(conn,cdim,new_t,tri,new_q,quad,tester);
          //create list of refinements on other procs
          //do for all for comm
          opquad[b][q] = tester;
          ntnew += new_t;
          nqnew += MAX(0,new_q-1);
        }
        
        //now, test markings
        //communicate markings to other procs using emaps, if in disagreement, mark changes++
        //if there is a disagreement, this will fix it (owning proc decides), rerun with proper mask, and then get elem count right
        if (num_procs > 1)
        {
          //first, set up mem 
          sendcnt = new int [num_procs];
          recvcnt = new int [num_procs];
          for (n = 0; n < num_procs; n++)
          {
            sendcnt[n] = 0;
            recvcnt[n] = 0;
          }
          srequest = new MPI_Request[num_procs];
          rrequest = new MPI_Request[num_procs];
          statuses = new MPI_Status[num_procs];

          //cycle thru elements of each type	
          for (n = 0; n < nt[b]; n++)
          {
            if ((p = tri_emap[b][n][0]) != my_rank)
              sendcnt[p]++;
          }
          for (n = 0; n < nq[b]; n++)
          {
            if ((p = quad_emap[b][n][0]) != my_rank)
              sendcnt[p]++;
          }

          #ifdef _DEBUG
          //for (n = 0 ; n < num_procs; n++)
            //fprintf(debug_f,"\n sendcnt[%d] = %d \n",n,sendcnt[n]);
          //fprintf(debug_f,"\nPrepared send counts.\n");
          //fflush(debug_f);
          #endif

          //now, send and recv count for all procs
          nreq_s = nreq_r = 0;
          for (p = 0; p < num_procs; p++)
          {
            if (p == my_rank)
              continue;

            MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
 
            MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }

          //now, wait for finish before unpacking
          MPI_Waitall(nreq_s,srequest,statuses);
          MPI_Waitall(nreq_r,rrequest,statuses);

          #ifdef _DEBUG
          //for (n = 0 ; n < num_procs; n++)
            //fprintf(debug_f,"\n recvcnt[%d] = %d \n",n,recvcnt[n]);
          //fprintf(debug_f,"\nSent counts.\n");
          //fflush(debug_f);
          #endif
 
          delete [] srequest;
          delete [] rrequest;
          delete [] statuses;
 
          srequest = new MPI_Request[num_procs];
          rrequest = new MPI_Request[num_procs];
          statuses = new MPI_Status[num_procs];

          //allocate buffers to send/recv node indices
          bdim2 = new int[num_procs];
          for (n = 0; n < num_procs; n++)
            bdim2[n] = MAX(sendcnt[n],recvcnt[n])*(3*sizeof(int) + 1*sizeof(long int)); //indices
          sbuff = new char*[num_procs];
          rbuff = new char*[num_procs];
          for (n = 0; n < num_procs; n++)
          {
            sbuff[n] = new char[bdim2[n]];
            rbuff[n] = new char[bdim2[n]];
          }
      
          //set up types
          int tritype = -100;
          int quadtype = -200;

          //now, pack node indices as per procs that will be sending them
          for (p = 0; p < num_procs; p++)
          {
            if (sendcnt[p] == 0 || p == my_rank)
              continue;

            //set position to 0
            sposition=0;

            //cycle thru elements of each type	
            for (c = 0; c < nt[b]; c++)
            {
              if (tri_emap[b][c][0] == p)
              {
                MPI_Pack(&(tritype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(tri_emap[b][c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(optri[b][c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
                MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }	
            for (c = 0; c < nq[b]; c++)
            {
              if (quad_emap[b][c][0] == p)
              {
                MPI_Pack(&(quadtype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(quad_emap[b][c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(opquad[b][c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
                MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }

          //now, send and recv packets from all procs
          nreq_s = nreq_r = 0;
          for (p = 0; p < num_procs; p++)
          {
            if (sendcnt[p] > 0 && p != my_rank) 
            {
              MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
              nreq_s++;
            }

            if (recvcnt[p] > 0 && p != my_rank) 
            {
              MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
              nreq_r++;
            }
          }

          //now, wait for finish before unpacking
          MPI_Waitall(nreq_s,srequest,statuses);
          MPI_Waitall(nreq_r,rrequest,statuses);

          #ifdef _DEBUG
          //fprintf(debug_f,"\nPacked and sent markings.\n");
          //fflush(debug_f);
          #endif
 
          delete [] srequest;
          delete [] rrequest;
          delete [] statuses;

          srequest = new MPI_Request[num_procs];
          rrequest = new MPI_Request[num_procs];
          statuses = new MPI_Status[num_procs];
 
          //now, delete sbuff and resize
          for (n = 0; n < num_procs; n++)
            if (sbuff[n] != 0) delete [] sbuff[n];
          delete [] sbuff;
 
          //resize
          bdim21 = new int[num_procs];
          for (n = 0; n < num_procs; n++)
            bdim21[n] = MAX(sendcnt[n],recvcnt[n])*(2*sizeof(int) + 1*sizeof(long int)); //marking and indices
          sbuff = new char*[num_procs];
          for (n = 0; n < num_procs; n++)
            sbuff[n] = new char[bdim21[n]];
 
          temp_node = temp_node2 = temp_node3 = 0;
          temp_node4 = temp_node5 = 0;
 
          //now, unpack type, index, mark, replace with type, op index, mark
          for (p = 0; p < num_procs; p++)
          {
            if (recvcnt[p] == 0 || p == my_rank)
              continue;
             
            sposition = rposition = 0; //reset pos
            for (n = 0; n < recvcnt[p]; n++)
            {
              //unpack type initially
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

              //switch
              switch(temp_node)
              {
                case -100:
                  //unpack index on proc
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack marking
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);
          
                  //now, find owning proc marking
                  //no change version
                  if (optri[b][temp_node2] == temp_node4)
                  {
                    //repack type
                    MPI_Pack(&(tritype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //finally, unpack other proc index into unused temp_node2
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                    //pack that immediately back for other proc
                    MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //pack back marking
                    MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  }
                  else if (optri[b][temp_node2] != temp_node4) //change version
                  {
                    //repack type
                    MPI_Pack(&(tritype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //finally, unpack other proc index into unused temp_node3
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                    //pack that immediately back for other proc
                    MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //pack new marking, as refinement of both sorts
                    temp_node5 = (optri[b][temp_node2]) | (temp_node4);
            
                    #ifdef _DEBUG
                    fprintf(debug_f,"\nTRI Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,optri[b][temp_node2],temp_node5);
                    fflush(debug_f);
                    #endif
                    //set on this proc too
                    optri[b][temp_node2] = temp_node5;
                    //pack back
                    MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //set up changes
                    changes++;
                  }
                  break;
                case -200:
                  //unpack index on proc
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack marking
                  MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);

                  //now, find owning proc marking
                  //no change version
                  if (opquad[b][temp_node2] == temp_node4)
                  {
                    //repack type
                    MPI_Pack(&(quadtype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //finally, unpack other proc index into unused temp_node2
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                    //pack that immediately back for other proc
                    MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //pack back marking
                    MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                  }
                  else if (opquad[b][temp_node2] != temp_node4) //change version
                  {
                    //repack type
                    MPI_Pack(&(quadtype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //finally, unpack other proc index into unused temp_node3
                    MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                    //pack that immediately back for other proc
                    MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //pack new marking, as refinement of both sorts
                    temp_node5 = (opquad[b][temp_node2]) | (temp_node4);
                    #ifdef _DEBUG
                    fprintf(debug_f,"\nQUAD Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,opquad[b][temp_node2],temp_node5);
                    fflush(debug_f);
                    #endif    
                    //set on this proc too
                    opquad[b][temp_node2] = temp_node5;
                    //pack back      
                    MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                    //set up changes
                    changes++;
                  }
                  break;
                default:
                  fprintf(stderr,"\nSomehow you have more than two boundary element types based on a bad index/misordering!  Exiting....\n");
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);  
                  break;
              }
            }
          }
          
          #ifdef _DEBUG
          //fprintf(debug_f,"\nFinished unpack, repack.\n");
          //fflush(debug_f);
          #endif

          //now, delete rbuff and bdim2, resize
          delete [] bdim2;
          for (n = 0; n < num_procs; n++)
            if (rbuff[n] != 0) delete [] rbuff[n];
          delete [] rbuff;

          //resize
          rbuff = new char*[num_procs];
          for (n = 0; n < num_procs; n++)
            rbuff[n] = new char[bdim21[n]];

          //now, send and recv packets from all procs
          nreq_s = nreq_r = 0;
          for (p = 0; p < num_procs; p++)
          {
            if (recvcnt[p] > 0 && p != my_rank) 
            {
              MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
              nreq_s++;
            }
  
            if (sendcnt[p] > 0 && p != my_rank) 
            {
              MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
              nreq_r++;
            }
          }
  
          //now, wait for finish before unpacking
          MPI_Waitall(nreq_s,srequest,statuses);
          MPI_Waitall(nreq_r,rrequest,statuses);

          #ifdef _DEBUG
          //fprintf(debug_f,"\nUnpacking bd elem.\n");
          //fflush(debug_f);
          #endif

          //finally, unpack elements and markings
          for (p = 0; p < num_procs; p++)
          {
            if (sendcnt[p] == 0 || p == my_rank)
              continue;
          
            rposition = 0; //reset pos
            for (n = 0; n < sendcnt[p]; n++)
            {
              //unpack type initially
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

              //switch
              switch(temp_node)
              {
                case -100:
                  //unpack index use
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack marking
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(optri[b][temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
                  //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,optet[temp_node2]);
                  //fflush(out_f);
                  break;
                case -200:
                  //unpack index use
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                  //unpack marking
                  MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(opquad[b][temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
                  //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,oppyr[temp_node2]);
                  //fflush(out_f);
                  break;
                default:
                  fprintf(stderr,"\nSomehow you have more than two boundary element types based on a bad index/misordering!  Exiting....\n");
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);  
                  break;
              }
            }
          }

          //finally, free MPI mem
          delete[] sendcnt;
          delete[] recvcnt;
          delete[] srequest;
          delete[] rrequest;
          delete[] statuses;
          for (n = 0; n < num_procs; n++)
          {
            if (sbuff[n] != 0) delete [] sbuff[n];
            if (rbuff[n] != 0) delete [] rbuff[n];
          }
          delete[] sbuff;
          delete[] rbuff;
          delete[] bdim21;
        }
      
        //now, assuming that we have made all changes, it will go back thru, agree, and move on (if all metrics observed, but it will eventually stop)  
        MPI_Allreduce(&changes,&total_changes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
        if (my_rank == 0)
        {
          fprintf(out_f,"\nPass %d, number of global changes = %d", pass, total_changes);
          fflush(out_f);
        }
  
      } while (total_changes > 0 && pass < 1000); //uses do loop since will affect element count
  
      if (pass == 1000)
      {
        fprintf(stderr,"\nProc %d has made too many passes at bd elem refinement codes.  Exiting....\n",my_rank);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
  
      #ifdef _DEBUG
      fprintf(debug_f,"\nFinished boundary loop for bd %d.\n",b);
      fflush(debug_f);
      #endif

      //save old elements first
      oldnt[b] = nt[b];
      oldnq[b] = nq[b];
      
      gmark=ntnew;
  
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&ntnew,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (gmark > 0 && my_rank == 0)
      {
        fprintf(out_f,"\nBoundary %d, # of new triangles = %d",b+1,gmark);
        fflush(out_f);
      }
      if (ntnew > 0)
      {
        if (num_procs > 1)
        {
          tri_map[b] = (int*)realloc((void*)tri_map[b],(nt[b]+ntnew)*sizeof(int));
          tri_emap[b] = (int**)realloc((void*)tri_emap[b],(nt[b]+ntnew)*sizeof(int*));
        }
        for (t=nt[b]; t < nt[b]+ntnew && num_procs > 1; t++)
        {
          tri_emap[b][t] = (int*)malloc(2*sizeof(int));
          tri_map[b][t] = -1;
          for (i = 0; i < 2; i++)
            tri_emap[b][t][i] = -1;
        }
        t_n[b] = (int**)realloc((void*)t_n[b],(nt[b]+ntnew)*sizeof(int*));
        for (t=nt[b]; t < nt[b]+ntnew; t++)
        {
          t_n[b][t] = (int*)malloc(3*sizeof(int));
          for (i=0; i < 3; i++)
            t_n[b][t][i] = -1;
        }
      }
      gmark=nqnew;
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&nqnew,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (gmark > 0 && my_rank == 0)
      {
        fprintf(out_f,"\nBoundary %d, # of new quadrilaterals = %d",b+1,gmark);
        fflush(out_f);
      }
      if (nqnew > 0)
      {
        if (num_procs > 1)
        {
          quad_map[b] = (int*)realloc((void*)quad_map[b],(nq[b]+nqnew)*sizeof(int));
          quad_emap[b] = (int**)realloc((void*)quad_emap[b],(nq[b]+nqnew)*sizeof(int*));
        }
        for (q=nq[b]; q < nq[b]+nqnew && num_procs > 1; q++)
        {
          quad_emap[b][q] = (int*)malloc(2*sizeof(int));
          quad_map[b][q] = -1;
          for (i = 0; i < 2; i++)
            quad_emap[b][q][i] = -1;
        }
        q_n[b] = (int**)realloc((void*)q_n[b],(nq[b]+nqnew)*sizeof(int*));
        for (q=nq[b]; q < nq[b]+nqnew; q++)
        {
          q_n[b][q] = (int*)malloc(4*sizeof(int));
          for (i=0; i < 4; i++)
            q_n[b][q][i] = -1;
        }
      }
      
      ntnew = nqnew = 0;
      for (t=0; t < nt[b]; t++)
      {
        n0 = t_n[b][t][0];
        n1 = t_n[b][t][1];
        n2 = t_n[b][t][2];
        n3 = -1;
        ltet=lpyr=lpri=lhex=-1;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (ltet >= 0) retrieve_face_nodes(n0,n1,n2,n3,ltet,-1,-1,-1,tet_n[ltet],conn);
        if (lpyr >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
        if (lpri >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
        if (num_procs > 1)
          tester = optri[b][t];
        else
          tester = -1;
        refine_tri(conn,cdim,new_t,tri,new_q,quad,tester);
        for (c=0; c < new_t; c++)
        {
          if (c == 0)
          {
            for (j=0; j < 3; j++)
              t_n[b][t][j] = tri[c][j];
            //ownership/map set based on old tri
          } else
          {
            for (j=0; j < 3; j++)
              t_n[b][nt[b]+ntnew][j] = tri[c][j];
            if (num_procs > 1)
            {
              //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
              tri_emap[b][nt[b]+ntnew][0] = tri_emap[b][t][0];
              if (tri_emap[b][nt[b]+ntnew][0] == my_rank)
              {
                tri_emap[b][nt[b]+ntnew][1] = nt[b]+ntnew;
              }
              else
              {
                tri_emap[b][nt[b]+ntnew][1] = -1;
              }
            }
            ntnew++;
          }
        }
        if (new_t == 0) // mark for deletion later
          for (j=0; j < 3; j++)
            t_n[b][t][j] = -1;
        for (c=0; c < new_q; c++)
        {
          for (j=0; j < 4; j++)
          {
            q_n[b][nq[b]+nqnew][j] = quad[c][j];
          }
          if (num_procs > 1)
          {
            //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
            quad_emap[b][nq[b]+nqnew][0] = tri_emap[b][t][0];
            if (quad_emap[b][nq[b]+nqnew][0] == my_rank)
            {
              quad_emap[b][nq[b]+nqnew][1] = nq[b]+nqnew;
            }
            else
            {
              quad_emap[b][nq[b]+nqnew][1] = -1;
            }
          }
          nqnew++;
        }
      }
      #ifdef _DEBUG
      //fprintf(debug_f,"\nntnew = %d, nqnew = %d after TRI refinement.\n",ntnew,nqnew);
      //fflush(debug_f);
      #endif

      for (q=0; q < nq[b]; q++)
      {
        n0 = q_n[b][q][0];
        n1 = q_n[b][q][1];
        n2 = q_n[b][q][2];
        n3 = q_n[b][q][3];
        ltet=lpyr=lpri=lhex=-1;
        neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
        if (lpyr >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
        if (lpri >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
        if (lhex >= 0) retrieve_face_nodes(n0,n1,n2,n3,-1,-1,-1,lhex,hex_n[lhex],conn);
        if (num_procs > 1)
          tester = opquad[b][q];
        else
          tester = -1;
        refine_quad(conn,cdim,new_t,tri,new_q,quad,tester);
        
        for (c=0; c < new_t; c++)
        {
          for (j=0; j < 3; j++)
            t_n[b][nt[b]+ntnew][j] = tri[c][j];
          if (num_procs > 1)
          {
            //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
            tri_emap[b][nt[b]+ntnew][0] = quad_emap[b][q][0];
            if (tri_emap[b][nt[b]+ntnew][0] == my_rank)
            {
              tri_emap[b][nt[b]+ntnew][1] = nt[b]+ntnew;
            }
            else
            {
              tri_emap[b][nt[b]+ntnew][1] = -1;
            }
          }
          ntnew++;
        }
        for (c=0; c < new_q; c++)
        {
          if (c == 0)
          {
            for (j=0; j < 4; j++)
              q_n[b][q][j] = quad[c][j];
            //ownership pre-set
          } else
          {
            for (j=0; j < 4; j++)
            {
              q_n[b][nq[b]+nqnew][j] = quad[c][j];
            }
            if (num_procs > 1)
            {
              //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
              quad_emap[b][nq[b]+nqnew][0] = quad_emap[b][q][0];
              if (quad_emap[b][nq[b]+nqnew][0] == my_rank)
              {
                quad_emap[b][nq[b]+nqnew][1] = nq[b]+nqnew;
              }
              else
              {
                quad_emap[b][nq[b]+nqnew][1] = -1;
              }
            }
            nqnew++;
          }
        }
        if (new_q == 0) // mark for deletion later
          for (j=0; j < 4; j++)
            q_n[b][q][j] = -1;
      }

      #ifdef _DEBUG
      //fprintf(debug_f,"\nntnew = %d, nqnew = %d after QUAD refinement.\n",ntnew,nqnew);
      //fflush(debug_f);
      #endif

      //set numbers anew
      nt[b] += ntnew;
      nq[b] += nqnew;
      
      gmark=nt[b];
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&nt[b],&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        fprintf(out_f,"\nBefore compression: Boundary %d, # of triangles = %d",b+1,gmark);
        fflush(out_f);
      }
      gmark=nq[b];
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&nq[b],&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        fprintf(out_f,"\nBefore compression: Boundary %d, # of quadrilaterals = %d",b+1,gmark);
        fflush(out_f);
      }
      
      // compress lists to eliminate holes
      int *map;
      if (nt[b] > 0)
      {
        map = new int[nt[b]];
        for (j=t=0; t < nt[b]; t++)
          if (t_n[b][t][0] >= 0)
            map[t] = j++;
          else
            map[t] = -1;
        for (t=0; t < nt[b]; t++)
          if (map[t] >= 0 && map[t] != t)
          {
            if (num_procs > 1)
            {
              tri_map[b][map[t]] = tri_map[b][t];
              tri_emap[b][map[t]][0] = tri_emap[b][t][0];
              if (tri_emap[b][map[t]][0] == my_rank)
                tri_emap[b][map[t]][1] = map[t];
              else
                tri_emap[b][map[t]][1] = -1;
            }
            for (i=0; i < 3; i++)
              t_n[b][map[t]][i] = t_n[b][t][i];
          }
        if (j < nt[b])
        {
          for (t=j; t < nt[b]; t++)
          {
            free(t_n[b][t]);
            if (num_procs > 1)
              free(tri_emap[b][t]);
          }
          nt[b] = j;
          t_n[b] = (int**)realloc((void*)t_n[b],nt[b]*sizeof(int*));
          if (num_procs > 1)
          {
            tri_map[b] = (int*)realloc((void*)tri_map[b],nt[b]*sizeof(int));
            tri_emap[b] = (int**)realloc((void*)tri_emap[b],nt[b]*sizeof(int*));  
          }          
        }
        delete[] map;
      }
      if (nq[b] > 0)
      {
        map = new int[nq[b]];
        for (j=q=0; q < nq[b]; q++)
          if (q_n[b][q][0] >= 0)
            map[q] = j++;
          else
            map[q] = -1;
        for (q=0; q < nq[b]; q++)
          if (map[q] >= 0 && map[q] != q)
          {
            if (num_procs > 1)
            {
              quad_map[b][map[q]] = quad_map[b][q];
              quad_emap[b][map[q]][0] = quad_emap[b][q][0];
              if (quad_emap[b][map[q]][0] == my_rank)
                quad_emap[b][map[q]][1] = map[q];
              else
                quad_emap[b][map[q]][1] = -1;
            }
            for (i=0; i < 4; i++)
             q_n[b][map[q]][i] = q_n[b][q][i];
          }
        if (j < nq[b])
        {
          for (q=j; q < nq[b]; q++)
          {
            free(q_n[b][q]);
            if (num_procs > 1)
              free(quad_emap[b][q]);
          }
          nq[b] = j;
          q_n[b] = (int**)realloc((void*)q_n[b],nq[b]*sizeof(int*));
          if (num_procs > 1)
          {
            quad_map[b] = (int*)realloc((void*)quad_map[b],nq[b]*sizeof(int));
            quad_emap[b] = (int**)realloc((void*)quad_emap[b],nq[b]*sizeof(int*));
          }
        }
        delete[] map;
      }

      gmark=nt[b];
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&nt[b],&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        fprintf(out_f,"\nAfter compression: Boundary %d, # of triangles = %d",b+1,gmark);
        fflush(out_f);
      }
      gmark=nq[b];
      if (num_procs > 1)
      {
        gmark=0;
        MPI_Allreduce(&nq[b],&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      }
      if (my_rank == 0)
      {
        fprintf(out_f,"\nAfter compression: Boundary %d, # of quadrilaterals = %d",b+1,gmark);
        fflush(out_f);
      }
    }
  }

  // delete node-edge hash table
  for (n=0; n < nn; n++)
    delete nhash[n];
  delete nhash;

  // delete current hash tables
  for (n=0; n < nnhash; n++)
  {
    delete tethash[n];
    delete pyrhash[n];
    delete prihash[n];
    delete hexhash[n];
  }
  delete[] tethash;
  delete[] pyrhash;
  delete[] prihash;
  delete[] hexhash;

  // delete edges
  delete[] edge;
  // delete element-edge connectivities
  if (ntet > 0)
    delete[] tet_edges;
  if (npyr > 0)
    delete[] pyr_edges;
  if (npri > 0)
    delete[] pri_edges;
  if (nhex > 0)
    delete[] hex_edges;

  gmark=newpts;
  if (num_procs > 1)
  {
    gmark=0;
    MPI_Allreduce(&newpts,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (gmark == 0)
    return(1);

  if (my_rank == 0)
  {
    fprintf(out_f,"\nCreating %d new nodes.",gmark);
    fflush(out_f);
  }

  //
  // expand the node array, as well as pmap
  //
  node = (double**)realloc((void*)node,(nn+newpts)*sizeof(double*));
  if (num_procs > 1)
    pmap = (int**)realloc((void*)pmap,(nn+newpts)*sizeof(int*));
  for (n=nn; n < nn+newpts; n++)
  {
    node[n] = (double*)malloc(3*sizeof(double));
    if (num_procs > 1)
    {
      pmap[n] = (int*)malloc(3*sizeof(int));
      for (i = 0 ; i < 3; i++)
        pmap[n][i] = -1;
    }
  }
  //increment to new nn, and save old for later
  int oldnn = nn;
  nn += newpts;
  
  // create new mid-edge, mid-face and centroid points
  //while we're at it, number new nodes, will get proper p_i in do loop, will comm to get off proc global numbers later
  //need to be sure ownership is consitent among procs
  mark = 0;
  for (c=0; c < ntet; c++)
  {
    for (i=4; i < 10; i++)
    {
      if ((n=tet_n[c][i]) < 0) continue;
      mark++;
      int i0, i1;
      switch(i)
      {
        case 4: i0 = 0; i1 = 1; break;
        case 5: i0 = 1; i1 = 2; break;
        case 6: i0 = 2; i1 = 0; break;
        case 7: i0 = 0; i1 = 3; break;
        case 8: i0 = 1; i1 = 3; break;
        case 9: i0 = 2; i1 = 3; break;
      }
      for (j=0; j < 3; j++)
        node[n][j] = (node[tet_n[c][i0]][j] + node[tet_n[c][i1]][j])*0.5;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
  }
  gmark = mark;
  if (num_procs > 1)
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of tetrahedra mid-edge nodes created = %d",gmark);
    fflush(out_f);
  }
  mark=0;
  int face=0;
  int cent=0;
  for (c=0; c < npyr; c++)
  {
    for (i=5; i < 13; i++)
    {
      if ((n=pyr_n[c][i]) < 0) continue;
      mark++;
      int i0, i1;
      switch(i)
      {
        case 5: i0 = 0; i1 = 1; break;
        case 6: i0 = 1; i1 = 2; break;
        case 7: i0 = 2; i1 = 3; break;
        case 8: i0 = 3; i1 = 0; break;
        case 9: i0 = 0; i1 = 4; break;
        case 10: i0 = 1; i1 = 4; break;
        case 11: i0 = 2; i1 = 4; break;
        case 12: i0 = 3; i1 = 4; break;
      }
      for (j=0; j < 3; j++)
        node[n][j] = (node[pyr_n[c][i0]][j] + node[pyr_n[c][i1]][j])*0.5;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pyr_n[c][13]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[pyr_n[c][0]][j] + node[pyr_n[c][1]][j] +
                      node[pyr_n[c][2]][j] + node[pyr_n[c][3]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pyr_n[c][14]) >= 0)
    {
      cent++;
      for (j=0; j < 3; j++)
        node[n][j] = ((node[pyr_n[c][0]][j] + node[pyr_n[c][1]][j] +
                       node[pyr_n[c][2]][j] + node[pyr_n[c][3]][j])*0.25+node[pyr_n[c][4]][j])*0.5;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
  }
  gmark = mark;
  gface = face;
  gcent = cent;
  if (num_procs > 1)
  {
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&face,&gface,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&cent,&gcent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of pyramid mid-edge nodes created = %d",gmark);
    fprintf(out_f,"\nNumber of pyramid mid-face nodes created = %d",gface);
    fprintf(out_f,"\nNumber of pyramid centroid nodes created = %d",gcent);
    fflush(out_f);
  }
  mark=cent=face=0;
  for (c=0; c < npri; c++)
  {
    for (i=6; i < 15; i++)
    {
      if ((n=pri_n[c][i]) < 0) continue;
      mark++;
      int i0, i1;
      switch(i)
      {
        case 6: i0 = 0; i1 = 1; break;
        case 7: i0 = 1; i1 = 2; break;
        case 8: i0 = 2; i1 = 0; break;
        case 9: i0 = 0; i1 = 3; break;
        case 10: i0 = 1; i1 = 4; break;
        case 11: i0 = 2; i1 = 5; break;
        case 12: i0 = 3; i1 = 4; break;
        case 13: i0 = 4; i1 = 5; break;
        case 14: i0 = 5; i1 = 3; break;
      }
      for (j=0; j < 3; j++)
        node[n][j] = (node[pri_n[c][i0]][j] + node[pri_n[c][i1]][j])*0.5;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pri_n[c][15]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[pri_n[c][0]][j] + node[pri_n[c][1]][j] +
                      node[pri_n[c][3]][j] + node[pri_n[c][4]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pri_n[c][16]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[pri_n[c][1]][j] + node[pri_n[c][2]][j] +
                      node[pri_n[c][4]][j] + node[pri_n[c][5]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pri_n[c][17]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[pri_n[c][2]][j] + node[pri_n[c][0]][j] +
                      node[pri_n[c][5]][j] + node[pri_n[c][3]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=pri_n[c][18]) >= 0)
    {
      cent++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[pri_n[c][0]][j] + node[pri_n[c][1]][j] +
                      node[pri_n[c][2]][j] + node[pri_n[c][3]][j] +
                      node[pri_n[c][4]][j] + node[pri_n[c][5]][j])/6.0;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
  }
  
  gmark = mark;
  gface = face;
  gcent = cent;
  if (num_procs > 1)
  {
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&face,&gface,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&cent,&gcent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of prism mid-edge nodes created = %d",gmark);
    fprintf(out_f,"\nNumber of prism mid-face nodes created = %d",gface);
    fprintf(out_f,"\nNumber of prism centroid nodes created = %d",gcent);
    fflush(out_f);
  }
  
  mark=cent=face=0;
  for (c=0; c < nhex; c++)
  {
    for (i=8; i < 20; i++)
    {
      if ((n=hex_n[c][i]) < 0) continue;
      mark++;
      int i0, i1;
      switch(i)
      {
        case 8: i0 = 0; i1 = 1; break;
        case 9: i0 = 1; i1 = 2; break;
        case 10: i0 = 2; i1 = 3; break;
        case 11: i0 = 3; i1 = 0; break;
        case 12: i0 = 0; i1 = 4; break;
        case 13: i0 = 1; i1 = 5; break;
        case 14: i0 = 2; i1 = 6; break;
        case 15: i0 = 3; i1 = 7; break;
        case 16: i0 = 4; i1 = 5; break;
        case 17: i0 = 5; i1 = 6; break;
        case 18: i0 = 6; i1 = 7; break;
        case 19: i0 = 7; i1 = 4; break;
      }
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][i0]][j] + node[hex_n[c][i1]][j])*0.5;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][20]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][0]][j] + node[hex_n[c][1]][j] +
                      node[hex_n[c][3]][j] + node[hex_n[c][4]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][21]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][0]][j] + node[hex_n[c][1]][j] +
                      node[hex_n[c][4]][j] + node[hex_n[c][5]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][22]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][1]][j] + node[hex_n[c][2]][j] +
                      node[hex_n[c][5]][j] + node[hex_n[c][6]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][23]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][2]][j] + node[hex_n[c][3]][j] +
                      node[hex_n[c][6]][j] + node[hex_n[c][7]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][24]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][3]][j] + node[hex_n[c][0]][j] +
                      node[hex_n[c][7]][j] + node[hex_n[c][4]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
    if ((n=hex_n[c][25]) >= 0)
    {
      face++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][4]][j] + node[hex_n[c][5]][j] +
                      node[hex_n[c][6]][j] + node[hex_n[c][7]][j])*0.25;
      if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }

    if ((n=hex_n[c][26]) >= 0)
    {
      cent++;
      for (j=0; j < 3; j++)
        node[n][j] = (node[hex_n[c][0]][j] + node[hex_n[c][1]][j] +
                      node[hex_n[c][2]][j] + node[hex_n[c][3]][j] +
                      node[hex_n[c][4]][j] + node[hex_n[c][5]][j] +
                      node[hex_n[c][6]][j] + node[hex_n[c][7]][j])*0.125;
     if (num_procs > 1 && pmap[n][1] < 0)
      {
        //no need to set global now, only do this once, although n and my_rank should be the same each time
        pmap[n][1] = my_rank;
        pmap[n][2] = n; //if ends up owned, will be right, else will be overwritten by proper proc p_i 
      }
    }
  }
  gmark = mark;
  gface = face;
  gcent = cent;
  if (num_procs > 1)
  {
    MPI_Allreduce(&mark,&gmark,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&face,&gface,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&cent,&gcent,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNumber of hex mid-edge nodes created = %d",gmark);
    fprintf(out_f,"\nNumber of hex mid-face nodes created = %d",gface);
    fprintf(out_f,"\nNumber of hex centroid nodes created = %d",gcent);
    fflush(out_f);
  }

  if (num_procs > 1)
  {
    //in order to assure that we walk around an edge to get the right p_i pair, we use a do loop
    total_changes = 1;
    pass = 0;
    do
    {
      changes = 0;
      total_changes = 0;
      pass++;  
      // get the ownership for nodes on elements owned by other procs in do loop since walking around edge
      sendcnt = new int [num_procs];
      recvcnt = new int [num_procs];
      for (n = 0; n < num_procs; n++)
      {
        sendcnt[n] = 0;
        recvcnt[n] = 0;
      }
      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];

      //cycle thru elements	
      for (n = 0; n < ntet; n++)
      { 
        if (tet_emap[n][0] != my_rank)
        {
          sendcnt[tet_emap[n][0]]++; 
        }
      }
      for (n = 0; n < npyr; n++)
      { 
        if (pyr_emap[n][0] != my_rank)
        {
          sendcnt[pyr_emap[n][0]]++; 
        }
      }
      for (n = 0; n < npri; n++)
      { 
        if (pri_emap[n][0] != my_rank)
        {
          sendcnt[pri_emap[n][0]]++; 
        }
      }
      for (n = 0; n < nhex; n++)
      { 
        if (hex_emap[n][0] != my_rank)
        {
          sendcnt[hex_emap[n][0]]++; 
        }
      }
      
      #ifdef _DEBUG
      /*for (n = 0 ; n < num_procs; n++)
        fprintf(debug_f,"\n sendcnt[%d] = %d \n",n,sendcnt[n]);
      fflush(debug_f);

      fprintf(debug_f,"\nPrepared send counts.\n");
      fflush(debug_f);*/    
      #endif
		
      //now, send and recv count for all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
  
        MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      #ifdef _DEBUG
      /*for (n = 0 ; n < num_procs; n++)
        fprintf(debug_f,"\n recvcnt[%d] = %d \n",n,recvcnt[n]);
      fflush(debug_f);

      fprintf(debug_f,"\nSent counts.\n");
      fflush(debug_f);*/     
      #endif
	  
      delete [] srequest;
      delete [] rrequest;
      delete [] statuses;
  
      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];

      //allocate buffers to send/recv
      bdim2 = new int[num_procs];
      for (n = 0; n < num_procs; n++)
        bdim2[n] = MAX(sendcnt[n],recvcnt[n])*41*sizeof(int); //type, elem number, return number, up to 19 nodes p_i
      sbuff = new char*[num_procs];
      rbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
      {
        sbuff[n] = new char[bdim2[n]];
        rbuff[n] = new char[bdim2[n]];
      }

      //now, pack indices
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;

        //set position to 0
        sposition=0;

        //cycle thru elements	
        for (n = 0; n < ntet; n++)
        { 
          if (tet_emap[n][0] == p)
          {
            MPI_Pack(&(tt),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(tet_emap[n][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            //now, go through nodes and pack
            for (j = 4; j < 10; j++)
            {
              if ((k = tet_n[n][j]) >= 0)
              {
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            } 
          }
        }
        for (n = 0; n < npyr; n++)
        { 
          if (pyr_emap[n][0] == p)
          {
            MPI_Pack(&(pyt),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(pyr_emap[n][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            for (j = 5; j < 15; j++)
            {
              if ((k = pyr_n[n][j]) >= 0)
              {
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
        }
        for (n = 0; n < npri; n++)
        { 
          if (pri_emap[n][0] == p)
          {
            MPI_Pack(&(prt),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(pri_emap[n][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            for (j = 6; j < 19; j++)
            {
              if ((k = pri_n[n][j]) >= 0)
              {
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
        }
        for (n = 0; n < nhex; n++)
        { 
          if (hex_emap[n][0] == p)
          {
            MPI_Pack(&(ht),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(hex_emap[n][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            for (j = 8; j < 27; j++)
            {
              if ((k = hex_n[n][j]) >= 0)
              {
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
        }
      }

      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }

        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      #ifdef _DEBUG
      //fprintf(debug_f,"\nSent nodes.\n");
      //fflush(debug_f);     
      #endif
 
      delete [] srequest;
      delete [] rrequest;
      delete [] statuses;

      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];
  
      //now, delete sbuff and resize
      for (n = 0; n < num_procs; n++)
        if (sbuff[n] != 0) delete [] sbuff[n];
      delete [] sbuff;
  
      //resize
      bdim21 = new int[num_procs];
      for (n = 0; n < num_procs; n++)
        bdim21[n] = MAX(sendcnt[n],recvcnt[n])*40*sizeof(int); //type, element, p_i for up to 19 nodes
      sbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
        sbuff[n] = new char[bdim21[n]];
  
      temp_node = 0; //to hold type
      temp_node2 = 0; //element number
      temp_node3 = 0; //to hold other proc index
  
      //now, unpack local indices and repack nodes with other index
      //we are assuming that A) indices are in the same order on both procs and B) if no refinement on one, then not on other either
      for (p = 0; p < num_procs; p++)
      {
        sposition = rposition = 0; //reset pos

        if (recvcnt[p] == 0 || p == my_rank)
          continue; //do not want it to go into do while loop unnecessarily

        for (n = 0; n < recvcnt[p]; n++)
        {
          //unpack type initially
          MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

          if (temp_node == -10)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack elem to be returned
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
            //pack back up
            MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            //now, go through nodes and repack
            for (j = 4; j < 10; j++)
            {
              if ((k = tet_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                if (temp_node < pmap[k][1])
                {
                  pmap[k][1] = temp_node;
                  pmap[k][2] = temp_node3;
                  changes++;
                }
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
          else if (temp_node == -30)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack elem to be returned
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
            //pack back up
            MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            //now, go through nodes and repack
            for (j = 5; j < 15; j++)
            {
              if ((k = pyr_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                if (temp_node < pmap[k][1])
                {
                  pmap[k][1] = temp_node;
                  pmap[k][2] = temp_node3;
                  changes++;
                }
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
          else if (temp_node == -20)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack elem to be returned
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
            //pack back up
            MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            //now, go through nodes and repack
            for (j = 6; j < 19; j++)
            {
              if ((k = pri_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                if (temp_node < pmap[k][1])
                {
                  pmap[k][1] = temp_node;
                  pmap[k][2] = temp_node3;
                  changes++;
                }
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
          else if (temp_node == -40)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack elem to be returned
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
            //pack back up
            MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
            //now, go through nodes and repack
            for (j = 8; j < 27; j++)
            {
              if ((k = hex_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                if (temp_node < pmap[k][1])
                {
                  pmap[k][1] = temp_node;
                  pmap[k][2] = temp_node3;
                  changes++;
                }
                MPI_Pack(&(pmap[k][1]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pmap[k][2]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
            }
          }
          else
          {
            fprintf(stderr,"\n Proc: %d to %d: You sent an element that cannot be categorized! Exiting....\n",p,my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
        }
      }

      #ifdef _DEBUG
      //fprintf(debug_f,"\nPack unpack done.\n");
      //fflush(debug_f);     
      #endif

      //now, delete rbuff and bdim2, resize
      delete [] bdim2;
      for (n = 0; n < num_procs; n++)
        if (rbuff[n] != 0) delete [] rbuff[n];
      delete [] rbuff;
 
      //resize
      rbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
        rbuff[n] = new char[bdim21[n]];
  	
      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }

        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
 
      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);
    
      temp_node = 0; //to hold type
      temp_node2 = 0; //element number
      temp_node3 = 0; //to hold other proc proc
      temp_node6 = 0; //to hold other proc index

      //finally, unpack nodes in proper position
      for (p = 0; p < num_procs; p++)
      {
        rposition = 0; //reset pos

        if (sendcnt[p] == 0 || p == my_rank)
          continue; //do not want it to go into do while loop unnecessarily

        for (n = 0; n < sendcnt[p]; n++)
        {
          //unpack type initially
          MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

          if (temp_node == -10)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //now, go through nodes and check
            for (j = 4; j < 10; j++)
            {
              if ((k = tet_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node6,1,MPI_INT,MPI_COMM_WORLD);
                //now compare
                if (temp_node3 < pmap[k][1])
                {
                  pmap[k][1] = temp_node3;
                  pmap[k][2] = temp_node6;
                  changes++;
                }
              }
            }
          }
          else if (temp_node == -30)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //now, go through nodes and check
            for (j = 5; j < 15; j++)
            {
              if ((k = pyr_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node6,1,MPI_INT,MPI_COMM_WORLD);
                //now compare
                if (temp_node3 < pmap[k][1])
                {
                  pmap[k][1] = temp_node3;
                  pmap[k][2] = temp_node6;
                  changes++;
                }
              }
            }
          }
          else if (temp_node == -20)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //now, go through nodes and check
            for (j = 6; j < 19; j++)
            {
              if ((k = pri_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node6,1,MPI_INT,MPI_COMM_WORLD);
                //now compare
                if (temp_node3 < pmap[k][1])
                {
                  pmap[k][1] = temp_node3;
                  pmap[k][2] = temp_node6;
                  changes++;
                }
              }
            }
          }
          else if (temp_node == -40)
          {
            //unpack elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //now, go through nodes and check
            for (j = 8; j < 27; j++)
            {
              if ((k = hex_n[temp_node2][j]) >= 0)
              {
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node6,1,MPI_INT,MPI_COMM_WORLD);
                //now compare
                if (temp_node3 < pmap[k][1])
                {
                  pmap[k][1] = temp_node3;
                  pmap[k][2] = temp_node6;
                  changes++;
                }
              }
            }
          }
          else
          {
            fprintf(stderr,"\n Proc: %d to %d: You sent an element that cannot be categorized! Exiting....\n",p,my_rank);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }
        }
      }	

      //finally, free MPI mem
      delete[] sendcnt;
      delete[] recvcnt;
      delete[] srequest;
      delete[] rrequest;
      delete[] statuses;
      for (n = 0; n < num_procs; n++)
      {
        if (sbuff[n] != 0) delete [] sbuff[n];
        if (rbuff[n] != 0) delete [] rbuff[n];
      }
      delete[] sbuff;
      delete[] rbuff;
      delete[] bdim21;

      MPI_Allreduce(&changes,&total_changes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
      if (my_rank == 0)
      {
        fprintf(out_f,"\nNode numbering pass %d, changes = %d",pass, total_changes);
        fflush(out_f);
      }

    } while (total_changes > 0 && pass < 1000);
  }
   
  if (pass == 1000)
  {
    fprintf(stderr,"\nProc %d has made too many passes proc/local node numbering.  Exiting....\n",my_rank);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
   
  #ifdef _DEBUG
  fprintf(debug_f,"\nP_I indices synced.\n");
  fflush(debug_f);     
  #endif
  
  if (num_procs > 1)
  {
    // assign ownership of new nodes, based on element ownership, must do now to take advantage of original element positioning in larger conn
    //find highest global nn on current proc
    int localmax = -2; //set low for debug
    //need a check to see if no nodes added, and then cycle through list again
    for (n = 0; n < oldnn; n++)
    {
      localmax = MAX(localmax,pmap[n][0]);
    }
          
    #ifdef _DEBUG
    //fprintf(debug_f,"\n max on proc %d = %d \n",my_rank,localmax);
    //fflush(debug_f);     
    #endif
        
    globalmax = -3; //set low to assure actually working
    //now, comm to find globalmax on all procs
    MPI_Allreduce(&localmax,&globalmax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        
    #ifdef _DEBUG
    //fprintf(debug_f,"\n globalmax = %d \n",globalmax);
    //fflush(debug_f);     
    #endif
        
    int *npp = 0;
    npp = (int*)calloc(num_procs,sizeof(int));
        
    int nnowned = 0;
        
    //we cannot just send b, since some nodes were created that are not unique (on other procs too)
    for (n = oldnn; n < nn; n++)
      if (pmap[n][1] == my_rank)
        nnowned++;
        	
    #ifdef _DEBUG
    //fprintf(debug_f,"\n nnowned = %d \n",nnowned);
    //fflush(debug_f);     
    #endif
            
    //now, we do an all gather so procs can determine range for new node numbers
    MPI_Allgather(&nnowned,1,MPI_INT,npp,1,MPI_INT,MPI_COMM_WORLD);
        
    #ifdef _DEBUG
    /*for (n = 0; n < num_procs; n++)
    {
      fprintf(debug_f,"\nnpp[%d] = %d\n",n,npp[n]);
      fflush(debug_f);
    }*/   
    #endif
        
    //need to increment globalmax by one to begin at next node, which will ripple through all procs
    //this gives us true next node, or boundary for nodes if no new created
    globalmax++;
        
    //now, we can easily create global node numbers for owned nodes
    //get proc starting number
    for (n = 0; n < my_rank; n++)
      globalmax += npp[n];
          
    #ifdef _DEBUG
    //fprintf(debug_f,"\n starting globalmax = %d \n",globalmax);
    //fflush(debug_f);     
    #endif
        
    //will free npp after debug check
    //also, now set local node number for owned nodes      
    for (n = oldnn; n < nn; n++)
    {
      if (pmap[n][1] == my_rank)
      { 
        pmap[n][0] = globalmax;
        globalmax++; //will actually leave us with one greater than max node
      }
    } 
        
    MPI_Allreduce(&globalmax,&z,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        
    #ifdef _DEBUG
    //fprintf(debug_f,"\n final globalmax, plus 1 = %d \n",z);
    //fflush(debug_f);     
    #endif

    globalmax = z; //for debug
        
    if (npp > 0)    
      free(npp);
		
    //finally, get the global number for nodes owned by other procs
    //no need for do loop as edges already walked around in getting p_i pair
    sendcnt = new int [num_procs];
    recvcnt = new int [num_procs];
    for (n = 0; n < num_procs; n++)
    {
      sendcnt[n] = 0;
      recvcnt[n] = 0;
    }
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];

    //cycle thru nodes	
    for (n = oldnn; n < nn; n++)
    { 
      if (pmap[n][1] != my_rank)
      {
        sendcnt[pmap[n][1]]++; 
      }
    }
      
    #ifdef _DEBUG
    /*for (n = 0 ; n < num_procs; n++)
      fprintf(debug_f,"\n sendcnt[%d] = %d \n",n,sendcnt[n]);
    fflush(debug_f);

    fprintf(debug_f,"\nPrepared send counts.\n");
    fflush(debug_f);*/    
    #endif
		
    //now, send and recv count for all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;
  
      MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
  
      MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }

    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    #ifdef _DEBUG
    /*for (n = 0 ; n < num_procs; n++)
      fprintf(debug_f,"\n recvcnt[%d] = %d \n",n,recvcnt[n]);
    fflush(debug_f);

    fprintf(debug_f,"\nSent counts.\n");
    fflush(debug_f); */    
    #endif
	  
    delete [] srequest;
    delete [] rrequest;
    delete [] statuses;
	  
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];
		
    //allocate buffers to send/recv node indices
    bdim2 = new int[num_procs];
    for (n = 0; n < num_procs; n++)
      bdim2[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //indices
    sbuff = new char*[num_procs];
    rbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
    {
      sbuff[n] = new char[bdim2[n]];
      rbuff[n] = new char[bdim2[n]];
    }

    //now, pack node indices as per procs that will be sending them
    for (p = 0; p < num_procs; p++)
    {
      if (sendcnt[p] == 0 || p == my_rank)
        continue;

      //set position to 0
      sposition=0;

      //pack node index on this proc last, preceeded by element number and spot/type
      for (n = oldnn; n < nn; n++)
      {
        if (pmap[n][1] == p)
        {
          MPI_Pack(&(pmap[n][2]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
          MPI_Pack(&(n),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
        }
      }
    }

    //now, send and recv packets from all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (sendcnt[p] > 0 && p != my_rank) 
      {
        MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }

      if (recvcnt[p] > 0 && p != my_rank) 
      {
        MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }

    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    #ifdef _DEBUG
    //fprintf(debug_f,"\nSent nodes.\n");
    //fflush(debug_f);     
    #endif
  
    delete [] srequest;
    delete [] rrequest;
    delete [] statuses;

    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];
  
    //now, delete sbuff and resize
    for (n = 0; n < num_procs; n++)
      if (sbuff[n] != 0) delete [] sbuff[n];
    delete [] sbuff;
  
    //resize
    bdim21 = new int[num_procs];
    for (n = 0; n < num_procs; n++)
      bdim21[n] = MAX(sendcnt[n],recvcnt[n])*2*sizeof(int); //nodes and indices
    sbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
      sbuff[n] = new char[bdim21[n]];
  
    temp_node = 0; //to hold spot/type
    temp_node2 = 0; //element number
  
    //now, unpack local indices and repack nodes with other index
    for (p = 0; p < num_procs; p++)
    {
      sposition = rposition = 0; //reset pos

      if (recvcnt[p] == 0 || p == my_rank)
        continue; //do not want it to go into do while loop unnecessarily

      for (n = 0; n < recvcnt[p]; n++)
      {
        //unpack index initially
        MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
        //unpack index to return
        MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
        //repack that
        MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
        //pack global
        MPI_Pack(&(pmap[temp_node][0]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
      }
    }

    #ifdef _DEBUG
    //fprintf(debug_f,"\nPack unpack done.\n");
    //fflush(debug_f);     
    #endif

    //now, delete rbuff and bdim2, resize
    delete [] bdim2;
    for (n = 0; n < num_procs; n++)
      if (rbuff[n] != 0) delete [] rbuff[n];
    delete [] rbuff;
 
    //resize
    rbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
      rbuff[n] = new char[bdim21[n]];
  	
    //now, send and recv packets from all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (recvcnt[p] > 0 && p != my_rank) 
      {
        MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }

      if (sendcnt[p] > 0 && p != my_rank) 
      {
        MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }
 
    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    //finally, unpack nodes in proper position
    for (p = 0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;
        
      rposition = 0; //reset pos
      for (n = 0; n < sendcnt[p]; n++)
      {
        //unpack index and node
        MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(pmap[temp_node][0]),1,MPI_INT,MPI_COMM_WORLD);
      }
    }	
 
    #ifdef _DEBUG
    fprintf(debug_f,"\nGlobal indices synced.\n");
    fflush(debug_f);     
    #endif

    //finally, free MPI mem
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] srequest;
    delete[] rrequest;
    delete[] statuses;
    for (n = 0; n < num_procs; n++)
    {
      if (sbuff[n] != 0) delete [] sbuff[n];
      if (rbuff[n] != 0) delete [] rbuff[n];
    }
    delete[] sbuff;
    delete[] rbuff;
    delete[] bdim21;
    
    //debug
    for (n = oldnn; n < nn; n++)
    {
      if (pmap[n][0] < 0 || pmap[n][0] >= globalmax)
      {
        fprintf(stderr,"\nCATASTROPHIC ERROR, proc %d:  GN for %d (%d - %d - %d) = %d, which is out of bounds!\n",my_rank,n,pmap[n][0],pmap[n][1],pmap[n][2],pmap[n][0]);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pmap[n][1] < 0 || pmap[n][1] >= num_procs)
      {
        fprintf(stderr,"\nCATASTROPHIC ERROR, proc %d:  Proc number for %d (%d - %d - %d) = %d, which is out of bounds!\n",my_rank,n,pmap[n][0],pmap[n][1],pmap[n][2],pmap[n][1]);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      if (pmap[n][1] == my_rank && pmap[n][2] != n)
      {
        fprintf(stderr,"\nCATASTROPHIC ERROR, proc %d:  Proc local index for %d (%d - %d - %d) = %d, which is incorrect!\n",my_rank,n,pmap[n][0],pmap[n][1],pmap[n][2],pmap[n][2]);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }

  #ifdef _DEBUG
  fprintf(debug_f,"\nAbout to remap bd elem.\n");
  fflush(debug_f);     
  #endif

  //before creating new elements, we need to retool the maps for bd elem, since we were unable to use the space before pmap was created with new nodes
  //also, it should be noted that while we fill in emaps for owned elements, the second index is effectively useless after refinement has been exchanged
  if (num_procs > 1)
  {
    fgnn = 0; //reset
    //to be sure we don't get a mix up, we set maps to be determined by -"gnn" or triangle node, allowing for having an element on a proc with no nodes owned by it
    MPI_Allreduce(&nn,&fgnn,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
    //need to reset all elements!
    for (b = 0; b < nb; b++)
    {
      for (n = 0; n < nt[b]; n++)
      {
        tri_map[b][n] = -1;
        if (tri_emap[b][n][0] == my_rank)
        {
          tri_map[b][n] = -fgnn;
        }
        else
        {
          flag = 0;
          for (j = 0; j < 3 && !flag; j++)
          {
            if (tri_emap[b][n][0] == pmap[t_n[b][n][j]][1])
            {
              tri_map[b][n] = -(pmap[t_n[b][n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
              flag = 1;
            }    
          }
        }
      }
      for (n = 0; n < nq[b]; n++)
      {
        quad_map[b][n] = -1;
        if (quad_emap[b][n][0] == my_rank)
        {
          quad_map[b][n] = -fgnn;
        }
        else
        {
          flag = 0;
          for (j = 0; j < 4 && !flag; j++)
          {
            if (quad_emap[b][n][0] == pmap[q_n[b][n][j]][1])
            {
              quad_map[b][n] = -(pmap[q_n[b][n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
              flag = 1;
            }
          }
        }
      }
    }
  }

  #ifdef _DEBUG
  fprintf(debug_f,"\nDone remapping bd elem.\n");
  fflush(debug_f);     
  #endif

  gtet = ntet;
  gpyr = npyr;
  gpri = npri;
  ghex = nhex;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntet,&gtet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npyr,&gpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npri,&gpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nhex,&ghex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nCurrent number of tetrahedra = %d",gtet);
    fprintf(out_f,"\nCurrent number of pyramids   = %d",gpyr);
    fprintf(out_f,"\nCurrent number of prisms     = %d",gpri);
    fprintf(out_f,"\nCurrent number of hexahedra  = %d",ghex);
    fflush(out_f);
  }  

  // save old number of elements first for init to -1
  int oldtet = ntet;
  int oldpyr = npyr;
  int oldpri = npri;
  int oldhex = nhex;

  total_changes = 1;
  pass = 0;
  do
  {
    changes = 0;
    total_changes = 0;
    pass++;
    // count the number of new elements
    lnn   = nn;
    ltet = ntet;
    lpyr = npyr;
    lpri = npri;
    lhex = nhex;
    for (c=0; c < ntet; c++)
    {
      tester = optet[c];
      refine_tet(tet_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
      optet[c] = tester;
      ltet += MAX(0,new_tet-1);
      lpyr += new_pyr;
      lpri += new_pri;
      lhex += new_hex;
    }
    for (c=0; c < npyr; c++)
    {
      tester = oppyr[c];
      refine_pyr(pyr_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
      oppyr[c] = tester;
      ltet += new_tet;
      lpyr += MAX(0,new_pyr-1);
      lpri += new_pri;
      lhex += new_hex;
    }
    for (c=0; c < npri; c++)
    {
      tester = oppri[c];
      refine_pri(pri_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
      oppri[c] = tester;
      ltet += new_tet;
      lpyr += new_pyr;
      lpri += MAX(0,new_pri-1);
      lhex += new_hex;
    }
    for (c=0; c < nhex; c++)
    {
      tester = ophex[c];
      refine_hex(hex_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
      ophex[c] = tester;
      ltet += new_tet;
      lpyr += new_pyr;
      lpri += new_pri;
      lhex += MAX(0,new_hex-1);
    }
  
    //now, test markings
    //communicate markings to other procs using emaps, if in disagreement, mark changes++
    //if there is a disagreement, this will fix it (owning proc decides), rerun with proper mask, and then get elem count right
    if (num_procs > 1)
    {
      //first, set up mem 
      sendcnt = new int [num_procs];
      recvcnt = new int [num_procs];
      for (n = 0; n < num_procs; n++)
      {
        sendcnt[n] = 0;
        recvcnt[n] = 0;
      }
      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];
  
      //cycle thru elements of each type	
      for (n = 0; n < ntet; n++)
      {
        if (tet_emap[n][0] != my_rank)
        {
          sendcnt[tet_emap[n][0]]++; 
        }
      }	
      for (n = 0; n < npri; n++)
      {
        if (pri_emap[n][0] != my_rank)
        {
          sendcnt[pri_emap[n][0]]++; 
        }
      }
      for (n = 0; n < npyr; n++)
      {
        if (pyr_emap[n][0] != my_rank)
        {
          sendcnt[pyr_emap[n][0]]++; 
        }
      }
      for (n = 0; n < nhex; n++)
      {
        if (hex_emap[n][0] != my_rank)
        {
          sendcnt[hex_emap[n][0]]++; 
        }
      }
      
      #ifdef _DEBUG
      //for (n = 0 ; n < num_procs; n++)
        //fprintf(debug_f,"\n sendcnt[%d] = %d \n",n,sendcnt[n]);
      //fprintf(debug_f,"\nPrepared send counts.\n");
      //fflush(debug_f);     
      #endif

      //now, send and recv count for all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;
  
        MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
  
        MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      #ifdef _DEBUG
      //for (n = 0 ; n < num_procs; n++)
        //fprintf(debug_f,"\n recvcnt[%d] = %d \n",n,recvcnt[n]);
      //fprintf(debug_f,"\nReceived counts.\n");
      //fflush(debug_f);     
      #endif
  
      delete [] srequest;
      delete [] rrequest;
      delete [] statuses;
 
      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];

      //allocate buffers to send/recv indices
      bdim2 = new int[num_procs];
      for (n = 0; n < num_procs; n++)
        bdim2[n] = MAX(sendcnt[n],recvcnt[n])*(3*sizeof(int) + 1*sizeof(long int)); //indices
      sbuff = new char*[num_procs];
      rbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
      {
        sbuff[n] = new char[bdim2[n]];
        rbuff[n] = new char[bdim2[n]];
      }
 
      //set up types
      int tettype = -100;
      int pyrtype = -200;
      int pritype = -300;
      int hextype = -400;

      //now, pack node indices as per procs that will be sending them
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;

        //set position to 0
        sposition=0;

        //cycle thru elements of each type	
        for (c = 0; c < ntet; c++)
        {
          if (tet_emap[c][0] == p)
          {
            MPI_Pack(&(tettype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(tet_emap[c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(optet[c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
          }
        }	
        for (c = 0; c < npri; c++)
        {
          if (pri_emap[c][0] == p)
          {
            MPI_Pack(&(pritype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(pri_emap[c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(oppri[c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
          }
        }
        for (c = 0; c < npyr; c++)
        {
          if (pyr_emap[c][0] == p)
          {
            MPI_Pack(&(pyrtype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(pyr_emap[c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(oppyr[c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
          }
        }
        for (c = 0; c < nhex; c++)
        {
          if (hex_emap[c][0] == p)
          {
            MPI_Pack(&(hextype),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(hex_emap[c][1]),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD);
            MPI_Pack(&(ophex[c]),1,MPI_LONG,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
            MPI_Pack(&(c),1,MPI_INT,sbuff[p],bdim2[p],&sposition,MPI_COMM_WORLD); 
          }
        }
      }

      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }

        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }

      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      #ifdef _DEBUG
      //fprintf(debug_f,"\nSent markings.\n");
      //fflush(debug_f);     
      #endif
            
      delete [] srequest;
      delete [] rrequest;
      delete [] statuses;

      srequest = new MPI_Request[num_procs];
      rrequest = new MPI_Request[num_procs];
      statuses = new MPI_Status[num_procs];
 
      //now, delete sbuff and resize
      for (n = 0; n < num_procs; n++)
        if (sbuff[n] != 0) delete [] sbuff[n];
      delete [] sbuff;
 
      //resize
      bdim21 = new int[num_procs];
      for (n = 0; n < num_procs; n++)
        bdim21[n] = MAX(sendcnt[n],recvcnt[n])*(2*sizeof(int) + 1*sizeof(long int)); //marking and indices
      sbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
        sbuff[n] = new char[bdim21[n]];
 
      temp_node = temp_node2 = temp_node3 = 0;
      temp_node4 = temp_node5 = 0;
 
      //now, unpack type, index, mark, replace with type, op index, mark
      for (p = 0; p < num_procs; p++)
      {
        if (recvcnt[p] == 0 || p == my_rank)
          continue;

        sposition = rposition = 0; //reset pos
        for (n = 0; n < recvcnt[p]; n++)
        {
          //unpack type initially
          MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

          #ifdef _DEBUG
          //fprintf(debug_f,"\n type = %d \n",temp_node);
          //fflush(debug_f);     
          #endif

          //switch
          switch(temp_node)
          {
            case -100:
              //unpack index on proc
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n type = %d, elem = %d \n",temp_node, temp_node2);
              //fflush(debug_f);     
              #endif
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n marking = %d \n",temp_node4);
              //fflush(debug_f);     
              #endif
          
              //now, find owning proc marking
              //no change version
              if (optet[temp_node2] == temp_node4)
              {
                //repack type
                MPI_Pack(&(tettype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node2
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack back marking
                MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              else if (optet[temp_node2] != temp_node4) //change version
              {
                //repack type
                MPI_Pack(&(tettype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node3
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack new marking, as refinement of both sorts
                temp_node5 = (optet[temp_node2]) | (temp_node4);
                #ifdef _DEBUG
                fprintf(debug_f,"\nTET Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,optet[temp_node2],temp_node5);
                fflush(debug_f);
                #endif

                //set on this proc too
                optet[temp_node2] = temp_node5;
                //pack back
                MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //set up changes
                changes++;
              }
              break;
            case -200:
              //unpack index on proc
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n type = %d, elem = %d \n",temp_node, temp_node2);
              //fflush(debug_f);     
              #endif
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n marking = %d \n",temp_node4);
              //fflush(debug_f);     
              #endif
          
              //now, find owning proc marking
              //no change version
              if (oppyr[temp_node2] == temp_node4)
              {
                //repack type
                MPI_Pack(&(pyrtype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node2
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack back marking
                MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              else if (oppyr[temp_node2] != temp_node4) //change version
              {
                //repack type
                MPI_Pack(&(pyrtype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node3
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack new marking, as refinement of both sorts
                temp_node5 = (oppyr[temp_node2]) | (temp_node4);
            
                #ifdef _DEBUG
                fprintf(debug_f,"\nPYR Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,oppyr[temp_node2],temp_node5);
                fflush(debug_f);
                #endif

                //set on this proc too
                oppyr[temp_node2] = temp_node5;
                //pack back
                MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //set up changes
                changes++;
              }
              break;
            case -300:
              //unpack index on proc
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n type = %d, elem = %d \n",temp_node, temp_node2);
              //fflush(debug_f);     
              #endif
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n marking = %d \n",temp_node4);
              //fflush(debug_f);     
              #endif
          
              //now, find owning proc marking
              //no change version
              if (oppri[temp_node2] == temp_node4)
              {
                //repack type
                MPI_Pack(&(pritype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node2
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack back marking
                MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              else if (oppri[temp_node2] != temp_node4) //change version
              {
                //repack type
                MPI_Pack(&(pritype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node3
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack new marking, as refinement of both sorts
                temp_node5 = (oppri[temp_node2]) | (temp_node4);
                #ifdef _DEBUG
                fprintf(debug_f,"\nPRI Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,oppri[temp_node2],temp_node5);
                fflush(debug_f);
                #endif

                //set on this proc too
                oppri[temp_node2] = temp_node5;
                //pack back
                MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //set up changes
                changes++;
              }
              break;
            case -400:
              //unpack index on proc
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n type = %d, elem = %d \n",temp_node, temp_node2);
              //fflush(debug_f);     
              #endif
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node4,1,MPI_LONG,MPI_COMM_WORLD);
              #ifdef _DEBUG
              //fprintf(debug_f,"\n marking = %d \n",temp_node4);
              //fflush(debug_f);     
              #endif
          
              //now, find owning proc marking
              //no change version
              if (ophex[temp_node2] == temp_node4)
              {
                //repack type
                MPI_Pack(&(hextype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node2
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack back marking
                MPI_Pack(&(temp_node4),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              else if (ophex[temp_node2] != temp_node4) //change version
              {
                //repack type
                MPI_Pack(&(hextype),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //finally, unpack other proc index into unused temp_node3
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node3,1,MPI_INT,MPI_COMM_WORLD);
                //pack that immediately back for other proc
                MPI_Pack(&(temp_node3),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //pack new marking, as refinement of both sorts
                temp_node5 = (ophex[temp_node2]) | (temp_node4);
            
                #ifdef _DEBUG
                fprintf(debug_f,"\nHEX Refinement was %ld on proc and %ld off proc, and is now %ld for both.\n",temp_node4,ophex[temp_node2],temp_node5);
                fflush(debug_f);
                #endif
                //set on this proc too
                ophex[temp_node2] = temp_node5;
                //pack back
                MPI_Pack(&(temp_node5),1,MPI_LONG,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                //set up changes
                changes++;
              }
              break;
            default:
              fprintf(stderr,"\nSomehow you have more than four basic element types based on a bad index/misordering!  Exiting....\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);  
              break;
          }
        }
      }
          
      //fprintf(out_f,"\nFinished unpack/repack.\n");
      //fflush(out_f);
          
      MPI_Barrier(MPI_COMM_WORLD);
          
      //fprintf(out_f,"\nAbout to resize rbuff.\n");
      //fflush(out_f);

      //now, delete rbuff and bdim2, resize
      delete [] bdim2;
      for (n = 0; n < num_procs; n++)
        if (rbuff[n] != 0) delete [] rbuff[n];
      delete [] rbuff;

      //resize
      rbuff = new char*[num_procs];
      for (n = 0; n < num_procs; n++)
        rbuff[n] = new char[bdim21[n]];
          
      //fprintf(out_f,"\nAbout to send new numbers.\n");
      //fflush(out_f);

      //now, send and recv packets from all procs
      nreq_s = nreq_r = 0;
      for (p = 0; p < num_procs; p++)
      {
        if (recvcnt[p] > 0 && p != my_rank) 
        {
          MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }
  
        if (sendcnt[p] > 0 && p != my_rank) 
        {
          MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
  
      //now, wait for finish before unpacking
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //fprintf(out_f,"\nAbout to unpack four basic types.\n");
      //fflush(out_f);

      //finally, unpack elements and markings
      for (p = 0; p < num_procs; p++)
      {
        if (sendcnt[p] == 0 || p == my_rank)
          continue;

        rposition = 0; //reset pos
        for (n = 0; n < sendcnt[p]; n++)
        {
          //unpack type initially
          MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

          //fprintf(out_f,"\nunpacking type = %d\n",temp_node);
          //fflush(out_f);

          //switch
          switch(temp_node)
          {
            case -100:
              //unpack index use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack marking
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(optet[temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,optet[temp_node2]);
              //fflush(out_f);
              break;
            case -200:
              //unpack index use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack marking
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(oppyr[temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,oppyr[temp_node2]);
              //fflush(out_f);
              break;
            case -300:
              //unpack index use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack marking
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(oppri[temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,oppri[temp_node2]);
              //fflush(out_f);
              break;
            case -400:
              //unpack index use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack marking
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(ophex[temp_node2]),1,MPI_LONG,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, mark = %d\n",temp_node2,ophex[temp_node2]);
              //fflush(out_f);
              break;
            default:
              fprintf(stderr,"\nSomehow you have more than four basic element types based on a bad index/misordering!  Exiting....\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);  
              break;
          }
        }
      }

      //finally, free MPI mem
      delete[] sendcnt;
      delete[] recvcnt;
      delete[] srequest;
      delete[] rrequest;
      delete[] statuses;
      for (n = 0; n < num_procs; n++)
      {
        if (sbuff[n] != 0) delete [] sbuff[n];
        if (rbuff[n] != 0) delete [] rbuff[n];
      }
      delete[] sbuff;
      delete[] rbuff;
      delete[] bdim21;
    }
    
    //now, assuming that we have made all changes, it will go back thru, agree, and move on (if all metrics observed, but it will eventually stop)  
    MPI_Allreduce(&changes,&total_changes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (my_rank == 0)
    {
      fprintf(out_f,"\nPass %d, number of global changes = %d", pass, total_changes);
      fflush(out_f);
    }
  
  } while (total_changes > 0 && pass < 1000); //uses do loop since will affect element count
  
  if (pass == 1000)
  {
    fprintf(stderr,"\nProc %d has made too many passes at elem refinement codes.  Exiting....\n",my_rank);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }

  gtet = ltet;
  gpyr = lpyr;
  gpri = lpri;
  ghex = lhex;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ltet,&gtet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&lpyr,&gpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&lpri,&gpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&lhex,&ghex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nNew number of tetrahedra = %d",gtet);
    fprintf(out_f,"\nNew number of pyramids   = %d",gpyr);
    fprintf(out_f,"\nNew number of prisms     = %d",gpri);
    fprintf(out_f,"\nNew number of hexahedra  = %d",ghex);
    fflush(out_f);
  }

  // increase cell dimensions
  if (ltet > ntet)
  {
    tet_n = (int**)realloc((void*)tet_n,ltet*sizeof(int*));
    if (num_procs > 1)
    {
      tet_map = (int*)realloc((void*)tet_map,ltet*sizeof(int));
      for (c = ntet; c < ltet; c++)
        tet_map[c] = -1;
      tet_emap = (int**)realloc((void*)tet_emap,ltet*sizeof(int*));
    }
    for (c=ntet; c < ltet; c++)
    {
      tet_n[c] = (int*)malloc(10*sizeof(int));
      if (num_procs > 1)
      {
        tet_emap[c] = (int*)malloc(2*sizeof(int));
        for (i=0; i < 2; i++)
          tet_emap[c][i] = -1;
      }
      for (i=0; i < 10; i++)
        tet_n[c][i] = -1;
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nReallocated tets.\n");
  //fflush(debug_f);     
  #endif

  if (lpyr > npyr)
  {
    pyr_n = (int**)realloc((void*)pyr_n,lpyr*sizeof(int*));
    if (num_procs > 1)
    {
      pyr_map = (int*)realloc((void*)pyr_map,lpyr*sizeof(int));
      for (c = npyr; c < lpyr; c++)
        pyr_map[c] = -1;
      pyr_emap = (int**)realloc((void*)pyr_emap,lpyr*sizeof(int*));
    }
    for (c=npyr; c < lpyr; c++)
    {
      pyr_n[c] = (int*)malloc(15*sizeof(int));
      if (num_procs > 1)
      {
        pyr_emap[c] = (int*)malloc(2*sizeof(int));
        for (i=0; i < 2; i++)
          pyr_emap[c][i] = -1;
      }
      for (i=0; i < 15; i++)
        pyr_n[c][i] = -1;
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nReallocated pyrs.\n");
  //fflush(debug_f);     
  #endif

  if (lpri > npri)
  {
    pri_n = (int**)realloc((void*)pri_n,lpri*sizeof(int*));
    if (num_procs > 1)
    {
      pri_map = (int*)realloc((void*)pri_map,lpri*sizeof(int));
      for (c = npri; c < lpri; c++)
        pri_map[c] = -1;
      pri_emap = (int**)realloc((void*)pri_emap,lpri*sizeof(int*));
    }
    for (c=npri; c < lpri; c++)
    {
      pri_n[c] = (int*)malloc(19*sizeof(int));
      if (num_procs > 1)
      {
        pri_emap[c] = (int*)malloc(2*sizeof(int));
        for (i=0; i < 2; i++)
          pri_emap[c][i] = -1;
      }
      for (i=0; i < 19; i++)
        pri_n[c][i] = -1;
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nReallocated pris.\n");
  //fflush(debug_f);     
  #endif

  if (lhex > nhex)
  {
    hex_n = (int**)realloc((void*)hex_n,lhex*sizeof(int*));
    if (num_procs > 1)
    {
      hex_map = (int*)realloc((void*)hex_map,lhex*sizeof(int));
      for (c = nhex; c < lhex; c++)
        hex_map[c] = -1;
      hex_emap = (int**)realloc((void*)hex_emap,lhex*sizeof(int*));
    }
    for (c=nhex; c < lhex; c++)
    {
      hex_n[c] = (int*)malloc(27*sizeof(int));
      if (num_procs > 1)
      {
        hex_emap[c] = (int*)malloc(2*sizeof(int));
        for (i=0; i < 2; i++)
          hex_emap[c][i] = -1;
      }
      for (i=0; i < 27; i++)
        hex_n[c][i] = -1;
    }
  }
  
  #ifdef _DEBUG
  //fprintf(debug_f,"\nReallocated hexes.\n");
  //fflush(debug_f);     
  #endif

  //
  // create new elements

  // reset local counters
  ltet = ntet;
  lpyr = npyr;
  lpri = npri;
  lhex = nhex;

  for (c=0; c < ntet; c++)
  {
    if (num_procs > 1)
      tester = optet[c];
    else
      tester = -1;
    refine_tet(tet_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);

    if (new_tet > 0)
    {
      // overwrite the first one
      tet_n[c][0] = tet[0][0];
      tet_n[c][1] = tet[0][1];
      tet_n[c][2] = tet[0][2];
      tet_n[c][3] = tet[0][3];
      //ownership pre-set
      for (i=4; i < 10; i++)
        tet_n[c][i] = -1;
      // store the new ones
      for (i=1; i < new_tet; i++)
      {
        tet_n[ltet][0] = tet[i][0];
        tet_n[ltet][1] = tet[i][1];
        tet_n[ltet][2] = tet[i][2];
        tet_n[ltet][3] = tet[i][3];
        if (num_procs > 1)
        { 
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          tet_emap[ltet][0] = tet_emap[c][0];
          if (tet_emap[ltet][0] == my_rank)
          {
            tet_emap[ltet][1] = ltet;
            //tet_map[ltet] = -fgnn;
          }
          else
          {
            tet_emap[ltet][1] = -1;
            /*flag = 0;
            for (j = 0; j < 4 && !flag; j++)
              if (pmap[tet_n[ltet][j]][1] == tet_emap[ltet][0])
              {
                tet_map[ltet] = -(pmap[tet_n[ltet][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        ltet++;
      }
    } else
      for (i=0; i < 10; i++)
        tet_n[c][i] = -1; // mark for deletion later
    if (new_pyr > 0)
    {
      for (i=0; i < new_pyr; i++)
      {
        pyr_n[lpyr][0] = pyr[i][0];
        pyr_n[lpyr][1] = pyr[i][1];
        pyr_n[lpyr][2] = pyr[i][2];
        pyr_n[lpyr][3] = pyr[i][3];
        pyr_n[lpyr][4] = pyr[i][4];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pyr_emap[lpyr][0] = tet_emap[c][0];
          if (pyr_emap[lpyr][0] == my_rank)
          {
            pyr_emap[lpyr][1] = lpyr;
            //pyr_map[lpyr] = -fgnn;
          }
          else
          {
            pyr_emap[lpyr][1] = -1;
            /*flag = 0;
            for (j = 0; j < 5 && !flag; j++)
              if (pmap[pyr_n[lpyr][j]][1] == pyr_emap[lpyr][0])
              {
                pyr_map[lpyr] = -(pmap[pyr_n[lpyr][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpyr++;
      }
    }
    if (new_pri > 0)
    {
      for (i=0; i < new_pri; i++)
      {
        pri_n[lpri][0] = pri[i][0];
        pri_n[lpri][1] = pri[i][1];
        pri_n[lpri][2] = pri[i][2];
        pri_n[lpri][3] = pri[i][3];
        pri_n[lpri][4] = pri[i][4];
        pri_n[lpri][5] = pri[i][5];
        if (num_procs > 1)
        {
        //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pri_emap[lpri][0] = tet_emap[c][0];
          if (pri_emap[lpri][0] == my_rank)
          {
            pri_emap[lpri][1] = lpri;
            //pri_map[lpri] = -fgnn;
          }
          else
          {
            pri_emap[lpri][1] = -1;
            /*flag = 0;
            for (j = 0; j < 6 && !flag; j++)
              if (pmap[pri_n[lpri][j]][1] == pri_emap[lpri][0])
              {
                pri_map[lpri] = -(pmap[pri_n[lpri][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpri++;
      }
    }
    if (new_hex > 0)
    {
      for (i=0; i < new_hex; i++)
      {
        hex_n[lhex][0] = hex[i][0];
        hex_n[lhex][1] = hex[i][1];
        hex_n[lhex][2] = hex[i][2];
        hex_n[lhex][3] = hex[i][3];
        hex_n[lhex][4] = hex[i][4];
        hex_n[lhex][5] = hex[i][5];
        hex_n[lhex][6] = hex[i][6];
        hex_n[lhex][7] = hex[i][7];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          hex_emap[lhex][0] = tet_emap[c][0];
          if (hex_emap[lhex][0] == my_rank)
          {
            hex_emap[lhex][1] = lhex;
            //hex_map[lhex] = -fgnn;
          }
          else
          {
            hex_emap[lhex][1] = -1;
            /*flag = 0;
            for (j = 0; j < 8 && !flag; j++)
              if (pmap[hex_n[lhex][j]][1] == hex_emap[lhex][0])
              {
                hex_map[lhex] = -(pmap[hex_n[lhex][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lhex++;
      }
    }
  }

  for (c=0; c < npyr; c++)
  {
    if (num_procs > 1)
      tester = oppyr[c];
    else
      tester = -1;
    refine_pyr(pyr_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
    if (new_pyr > 0)
    {
      pyr_n[c][0] = pyr[0][0];
      pyr_n[c][1] = pyr[0][1];
      pyr_n[c][2] = pyr[0][2];
      pyr_n[c][3] = pyr[0][3];
      pyr_n[c][4] = pyr[0][4];
      //ownership pre-set
      for (i=5; i < 15; i++)
        pyr_n[c][i] = -1;
      for (i=1; i < new_pyr; i++)
      {
        pyr_n[lpyr][0] = pyr[i][0];
        pyr_n[lpyr][1] = pyr[i][1];
        pyr_n[lpyr][2] = pyr[i][2];
        pyr_n[lpyr][3] = pyr[i][3];
        pyr_n[lpyr][4] = pyr[i][4];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pyr_emap[lpyr][0] = pyr_emap[c][0];
          if (pyr_emap[lpyr][0] == my_rank)
          {
            pyr_emap[lpyr][1] = lpyr;
            //pyr_map[lpyr] = -fgnn;
          }
          else
          {
            pyr_emap[lpyr][1] = -1;
            /*flag = 0;
            for (j = 0; j < 5 && !flag; j++)
              if (pmap[pyr_n[lpyr][j]][1] == pyr_emap[lpyr][0])
              {
                pyr_map[lpyr] = -(pmap[pyr_n[lpyr][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpyr++;
      }
    } else
      for (i=0; i < 15; i++)
        pyr_n[c][i] = -1; // mark for deletion later
    if (new_tet > 0)
    {
      for (i=0; i < new_tet; i++)
      {
        tet_n[ltet][0] = tet[i][0];
        tet_n[ltet][1] = tet[i][1];
        tet_n[ltet][2] = tet[i][2];
        tet_n[ltet][3] = tet[i][3];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          tet_emap[ltet][0] = pyr_emap[c][0];
          if (tet_emap[ltet][0] == my_rank)
          {
            tet_emap[ltet][1] = ltet;
            //tet_map[ltet] = -fgnn;
          }
          else
          {
            tet_emap[ltet][1] = -1;
            /*flag = 0;
            for (j = 0; j < 4 && !flag; j++)
              if (pmap[tet_n[ltet][j]][1] == tet_emap[ltet][0])
              {
                tet_map[ltet] = -(pmap[tet_n[ltet][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        ltet++;
      }
    }
    if (new_pri > 0)
    {
      for (i=0; i < new_pri; i++)
      {
        pri_n[lpri][0] = pri[i][0];
        pri_n[lpri][1] = pri[i][1];
        pri_n[lpri][2] = pri[i][2];
        pri_n[lpri][3] = pri[i][3];
        pri_n[lpri][4] = pri[i][4];
        pri_n[lpri][5] = pri[i][5];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pri_emap[lpri][0] = pyr_emap[c][0];
          if (pri_emap[lpri][0] == my_rank)
          {
            pri_emap[lpri][1] = lpri;
            //pri_map[lpri] = -fgnn;
          }
          else
          {
            pri_emap[lpri][1] = -1;
            /*flag = 0;
            for (j = 0; j < 6 && !flag; j++)
              if (pmap[pri_n[lpri][j]][1] == pri_emap[lpri][0])
              {
                pri_map[lpri] = -(pmap[pri_n[lpri][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpri++;
      }
    }
    if (new_hex > 0)
    {
      for (i=0; i < new_hex; i++)
      {
        hex_n[lhex][0] = hex[i][0];
        hex_n[lhex][1] = hex[i][1];
        hex_n[lhex][2] = hex[i][2];
        hex_n[lhex][3] = hex[i][3];
        hex_n[lhex][4] = hex[i][4];
        hex_n[lhex][5] = hex[i][5];
        hex_n[lhex][6] = hex[i][6];
        hex_n[lhex][7] = hex[i][7];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          hex_emap[lhex][0] = pyr_emap[c][0];
          if (hex_emap[lhex][0] == my_rank)
          {
            hex_emap[lhex][1] = lhex;
            //hex_map[lhex] = -fgnn;
          }
          else
          {
            hex_emap[lhex][1] = -1;
            /*flag = 0;
            for (j = 0; j < 8 && !flag; j++)
              if (pmap[hex_n[lhex][j]][1] == hex_emap[lhex][0])
              {
                hex_map[lhex] = -(pmap[hex_n[lhex][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lhex++;
      }
    }
  }
  for (c=0; c < npri; c++)
  {
    if (num_procs > 1)
      tester = oppri[c];
    else
      tester = -1;
    refine_pri(pri_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
    if (new_pri > 0)
    {
      pri_n[c][0] = pri[0][0];
      pri_n[c][1] = pri[0][1];
      pri_n[c][2] = pri[0][2];
      pri_n[c][3] = pri[0][3];
      pri_n[c][4] = pri[0][4];
      pri_n[c][5] = pri[0][5];
      //ownership pre-set
      for (i=6; i < 19; i++)
        pri_n[c][i] = -1;
      for (i=1; i < new_pri; i++)
      {
        pri_n[lpri][0] = pri[i][0];
        pri_n[lpri][1] = pri[i][1];
        pri_n[lpri][2] = pri[i][2];
        pri_n[lpri][3] = pri[i][3];
        pri_n[lpri][4] = pri[i][4];
        pri_n[lpri][5] = pri[i][5];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pri_emap[lpri][0] = pri_emap[c][0];
          if (pri_emap[lpri][0] == my_rank)
          {
            pri_emap[lpri][1] = lpri;
            //pri_map[lpri] = -fgnn;
          }
          else
          {
            pri_emap[lpri][1] = -1;
            /*flag = 0;
            for (j = 0; j < 6 && !flag; j++)
              if (pmap[pri_n[lpri][j]][1] == pri_emap[lpri][0])
              {
                pri_map[lpri] = -(pmap[pri_n[lpri][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpri++;
      }
    } else
      for (i=0; i < 19; i++)
        pri_n[c][i] = -1; // mark for deletion later
    if (new_tet > 0)
    {
      for (i=0; i < new_tet; i++)
      {
        tet_n[ltet][0] = tet[i][0];
        tet_n[ltet][1] = tet[i][1];
        tet_n[ltet][2] = tet[i][2];
        tet_n[ltet][3] = tet[i][3];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          tet_emap[ltet][0] = pri_emap[c][0];
          if (tet_emap[ltet][0] == my_rank)
          {
            tet_emap[ltet][1] = ltet;
            //tet_map[ltet] = -fgnn;
          }
          else
          {
            tet_emap[ltet][1] = -1;
            /*flag = 0;
            for (j = 0; j < 4 && !flag; j++)
              if (pmap[tet_n[ltet][j]][1] == tet_emap[ltet][0])
              {
                tet_map[ltet] = -(pmap[tet_n[ltet][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        ltet++;
      }
    }
    if (new_pyr > 0)
    {
      for (i=0; i < new_pyr; i++)
      {
        pyr_n[lpyr][0] = pyr[i][0];
        pyr_n[lpyr][1] = pyr[i][1];
        pyr_n[lpyr][2] = pyr[i][2];
        pyr_n[lpyr][3] = pyr[i][3];
        pyr_n[lpyr][4] = pyr[i][4];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pyr_emap[lpyr][0] = pri_emap[c][0];
          if (pyr_emap[lpyr][0] == my_rank)
          {
            pyr_emap[lpyr][1] = lpyr;
            //pyr_map[lpyr] = -fgnn;
          }
          else
          {
            pyr_emap[lpyr][1] = -1;
            /*flag = 0;
            for (j = 0; j < 5 && !flag; j++)
              if (pmap[pyr_n[lpyr][j]][1] == pyr_emap[lpyr][0])
              {
                pyr_map[lpyr] = -(pmap[pyr_n[lpyr][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpyr++;
      }
    }
    if (new_hex > 0)
    {
      for (i=0; i < new_hex; i++)
      {
        hex_n[lhex][0] = hex[i][0];
        hex_n[lhex][1] = hex[i][1];
        hex_n[lhex][2] = hex[i][2];
        hex_n[lhex][3] = hex[i][3];
        hex_n[lhex][4] = hex[i][4];
        hex_n[lhex][5] = hex[i][5];
        hex_n[lhex][6] = hex[i][6];
        hex_n[lhex][7] = hex[i][7];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          hex_emap[lhex][0] = pri_emap[c][0];
          if (hex_emap[lhex][0] == my_rank)
          {
            hex_emap[lhex][1] = lhex;
            //hex_map[lhex] = -fgnn;
          }
          else
          {
            hex_emap[lhex][1] = -1;
            /*flag = 0;
            for (j = 0; j < 8 && !flag; j++)
              if (pmap[hex_n[lhex][j]][1] == hex_emap[lhex][0])
              {
                hex_map[lhex] = -(pmap[hex_n[lhex][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lhex++;
      }
    }
  }
  for (c=0; c < nhex; c++)
  {
    if (num_procs > 1)
      tester = ophex[c];
    else
      tester = -1;
    refine_hex(hex_n[c],cdim,new_tet,tet,new_pyr,pyr,new_pri,pri,new_hex,hex,tester);
    if (new_hex > 0)
    {
      hex_n[c][0] = hex[0][0];
      hex_n[c][1] = hex[0][1];
      hex_n[c][2] = hex[0][2];
      hex_n[c][3] = hex[0][3];
      hex_n[c][4] = hex[0][4];
      hex_n[c][5] = hex[0][5];
      hex_n[c][6] = hex[0][6];
      hex_n[c][7] = hex[0][7];
      //ownership pre-set
      for (i=8; i < 27; i++)
        hex_n[c][i] = -1;
      for (i=1; i < new_hex; i++)
      {
        hex_n[lhex][0] = hex[i][0];
        hex_n[lhex][1] = hex[i][1];
        hex_n[lhex][2] = hex[i][2];
        hex_n[lhex][3] = hex[i][3];
        hex_n[lhex][4] = hex[i][4];
        hex_n[lhex][5] = hex[i][5];
        hex_n[lhex][6] = hex[i][6];
        hex_n[lhex][7] = hex[i][7];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          hex_emap[lhex][0] = hex_emap[c][0];
          if (hex_emap[lhex][0] == my_rank)
          {
            hex_emap[lhex][1] = lhex;
            //hex_map[lhex] = -fgnn;
          }
          else
          {
            hex_emap[lhex][1] = -1;
            /*flag = 0;
            for (j = 0; j < 8 && !flag; j++)
              if (pmap[hex_n[lhex][j]][1] == hex_emap[lhex][0])
              {
                hex_map[lhex] = -(pmap[hex_n[lhex][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lhex++;
      }
    } else
      for (i=0; i < 27; i++)
        hex_n[c][i] = -1; // mark for deletion later
    if (new_tet > 0)
    {
      for (i=0; i < new_tet; i++)
      {
        tet_n[ltet][0] = tet[i][0];
        tet_n[ltet][1] = tet[i][1];
        tet_n[ltet][2] = tet[i][2];
        tet_n[ltet][3] = tet[i][3];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          tet_emap[ltet][0] = hex_emap[c][0];
          if (tet_emap[ltet][0] == my_rank)
          {
            tet_emap[ltet][1] = ltet;
            //tet_map[ltet] = -fgnn;
          }
          else
          {
            tet_emap[ltet][1] = -1;
            /*flag = 0;
            for (j = 0; j < 4 && !flag; j++)
              if (pmap[tet_n[ltet][j]][1] == tet_emap[ltet][0])
              {
                tet_map[ltet] = -(pmap[tet_n[ltet][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        ltet++;
      }
    }
    if (new_pyr > 0)
    {
      for (i=0; i < new_pyr; i++)
      {
        pyr_n[lpyr][0] = pyr[i][0];
        pyr_n[lpyr][1] = pyr[i][1];
        pyr_n[lpyr][2] = pyr[i][2];
        pyr_n[lpyr][3] = pyr[i][3];
        pyr_n[lpyr][4] = pyr[i][4];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pyr_emap[lpyr][0] = hex_emap[c][0];
          if (pyr_emap[lpyr][0] == my_rank)
          {
            pyr_emap[lpyr][1] = lpyr;
            //pyr_map[lpyr] = -fgnn;
          }
          else
          {
            pyr_emap[lpyr][1] = -1;
            /*flag = 0;
            for (j = 0; j < 5 && !flag; j++)
              if (pmap[pyr_n[lpyr][j]][1] == pyr_emap[lpyr][0])
              {
                pyr_map[lpyr] = -(pmap[pyr_n[lpyr][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpyr++;
      }
    }
    if (new_pri > 0)
    {
      for (i=0; i < new_pri; i++)
      {
        pri_n[lpri][0] = pri[i][0];
        pri_n[lpri][1] = pri[i][1];
        pri_n[lpri][2] = pri[i][2];
        pri_n[lpri][3] = pri[i][3];
        pri_n[lpri][4] = pri[i][4];
        pri_n[lpri][5] = pri[i][5];
        if (num_procs > 1)
        {
          //set ownership (for bd elem, will have to set map to p-i of node to avoid global search after pmap is set)
          pri_emap[lpri][0] = hex_emap[c][0];
          if (pri_emap[lpri][0] == my_rank)
          {
            pri_emap[lpri][1] = lpri;
            //pri_map[lpri] = -fgnn;
          }
          else
          {
            pri_emap[lpri][1] = -1;
            /*flag = 0;
            for (j = 0; j < 6 && !flag; j++)
              if (pmap[pri_n[lpri][j]][1] == pri_emap[lpri][0])
              {
                pri_map[lpri] = -(pmap[pri_n[lpri][j]][2] + 2);
                flag = 1; 
              }*/
          }
        }
        lpri++;
      }
    }
  }
  
  #ifdef _DEBUG
  fprintf(debug_f,"\nCreated new elements.\n");
  fflush(debug_f);     
  #endif

  MPI_Barrier(MPI_COMM_WORLD);     

  //set new number of elements
  ntet = ltet;
  npyr = lpyr;
  npri = lpri;
  nhex = lhex;

  // compress the element lists (eliminate holes)

  gtet = ntet;
  gpyr = npyr;
  gpri = npri;
  ghex = nhex;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntet,&gtet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npyr,&gpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npri,&gpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nhex,&ghex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nBefore compressions, number of tetrahedra = %d",gtet);
    fprintf(out_f,"\nBefore compressions, number of pyramids   = %d",gpyr);
    fprintf(out_f,"\nBefore compressions, number of prisms     = %d",gpri);
    fprintf(out_f,"\nBefore compressions, number of hexahedra  = %d",ghex);
    fflush(out_f);
  }

  for (c=0; c < ntet; c++)
  {
    if (tet_n[c][0] < 0)
    {
      while (tet_n[ntet-1][0] < 0 && ntet > c)
        ntet--;
      if (ntet > c)
      {
        if (num_procs > 1)
          tet_map[c] = tet_map[ntet-1];
        for (i=0; i < 4; i++)
          tet_n[c][i] = tet_n[ntet-1][i];
        if (num_procs > 1)
        {
          tet_emap[c][0] = tet_emap[ntet-1][0];
          if (tet_emap[c][0] == my_rank)
            tet_emap[c][1] = c;
          else
            tet_emap[c][1] = -1;
        }
        ntet--;
      }
    }
  }

  for (c=0; c < npyr; c++)
  {
    if (pyr_n[c][0] < 0)
    {
      while (pyr_n[npyr-1][0] < 0 && npyr > c)
        npyr--;
      if (npyr > c)
      {
        if (num_procs > 1)
          pyr_map[c] = pyr_map[npyr-1];
        for (i=0; i < 5; i++)
          pyr_n[c][i] = pyr_n[npyr-1][i];
        if (num_procs > 1)
        {
          pyr_emap[c][0] = pyr_emap[npyr-1][0];
          if (pyr_emap[c][0] == my_rank)
            pyr_emap[c][1] = c;
          else
            pyr_emap[c][1] = -1;
        }
        npyr--;
      }
    }
  }

  for (c=0; c < npri; c++)
  {
    if (pri_n[c][0] < 0)
    {
      while (pri_n[npri-1][0] < 0 && npri > c)
        npri--;
      if (npri > c)
      {
        if (num_procs > 1)
          pri_map[c] = pri_map[npri-1];
        for (i=0; i < 6; i++)
          pri_n[c][i] = pri_n[npri-1][i];
        if (num_procs > 1)
        {
          pri_emap[c][0] = pri_emap[npri-1][0];
          if (pri_emap[c][0] == my_rank)
            pri_emap[c][1] = c;
          else
            pri_emap[c][1] = -1;
        }
        npri--;
      }
    }
  }

  for (c=0; c < nhex; c++)
  {
    if (hex_n[c][0] < 0)
    {
      while (hex_n[nhex-1][0] < 0 && nhex > c)
        nhex--;
      if (nhex > c)
      {
        if (num_procs > 1)
          hex_map[c] = hex_map[nhex-1];
        for (i=0; i < 8; i++)
          hex_n[c][i] = hex_n[nhex-1][i];
        if (num_procs > 1)
        {
          hex_emap[c][0] = hex_emap[nhex-1][0];
          if (hex_emap[c][0] == my_rank)
            hex_emap[c][1] = c;
          else
            hex_emap[c][1] = -1;
        }
        nhex--;
      }
    }
  }

  //Resizing to actual count and reduce to linear elements
  if (ltet > ntet)
  {
    for (c=ntet; c < ltet; c++)
    {
      free(tet_n[c]);
      if (num_procs > 1)
        free(tet_emap[c]);
    }
    tet_n = (int**)realloc((void*)tet_n,ntet*sizeof(int*));
    if (num_procs > 1)
    {
      tet_emap = (int**)realloc((void*)tet_emap,ntet*sizeof(int*));
      tet_map = (int*)realloc((void*)tet_map,ntet*sizeof(int));
    }
    for (c=0; c < ntet; c++)
    {
      tet_n[c] = (int*)realloc((void*)tet_n[c],4*sizeof(int));
      if (num_procs > 1)
        tet_emap[c] = (int*)realloc((void*)tet_emap[c],2*sizeof(int));
    }
  }
  if (lpyr > npyr)
  {
    for (c=npyr; c < lpyr; c++)
    {
      free(pyr_n[c]);
      if (num_procs > 1)
        free(pyr_emap[c]);
    }
    pyr_n = (int**)realloc((void*)pyr_n,npyr*sizeof(int*));
    if (num_procs > 1)
    {
      pyr_emap = (int**)realloc((void*)pyr_emap,npyr*sizeof(int*));
      pyr_map = (int*)realloc((void*)pyr_map,npyr*sizeof(int));
    }
    for (c=0; c < npyr; c++)
    {
      pyr_n[c] = (int*)realloc((void*)pyr_n[c],5*sizeof(int));
      if (num_procs > 1)
        pyr_emap[c] = (int*)realloc((void*)pyr_emap[c],2*sizeof(int));
    }
  }
  if (lpri > npri)
  {
    for (c=npri; c < lpri; c++)
    {
      free(pri_n[c]);
      if (num_procs > 1)
        free(pri_emap[c]);
    }
    if (num_procs > 1)
    {
      pri_emap = (int**)realloc((void*)pri_emap,npri*sizeof(int*));
      pri_map = (int*)realloc((void*)pri_map,npri*sizeof(int));
    }
    pri_n = (int**)realloc((void*)pri_n,npri*sizeof(int*));
    for (c=0; c < npri; c++)
    {
      pri_n[c] = (int*)realloc((void*)pri_n[c],6*sizeof(int));
      if (num_procs > 1)
        pri_emap[c] = (int*)realloc((void*)pri_emap[c],2*sizeof(int));
    }
  }
  if (lhex > nhex)
  {
    for (c=nhex; c < lhex; c++)
    {
      free(hex_n[c]);
      if (num_procs > 1)
        free(hex_emap[c]);
    }
    if (num_procs > 1)
    {
      hex_emap = (int**)realloc((void*)hex_emap,nhex*sizeof(int*));
      hex_map = (int*)realloc((void*)hex_map,nhex*sizeof(int));
    }
    hex_n = (int**)realloc((void*)hex_n,nhex*sizeof(int*));
    for (c=0; c < nhex; c++)
    {
      hex_n[c] = (int*)realloc((void*)hex_n[c],8*sizeof(int));
      if (num_procs > 1)
        hex_emap[c] = (int*)realloc((void*)hex_emap[c],2*sizeof(int));
    }
  }

  gtet = ntet;
  gpyr = npyr;
  gpri = npri;
  ghex = nhex;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntet,&gtet,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npyr,&gpyr,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&npri,&gpri,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nhex,&ghex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
  if (my_rank == 0)
  {
    fprintf(out_f,"\nAfter compressions, number of tetrahedra = %d",gtet);
    fprintf(out_f,"\nAfter compressions, number of pyramids   = %d",gpyr);
    fprintf(out_f,"\nAfter compressions, number of prisms     = %d",gpri);
    fprintf(out_f,"\nAfter compressions, number of hexahedra  = %d",ghex);
    fflush(out_f);
  }
  
  //delete memory created for opelem
  if (optet > 0)
    free(optet);
  if (ophex > 0)
    free(ophex);
  if (oppri > 0)
    free(oppri);
  if (oppyr > 0)
    free(oppyr);
  for (n = 0; n < nb; n++)
    if (optri[n] > 0)
      free(optri[n]);
  if (optri > 0)
    free(optri);
  for (n = 0; n < nb; n++)
    if (opquad[n] > 0)
      free(opquad[n]);
  if (opquad > 0)
    free(opquad);

  if (num_procs > 1)
  {
    //now, go through elements and reset maps
    for (n = 0; n < ntet; n++)
    {
      tet_map[n] = -1; //reset in case off proc tet
      if (tet_emap[n][0] == my_rank)
      {
        tet_map[n] = -fgnn;
      }
      else
      {
        flag = 0;
        for (j = 0; j < 4 && !flag; j++)
        {
          if (tet_emap[n][0] == pmap[tet_n[n][j]][1])
          {
            tet_map[n] = -(pmap[tet_n[n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
            flag = 1;
          }    
        }
      }
    }
    
    for (n = 0; n < npyr; n++)
    {
      pyr_map[n] = -1; //reset in case off proc tet
      if (pyr_emap[n][0] == my_rank)
      {
        pyr_map[n] = -fgnn;
      }
      else
      {
        flag = 0;
        for (j = 0; j < 5 && !flag; j++)
        {
          if (pyr_emap[n][0] == pmap[pyr_n[n][j]][1])
          {
            pyr_map[n] = -(pmap[pyr_n[n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
            flag = 1;
          }    
        }
      }
    }
      
    for (n = 0; n < npri; n++)
    {
      pri_map[n] = -1; //reset in case off proc tet
      if (pri_emap[n][0] == my_rank)
      {
        pri_map[n] = -fgnn;
      }
      else
      {
        flag = 0;
        for (j = 0; j < 6 && !flag; j++)
        {
          if (pri_emap[n][0] == pmap[pri_n[n][j]][1])
          {
            pri_map[n] = -(pmap[pri_n[n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
            flag = 1;
          }    
        }
      }
    }
      
    for (n = 0; n < nhex; n++)
    {
      hex_map[n] = -1; //reset in case off proc tet
      if (hex_emap[n][0] == my_rank)
      {
        hex_map[n] = -fgnn;
      }
      else
      {
        flag = 0;
        for (j = 0; j < 8 && !flag; j++)
        {
          if (hex_emap[n][0] == pmap[hex_n[n][j]][1])
          {
            hex_map[n] = -(pmap[hex_n[n][j]][2] + 2); //simply use node from triangle, plus two to avoid -0 issues...if no node exists, will be detected (-1) and do a more thorough search on other end
            flag = 1;
          }    
        }
      }
    }
  }

  //
  // create hash tables for elements
  //
  nnhash = nn;
  tethash = new List*[nnhash];
  pyrhash = new List*[nnhash];
  prihash = new List*[nnhash];
  hexhash = new List*[nnhash];
  for (n=0; n < nnhash; n++)
  {
    tethash[n] = new List();
    pyrhash[n] = new List();
    prihash[n] = new List();
    hexhash[n] = new List();
  }
  
  //
  // create node-to-tetrahedra hash list
  //
  for (c=0; c < ntet; c++)
  {
    for (i=0; i < 4; i++)
    {
      n = tet_n[c][i];
      tethash[n]->Add_To_List(c);
    }
  }
  //
  // create node-to-pyramid hash list
  //
  for (c=0; c < npyr; c++)
  {
    for (i=0; i < 5; i++)
    {
      n = pyr_n[c][i];
      pyrhash[n]->Add_To_List(c);
    }
  }
  //
  // create node-to-prism hash list
  //
  for (c=0; c < npri; c++)
  {
    for (i=0; i < 6; i++)
    {
      n = pri_n[c][i];
      prihash[n]->Add_To_List(c);
    }
  }
  //
  // create node-to-hexahedra hash list
  //
  for (c=0; c < nhex; c++)
  {
    for (i=0; i < 8; i++)
    {
      n = hex_n[c][i];
      hexhash[n]->Add_To_List(c);
    }
  }

  if (my_rank == 0)
  {
    fprintf(out_f,"\nNode-to-element hash tables created.");
    fflush(out_f);
  }
  
  // set up maps
  // once this is done, DO NOT run update_emap for bd and reg elements, since there is no need for the index part of the p_i and the proc has been set properly
  // now that elements are created and mem allocated post compression for maps, create maps
  // note that the emaps will still be incomplete, but they are not necessary at this point
  // also, unlike in mapping nodes, since some elements are deleted, we need to begin anew with maps
  if (num_procs > 1)
  {
    //set up array to hold each range for each proc on each proc
    int **new_elem;
    new_elem = (int**)calloc((4 + 2*nb),sizeof(int*));
    for (n = 0; n < (4 + 2*nb); n++)
      new_elem[n] = (int*)calloc(num_procs,sizeof(int));
 
    for (n = 0; n < (4 + 2*nb); n++)
      for (p = 0; p < num_procs; p++)
        new_elem[n][p] = -1; //init for debug, since overwritten in gather
		
    //hold list of number of owned elements, based on emaps
    int *nelem_owned;
      nelem_owned = (int*)calloc((4 + 2*nb),sizeof(int));

    //fprintf(out_f,"\nInitialized arrays.\n");
    //fflush(out_f);

    //need to sift through and determine who will name which element (as long as unique, no matter who names)
    for (n = 0; n < 4; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (p = 0; p < ntet; p++)
          {
            if (tet_emap[p][0] == my_rank)
              nelem_owned[n]++;
          }
          break;
        case 1:
          for (p = 0; p < npyr; p++)
          {
            if (pyr_emap[p][0] == my_rank)
              nelem_owned[n]++;
          }
          break;
        case 2:
          for (p = 0; p < npri; p++)
          {
            if (pri_emap[p][0] == my_rank)
              nelem_owned[n]++;
          }
          break;
        case 3:
          for (p = 0; p < nhex; p++)
          {
            if (hex_emap[p][0] == my_rank)
              nelem_owned[n]++;
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than four basic element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    #ifdef _DEBUG
    fprintf(debug_f,"\nDetermined ownership of four basic types.\n");
    fflush(debug_f);
    #endif
        
    for (n = 0; n < 2; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nt[z]; p++)
            {
              if (tri_emap[z][p][0] == my_rank)                   
              {
                nelem_owned[4 + z]++;
              }
            }
          }
          break;
        case 1:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nq[z]; p++)
            {
              if (quad_emap[z][p][0] == my_rank)                   
              {
                nelem_owned[4 + nb + z]++;
              }
            }
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than two boundary element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }

    #ifdef _DEBUG
    fprintf(debug_f,"\nDetermined ownership of bd elem types.\n");
    fflush(debug_f);
    #endif 
  
    //now, we do an all gather so procs can determine range for new elem numbers
    for (n = 0; n < (4 + 2*nb); n++)
      MPI_Allgather(&(nelem_owned[n]),1,MPI_INT,new_elem[n],1,MPI_INT,MPI_COMM_WORLD);
      
    #ifdef _DEBUG
    /*for (n = 0; n < (4 + 2*nb); n++)
    {
      for (p = 0; p < num_procs; p++)
      {
        fprintf(debug_f,"new_elem[%d][%d] = %d\n",n,p,new_elem[n][p]);
      }
      fprintf(debug_f,"\n");
    }
    fflush(debug_f);*/
    #endif
    
    //set up for global max of each type for debug and ease
    int *glmax;
    glmax = (int*)calloc((4 + 2*nb),sizeof(int));

    //now, we can easily create global elem numbers for owned elem
    //get proc starting number
    for (n = 0; n < my_rank; n++)
      for (p = 0; p < (4 + 2*nb); p++)
        glmax[p] += new_elem[p][n]; 
  
    //need to give owned elem global numbers
    for (n = 0; n < 4; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (p = 0; p < ntet; p++)
          {
            if (tet_emap[p][0] == my_rank)
            {
              tet_map[p] = glmax[n];
              glmax[n]++;
            }
          }
          break;
        case 1:
          for (p = 0; p < npyr; p++)
          {
            if (pyr_emap[p][0] == my_rank)
            {
              pyr_map[p] = glmax[n];
              glmax[n]++;
            }
          }
          break;
        case 2:
          for (p = 0; p < npri; p++)
          {
            if (pri_emap[p][0] == my_rank)
            {
              pri_map[p] = glmax[n];
              glmax[n]++;
            }
          }
          break;
        case 3:
          for (p = 0; p < nhex; p++)
          {
            if (hex_emap[p][0] == my_rank)
            {
              hex_map[p] = glmax[n];
              glmax[n]++;
            }
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than four basic element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    //fprintf(out_f,"\nSet map for four basic types.\n");
    //fflush(out_f);
        
    //need to give owned elem global numbers
    for (n = 0; n < 2; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nt[z]; p++)
            {
              if (tri_emap[z][p][0] == my_rank)
              {
                tri_map[z][p] = glmax[4 + z];
                glmax[4 + z]++;
              }
            }
          }
          break;
        case 1:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nq[z]; p++)
            {
              if (quad_emap[z][p][0] == my_rank)
              {
                quad_map[z][p] = glmax[4 + nb + z];
                glmax[4 + nb + z]++;
              }
            }
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than two boundary element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    //fprintf(out_f,"\nSet map for boundary elements.\n");
    //fflush(out_f);
    
    //reset glmax first
    for (p = 0; p < (4 + 2*nb); p++)
      glmax[p] = 0;
        
    //now, use new_elem to find globalmax on all procs for debug    
    for (n = 0; n < num_procs; n++)
      for (p = 0; p < (4 + 2*nb); p++)
        glmax[p] += new_elem[p][n];  
      
    #ifdef _DEBUG
    /*for (n = 0; n < (4 + 2*nb); n++)
    {
      fprintf(debug_f,"\nglmax[%d] = %d, which num_elem, not index!\n",n,glmax[n]);
    }
    fflush(debug_f);*/
    #endif           

    if (new_elem > 0)
    {
      for (n = 0; n < (4 + 2*nb); n++)
        free(new_elem[n]);
      free(new_elem);
    }
    //will free glmax after check
    if (nelem_owned > 0)
      free(nelem_owned);

    //fprintf(out_f,"\nAbout to communicate new global element numbers.\n");
    //fflush(out_f);
        
    MPI_Barrier(MPI_COMM_WORLD);

    //finally, get the global number for elem owned by other procs
    sendcnt = new int [num_procs];
    recvcnt = new int [num_procs];
    for (n = 0; n < num_procs; n++)
    {
      sendcnt[n] = 0;
      recvcnt[n] = 0;
    }
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];

    //cycle thru elem	
    for (n = 0; n < 4; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (p = 0; p < ntet; p++)
          {
            if (tet_emap[p][0] != my_rank)
              sendcnt[tet_emap[p][0]]++;
          }
          break;
        case 1:
          for (p = 0; p < npyr; p++)
          {
            if (pyr_emap[p][0] != my_rank)
              sendcnt[pyr_emap[p][0]]++;
          }
          break;
        case 2:
          for (p = 0; p < npri; p++)
          {
            if (pri_emap[p][0] != my_rank)
              sendcnt[pri_emap[p][0]]++;
          }
          break;
        case 3:
          for (p = 0; p < nhex; p++)
          {
            if (hex_emap[p][0] != my_rank)
              sendcnt[hex_emap[p][0]]++;
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than four basic element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    //cycle thru elem	
    for (n = 0; n < 2; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (z = 0; z < nb; z++)
            for (p = 0; p < nt[z]; p++)
              if (tri_emap[z][p][0] != my_rank)
                sendcnt[tri_emap[z][p][0]]++;
          break;
        case 1:
          for (z = 0; z < nb; z++)
            for (p = 0; p < nq[z]; p++)
              if (quad_emap[z][p][0] != my_rank)
                sendcnt[quad_emap[z][p][0]]++;
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than two boundary element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }

    //fprintf(out_f,"\nAbout to send counts.\n");
    //fflush(out_f);

    //now, send and recv count for all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;
 
      MPI_Isend(&(sendcnt[p]),1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
      nreq_s++;
 
      MPI_Irecv(&(recvcnt[p]),1,MPI_INT,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
      nreq_r++;
    }
 
    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    delete [] srequest;
    delete [] rrequest;
    delete [] statuses;
 
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];
    
    #ifdef _DEBUG
    /*for (n = 0 ; n < num_procs; n++)
      fprintf(debug_f,"\n recvcnt[%d] = %d \n",n,recvcnt[n]);
    fflush(debug_f);
        
    for (n = 0 ; n < num_procs; n++)
      fprintf(debug_f,"\n sendcnt[%d] = %d \n",n,sendcnt[n]);
    fflush(debug_f);*/
    #endif

    //allocate buffers to send/recv elem/node indices
    bdim2 = new int[num_procs];
    for (n = 0; n < num_procs; n++)
      bdim2[n] = MAX(sendcnt[n],recvcnt[n])*11*sizeof(int); //indices (3 qualifiers + up to 8 nodes)
    sbuff = new char*[num_procs];
    rbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
    {
      sbuff[n] = new char[bdim2[n]];
      rbuff[n] = new char[bdim2[n]];
    }

    //another benefit here is they can be used on all procs the same
    int *trit = 0;
    int *quadt = 0;
        
    trit = (int*)calloc(nb,sizeof(int));
    quadt = (int*)calloc(nb,sizeof(int));
        
    q = -50; //start
    //init to higher than four basic
    for (n = 0; n < nb; n++)
    {
      trit[n] = q;
      q-=10;
    }
    for (n = 0; n < nb; n++)
    {
      quadt[n] = q;
      q-=10;
    }

    #ifdef _DEBUG
    fprintf(debug_f,"\nAbout to pack elem info.\n");
    fflush(debug_f);
    #endif

    //now, pack elem indices as per procs that will be sending them
    for (q = 0; q < num_procs; q++)
    {
        
      if (sendcnt[q] == 0 || q == my_rank)
        continue;

      //set position to 0
      sposition=0;
 
      //pack elem type, index (for return) and nodes
      for (n = 0; n < 4; n++)
      {
          
        //loop over elements
        switch(n)
        {
          case 0:
            for (p = 0; p < ntet; p++)
            {
              //fprintf(out_f,"\nPacking tets %d of %d - %d.\n",p,0,ntet);
              //fflush(out_f);
              if (tet_emap[p][0] == q && tet_map[p] < 0)
              {
                MPI_Pack(&(tt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(tet_map[p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);        
                for (j = 0; j < 4; j++)
                  MPI_Pack(&(pmap[tet_n[p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
              }
              //fprintf(out_f,"\nPACK:Type = %d, tet = %d, GN = %d %d %d %d\n",tt,p,pmap[tet_n[p][0]][0],pmap[tet_n[p][1]][0],pmap[tet_n[p][2]][0],pmap[tet_n[p][3]][0]);
              //fflush(out_f);
            }
            break;
          case 1:
            for (p = 0; p < npyr; p++)
            {
              //fprintf(out_f,"\nPacking pyrs %d of %d - %d.\n",p,0,npyr);
              //fflush(out_f);
              if (pyr_emap[p][0] == q && pyr_map[p] < 0)
              {
                MPI_Pack(&(pyt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pyr_map[p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                for (j = 0; j < 5; j++)
                  MPI_Pack(&(pmap[pyr_n[p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
              }
              //fprintf(out_f,"\nPACK:Type = %d, pyr = %d, GN = %d %d %d %d\n",pyt,p,pmap[pyr_n[p][0]][0],pmap[pyr_n[p][1]][0],pmap[pyr_n[p][2]][0],pmap[pyr_n[p][3]][0]);
              //fflush(out_f);
            }
            break;
          case 2:
            for (p = 0; p < npri; p++)
            {
              //fprintf(out_f,"\nPacking pris %d of %d - %d.\n",p,0,npri);
              //fflush(out_f);
              if (pri_emap[p][0] == q && pri_map[p] < 0)
              {
                MPI_Pack(&(prt),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pri_map[p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                for (j = 0; j < 6; j++)
                  MPI_Pack(&(pmap[pri_n[p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
              }
              //fprintf(out_f,"\nPACK:Type = %d, pri = %d, GN = %d %d %d %d\n",prt,p,pmap[pri_n[p][0]][0],pmap[pri_n[p][1]][0],pmap[pri_n[p][2]][0],pmap[pri_n[p][3]][0]);
              //fflush(out_f);
            }
            break;
          case 3:
            for (p = 0; p < nhex; p++)
            {
              //fprintf(out_f,"\nPacking hexs %d of %d - %d.\n",p,0,nhex);
              //fflush(out_f);
              if (hex_emap[p][0] == q && hex_map[p] < 0)
              {
                MPI_Pack(&(ht),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(hex_map[p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                for (j = 0; j < 8; j++)
                  MPI_Pack(&(pmap[hex_n[p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
              }
              //fprintf(out_f,"\nPACK:Type = %d, hex = %d, GN = %d %d %d %d\n",ht,p,pmap[hex_n[p][0]][0],pmap[hex_n[p][1]][0],pmap[hex_n[p][2]][0],pmap[hex_n[p][3]][0]);
              //fflush(out_f);
            }
            break;
          default:
            fprintf(stderr,"\nSomehow you have more than four basic element types!  Exiting....\n");
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
            break;
        }
      }
          
      //fprintf(out_f,"\nAbout to pack boundary elements.\n");
      //fflush(out_f);
          
      //pack elem type, index (for return) and nodes
      for (n = 0; n < 2; n++)
      {
        //loop over elements
        switch(n)
        {
          case 0:
            for (z = 0; z < nb; z++)
            {
              for (p = 0; p < nt[z]; p++)
              {
                //fprintf(out_f,"\nPacking tri %d of %d - %d, bnd %d.\n",p,0,nt[z],z);
                //fflush(out_f);
                if (tri_emap[z][p][0] == q)
                {
                  MPI_Pack(&(trit[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  //MPI_Pack(&(tri_map[z][p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  for (j = 0; j < 3; j++)
                    MPI_Pack(&(pmap[t_n[z][p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                }
                //fprintf(out_f,"\nPACK:Type = %d, tri/bd = %d/%d, GN = %d %d %d\n",trit[z],p,z,pmap[t_n[z][p][0]][0],pmap[t_n[z][p][1]][0],pmap[t_n[z][p][2]][0]);
                //fflush(out_f);
              }
            }
            break;
          case 1:
            for (z = 0; z < nb; z++)
            {
              for (p = 0; p < nq[z]; p++)
              {
                //fprintf(out_f,"\nPacking quad %d of %d - %d, bnd %d.\n",p,0,nq[z],z);
                //fflush(out_f);
                if (quad_emap[z][p][0] == q)
                {
                  MPI_Pack(&(quadt[z]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(p),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  //MPI_Pack(&(quad_map[z][p]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                  for (j = 0; j < 4; j++)
                    MPI_Pack(&(pmap[q_n[z][p][j]][0]),1,MPI_INT,sbuff[q],bdim2[q],&sposition,MPI_COMM_WORLD);
                }
                //fprintf(out_f,"\nPACK:Type = %d, quad/bd = %d/%d, GN = %d %d %d %d\n",quadt[z],p,z,pmap[q_n[z][p][0]][0],pmap[q_n[z][p][1]][0],pmap[q_n[z][p][2]][0],pmap[q_n[z][p][3]][0]);
                //fflush(out_f);  
              }
            }
            break;
          default:
            fprintf(stderr,"\nSomehow you have more than two boundary element types!  Exiting....\n");
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
            break;
        }
      }  
    }

    free(trit);
    free(quadt);

    MPI_Barrier(MPI_COMM_WORLD);
 
    #ifdef _DEBUG
    fprintf(debug_f,"\nAbout to send elem info.\n");
    fflush(debug_f);
    #endif

    //now, send and recv packets from all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (sendcnt[p] > 0 && p != my_rank) 
      {
        MPI_Isend(sbuff[p],bdim2[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }

      if (recvcnt[p] > 0 && p != my_rank) 
      {
        MPI_Irecv(rbuff[p],bdim2[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }

    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);
 
    delete [] srequest;
    delete [] rrequest;
    delete [] statuses;
 
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses = new MPI_Status[num_procs];
 
    //now, delete sbuff and resize
    for (n = 0; n < num_procs; n++)
      if (sbuff[n] != 0) delete [] sbuff[n];
    delete [] sbuff;
 
    //resize
    bdim21 = new int[num_procs];
    for (n = 0; n < num_procs; n++)
      bdim21[n] = MAX(sendcnt[n],recvcnt[n])*3*sizeof(int); //nodes and indices
    sbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
      sbuff[n] = new char[bdim21[n]];
 
    temp_node = temp_node2 = 0;

    int el = 0; //will hold matching element
       
    //will hold nodes needed to check against
    List check;
      
    check.Redimension(0);
 
    MPI_Barrier(MPI_COMM_WORLD);

    #ifdef _DEBUG
    fprintf(debug_f,"\nAbout to unpack four basic types.\n");
    fflush(debug_f);
    #endif
 
    //now, unpack type, index, nodes, replace with type, index, global index
    for (p = 0; p < num_procs; p++)
    {
      if (recvcnt[p] == 0 || p == my_rank)
        continue;                      
      sposition = rposition = 0; //reset pos
      for (n = 0; n < recvcnt[p]; n++)
      {
        //unpack type initially
        MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);            
        //switch
        switch(temp_node)
        {
          case -10:

            //unpack index to send back
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);

            //fprintf(debug_f,"\nReturn index: %d\n",temp_node2);
            //fflush(debug_f);

            //unpack node in element
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);

            //fprintf(debug_f,"\nNode in element: %d\n",el);
            //fflush(debug_f);
                                
            //unpack other nodes into list
            for (q = 0; q < 4; q++)
            {
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
              check.Add_To_List(i);

              //fprintf(debug_f,"\nOther nodes: %d\n",i);
              //fflush(debug_f);
            }
            //need to deal with fact that occasionally we will have no nodes owned by owning elem proc  
            if (el == -1)
            {
              for (q = 0; q < nn; q++)
                if (pmap[q][0] == check.list[0])
                  el = q;

              //fprintf(debug_f,"\nOff-proc el/ node: %d\n",el);
              //fflush(debug_f);
            }
            else
            {
              el = -(el + 2); //change sign, as no longer needed to parse new from old 

              //fprintf(debug_f,"\nNode on proc: %d\n",el);
              //fflush(debug_f);
            }
            //set el as hash
            flag = 0;
            for (q = 0; q < tethash[el]->max && !flag; q++)
            {
              i = 0; //reset
              for (j = 0; j < 4; j++)
              {
                if (check.Is_In_List(pmap[tet_n[tethash[el]->list[q]][j]][0]))
                  i++;
                //fprintf(debug_f,"\ntet %d, node %d-%d (%d)\n",tethash[el]->list[q],pmap[tet_n[tethash[el]->list[q]][j]][0],pmap[tet_n[tethash[el]->list[q]][j]][1],pmap[tet_n[tethash[el]->list[q]][j]][2]);
                //fflush(debug_f);
              }
              if (i == 4)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",tt,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(tt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(tet_map[tethash[el]->list[q]]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                flag = 1;
              }
            }
                
            if (i < 4)
            {
              fprintf(stderr,"\nWARNING: You have a bad tet element sent that cannot be found!\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            check.Redimension(0);
            break;
          case -20:

            //unpack index to send back
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack node in element
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);

            //unpack other nodes into list
            for (q = 0; q < 6; q++)
            {
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
              check.Add_To_List(i);
            }

            //need to deal with fact that occasionally we will have no nodes owned by owning elem proc  
            if (el == -1)
            {
              for (q = 0; q < nn; q++)
                if (pmap[q][0] == check.list[0])
                  el = q;
            }
            else
            {
              el = -(el + 2); //change sign, as no longer needed to parse new from old 
            }
                  
            //set el as hash
            flag = 0;
            for (q = 0; q < prihash[el]->max && !flag; q++)
            {
              i = 0; //reset
              for (j = 0; j < 6; j++)
              {
                if (check.Is_In_List(pmap[pri_n[prihash[el]->list[q]][j]][0]))
                  i++;
              }
              if (i == 6)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",tt,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(prt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pri_map[prihash[el]->list[q]]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                flag = 1;
              }
            }
                
            if (i < 6)
            {
              fprintf(stderr,"\nWARNING: You have a bad pri element sent that cannot be found!\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            check.Redimension(0);
            break;
          case -30:

            //unpack index to send back
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack node in element
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);

            //unpack other nodes into list
            for (q = 0; q < 5; q++)
            {
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
              check.Add_To_List(i);
            }
                  
            //need to deal with fact that occasionally we will have no nodes owned by owning elem proc  
            if (el == -1)
            {
              for (q = 0; q < nn; q++)
                if (pmap[q][0] == check.list[0])
                  el = q;
            }
            else
            {
              el = -(el + 2); //change sign, as no longer needed to parse new from old 
            }
                  
            //set el as hash
            flag = 0;
            for (q = 0; q < pyrhash[el]->max && !flag; q++)
            {
              i = 0; //reset
              for (j = 0; j < 5; j++)
              {
                if (check.Is_In_List(pmap[pyr_n[pyrhash[el]->list[q]][j]][0]))
                  i++;
              }
              if (i == 5)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",tt,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(pyt),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(pyr_map[pyrhash[el]->list[q]]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                flag = 1;
              }
            }
               
            if (i < 5)
            {
              fprintf(stderr,"\nWARNING: You have a bad pyr element sent that cannot be found!\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            check.Redimension(0);
            break;
          case -40:

            //unpack index to send back
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack node in element
            MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);

            //unpack other nodes into list
            for (q = 0; q < 8; q++)
            {
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
              check.Add_To_List(i);
            }

            //need to deal with fact that occasionally we will have no nodes owned by owning elem proc  
            if (el == -1)
            {
              for (q = 0; q < nn; q++)
                if (pmap[q][0] == check.list[0])
                  el = q;
            }
            else
            {
              el = -(el + 2); //change sign, as no longer needed to parse new from old 
            }
                 
            //set el as hash
            flag = 0;
            for (q = 0; q < hexhash[el]->max && !flag; q++)
            {
              i = 0; //reset
              for (j = 0; j < 8; j++)
              {
                if (check.Is_In_List(pmap[hex_n[hexhash[el]->list[q]][j]][0]))
                  i++;
              }
              if (i == 8)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",tt,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(ht),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(hex_map[hexhash[el]->list[q]]),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                flag = 1;
              }
            }

            if (i < 8)
            {
              fprintf(stderr,"\nWARNING: You have a bad hex element sent that cannot be found!\n");
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            check.Redimension(0);
            break;
          //will use this for boundary element types
          default:
            //NOTE:  we are doing a global search since a bd hash would use lots of mem and bd are split and small
            //also, we already needed the elem hashes for make_nabor
            //if this is too slow, we can make a bd hash per nn (letting interior nodes be wasted) with elem indices incremented by nt[b..] then take the results and cut off nt[b..]
            //fprintf(out_f,"\nUnpacking bd elem.\n");
            //fflush(out_f);
            if ((temp_node <= -50) && (temp_node > (-50 - 10*nb)))
            {
              //unpack index to send back
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //fprintf(out_f,"\ntype = %d, elem = %d\n",temp_node,temp_node2);
              //fflush(out_f);
              for (q = 0; q < 3; q++)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                //fprintf(out_f,"\nnode = %d\n",el);
                //fflush(out_f);
                check.Add_To_List(el);
              }
              //check.print(out_f);
              flag = 0;
              for (q = 0; q < nt[(-50 - temp_node)/10] && !flag; q++)
              {
                el = 0; //reset
                for (z = 0; z < 3; z++)
                {
                  if (check.Is_In_List(pmap[t_n[(-50 - temp_node)/10][q][z]][0]))
                  {
                    el++;
                  }
                }
                if (el == 3)
                {
                  //fprintf(out_f,"\nel = %d, q = %d\n",el,q);
                  //fflush(out_f);
                  el = tri_map[(-50 - temp_node)/10][q];
                  flag = 1;
                  //fprintf(out_f,"\nfound element %d\n",el);
                  //fflush(out_f);
                } 
                else if (el != 3 && q == nt[(-50 - temp_node)/10]-1)
                {  
                  fprintf(stderr,"\nWARNING: You have a bad tri bd element sent that cannot be found!\n");
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
              }
              if (flag)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",temp_node,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              check.Redimension(0); //reset list
            }
            else if ((temp_node <= (-50 - 10*nb)) && (temp_node > (-50 - 20*nb)))
            {
              //unpack index to send back
              MPI_Unpack(rbuff[p],bdim2[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              for (q = 0; q < 4; q++)
              {
                MPI_Unpack(rbuff[p],bdim2[p],&rposition,&el,1,MPI_INT,MPI_COMM_WORLD);
                check.Add_To_List(el);
              }
              flag = 0;
              for (q = 0; q < nq[(-50 - 10*nb - temp_node)/10] && !flag; q++)
              {
                el = 0; //reset
                for (z = 0; z < 4; z++)
                {
                  if (check.Is_In_List(pmap[q_n[(-50 - 10*nb - temp_node)/10][q][z]][0]))
                  {
                    el++;
                  }
                }
                if (el == 4)
                {
                  el = quad_map[(-50 - 10*nb - temp_node)/10][q];
                  flag = 1;
                } 
                else if (el != 4 && q == nq[(-50 - 10*nb - temp_node)/10]-1)
                {  
                  fprintf(stderr,"\nWARNING: You have a bad quad bd element sent that cannot be found!\n");
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
              }
              if (flag)
              {
                //fprintf(out_f,"\npacking %d %d %d\n",temp_node,temp_node2,el);
                //fflush(out_f);
                MPI_Pack(&(temp_node),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(temp_node2),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(el),1,MPI_INT,sbuff[p],bdim21[p],&sposition,MPI_COMM_WORLD);
              }
              check.Redimension(0); //reset list
            }
            break;
        }
      }
    }
        
    #ifdef _DEBUG
    fprintf(debug_f,"\nFinished unpack/repack of elem.\n");
    fflush(debug_f);
    #endif
        
    MPI_Barrier(MPI_COMM_WORLD);
        
    //fprintf(out_f,"\nAbout to resize rbuff.\n");
    //fflush(out_f);

    //now, delete rbuff and bdim2, resize
    delete [] bdim2;
    for (n = 0; n < num_procs; n++)
      if (rbuff[n] != 0) delete [] rbuff[n];
    delete [] rbuff;

    //resize
    rbuff = new char*[num_procs];
    for (n = 0; n < num_procs; n++)
      rbuff[n] = new char[bdim21[n]];
  	
    //fprintf(out_f,"\nAbout to send new numbers.\n");
    //fflush(out_f);
    
    //now, send and recv packets from all procs
    nreq_s = nreq_r = 0;
    for (p = 0; p < num_procs; p++)
    {
      if (recvcnt[p] > 0 && p != my_rank) 
      {
        MPI_Isend(sbuff[p],bdim21[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        nreq_s++;
      }

      if (sendcnt[p] > 0 && p != my_rank) 
      {
        MPI_Irecv(rbuff[p],bdim21[p],MPI_CHAR,p,p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_r++;
      }
    }
  
    //now, wait for finish before unpacking
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    //fprintf(out_f,"\nAbout to unpack four basic types.\n");
    //fflush(out_f);

    //finally, unpack nodes in proper position
    for (p = 0; p < num_procs; p++)
    {
      if (sendcnt[p] == 0 || p == my_rank)
        continue;
        
      rposition = 0; //reset pos
      for (n = 0; n < sendcnt[p]; n++)
      {
        //unpack type initially
        MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node,1,MPI_INT,MPI_COMM_WORLD);

        //fprintf(out_f,"\nunpacking type = %d\n",temp_node);
        //fflush(out_f);

        //switch
        switch(temp_node)
        {
          case -10:
            //unpack index use
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack global elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(tet_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
            //fprintf(out_f,"\nunpacking LEN = %d, GEN = %d\n",temp_node2,tet_map[temp_node2]);
            //fflush(out_f);
            break;
          case -20:
            //unpack index use
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack global elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(pri_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
            //fprintf(out_f,"\nunpacking LEN = %d, GEN = %d\n",temp_node2,pri_map[temp_node2]);
            //fflush(out_f);
            break;
          case -30:
            //unpack index use
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack global elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(pyr_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
            //fprintf(out_f,"\nunpacking LEN = %d, GEN = %d\n",temp_node2,pyr_map[temp_node2]);
            //fflush(out_f);
            break;
          case -40:
            //unpack index use
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
            //unpack global elem number
            MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(hex_map[temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
            //fprintf(out_f,"\nunpacking LEN = %d, GEN = %d\n",temp_node2,hex_map[temp_node2]);
            //fflush(out_f);
            break;
          //will use this for boundary element types
          default:
            //fprintf(out_f,"\nUnpacking a boundary element.\n");
            //fflush(out_f);
            if ((temp_node <= -50) && (temp_node > (-50 - 10*nb)))
            {
              //unpack index to use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack global elem number
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(tri_map[(-50 - temp_node)/10][temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, b = %d, GEN = %d\n",temp_node2,(-50 - temp_node)/10,tri_map[(-50 - temp_node)/10][temp_node2]);
              //fflush(out_f);
            }
            if ((temp_node <= (-50 - 10*nb)) && (temp_node > (-50 - 20*nb)))
            {
              //unpack index to use
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&temp_node2,1,MPI_INT,MPI_COMM_WORLD);
              //unpack global elem number
              MPI_Unpack(rbuff[p],bdim21[p],&rposition,&(quad_map[(-50 - 10*nb - temp_node)/10][temp_node2]),1,MPI_INT,MPI_COMM_WORLD);
              //fprintf(out_f,"\nunpacking LEN = %d, b = %d, GEN = %d\n",temp_node2,(-50 - 10*nb - temp_node)/10,quad_map[(-50 - 10*nb - temp_node)/10][temp_node2]);
              //fflush(out_f);
            }
            break;
        }
      }
    }	

    //finally, free MPI mem
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] srequest;
    delete[] rrequest;
    delete[] statuses;
    for (n = 0; n < num_procs; n++)
    {
      if (sbuff[n] != 0) delete [] sbuff[n];
      if (rbuff[n] != 0) delete [] rbuff[n];
    }
    delete[] sbuff;
    delete[] rbuff;
    delete[] bdim21;

    //check for negatives and numbers too large
    for (n = 0; n < 4; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (p = 0; p < ntet; p++)
          {
            if (tet_map[p] < 0 || tet_map[p] >= glmax[n])
            {
              fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has map %d which is out of range.\n",p,tet_map[p]);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            for (q = 0; q < 4; q++)
            {
              //already checked all nodes, so looking for negative indices or nodes that don't match up
              if (tet_n[p][q] < 0)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has node %d(L = %d) which is -1.\n",p,q,tet_n[p][q]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              if (pmap[tet_n[p][q]][0] < 0 || pmap[tet_n[p][q]][0] >= globalmax)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: tet %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,tet_n[p][q],pmap[tet_n[p][q]][0]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
            }
          }
          break;
        case 1:
          for (p = 0; p < npyr; p++)
          {
            if (pyr_map[p] < 0 || pyr_map[p] >= glmax[n])
            {
              fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has map %d which is out of range.\n",p,pyr_map[p]);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            for (q = 0; q < 5; q++)
            {
              //already checked all nodes, so looking for negative indices or nodes that don't match up
              if (pyr_n[p][q] < 0)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has node %d(L = %d) which is -1.\n",p,q,pyr_n[p][q]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              if (pmap[pyr_n[p][q]][0] < 0 || pmap[pyr_n[p][q]][0] >= globalmax)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: pyr %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,pyr_n[p][q],pmap[pyr_n[p][q]][0]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
            }
          }
          break;
        case 2:
          for (p = 0; p < npri; p++)
          {
            if (pri_map[p] < 0 || pri_map[p] >= glmax[n])
            {
              fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has map %d which is out of range.\n",p,pri_map[p]);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            for (q = 0; q < 6; q++)
            {
              //already checked all nodes, so looking for negative indices or nodes that don't match up
              if (pri_n[p][q] < 0)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has node %d(L = %d) which is -1.\n",p,q,pri_n[p][q]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              if (pmap[pri_n[p][q]][0] < 0 || pmap[pri_n[p][q]][0] >= globalmax)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: pri %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,pri_n[p][q],pmap[pri_n[p][q]][0]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
            }
          }
          break;
        case 3:
          for (p = 0; p < nhex; p++)
          {
            if (hex_map[p] < 0 || hex_map[p] >= glmax[n])
            {
              fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has map %d which is out of range.\n",p,hex_map[p]);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            for (q = 0; q < 8; q++)
            {
              //already checked all nodes, so looking for negative indices or nodes that don't match up
              if (hex_n[p][q] < 0)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has node %d(L = %d) which is -1.\n",p,q,hex_n[p][q]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              if (pmap[hex_n[p][q]][0] < 0 || pmap[hex_n[p][q]][0] >= globalmax)
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: hex %d has node %d(L = %d) which is out of range globally @ %d.\n",p,q,hex_n[p][q],pmap[hex_n[p][q]][0]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
            }
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than four basic element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    //bd elem
    for (n = 0; n < 2; n++)
    {
      //loop over elements
      switch(n)
      {
        case 0:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nt[z]; p++)
            {
              if (tri_map[z][p] < 0 || tri_map[z][p] >= glmax[4 + z])
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has map %d which is out of range.\n",p,z,tri_map[z][p]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              for (q = 0; q < 3; q++)
              {
                //already checked all nodes, so looking for negative indices or nodes that don't match up
                if (t_n[z][p][q] < 0)
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has node %d(L = %d) which is -1.\n",p,z,q,t_n[z][p][q]);
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
                if (pmap[t_n[z][p][q]][0] < 0 || pmap[t_n[z][p][q]][0] >= globalmax)
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: tri %d on %d has node %d(L = %d) which is out of range globally @ %d.\n",p,z,q,t_n[z][p][q],pmap[t_n[z][p][q]][0]);
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
              }
            }
          }
          break;
        case 1:
          for (z = 0; z < nb; z++)
          {
            for (p = 0; p < nq[z]; p++)
            {
              if (quad_map[z][p] < 0 || quad_map[z][p] >= glmax[4 + nb + z])
              {
                fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has map %d which is out of range.\n",p,z,quad_map[z][p]);
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,0);
                exit(0);
              }
              for (q = 0; q < 4; q++)
              {
                //already checked all nodes, so looking for negative indices or nodes that don't match up
                if (q_n[z][p][q] < 0)
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has node %d(L = %d) which is -1.\n",p,z,q,q_n[z][p][q]);
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
                if (pmap[q_n[z][p][q]][0] < 0 || pmap[q_n[z][p][q]][0] >= globalmax)
                {
                  fprintf(stderr,"\nCATASTROPHIC ERROR: quad %d on %d has node %d(L = %d) which is out of range globally @ %d.\n",p,z,q,q_n[z][p][q],pmap[q_n[z][p][q]][0]);
                  fflush(stderr);
                  MPI_Abort(MPI_COMM_WORLD,0);
                  exit(0);
                }
              }
            }
          }
          break;
        default:
          fprintf(stderr,"\nSomehow you have more than two boundary element types!  Exiting....\n");
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
          break;
      }
    }
        
    //since used for debug, freed here
    if (glmax > 0)
      free(glmax);
    if (oldnt > 0)
      free(oldnt);
    if (oldnq > 0)
      free(oldnq);

    #ifdef _DEBUG
    /*fprintf(debug_f,"\nLooking for duplicate bd elements1!\n");
    fflush(debug_f);
    List trifind;
    for (b = 0; b < nb; b++)
    {
      for (n = 0; n < nq[b]; n++)
      {
        trifind.Redimension(0);
        trifind.Add_To_List(pmap[q_n[b][n][0]][0]);
        trifind.Add_To_List(pmap[q_n[b][n][1]][0]);
        trifind.Add_To_List(pmap[q_n[b][n][2]][0]);
        trifind.Add_To_List(pmap[q_n[b][n][3]][0]);
        m = 0;
        for (i = 0; i < nt[b] && m < 2; i++)
        {
          k = 0; //cntr
          for (j = 0; j < 3; j++)
          {
            if (trifind.Is_In_List(pmap[t_n[b][i][j]][0]))
              k++;
          }
          if (k == 3)
          {
            fprintf(debug_f,"\nNODALLY FOUND ONE TRIANGLE %d (%d - %d) WITH THREE CONCURRENT NODES WITH QUADRILATERAL %d (%d - %d)!",i,tri_emap[b][i][0],tri_emap[b][i][1],n,quad_emap[b][n][0],quad_emap[b][n][1]);
            fflush(debug_f);
            m++;
          }
        }
        if (m == 2)
        {
          fprintf(debug_f,"\nNODALLY FOUND TWO TRIANGLES COVERING QUADRILATERAL %d (%d-%d)!",n,quad_emap[b][n][0],quad_emap[b][n][1]);
          fflush(debug_f);
        }
      }
    }
    List quadfind;
    for (b = 0; b < nb; b++)
    {
      for (n = 0; n < nq[b]; n++)
      {
        quadfind.Redimension(0);
        quadfind.Add_To_List(pmap[q_n[b][n][0]][0]);
        quadfind.Add_To_List(pmap[q_n[b][n][1]][0]);
        quadfind.Add_To_List(pmap[q_n[b][n][2]][0]);
        quadfind.Add_To_List(pmap[q_n[b][n][3]][0]);
        for (i = 0; i < nq[b]; i++)
        {
          if (i == n)
            continue;
          k = 0; //cntr
          for (j = 0; j < 4; j++)
          {
            if (quadfind.Is_In_List(pmap[q_n[b][i][j]][0]))
              k++;
          }
          if (k == 4)
          {
            fprintf(debug_f,"\nCOORDINATELY FOUND QUAD %d (%d - %d) WITH SAME NODES AS QUAD %d (%d - %d)!",i,quad_emap[b][i][0],quad_emap[b][i][1],n,quad_emap[b][n][0],quad_emap[b][n][1]);
            fflush(debug_f);
          }
          k = 0; //cntr
          for (j = 0; j < 4; j++)
          {
            if (node[q_n[b][i][j]][0] == node[q_n[b][n][j]][0] && node[q_n[b][i][j]][1] == node[q_n[b][n][j]][1] && node[q_n[b][i][j]][2] == node[q_n[b][n][j]][2] && node[q_n[b][i][j]][3] == node[q_n[b][n][j]][3])
              k++;
          }
          if (k == 4)
          {
            fprintf(debug_f,"\nNODALLY FOUND QUAD %d (%d - %d) WITH SAME NODES AS QUAD %d (%d - %d)!",i,quad_emap[b][i][0],quad_emap[b][i][1],n,quad_emap[b][n][0],quad_emap[b][n][1]);
            fflush(debug_f);
          }
        }
      }
    }*/
    #endif
    
    //DON'T re-run, to avoid heredity problems in neighbor check
    //element_pi();
  }

  //need emaps and maps for this as a DEBUG CHECK...hashes created prior to maps for use there.
  if (num_procs == 1)
    make_nabors(tethash,pyrhash,prihash,hexhash,-1);

  // delete current hash tables
  for (n=0; n < nnhash; n++)
  {
    delete tethash[n];
    delete pyrhash[n];
    delete prihash[n];
    delete hexhash[n];
  }
  delete[] tethash;
  delete[] pyrhash;
  delete[] prihash;
  delete[] hexhash;

  return(1);
}

//now works in parallel
int mesh_obj::make_nabors( List **tethash, List **pyrhash, List **prihash, List **hexhash, int timing)
{
  int b, c, i, j, k, m, n, n0, n1, n2, n3, n4, n5, n6, n7, t, q, s, glob, glob2;
  int ltet, lpyr, lpri, lhex, nnhash;
  int conn[9];
  
  if (my_rank == 0)
  {
    fprintf(out_f,"\nCreating neighbor connectivity to check mesh.");
    fflush(out_f);
  }

  // create neighbor connectivity to check mesh
  int (*tet_nbr)[4], (*pyr_nbr)[5], (*pri_nbr)[5], (*hex_nbr)[6];
  tet_nbr = new int[ntet][4];
  pyr_nbr = new int[npyr][5];
  pri_nbr = new int[npri][5];
  hex_nbr = new int[nhex][6];

  int ntri = 0;
  int nquad= 0;

  for (c=0; c < ntet; c++)
  {
    n3 = -1;
    for (i=0; i < 4; i++)
    {
      tet_nbr[c][i] = -1;
      switch (i)
      {
        case 0: n0 = tet_n[c][0]; n1 = tet_n[c][2]; n2 = tet_n[c][1]; break;
        case 1: n0 = tet_n[c][0]; n1 = tet_n[c][1]; n2 = tet_n[c][3]; break;
        case 2: n0 = tet_n[c][1]; n1 = tet_n[c][2]; n2 = tet_n[c][3]; break;
        case 3: n0 = tet_n[c][2]; n1 = tet_n[c][0]; n2 = tet_n[c][3]; break;
      }

      ltet=lpyr=lpri=lhex=-1;
      ltet=c;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (ltet >= 0)
      {
        if ((s = tet_side(tet_n[ltet],n0,n2,n1)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for tet %d, side %d to tet %d",c,i,ltet);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          tet_nbr[c][i] = ltet;
      } else if (lpyr >= 0)
      {
        if ((s = pyramid_side(pyr_n[lpyr],n0,n2,n1,n3)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for tet %d, side %d to pyramid %d",c,i,lpyr);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          tet_nbr[c][i] = ntet+lpyr;
      } else if (lpri >= 0)
      {
        if ((s = prism_side(pri_n[lpri],n0,n2,n1,n3)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for tet %d, side %d to prism %d",c,i,lpri);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          tet_nbr[c][i] = ntet+npyr+lpri;
      } else
        ntri++;
    }
  }
  for (c=0; c < npyr; c++)
  {
    for (i=0; i < 5; i++)
    {
      pyr_nbr[c][i] = -1;
      switch (i)
      {
        case 0: n0 = pyr_n[c][0]; n1 = pyr_n[c][3]; n2 = pyr_n[c][2]; n3 = pyr_n[c][1]; break;
        case 1: n0 = pyr_n[c][0]; n1 = pyr_n[c][1]; n2 = pyr_n[c][4]; n3 = -1; break;
        case 2: n0 = pyr_n[c][1]; n1 = pyr_n[c][2]; n2 = pyr_n[c][4]; n3 = -1; break;
        case 3: n0 = pyr_n[c][2]; n1 = pyr_n[c][3]; n2 = pyr_n[c][4]; n3 = -1; break;
        case 4: n0 = pyr_n[c][3]; n1 = pyr_n[c][0]; n2 = pyr_n[c][4]; n3 = -1; break;
      }

      ltet=lpyr=lpri=lhex=-1;
      lpyr=c;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (n3 >= 0)
      {
        if (lpyr >= 0)
        {
          if ((s = pyramid_side(pyr_n[lpyr],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to pyramid %d",c,i,lpyr);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ntet+lpyr;
        } else if (lpri >= 0)
        {
          if ((s = prism_side(pri_n[lpri],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to prism %d",c,i,lpri);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ntet+npyr+lpri;
        } else if (lhex >= 0)
        {
          if ((s = hex_side(hex_n[lhex],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to hex %d",c,i,lhex);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ntet+npyr+npri+lhex;
        } else
          nquad++;
      } else
      {
        if (ltet >= 0)
        {
          if ((s = tet_side(tet_n[ltet],n0,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to tet %d",c,i,ltet);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ltet;
        } else if (lpyr >= 0)
        {
          if ((s = pyramid_side(pyr_n[lpyr],n0,n2,n1,n3)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to pyramid %d",c,i,lpyr);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ntet+lpyr;
        } else if (lpri >= 0)
        {
          if ((s = prism_side(pri_n[lpri],n0,n2,n1,n3)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for pyramid %d, side %d to prism %d",c,i,lpri);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pyr_nbr[c][i] = ntet+npyr+lpri;
        } else
          ntri++;
      }
    }
  }
  for (c=0; c < npri; c++)
  {
    for (i=0; i < 5; i++)
    {
      pri_nbr[c][i] = -1;
      switch (i)
      {
        case 0: n0 = pri_n[c][0]; n1 = pri_n[c][1]; n2 = pri_n[c][4]; n3 = pri_n[c][3]; break;
        case 1: n0 = pri_n[c][1]; n1 = pri_n[c][2]; n2 = pri_n[c][5]; n3 = pri_n[c][4]; break;
        case 2: n0 = pri_n[c][2]; n1 = pri_n[c][0]; n2 = pri_n[c][3]; n3 = pri_n[c][5]; break;
        case 3: n0 = pri_n[c][0]; n1 = pri_n[c][2]; n2 = pri_n[c][1]; n3 = -1; break;
        case 4: n0 = pri_n[c][3]; n1 = pri_n[c][4]; n2 = pri_n[c][5]; n3 = -1; break;
      }

      ltet=lpyr=lpri=lhex=-1;
      lpri=c;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (n3 >= 0)
      {
        if (lpyr >= 0)
        {
          if ((s = pyramid_side(pyr_n[lpyr],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to pyramid %d",c,i,lpyr);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ntet+lpyr;
        } else if (lpri >= 0)
        {
          if ((s = prism_side(pri_n[lpri],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to prism %d",c,i,lpri);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ntet+npyr+lpri;
        } else if (lhex >= 0)
        {
          if ((s = hex_side(hex_n[lhex],n0,n3,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to hex %d",c,i,lhex);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ntet+npyr+npri+lhex;
        } else
          nquad++;
      } else
      {
        if (ltet >= 0)
        {
          if ((s = tet_side(tet_n[ltet],n0,n2,n1)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to tet %d",c,i,ltet);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ltet;
        } else if (lpyr >= 0)
        {
          if ((s = pyramid_side(pyr_n[lpyr],n0,n2,n1,n3)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to pyramid %d",c,i,lpyr);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ntet+lpyr;
        } else if (lpri >= 0)
        {
          if ((s = prism_side(pri_n[lpri],n0,n2,n1,n3)) < 0)
          {
            fprintf(stderr,"\nIncorrect face winding for prism %d, side %d to prism %d",c,i,lpri);
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          } else
            pri_nbr[c][i] = ntet+npyr+lpri;
        } else
          ntri++;
      }
    }
  }
  for (c=0; c < nhex; c++)
  {
    for (i=0; i < 6; i++)
    {
      hex_nbr[c][i] = -1;
      switch (i)
      {
        case 0: n0 = hex_n[c][0]; n1 = hex_n[c][3]; n2 = hex_n[c][2]; n3 = hex_n[c][1]; break;
        case 1: n0 = hex_n[c][0]; n1 = hex_n[c][1]; n2 = hex_n[c][5]; n3 = hex_n[c][4]; break;
        case 2: n0 = hex_n[c][1]; n1 = hex_n[c][2]; n2 = hex_n[c][6]; n3 = hex_n[c][5]; break;
        case 3: n0 = hex_n[c][2]; n1 = hex_n[c][3]; n2 = hex_n[c][7]; n3 = hex_n[c][6]; break;
        case 4: n0 = hex_n[c][0]; n1 = hex_n[c][4]; n2 = hex_n[c][7]; n3 = hex_n[c][3]; break;
        case 5: n0 = hex_n[c][4]; n1 = hex_n[c][5]; n2 = hex_n[c][6]; n3 = hex_n[c][7]; break;
      }

      ltet=lpyr=lpri=lhex=-1;
      lhex=c;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (lpyr >= 0)
      {
        if ((s = pyramid_side(pyr_n[lpyr],n0,n3,n2,n1)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for hex %d, side %d to pyramid %d",c,i,lpyr);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          hex_nbr[c][i] = ntet+lpyr;
      } else if (lpri >= 0)
      {
        if ((s = prism_side(pri_n[lpri],n0,n3,n2,n1)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for hex %d, side %d to prism %d",c,i,lpri);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          hex_nbr[c][i] = ntet+npyr+lpri;
      } else if (lhex >= 0)
      {
        if ((s = hex_side(hex_n[lhex],n0,n3,n2,n1)) < 0)
        {
          fprintf(stderr,"\nIncorrect face winding for hex %d, side %d to hex %d",c,i,lhex);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
          hex_nbr[c][i] = ntet+npyr+npri+lhex;
      } else
        nquad++;
    }
  }

  glob = ntri;
  glob2 = nquad;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntri,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nquad,&glob2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  if (my_rank == 0)
  {
    fprintf(out_f,"\nAfter element-to-element check, # of unknown neighbors = %d tri-faces, %d quad-faces", glob, glob2);
    fflush(out_f);
  }

  int belems=ntet+npyr+npri+nhex;
  for (b=0; b < nb; b++)
  {
    for (t=0; t < nt[b]; t++)
    {
      n0 = t_n[b][t][0];
      n1 = t_n[b][t][1];
      n2 = t_n[b][t][2];
      n3 = -1;
      ltet=lpyr=lpri=lhex=-1;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (num_procs > 1)
      {
        if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
        {
          double cg[3];
          for (i=0; i < 3; i++)
            cg[i] = (node[n0][i]+node[n1][i]+node[n2][i])/3.0;
          //VCB: debug only to see all triangles with issues          
          //fprintf(debug_f,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR TRIANGLE %d (%d-%d)!",my_rank,b,t,tri_emap[b][t][0],tri_emap[b][t][1]);
          //fprintf(debug_f,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
          //for (i=0; i < 3; i++)
            //fprintf(debug_f,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",t_n[b][t][i],pmap[t_n[b][t][i]][0],pmap[t_n[b][t][i]][1],pmap[t_n[b][t][i]][2],node[t_n[b][t][i]][0],node[t_n[b][t][i]][1],node[t_n[b][t][i]][2]);
          //fflush(debug_f);
          fprintf(stderr,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR TRIANGLE %d (%d-%d)!",my_rank,b,t,tri_emap[b][t][0],tri_emap[b][t][1]);
          fprintf(stderr,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
          for (i=0; i < 3; i++)
            fprintf(stderr,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",t_n[b][t][i],pmap[t_n[b][t][i]][0],pmap[t_n[b][t][i]][1],pmap[t_n[b][t][i]][2],node[t_n[b][t][i]][0],node[t_n[b][t][i]][1],node[t_n[b][t][i]][2]);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } 
        else
        {
          if (ltet >= 0)
          {
            if ((j = tet_side(tet_n[ltet],n0,n2,n1)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to tet %d",b,t,ltet);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,ltet,-1,-1,-1,tet_n[ltet],conn);
            tet_nbr[ltet][j] = belems + t;
            ntri--;
          }
          if (lpyr >= 0)
          {
            if ((j = pyramid_side(pyr_n[lpyr],n0,n2,n1,n3)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to pyramid %d",b,t,lpyr);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
            pyr_nbr[lpyr][j] = belems + t;
            ntri--;
          }
          if (lpri >= 0)
          {
            if ((j = prism_side(pri_n[lpri],n0,n2,n1,n3)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to prism %d",b,t,lpri);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
            pri_nbr[lpri][j] = belems + t;
            ntri--;
          }
        }
      }
      if (num_procs == 1)
      {
        if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
        {
          double cg[3];
          for (i=0; i < 3; i++)
            cg[i] = (node[n0][i]+node[n1][i]+node[n2][i])/3.0;
          fprintf(stderr,"\nNO MATCHING VOLUME ELEMENT FOUND FOR BOUNDARY %d, TRIANGLE %d, (%lg, %lg, %lg)",b,t,cg[0],cg[1],cg[2]);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } 
        else
        {
          if (ltet >= 0)
          {
            if ((j = tet_side(tet_n[ltet],n0,n2,n1)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to tet %d",b,t,ltet);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,ltet,-1,-1,-1,tet_n[ltet],conn);
            tet_nbr[ltet][j] = belems + t;
            ntri--;
          }
          if (lpyr >= 0)
          {
            if ((j = pyramid_side(pyr_n[lpyr],n0,n2,n1,n3)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to pyramid %d",b,t,lpyr);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
            pyr_nbr[lpyr][j] = belems + t;
            ntri--;
          }
          if (lpri >= 0)
          {
            if ((j = prism_side(pri_n[lpri],n0,n2,n1,n3)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, triangle %d to prism %d",b,t,lpri);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
            pri_nbr[lpri][j] = belems + t;
            ntri--;
          }
        }
      }
    }
    belems += nt[b];
    for (q=0; q < nq[b]; q++)
    {
      n0 = q_n[b][q][0];
      n1 = q_n[b][q][1];
      n2 = q_n[b][q][2];
      n3 = q_n[b][q][3];
      ltet=lpyr=lpri=lhex=-1;
      neighbor(n0,n1,n2,n3,ltet,lpyr,lpri,lhex,tethash,pyrhash,prihash,hexhash);
      if (num_procs > 1)
      {
        if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
        {
          double cg[3];
          for (i=0; i < 3; i++)
            cg[i] = (node[n0][i]+node[n1][i]+node[n2][i]+node[n3][i])/4.0;
          //VCB: debug only to see all quads with issues          
          //fprintf(debug_f,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR QUADRILATERAL %d (%d-%d)!",my_rank,b,q,quad_emap[b][q][0],quad_emap[b][q][1]);
          //fprintf(debug_f,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
          //for (i=0; i < 4; i++)
            //fprintf(debug_f,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",q_n[b][q][i],pmap[q_n[b][q][i]][0],pmap[q_n[b][q][i]][1],pmap[q_n[b][q][i]][2],node[q_n[b][q][i]][0],node[q_n[b][q][i]][1],node[q_n[b][q][i]][2]);
          //fflush(debug_f);
          fprintf(stderr,"\nProcess %d, BD %d: NO MATCHING VOLUME ELEMENT FOUND FOR QUADRILATERAL %d (%d-%d)!",my_rank,b,q,quad_emap[b][q][0],quad_emap[b][q][1]);
          fprintf(stderr,"\ncentered at: %16.10e %16.10e %16.10e",cg[0],cg[1],cg[2]);
          for (i=0; i < 4; i++)
            fprintf(stderr,"\nnode %d ( %d - %d - %d): %16.10e %16.10e %16.10e",q_n[b][q][i],pmap[q_n[b][q][i]][0],pmap[q_n[b][q][i]][1],pmap[q_n[b][q][i]][2],node[q_n[b][q][i]][0],node[q_n[b][q][i]][1],node[q_n[b][q][i]][2]);
          fflush(stderr);
          //look to see if two triangles by numbers
          List trifind;
          trifind.Redimension(0);
          trifind.Add_To_List(n0);
          trifind.Add_To_List(n1);
          trifind.Add_To_List(n2);
          trifind.Add_To_List(n3);
          m = 0;
          for (i = 0; i < nt[b] && m < 2; i++)
          {
            k = 0; //cntr
            for (j = 0; j < 3; j++)
            {
              if (trifind.Is_In_List(t_n[b][i][j]))
                k++;
            }
            if (k == 3)
            {
              m++;
            }
          }
          if (m == 2)
          {
            fprintf(stderr,"\nNODALLY FOUND TWO TRIANGLES COVERING QUADRILATERAL %d (%d-%d)!",q,quad_emap[b][q][0],quad_emap[b][q][1]);
            fflush(stderr);
          }
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
        {
          if (lpyr >= 0)
          {
            if ((j = pyramid_side(pyr_n[lpyr],n0,n3,n2,n1)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to pyramid %d",b,q,lpyr);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
            pyr_nbr[lpyr][j] = belems + q;
            nquad--;
          }
          if (lpri >= 0)
          {
            if ((j = prism_side(pri_n[lpri],n0,n3,n2,n1)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to prism %d",b,q,lpri);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
            pri_nbr[lpri][j] = belems + q;
            nquad--;
          }
          if (lhex >= 0)
          {
            if ((j = hex_side(hex_n[lhex],n0,n3,n2,n1)) < 0)
            {
              //if (my_rank == 0)
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to hex %d",b,q,lhex);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,-1,lhex,hex_n[lhex],conn);
            hex_nbr[lhex][j] = belems + q;
            nquad--;
          }
        }
      }
      if (num_procs == 1)
      {
        if (ltet < 0 && lpyr < 0 && lpri < 0 && lhex < 0)
        {
          double cg[3];
          for (i=0; i < 3; i++)
            cg[i] = (node[n0][i]+node[n1][i]+node[n2][i]+node[n3][i])/4.0;
          fprintf(stderr,"\nNO MATCHING VOLUME ELEMENT FOUND FOR BOUNDARY %d, QUADRILATERAL %d, (%lg, %lg, %lg)",b,q,cg[0],cg[1],cg[2]);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        } else
        {
          if (lpyr >= 0)
          {
            if ((j = pyramid_side(pyr_n[lpyr],n0,n3,n2,n1)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to pyramid %d",b,q,lpyr);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,lpyr,-1,-1,pyr_n[lpyr],conn);
            pyr_nbr[lpyr][j] = belems + q;
            nquad--;
          }
          if (lpri >= 0)
          {
            if ((j = prism_side(pri_n[lpri],n0,n3,n2,n1)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to prism %d",b,q,lpri);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,lpri,-1,pri_n[lpri],conn);
            pri_nbr[lpri][j] = belems + q;
            nquad--;
          }
          if (lhex >= 0)
          {
            if ((j = hex_side(hex_n[lhex],n0,n3,n2,n1)) < 0)
            {
              fprintf(stderr,"\nIncorrect face winding for boundary %d, quadrilateral %d to hex %d",b,q,lhex);
              fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            //j=retrieve_face_nodes(n0,n1,n2,n3,-1,-1,-1,lhex,hex_n[lhex],conn);
            hex_nbr[lhex][j] = belems + q;
            nquad--;
          }
        }
      }
    }
    belems += nq[b];
  }

  glob = ntri;
  glob2 = nquad;
  if (num_procs > 1)
  {
    MPI_Allreduce(&ntri,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nquad,&glob2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  if (my_rank == 0)
  {
    fprintf(out_f,"\nAfter boundary-to-element check, # of unknown neighbors = %d tri-faces, %d quad-faces", glob, glob2);
    fprintf(out_f,"\nThis may not be zero in parallel, due to ghost elements, so check for ***** comments.\n");
    fflush(out_f);
  }
  
  //rest is no longer valid based on heredity
  if (timing == -1)
  {
    delete[] tet_nbr;
    delete[] pyr_nbr;
    delete[] pri_nbr;
    delete[] hex_nbr;

    return(1);
  }
  
  if (my_rank == 0) 
  {
    fprintf(out_f,"\nChecking neighbor connectivity.");
    fflush(out_f);
  }

  double x, y, z;
  t = q = 0;
  for (c=0; c < ntet; c++)
  {
    n0 = tet_n[c][0];
    n1 = tet_n[c][1];
    n2 = tet_n[c][2];
    n3 = tet_n[c][3];
    x = (node[n0][0]+node[n1][0]+node[n2][0]+node[n3][0])*0.25;
    y = (node[n0][1]+node[n1][1]+node[n2][1]+node[n3][1])*0.25;
    z = (node[n0][2]+node[n1][2]+node[n2][2]+node[n3][2])*0.25;
    if (num_procs > 1)
    {
      for (i=0; i < 4; i++)
        if (tet_nbr[c][i] < 0 && tet_emap[c][0] == my_rank && (pmap[n0][1] == my_rank && pmap[n1][1] == my_rank && pmap[n2][1] == my_rank && pmap[n3][1] == my_rank))
        {
          t++;
          #ifdef _DEBUG
          fprintf(debug_f,"\nTet %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(debug_f);
          #endif
        }
    }
    if (num_procs == 1)
    {
      for (i=0; i < 4; i++)
        if (tet_nbr[c][i] < 0)
        {
          t++;
          fprintf(out_f,"\nTet %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(out_f);
        }
    }
  }
    
  glob = t;
  if (num_procs > 1)
    MPI_Allreduce(&t,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET TET NEIGHBORS = %d",glob);
    fflush(out_f);
  }  
    
  t = q = 0;
  for (c=0; c < npyr; c++)
  {
    n0 = pyr_n[c][0];
    n1 = pyr_n[c][1];
    n2 = pyr_n[c][2];
    n3 = pyr_n[c][3];
    n4 = pyr_n[c][4];
    x = (node[n0][0]+node[n1][0]+node[n2][0]+node[n3][0]+node[n4][0])/5.0;
    y = (node[n0][1]+node[n1][1]+node[n2][1]+node[n3][1]+node[n4][1])/5.0;
    z = (node[n0][2]+node[n1][2]+node[n2][2]+node[n3][2]+node[n4][2])/5.0;
    if (num_procs > 1)
    {
      for (i=0; i < 5; i++)
      {
        if (pyr_nbr[c][i] < 0 && pyr_emap[c][0] == my_rank && (pmap[n0][1] == my_rank && pmap[n1][1] == my_rank && pmap[n2][1] == my_rank && pmap[n3][1] == my_rank && pmap[n4][1] == my_rank))
        {
          switch(i)
          {
            case (0):
              q++;
              break;
            case (1):
            case (2):
            case (3):
            case (4):
              t++;
              break;
          }
          #ifdef _DEBUG
          fprintf(debug_f,"\nPyramid %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(debug_f);
          #endif
        }
      }
    }
    if (num_procs == 1)
    {
      for (i=0; i < 5; i++)
      {
        if (pyr_nbr[c][i] < 0)
        {
          switch(i)
          {
            case (0):
              q++;
              break;
            case (1):
            case (2):
            case (3):
            case (4):
              t++;
              break;
          }
          fprintf(out_f,"\nPyramid %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(out_f);
        }
      }
    }
  }
  
  glob = t;
  if (num_procs > 1)
    MPI_Allreduce(&t,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET PYRAMID Tri-face NEIGHBORS = %d",glob);
    fflush(out_f);
  }
  glob = q;
  if (num_procs > 1)
    MPI_Allreduce(&q,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET PYRAMID Quad-face NEIGHBORS = %d",glob);
    fflush(out_f);
  }

  t = q = 0;
  for (c=0; c < npri; c++)
  {
    n0 = pri_n[c][0];
    n1 = pri_n[c][1];
    n2 = pri_n[c][2];
    n3 = pri_n[c][3];
    n4 = pri_n[c][4];
    n5 = pri_n[c][5];
    x = (node[n0][0]+node[n1][0]+node[n2][0]+node[n3][0]+node[n4][0]+node[n5][0])/6.0;
    y = (node[n0][1]+node[n1][1]+node[n2][1]+node[n3][1]+node[n4][1]+node[n5][1])/6.0;
    z = (node[n0][2]+node[n1][2]+node[n2][2]+node[n3][2]+node[n4][2]+node[n5][2])/6.0;
    if (num_procs > 1)
    {
      for (i=0; i < 5; i++)
      {
        if (pri_nbr[c][i] < 0 && pri_emap[c][0] == my_rank && (pmap[n0][1] == my_rank && pmap[n1][1] == my_rank && pmap[n2][1] == my_rank && pmap[n3][1] == my_rank && pmap[n4][1] == my_rank && pmap[n5][1] == my_rank))
        {
          switch(i)
          {
            case (0):
            case (1):
            case (2):
              q++;
              break;
            case (3):
            case (4):
              t++;
              break;
          }
          #ifdef _DEBUG
          fprintf(debug_f,"\nPrism %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(debug_f);
          #endif
        }
      }
    }
    if (num_procs == 1)
    {
      for (i=0; i < 5; i++)
      {
        if (pri_nbr[c][i] < 0)
        {
          switch(i)
          {
            case (0):
            case (1):
            case (2):
              q++;
              break;
            case (3):
            case (4):
              t++;
              break;
          }
          fprintf(out_f,"\nPrism %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(out_f);
        }
      }
    }
  }
  
  glob = t;
  if (num_procs > 1)
    MPI_Allreduce(&t,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET PRISM Tri-face NEIGHBORS = %d",glob);
    fflush(out_f);
  }
  glob = q;
  if (num_procs > 1)
    MPI_Allreduce(&q,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET PRISM Quad-face NEIGHBORS = %d",glob);
    fflush(out_f);
  }

  t = q = 0;
  for (c=0; c < nhex; c++)
  {
    n0 = hex_n[c][0];
    n1 = hex_n[c][1];
    n2 = hex_n[c][2];
    n3 = hex_n[c][3];
    n4 = hex_n[c][4];
    n5 = hex_n[c][5];
    n6 = hex_n[c][6];
    n7 = hex_n[c][7];
    x = (node[n0][0]+node[n1][0]+node[n2][0]+node[n3][0]+node[n4][0]+node[n5][0]+node[n6][0]+node[n7][0])/8.0;
    y = (node[n0][1]+node[n1][1]+node[n2][1]+node[n3][1]+node[n4][1]+node[n5][1]+node[n6][1]+node[n7][1])/8.0;
    z = (node[n0][2]+node[n1][2]+node[n2][2]+node[n3][2]+node[n4][2]+node[n5][2]+node[n6][2]+node[n7][2])/8.0;
    if (num_procs > 1)
    {
      for (i=0; i < 6; i++)
      {
        if (hex_nbr[c][i] < 0 && hex_emap[c][0] == my_rank && (pmap[n0][1] == my_rank && pmap[n1][1] == my_rank && pmap[n2][1] == my_rank && pmap[n3][1] == my_rank && pmap[n4][1] == my_rank && pmap[n5][1] == my_rank && pmap[n6][1] == my_rank && pmap[n7][1] == my_rank))
        {
          q++;
          #ifdef _DEBUG
          fprintf(debug_f,"\nHex %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(debug_f);
          #endif
        }
      }
    }
    if (num_procs == 1)
    {
      for (i=0; i < 6; i++)
      {
        if (hex_nbr[c][i] < 0)
        {
          q++;
          fprintf(out_f,"\nHex %d, side %d, no neighbor! ( %lf, %lf, %lf)",c,i,x,y,z);
          fflush(out_f);
        }
      }
    }
  }
  glob = q;
  if (num_procs > 1)
  MPI_Allreduce(&q,&glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (glob > 0 && my_rank == 0)
  {
    fprintf(out_f,"\n***** NUMBER OF UNSET HEX NEIGHBORS = %d",glob);
    fflush(out_f);
  }
    
  fflush(out_f); //just to be sure

  delete[] tet_nbr;
  delete[] pyr_nbr;
  delete[] pri_nbr;
  delete[] hex_nbr;

  return(1);
}

