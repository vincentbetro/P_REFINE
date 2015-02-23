#include <stdio.h>
#include <stdlib.h>

void refine_tri(int *conn, int cdim, int &nt, int tri[][3], int &nq, int quad[][4], long int &refinetype);
void refine_quad(int *conn, int cdim, int &nt, int tri[][3], int &nq, int quad[][4], long int &refinetype);
void quad_fill(int cdim, int &ntet, int tet[][4],
                         int &npyr, int pyr[][5],
                         int &npri, int pri[][6],
                         int &nhex, int hex[][8], int mid, int fn[9]);
int refine_tet(int conn[10], int cdim, int &ntet, int tet[][4], int &npyr, int pyr[][5],
                                       int &npri, int pri[][6], int &nhex, int hex[][8], long int &refinetype);
int refine_pyr(int conn[15], int cdim, int &ntet, int tet[][4], int &npyr, int pyr[][5],
                                       int &npri, int pri[][6], int &nhex, int hex[][8], long int &refinetype);
int refine_pri(int conn[19], int cdim, int &ntet, int tet[][4], int &npyr, int pyr[][5],
                                       int &npri, int pri[][6], int &nhex, int hex[][8], long int &refinetype);
int refine_hex(int conn[27], int cdim, int &ntet, int tet[][4], int &npyr, int pyr[][5],
                                       int &npri, int pri[][6], int &nhex, int hex[][8], long int &refinetype);

