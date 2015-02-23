int CGNS_read(char *name, int &nn, double ***p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn);

int CGNS_write(char *name, int nn, double **p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn);

int P_CGNS_read(int parallel, char *fname, int &nn, double ***p, int &nb, char ***b_name,
              int **nt, int ****tri_conn,
              int **nq, int ****quad_conn,
              int &ntet, int ***tet_conn,
              int &npyr, int ***pyr_conn,
              int &npri, int ***pri_conn,
              int &nhex, int ***hex_conn,
              int ***nmap, int ***tri_map, int ***quad_map,
              int **tet_map, int **pyr_map, int **pri_map, int **hex_map);

int P_CGNS_write(int parallel, char *fname, int nn, double **p, int nb, char **b_name,
               int *nt, int ***tri_conn,
               int *nq, int ***quad_conn,
               int ntet, int **tet_conn,
               int npyr, int **pyr_conn,
               int npri, int **pri_conn,
               int nhex, int **hex_conn,
               int **nmap, int **tri_map, int **quad_map,
               int *tet_map, int *pyr_map, int *pri_map, int *hex_map);
