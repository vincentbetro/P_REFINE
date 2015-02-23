#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "P_Adapt.h"
#include "mesh.h"
#include "CGNS.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

void adapt_mesh(char aname[], char sname[], int target, int tri_flag, int pri_stack, double threshold, int mltype, int bnd_flag)
{
  const int bdim = 80;
  char buff[bdim];
  mesh_obj *mesh = new mesh_obj();

  if (num_procs > 1)
  {
    sprintf(buff,"%s",sname);
    char *ptr = strstr(buff,".cgns");
    if (ptr == NULL)
    {
      fprintf(out_f,"\nCGNS suffix <.cgns> not found in file name!");
      fflush(out_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    } else
    {
      //reset cursor to overwrite extension with new extension and proc num
      *ptr = '\0';
    }
    sprintf(sname,"%s_%d.cgns",buff,my_rank);
  }
  
  #ifdef _DEBUG
  fprintf(debug_f,"\nAbout to read in mesh.\n");
  fflush(debug_f);
  #endif

  // read mesh data
  mesh->refine_io(-1,sname);
  
  #ifdef _DEBUG
  fprintf(debug_f,"\nDone reading in mesh.\n");
  fflush(debug_f);
  #endif
  
  if (num_procs > 1)
  {
    #ifdef _DEBUG
    fprintf(debug_f,"\nAbout to create emaps.\n");
    fflush(debug_f);
    #endif
    //create hash, determine element p_i
    mesh->element_pi();
    #ifdef _DEBUG
    fprintf(debug_f,"\nDone creating emaps.\n");
    fflush(debug_f);
    #endif
  }

  #ifdef _DEBUG
  fprintf(debug_f,"\nBeginning refinement.\n");
  fflush(debug_f);
  #endif  
  // call refinement function
  mesh->refine(target,tri_flag,pri_stack,threshold,aname,mltype,bnd_flag);
  #ifdef _DEBUG
  fprintf(debug_f,"\nFinished refinement.\n");
  fflush(debug_f);
  #endif 

  #ifdef _DEBUG
  fprintf(debug_f,"\nWriting out mesh.\n");
  fflush(debug_f);
  #endif 
  // overwrite mesh data
  mesh->refine_io(1,sname);
  #ifdef _DEBUG
  fprintf(debug_f,"\nFinished writing out mesh.\n");
  fflush(debug_f);
  #endif

  delete mesh;

}
