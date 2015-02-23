#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "P_Adapt.h"
#include "journal.h"
#include "Spacing_Field.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

//global variables
int my_rank; /*rank of process*/
int num_procs; /*number of processes*/
FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

int main(int argcs, char* pArgs[])
{
  int source; //rank of sender
  int tag = 0; //tag for messages
  int i, k, nf, digits, target, tri_flag, pri_stack, mltype, bnd_flag;
  double threshold;
  const int bdim = 132; //buffer dim
  char extension[bdim]; //file extension to allow padded digits
  char adapt_file[bdim];
  char buff[bdim]; //buffer
  char sname[bdim]; //mesh name storage
  time_t tm;
  clock_t time_0, time_1;
  char *t_char;
 
  my_rank = 0;
  num_procs = 1;

  int dest; //rank of receiver
  MPI_Status status; //return status for receive
  //Start up MPI
  MPI_Init(&argcs, &pArgs);
   
  //Find out process rank of current instance
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
  //Find out number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  in_f = stdin;
  out_f = stdout;

  if (my_rank == 0)
  {
    //create journal file
    if ((jou_f=fopen("P_REFINE.jou","w")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file journal file");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  
    //check for standard input
    if (--argcs < 1)
    {
      fprintf(stdout,"\nNo input file specified!");
      fprintf(stdout,"\nUsing standard input!");
      fflush(stdout);
    } else
    {
      if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file <%s>\n",pArgs[argcs]);
        fflush(stderr);
        fprintf(jou_f,"\nCouldn't open file <%s>\n",pArgs[argcs]);
        fflush(jou_f);
        MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }

  if (my_rank == 0 && in_f != stdin)
  {
    sprintf(buff,"P_REFINE.out");

    if ((out_f=fopen(buff,"w")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file output file %s",buff);
      fflush(stderr);
      fclose(jou_f);
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }

  #ifdef _DEBUG
  buff[0] = '\0';
  
  sprintf(buff,"debug.%d",my_rank);

  if ((debug_f=fopen(buff,"w")) == NULL)
  {
    fprintf(stderr,"\nCouldn't open debug output file %s",buff);
    fflush(stderr);
    fclose(jou_f);
    fclose(out_f);
    MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
  #endif

  MPI_Barrier(MPI_COMM_WORLD);

  // output the beginning wall time
  time(&tm);
  t_char = ctime(&tm);
  time_0 = clock();
  if (my_rank == 0)
  {
    fprintf(out_f,"\nP_REFINE run started at %s",t_char);

    //print program info 
    fprintf(out_f,"\n======================================================================");
    fprintf(out_f,"\n| COPYRIGHT 2003-2012 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA     |");
    fprintf(out_f,"\n|                                                                    |");
    fprintf(out_f,"\n|                    RIGHTS IN DATA                                  |");
    fprintf(out_f,"\n|                                                                    |");
    fprintf(out_f,"\n| THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT |");
    fprintf(out_f,"\n|    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO       |");
    fprintf(out_f,"\n|         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL      |");
    fprintf(out_f,"\n|              RIGHTS IN DATA CONTRACT CLAUSES.                      |");
    fprintf(out_f,"\n|       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY   |");
    fprintf(out_f,"\n|              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA            |");
    fprintf(out_f,"\n|                                                                    |");
    fprintf(out_f,"\n| Parallel Hybrid Mesh Refinement Code (P_REFINE)                    |");
    fprintf(out_f,"\n| NOTE: This data includes the SimCenter at Chattanooga P_REFINE     |");
    fprintf(out_f,"\n| code, which was developed under private non-government funding.    |");
    fprintf(out_f,"\n| This software is submitted with limited rights to use, reproduce,  |");
    fprintf(out_f,"\n| and disclose this data for Government Purposes only.               |");
    fprintf(out_f,"\n| Requests for access to the software for non-governmental purposes  |");
    fprintf(out_f,"\n| should be referred to                                              |"); 
    fprintf(out_f,"\n|                                                                    |");
    fprintf(out_f,"\n|    Dr. Steve Karman                                                |"); 
    fprintf(out_f,"\n|    Steve-Karman@utc.edu                                            |"); 
    fprintf(out_f,"\n|    423-425-5492  or  423-425-5470                                  |"); 
    fprintf(out_f,"\n|                                                                    |");
    fprintf(out_f,"\n|    SimCenter: National Center for Computational Engineering        |"); 
    fprintf(out_f,"\n|    701 East M. L. King Boulevard                                   |"); 
    fprintf(out_f,"\n|    Chattanooga, TN 37403                                           |"); 
    fprintf(out_f,"\n======================================================================\n");
    fflush(out_f);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  adapt_file[0] = '\0';

  if (my_rank == 0)
  {
    journal(in_f, out_f, jou_f,"#Adaption file name ->",adapt_file);
    journal(in_f, out_f, jou_f,"#Enter physical grid file name (serial version) ->",sname);

    if (num_procs > 1)
    {
      // inform slave processes of adaptation file name
      tag = my_rank;
      for (dest = 1; dest < num_procs; dest++)
        MPI_Send(adapt_file, strlen(adapt_file)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);

      // inform slave processes of mesh file name
      tag = my_rank;
      for (dest = 1; dest < num_procs; dest++)
        MPI_Send(sname, strlen(sname)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);	
    }
  } else
  {
    source = 0;
    tag = 0;
    MPI_Recv(adapt_file, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);

    MPI_Recv(sname, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    journal(in_f, out_f, jou_f,"#Enter target number of new points ->",target);
    journal(in_f, out_f, jou_f,"#Prohibit 2-node triangle marking [0 1] ->",tri_flag);
    journal(in_f, out_f, jou_f,"#Force prism stack refinement [0 1] ->",pri_stack);
    journal(in_f, out_f, jou_f,"#Enter metric length refinement threshold [ >= 1.0] ->",threshold);
    journal(in_f, out_f, jou_f,"#Enter edge refinement type [ 0-Max ML, 1-Weighted Avg ML] ->",mltype);
    journal(in_f, out_f, jou_f,"#Refine boundary edges only [ 0 1] ->",bnd_flag);
  
    if (in_f != stdin) fclose(in_f);

    fclose(jou_f);

    fprintf(out_f,"\nTarget number of new nodes = %i",target);
    fprintf(out_f,"\nProhibit 2-node triangle marking = %i",tri_flag);
    fprintf(out_f,"\nForce prism stack refinement = %i",pri_stack);
    fprintf(out_f,"\nMetric length refinement threshold = %lf",threshold);
    fprintf(out_f,"\nEdge refinement type [ 0-Max ML, 1-Weighted Avg ML] = %i",mltype);
    fprintf(out_f,"\nRefine boundary edges only [ 0 1] = %i",bnd_flag);
    fflush(out_f);
  }

  if (num_procs > 1)
  {
    MPI_Bcast(&target, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tri_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pri_stack, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&threshold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mltype, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bnd_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  adapt_mesh(adapt_file, sname, target, tri_flag, pri_stack, threshold, mltype, bnd_flag);

  if (my_rank == 0)
  {
    fprintf(out_f,"\nAdapted mesh.\n");
    fflush(out_f);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //perform end time check
  time(&tm);
  t_char = ctime(&tm);
  time_1 = clock();
  if (my_rank == 0)
  {
    fprintf(out_f,"\nTotal adaptation time in seconds = %14.7e\n",(float)(time_1-time_0)/CLOCKS_PER_SEC);
    fflush(out_f);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0)
  {
    fprintf(out_f,"\nP_REFINE run completed at %s",t_char);
    fflush(out_f);
  }

  if (out_f != stdout) fclose(out_f);

  #ifdef _DEBUG
  if (in_f != stdin)
    fclose(debug_f);
  #endif
  
  MPI_Finalize(); 

  return(0);
}
