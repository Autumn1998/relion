

/**********************************************************************************main.c*/

#include "mrcfile_atom.h"
  

int  ATOM(char *inf, char *outf, char *coef, int myid,int mypro);
/**********************************************************************************/


int main(int argc, char *argv[])
{
     //double elapsed_time; /*parallel execution time*/
   int id;       //process id number
   int p;       //number of process
   //double start_time,end_time;
   int namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   
   MPI_Init(&argc,&argv); //parallel init
   MPI_Comm_rank(MPI_COMM_WORLD,&id);
   MPI_Comm_size(MPI_COMM_WORLD,&p);
   MPI_Get_processor_name(processor_name,&namelen); 

   char *inf=argv[1];
   char *outf=argv[2];
//   char *filter_prj=argv[3];
   char *coef=argv[3];
   ATOM(inf, outf,coef, id,p);
   


   
    MPI_Finalize();                 //parallel finish
  
      return 0;
  
}




