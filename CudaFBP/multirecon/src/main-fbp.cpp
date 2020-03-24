

/**********************************************************************************main.c*/

#include "mrcfile_atom.h"
  

int  ATOM(char *inf, char *outf, char *coef, int myid,int mypro);
//int  ATOM(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro);
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
   
   // eel-tomo.st  doublerecon.mrc eel-tomo4.txbr

   
    MPI_Finalize();                 //parallel finish
  
      return 0;
  
}
// eel-tomo4a.st  singlebpt.mrc eel-tomo4a.txbr
// eel-tomo.st  doublebpt.mrc cut_eel-tomo4.txbr


