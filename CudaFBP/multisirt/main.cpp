

/**********************************************************************************main.c*/

#include "mrcfile_atom.h"
  

int  ATOM(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag) ;//Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
int  ATOM_NCC(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag);
int  ATOM_SIMPLE(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag);
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

   float ATOM_ITR_STEP=0.3;
   int ATOM_ITR_NUM=10;

   int initflag=1;
   int i;
   for(i=6;i<argc;i=i+2)
   {
      if(argv[i][0]!='-')
      {
         printf("\"%s\" error! No such option!\n",argv[i]);
         return 0;
      }

      switch(argv[i][1])
      {
         case 'n':
           ATOM_ITR_NUM=atoi(argv[i+1]);
           break;
         case 't':
           ATOM_ITR_STEP=atof(argv[i+1]);
           break;
         case 'i':
           initflag=atof(argv[i+1]);
           break;
      }
   }


   char *inf=argv[1];
   char *outf=argv[2];
   char *fbpf=argv[3];
   char *coef=argv[4];

   if(!id)
   {
	   for(int i=1;i<argc;i++)
		   printf("%s ",argv[i]);
	   printf("\n");
   }

   if(!strcmp(argv[5],"SIRT"))
	   ATOM(inf, outf,  fbpf, coef, ATOM_ITR_STEP, ATOM_ITR_NUM,"SIRT",id,p,initflag);
   if(!strcmp(argv[5],"SIRTNCC"))
	   ATOM_NCC(inf, outf,  fbpf, coef, ATOM_ITR_STEP, ATOM_ITR_NUM,"SIRTNCC",id,p,initflag);
   if(!strcmp(argv[5],"SIRTSIMPLE"))
	   ATOM_SIMPLE(inf, outf,  fbpf, coef, ATOM_ITR_STEP, ATOM_ITR_NUM,"SIRTSIMPLE",id,p,initflag);
   
   
    MPI_Finalize();                 //parallel finish
  
      return 0;
  
}

//eel-tomo.st doublerecon_SIRT1-0.2.mrc doublereconbpt.mrc cut_eel-tomo4.txbr SIRT -n 10 -t 0.2
//single
//eel-tomo4a.st singlerecon_SIRT1-0.2.mrc  singlebpt.mrc eel-tomo4a.txbr SIRT -n 10 -t 0.2



