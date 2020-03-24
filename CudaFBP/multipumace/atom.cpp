//gcc -o atom atom.c mrcfile.c -lm
#include <time.h>
#include "mrcfile_atom.h"

#include "atom.h"
#ifndef TEXT_LINE_MAX
#define TEXT_LINE_MAX 500
#endif

#define MILLION 1000000

//float ANG[ANG_MAX];
//float BNG[ANG_MAX];
/**********************************************************************************/

/**********************************************************************************/
int read_coef(double *x_coef, double *y_coef, FILE *f_coef)
{
  char *lstr;
  lstr=(char *)malloc(TEXT_LINE_MAX); 
  char *tmp;
  int i=0;
  int j=0;
  int ang_num=0;

  while((fgets(lstr,TEXT_LINE_MAX,f_coef)) != NULL)
     {
       if(lstr[0]=='l') ang_num++;
       else if(lstr[0]=='x')
            { 
             /*tmp=strtok(lstr," ");
               
             for(tmp=strtok(lstr," ");tmp;tmp=strtok(NULL," "))
                  x_coef[i++] = strtod(tmp,NULL);*/
              tmp=strtok(lstr," ");
                 //  printf("%s\n",tmp);
              while(tmp!=NULL)
                {
                tmp=strtok(NULL," ");
                if(tmp!=NULL) 
                  {               
                   x_coef[i++]=strtod(tmp,NULL);              
                 //  printf("%lf\n",x);
                  }
                }
            }
            else if(lstr[0]=='y')
              { 
             /*tmp=strtok(lstr," ");
             for(tmp=strtok(lstr," ");tmp;tmp=strtok(NULL," "))
                  y_coef[j++] = strtod(tmp,NULL);*/
              tmp=strtok(lstr," ");
                 //  printf("%s\n",tmp);
              while(tmp!=NULL)
                 {
                tmp=strtok(NULL," ");
                if(tmp!=NULL) 
                   {               
                   y_coef[j++]=strtod(tmp,NULL);              
                 //  printf("%lf\n",x);
                   }
                 }
             } 
      }//end while

 //  printf("hello world/n");
  return TRUE;
}
  

/**********************************************************************************/
int val_coef(Pixel pixel, int angle,double *x_coef, double *y_coef, Weight *comp_prj)
{
  double x,y;
 
  int index;
 //order=2;
  /*
  index=10*angle;
  x=x_coef[index]+x_coef[index+1]*pixel.X+x_coef[index+2]*pixel.Y+x_coef[index+3]*pixel.Z+x_coef[index+4]*pixel.X*pixel.X+x_coef[index+5]*pixel.X*pixel.Y+x_coef[index+6]*pixel.X*pixel.Z+x_coef[index+7]*pixel.Y*pixel.Y+x_coef[index+8]*pixel.Y*pixel.Z+x_coef[index+9]*pixel.Z*pixel.Z;
  y=y_coef[index]+y_coef[index+1]*pixel.X+y_coef[index+2]*pixel.Y+y_coef[index+3]*pixel.Z+y_coef[index+4]*pixel.X*pixel.X+y_coef[index+5]*pixel.X*pixel.Y+y_coef[index+6]*pixel.X*pixel.Z+y_coef[index+7]*pixel.Y*pixel.Y+y_coef[index+8]*pixel.Y*pixel.Z+y_coef[index+9]*pixel.Z*pixel.Z;
*/

   
//  order=1
   index=4*angle;
   x=x_coef[index]+x_coef[index+1]*pixel.X+x_coef[index+2]*pixel.Y+x_coef[index+3]*pixel.Z;
   y=y_coef[index]+y_coef[index+1]*pixel.X+y_coef[index+2]*pixel.Y+y_coef[index+3]*pixel.Z;
   
   x--;
   y--;

   comp_prj->x_min=floor(x);
   comp_prj->y_min=floor(y);
   
   comp_prj->x_min_del=x-comp_prj->x_min;      
   comp_prj->y_min_del=y-comp_prj->y_min;

   
   
  return TRUE;

}


/**********************************************************************************/  
double compX_ini(int index, double *x_coef,int X,int Y,int Z)
{
     return x_coef[index]+x_coef[index+1]*X+x_coef[index+2]*Y+x_coef[index+3]*Z+x_coef[index+4]*X*X+x_coef[index+5]*X*Y+x_coef[index+6]*X*Z+x_coef[index+7]*Y*Y+x_coef[index+8]*Y*Z+x_coef[index+9]*Z*Z;
}
 

/**********************************************************************************/

 
double compY_ini(int index, double *y_coef,int X,int Y,int Z)
{
        return y_coef[index]+y_coef[index+1]*X+y_coef[index+2]*Y+y_coef[index+3]*Z+y_coef[index+4]*X*X+y_coef[index+5]*X*Y+y_coef[index+6]*X*Z+y_coef[index+7]*Y*Y+y_coef[index+8]*Y*Z+y_coef[index+9]*Z*Z;
}

/**********************************************************************************/

double q_val(int i, int m, int n, Pixel pixel, int index, double *x_coef,double *x_val)
{
  if (0 <= i <= n && m==0)
      return x_val[n];
    
  else if (i==0 && 0<m<=n)
         return (q_val(i+1,m-1,n,pixel,index,x_coef,x_val)-q_val(i,m,n,pixel,index,x_coef,x_val));
  

        else {
            if (m==n) return q_val(0,n,n,pixel,index,x_coef,x_val);
            else      return (q_val(i-1,m,n,pixel,index,x_coef,x_val)+q_val(i-1,m+1,n,pixel,index,x_coef,x_val));
             }

}








/**********************************************************************************/
int SIRT_Slice_reproj_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_calc_num_all,float *slc_data, int Z_start, int Z_end)
{
   
   int n=1;//order
   int j;


  // printf("hello\n");
   int pixel_index,prj_index;
           

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
           pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
           
            for(j=0;j<prj.AngN;j++)
            //for(j=0;j<1;j++)

            {
              val_coef(pixel, j, x_coef, y_coef, comp_prj); 
            // val_coef_recur(n, pixel, vol, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);
            /* if (pixel.X==-62 && pixel.Y==82 && pixel.Z==-12)
                printf("comp_prj->x is %d, comp_prj->y is %d in angle=%d\n", comp_prj->x_min, comp_prj->y_min,j);*/
            


              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                 {
                
                   prj_index=comp_prj->x_min+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj index
                   prj_calc_all[prj_index]+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*slc_data[pixel_index];
                   prj_calc_num_all[prj_index] += (1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);
    
                 }//end if
              if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min+1, y_min)
                 {
                
                   prj_index=comp_prj->x_min+1+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj index
                   prj_calc_all[prj_index]+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*slc_data[pixel_index];
                   prj_calc_num_all[prj_index] += comp_prj->x_min_del*(1-comp_prj->y_min_del);
    
                 }//end if
              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                 {
                    prj_index=comp_prj->x_min+(comp_prj->y_min+1)*prj.X+j*prj.X*prj.Y; //prj index
                    prj_calc_all[prj_index]+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*slc_data[pixel_index];
                    prj_calc_num_all[prj_index] += (1-comp_prj->x_min_del)*comp_prj->y_min_del;           
                 }
              if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min+1, y_min+1)
                 {
                    prj_index=(comp_prj->x_min+1)+(comp_prj->y_min+1)*prj.X+j*prj.X*prj.Y; //prj index
                    prj_calc_all[prj_index]+=comp_prj->x_min_del*comp_prj->y_min_del*slc_data[pixel_index];
                    prj_calc_num_all[prj_index] += comp_prj->x_min_del*comp_prj->y_min_del;   
                 }

            }//end for prj.AngN                
         }//end for pixel

     //printf("slice %d is finished to reproj.\n",pixel.Z);
     }//end for pixel.Z
   
   return TRUE;

}



/**********************************************************************************/
int SIRT_one_update_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *prj_real,float *slc_data, float ASART_ITR_STEP, int Z_start, int Z_end)
{

   int n=1;//order
   int j;

   int pixel_index,prj_cal_index;
   double s,c;


   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
              pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
              s=0;
              c=0;
              //for(j=0;j<1;j++)
             for(j=0;j<prj.AngN;j++)
              {
               val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //  val_coef_recur(n, pixel, vol,j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);

            /*   if (pixel.X==-62 && pixel.Y==82 && pixel.Z==-12)
                  printf("comp_prj->x is %d, comp_prj->y is %d in angle=%d\n", comp_prj->x_min, comp_prj->y_min,j);*/

            
              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                {
                   prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);                    
                 } //end if
                
              if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)//(x_min+1, y_min)
                {
                   prj_cal_index=comp_prj->x_min+1+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index]; 
                   c+=comp_prj->x_min_del*(1-comp_prj->y_min_del);                    
                 } //end if
              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+(1+comp_prj->y_min)*prj.X+j*prj.X*prj.Y; //prj_calc index
               //    prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*comp_prj->y_min_del;                    
                 } //end if
               if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y)   //(x_min+1, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+1+(comp_prj->y_min+1)*prj.X+j*prj.X*prj.Y; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   c+=comp_prj->x_min_del*comp_prj->y_min_del;                    
                 } //end if
              
              }//end for prj.AngN

           if(c!=0)
                slc_data[pixel_index] += (float)(s/c)*ASART_ITR_STEP;

                
         }// end pixel
       //printf("slice %d is completed to update.\n",pixel.Z);
                
     }// end pixel.Z

  // printf("block in Z_start %d is completed\n",Z_start);

   return TRUE;

}

/**********************************************************************************/
int SIRT_update_slice(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_calc_num_all, float *prj_real,float *slc_data, float *global_prj_calc_all, float *global_prj_calc_num_all, float ATOM_ITR_STEP,int Z_start, int Z_end)
{

   int n;
   int line_allnum=prj.X*prj.Y*prj.AngN;

   int slcN_per=0;
   int global_slcN=0;


 
   SIRT_Slice_reproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all, slc_data, Z_start, Z_end);



    MPI_Allreduce(prj_calc_all,global_prj_calc_all,line_allnum,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(prj_calc_num_all,global_prj_calc_num_all,line_allnum,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

    for(n=0;n<line_allnum;n++)
    {
       if(global_prj_calc_num_all[n] != 0)
    	   global_prj_calc_all[n]/=global_prj_calc_num_all[n];
       global_prj_calc_all[n]=prj_real[n]-global_prj_calc_all[n];  //prj error
    }

     SIRT_one_update_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc_all, prj_real, slc_data, ATOM_ITR_STEP,Z_start,Z_end);

   MPI_Barrier(MPI_COMM_WORLD);


  

   return TRUE;

}






/**********************************************************************************/
int  ATOM_MPI_SIRT(char *outf,char *fbpf, float ATOM_ITR_STEP, int ATOM_ITR_NUM,int myid,int mypro,int flagpart)
//Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
{
  
   Volume vol;
   Projection prj;
   char *inf;
   char *coef;
/************************File and file head processing  ****************************************/
   MPI_File fout, f_fbp;

   if(flagpart==0)
   {
	   inf="eel-tomo4a.st";
   }
   else
   {
	   inf="eel-tomo4b.st";
   }
   if(flagpart==0)
   {
	   coef="eel-tomo4a.txbr";
   }
   else
   {
	   coef="eel-tomo4b.txbr";
   }

   MPI_File_open(MPI_COMM_WORLD, outf, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fout);
   MPI_File_open(MPI_COMM_WORLD, fbpf, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &f_fbp);
   MrcHeader  *inhead,  *outhead, *fbphead;
   inhead=(MrcHeader *)malloc(sizeof(MrcHeader));
   outhead=(MrcHeader *)malloc(sizeof(MrcHeader));
   fbphead=(MrcHeader *)malloc(sizeof(MrcHeader));
   fbphead=(MrcHeader *)malloc(sizeof(MrcHeader));
   if (myid==0)
   {
     FILE *fin;
     fin=fopen(inf,"r");
     if(!fin)
      {
       printf("Cannot open file inf strike any key exit!");
       return 1;
       }
     mrc_read_head(fin,inhead);
     fclose(fin);

     mrc_read_head_MPI(f_fbp,fbphead);

    } 
    MPI_Bcast(inhead, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(fbphead, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);

//   printf("head->cmap[0] is %c\n",inhead->cmap[0]);

   prj.X=inhead->nx;
   prj.Y=inhead->ny;
   prj.AngN=inhead->nz;

   vol.X=fbphead->nx;
   vol.Y=fbphead->ny;
   vol.Z=fbphead->nz;

   vol.Xstart=-37;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=-61;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-33;
   vol.Zend=vol.Zstart+vol.Z;
  
   mrc_init_head(outhead);



   outhead->nx=vol.X;
   outhead->ny=vol.Y;
   outhead->nz=vol.Z;

   if(myid==0)
   {
	   mrc_write_head(fout,outhead);

   }


 
/************************malloc memory for prj, slc, slc->data****************************************/
   double *x_coef;
   x_coef=(double *)malloc(sizeof(double)*prj.AngN*10);
   memset(x_coef, 0 , sizeof(double)*prj.AngN*10);
   
   double *y_coef;
   y_coef=(double *)malloc(sizeof(double)*prj.AngN*10);
   memset(y_coef, 0 , sizeof(double)*prj.AngN*10);

 
   if(myid==0)
   {   
       FILE *f_coef;
      f_coef=fopen(coef,"r");
      if(!f_coef)
      {
       printf("Cannot open file coef strike any key exit!");
       return 1;
       }
     
      read_coef(x_coef, y_coef, f_coef);

      fclose(f_coef);
   }
    MPI_Bcast(x_coef, prj.AngN*10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_coef, prj.AngN*10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);  

  //  printf("x_coef[35] is %lf, y_coef[35] is %lf in process %d\n", x_coef[35], y_coef[35], myid);
   int line_allnum=prj.X*prj.Y*prj.AngN;
  // printf("line_allnum is %d in process %d\n",line_allnum,myid);

   int line_num=prj.X*prj.Y;
   int pixel_num=vol.X*vol.Y*vol.Z;
   int volZ_per;
   int volZ_add=vol.Z%mypro;
   int Z_start;   //the start slice of reproject per process
   int Z_end  ;   //the end slice of reproject per process
   if(myid<volZ_add)
   {
     Z_start=vol.Zstart+(vol.Z/mypro+1)*myid;
     Z_end = vol.Zstart+(vol.Z/mypro+1)*(myid+1);
    /* Z_start=-43+(vol.Z/mypro+1)*myid;
     Z_end = -43+(vol.Z/mypro+1)*(myid+1);*/
     volZ_per=vol.Z/mypro+1;
   }
   else
   {
     Z_start=vol.Zstart+(vol.Z/mypro)*myid+volZ_add;
     Z_end = vol.Zstart+(vol.Z/mypro)*(myid+1)+volZ_add;
    /* Z_start=-43+(vol.Z/mypro)*myid+volZ_add;
     Z_end = -43+(vol.Z/mypro)*(myid+1)+volZ_add;*/
     volZ_per=vol.Z/mypro;
   }
   printf("Z_start is %d, Z_end is %d in process %d\n",Z_start,Z_end,myid);

     
   int pixel_num_per=vol.X*vol.Y*volZ_per;
  // printf("pixel_num_per is %d in process %d\n",pixel_num_per,myid);

   float *prj_real, *slc_data;

/*for input file*/
   if((prj_real=(float *)malloc(sizeof(float)*line_allnum))==NULL)
   {
       printf("Error with Function 'ATOM()'!Can't malloc memery for 'prj_real'!");
       return FALSE;
   }
   memset(prj_real, 0 , sizeof(float)*line_allnum);
 //  printf("hello1\n");

/*for output file*/
   if((slc_data=(float *)malloc(sizeof(float)*pixel_num_per))==NULL)
   {
       printf("Error with Function 'ATOM()'!Can't malloc memery for 'slc_data'!");
       return FALSE;
   }
   memset(slc_data, 0 , sizeof(float)*pixel_num_per);
//printf("hello2\n");
   Pixel pixel;
   Weight *comp_prj;
   comp_prj=(Weight *)malloc(sizeof(Weight));
   if((comp_prj=(Weight *)malloc(sizeof(Weight)))==NULL)
   {
       printf("Error with Function 'ATOM()'!Can't malloc memery for 'comp_prj'!");
       return FALSE;
   }


   if(myid==0)
   {
     FILE *fin;
     fin=fopen(inf,"r");
     if(!fin)
      {
       printf("Cannot open file inf strike any key exit!");
       return 1;
       }
     mrc_read_all(fin, inhead, prj_real);
     fclose(fin);
   }

   MPI_Bcast(prj_real, line_allnum, MPI_FLOAT, 0, MPI_COMM_WORLD);

    int start=Z_start-vol.Zstart;
    int end=Z_end-vol.Zstart;
    mrc_read_block(f_fbp, fbphead, start, end, 'Z', slc_data); 


    MPI_Barrier(MPI_COMM_WORLD);

    int i;


      {

       float  *prj_calc_all, *prj_calc_num_all, *global_prj_calc_all, *global_prj_calc_num_all; 
      
       if((prj_calc_all=(float *)malloc(sizeof(float)*line_allnum))==NULL)
        {
          printf("Error with Function 'ATOM()'!Can't malloc memery for 'prj_calc_all'!");
          return FALSE;
        }
       memset(prj_calc_all, 0 , sizeof(float)*line_allnum);

       if((prj_calc_num_all=(float *)malloc(sizeof(float)*line_allnum))==NULL)
        {
          printf("Error with Function 'ATOM()'!Can't malloc memery for 'prj_calc_num_all'!");
          return FALSE;
        }
       memset(prj_calc_num_all, 0 , sizeof(float)*line_allnum);

       if((global_prj_calc_all=(float *)malloc(sizeof(float)*line_allnum))==NULL)
       {
          printf("Error with Function 'ATOM()'!Can't malloc memery for 'global_prj_calc_all'!");
          return FALSE;
       }
       memset(global_prj_calc_all, 0 , sizeof(float)*line_allnum);

       if((global_prj_calc_num_all=(float *)malloc(sizeof(float)*line_allnum))==NULL)
       {
          printf("Error with Function 'ATOM()'!Can't malloc memery for 'global_prj_calc_num_all'!");
          return FALSE;
       }
       memset(global_prj_calc_num_all, 0 , sizeof(float)*line_allnum);



        for(i=0;i<ATOM_ITR_NUM;i++)
          {  
              SIRT_update_slice(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all, prj_real,slc_data, global_prj_calc_all, global_prj_calc_num_all, ATOM_ITR_STEP, Z_start,Z_end);

              MPI_Barrier(MPI_COMM_WORLD);


              memset(prj_calc_all, 0 , sizeof(float)*line_allnum);
              memset(prj_calc_num_all, 0 , sizeof(float)*line_allnum);

              memset(global_prj_calc_all, 0 , sizeof(float)*line_allnum);  
              memset(global_prj_calc_num_all, 0 , sizeof(float)*line_allnum); 
               
          }//end for i
        free(prj_calc_num_all);
        free(prj_calc_all);
        free(global_prj_calc_all);
        free(global_prj_calc_num_all);


        }//end if SIRT


   

   MPI_Barrier(MPI_COMM_WORLD);
   Z_start -= vol.Zstart;
   Z_end -= vol.Zstart;
     
   mrc_write_all(fout, outhead, Z_start,Z_end,slc_data);



   MPI_Barrier(MPI_COMM_WORLD);
    if(!myid) {
   // mrc_flipyz(out_temp,outf,id,p);
        mrc_update_head(fout);
   // mrc_update_head(inf);
          }     
 
   free(comp_prj);
   free(slc_data);
    
   free(prj_real);
 
   free(y_coef);
   free(x_coef); 

   free(inhead);
   inhead=NULL;

   free(fbphead);
   fbphead=NULL;

   free(outhead); 
   outhead=NULL;
   
   MPI_File_close(&f_fbp);
   MPI_File_close(&fout);
   return TRUE;
}







