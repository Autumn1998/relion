//gcc -o atom atom.c mrcfile.c -lm
#include <time.h>
#include "mrcfile_atom.h"

#include "atom.h"
#include "filter-prj.h"

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
/* //order=2;
  index=10*angle;
  x=x_coef[index]+x_coef[index+1]*pixel.X+x_coef[index+2]*pixel.Y+x_coef[index+3]*pixel.Z+x_coef[index+4]*pixel.X*pixel.X+x_coef[index+5]*pixel.X*pixel.Y+x_coef[index+6]*pixel.X*pixel.Z+x_coef[index+7]*pixel.Y*pixel.Y+x_coef[index+8]*pixel.Y*pixel.Z+x_coef[index+9]*pixel.Z*pixel.Z;
  y=y_coef[index]+y_coef[index+1]*pixel.X+y_coef[index+2]*pixel.Y+y_coef[index+3]*pixel.Z+y_coef[index+4]*pixel.X*pixel.X+y_coef[index+5]*pixel.X*pixel.Y+y_coef[index+6]*pixel.X*pixel.Z+y_coef[index+7]*pixel.Y*pixel.Y+y_coef[index+8]*pixel.Y*pixel.Z+y_coef[index+9]*pixel.Z*pixel.Z;
*/

   
//  order=1
   index=4*angle;
   x=x_coef[index]+x_coef[index+1]*pixel.X+x_coef[index+2]*pixel.Y+x_coef[index+3]*pixel.Z;
   y=y_coef[index]+y_coef[index+1]*pixel.X+y_coef[index+2]*pixel.Y+y_coef[index+3]*pixel.Z;
   
 /*  x--;
   y--;*/

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
int val_coef_recur(int n, Pixel pixel, Volume vol, int angle,double *x_coef, double *y_coef,  Weight *comp_prj, double *x_val, double *qx_pre_val, double *qx_val, double *y_val, double *qy_pre_val, double *qy_val)
{

  double x,y;
/*  int x_lef,y_low;
  double x_del,y_del;*/
  int index;
  index=10*angle; 
  int i,m;



 // if (pixel.X==-87)
   if (pixel.X==vol.Xstart)
   {
     for(i=0;i<=n;i++)
      {
       x_val[i]=compX_ini(index,x_coef,pixel.X+i,pixel.Y,pixel.Z);
       y_val[i]=compY_ini(index,y_coef,pixel.X+i,pixel.Y,pixel.Z);
      }

     qx_val[0]=x_val[0];
     qx_val[1]=x_val[1]-x_val[0];
     qx_val[2]=x_val[2]-x_val[1]-qx_val[1];       

     qy_val[0]=y_val[0];
     qy_val[1]=y_val[1]-y_val[0];
     qy_val[2]=y_val[2]-y_val[1]-qy_val[1]; 
/*if (pixel.Y==0 && pixel.Z==0 && index==0)
        printf("qx_val[0] is %f, qx_val[1] is %f, qx_val[2] is %f, qy_val[0] is %f, qy_val[1] is %f, qy_val[2] is %f in pixel(%d,%d,%d)\n", qx_val[0],qx_val[1],qx_val[2],qy_val[0],qy_val[1],qy_val[2],pixel.X,pixel.Y,pixel.Z);*/
    }

   else 
    {
     for(m=0;m<=n;m++)
      {
       qx_pre_val[m]=qx_val[m];
       qy_pre_val[m]=qy_val[m];
      }
     
     for(m=n-1;m>=0;m--)
      {
       qx_val[m]=qx_pre_val[m]+qx_pre_val[m+1];
       qy_val[m]=qy_pre_val[m]+qy_pre_val[m+1];
      }
     }
    x=qx_val[0];
    y=qy_val[0];


    x--;
    y--;


    comp_prj->x_min=floor(x);
    comp_prj->y_min=floor(y);
   
    comp_prj->x_min_del=x-comp_prj->x_min;      
    comp_prj->y_min_del=y-comp_prj->y_min;
   

   
  return TRUE;

}

/**********************************************************************************/
int Slice_backproj_nearest(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_real,float *slc_data,int Z_start, int Z_end, int slcN_per)
{
  double s;
  double c;

  int angle,n,index;
  int slc=0;

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
   //   printf("slice %d is beginning.\n",pixel.Z);
      
      for (pixel.Y=-31;pixel.Y<653;pixel.Y++)
     // for (pixel.Y=4;pixel.Y<495;pixel.Y++)
       {
       for (pixel.X=-87;pixel.X<564;pixel.X++)
      // for (pixel.X=-51;pixel.X<547;pixel.X++)
         {
         s=0;
         c=0;
         for(angle=0;angle<prj.AngN;angle++)
             {
             val_coef(pixel, angle,x_coef, y_coef, comp_prj); 

            /* if (pixel.X==-62 && pixel.Z==82 && pixel.Z==-12)
                printf("comp_prj->x is %d, comp_prj->y is %d in angle=%d\n", comp_prj->x, comp_prj->y,angle);*/
             if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
             if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;

             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                 {
                n=comp_prj->x_min+comp_prj->y_min*prj.X+angle*prj.X*prj.Y; //prj index
                //s+=comp_prj->del*prj_real[n];
                  s+=prj_real[n];
                //c+=comp_prj->del;
                   c++;
                 }
            }//end angle              
            index=(pixel.X+87)+(pixel.Y+31)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
            //index=(pixel.X+51)+(pixel.Y-4)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y; 
            if(c!=0.0f)
                slc_data[index]=(float)(s/c);
            // slc_data[index]=(float)(s/prj.AngN);
          }      
         }// end pixel
    //   printf("slice %d is finished to backprojection.\n",pixel.Z);  
       slcN_per++;     
                
     }// end pixel.Z

   printf("%d slices has been finished in Z_start%d\n",slcN_per,Z_start);

 /*  MPI_Reduce(&slcN_per,&global_slcN,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

   if(!myid) printf("%d slice is reconstructed by BPT!\n",global_slcN);*/

   return TRUE;

}
/**********************************************************************************/
int Slice_backproj_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_real,float *slc_data,int Z_start, int Z_end, int slcN_per)
{
  double s;
  double c;
  

  int angle,n,index;
  int slc=0;

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
      printf("slice %d is beginning.\n",pixel.Z);
      
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)

      {
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)

         {
         index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index

         s=0;
         c=0;
         for(angle=0;angle<prj.AngN;angle++)
           {
             val_coef(pixel, angle,x_coef, y_coef, comp_prj); 
        /* if (pixel.X==-62 && pixel.Y==82 && pixel.Z==-12)
                printf("comp_prj->x_min is %d, comp_prj->y_min is %d in angle=%d\n", comp_prj->x_min, comp_prj->y_min,angle);*/
            /* comp_prj->x_min=pixel.X;
              comp_prj->x_min_del=0;
              comp_prj->y_min=pixel.Y;
              comp_prj->y_min_del=0;*/
             
             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                 {
                    n=comp_prj->x_min+comp_prj->y_min*prj.X+angle*prj.X*prj.Y; //prj index
                    s+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*prj_real[n];
                    c+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);
                  }
             if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min+1, y_min)
                 {
                    n=comp_prj->x_min+1+comp_prj->y_min*prj.X+angle*prj.X*prj.Y; //prj index
                    s+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*prj_real[n];
                    c+=comp_prj->x_min_del*(1-comp_prj->y_min_del);
                 }
             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                 {
                    n=comp_prj->x_min+(comp_prj->y_min+1)*prj.X+angle*prj.X*prj.Y; //prj index
                    s+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*prj_real[n];
                    c+=(1-comp_prj->x_min_del)*comp_prj->y_min_del;            
                 }
             if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min+1, y_min+1)
                 {
                    n=comp_prj->x_min+1+(comp_prj->y_min+1)*prj.X+angle*prj.X*prj.Y; //prj index
                    s+=comp_prj->x_min_del*comp_prj->y_min_del*prj_real[n];
                    c+=comp_prj->x_min_del*comp_prj->y_min_del;                  
                 }
             // printf("hello world in angle %d (%d,%d,%d)\n",angle,pixel.X,pixel.Y,pixel.Z);
             
            }//end angle              
            
            if(c!=0.0f)
                slc_data[index]=(float)(s/c);
            // slc_data[index]=(float)(s/prj.AngN);
          }      
         }// end pixel
    //   printf("slice %d is finished to backprojection.\n",pixel.Z);  
       slcN_per++;     
                
     }// end pixel.Z

   printf("%d slices has been finished in Z_start%d\n",slcN_per,Z_start);

 /*  MPI_Reduce(&slcN_per,&global_slcN,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

   if(!myid) printf("%d slice is reconstructed by BPT!\n",global_slcN);*/

   return TRUE;

}
/**********************************************************************************/
int  ATOM(char *inf, char *outf, char *coef, int myid, int mypro) //Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
{
  
   Volume vol;
   Projection prj;

/************************File and file head processing  ****************************************/
   MPI_File fout;
   
  // MPI_File_open(MPI_COMM_WORLD, inf, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fin);
   MPI_File_open(MPI_COMM_WORLD, outf, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fout);

 /*  MPI_File f_filter;
   MPI_File_open(MPI_COMM_WORLD, filter_prj, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &f_filter);
*/
   MrcHeader  *inhead,  *outhead;
   inhead=(MrcHeader *)malloc(sizeof(MrcHeader));
   outhead=(MrcHeader *)malloc(sizeof(MrcHeader));
 //  filter_head=(MrcHeader *)malloc(sizeof(MrcHeader));
   

   
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
    } 
    MPI_Bcast(inhead, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);

//   printf("head->cmap[0] is %c\n",inhead->cmap[0]);

   prj.X=inhead->nx;
   prj.Y=inhead->ny;
   prj.AngN=inhead->nz;

   vol.X=1467;//x->-87.000000:1.000000:563.000000
   vol.Y=1521; //y->-31.000000:1.000000:652.000000
   vol.Z=58; //z->-62.000000:1.000000:38.000000

   vol.Xstart=-259;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=-209;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-32;
   vol.Zend=vol.Zstart+vol.Z;
  
   mrc_init_head(outhead);
   outhead->nx=vol.X;
   outhead->ny=vol.Y;
   outhead->nz=vol.Z;

   outhead->nxstart=vol.Xstart;
   outhead->nystart=vol.Ystart;
   outhead->nzstart=vol.Zstart;
/*
   outhead->mx=651;
   outhead->my=684;
   outhead->mz=101;

   outhead->xlen=13020;
   outhead->ylen=13680;
   outhead->zlen=2020;*/

   if(myid==0)
   mrc_write_head(fout,outhead);

/*   mrc_init_head(filter_head);
   filter_head->nx=prj.X;
   filter_head->ny=prj.Y;
   filter_head->nz=prj.AngN;
   if(myid==0)
       mrc_write_head(f_filter,filter_head);*/



 
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

 // printf("hello3\n");

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
   
  // printf("prj_real[345] is %f in process %d\n",prj_real[345],myid);


  
   int slcN_per=0;
   int global_slcN=0;

   int filtlength=251;
   //char *filter="RamLak";
   char *filter="SheppLogan";
   int row_pad = 2;  //x-padded 
   int col_pad = 2;  //y-padded

   int symmetrize_2D_flag=1;
 
   int angle=0;



   for(angle=0;angle<prj.AngN;angle++)
       filter_prj_sym(prj_real, filter, filtlength, prj, row_pad, col_pad, symmetrize_2D_flag, angle, myid);
/*
   if(!myid)
       mrc_write_all(f_filter,filter_head,0,prj.AngN,prj_real);*/
       
     // Slice_backproj_nearest(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_real,slc_data,Z_start,Z_end,slcN_per); //nearest
   Slice_backproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_real,slc_data,Z_start,Z_end,slcN_per); //four-weight
  
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


   free(outhead); 
   outhead=NULL;
   
  // MPI_File_close(&fin);
   MPI_File_close(&fout);
   return TRUE;
}







