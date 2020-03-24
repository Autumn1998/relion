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
   /* x_lef=(int)floor(x);
    y_low=(int)floor(y);

   if((x-x_lef)<=0.5)
     {
       comp_prj->x=x_lef;
       x_del=x-x_lef;
     }
   else
     {  
       comp_prj->x=x_lef+1;
       x_del=x_lef+1-x;
     }
   if((y-y_low)<=0.5)
     {
       comp_prj->y=y_low;
       y_del=y-y_low;
     }
   else
     {  
       comp_prj->y=y_low+1;
       y_del=y_low+1-y;
     }*/

   
  return TRUE;

}


//SIRT
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
      
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       {
       for (pixel.X=vol.Xstart;pixel.X<vol.Xstart;pixel.X++)
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
            index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
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
    //  printf("slice %d is beginning.\n",pixel.Z);
      
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

  // printf("%d slices has been finished in Z_start%d\n",slcN_per,Z_start);

 /*  MPI_Reduce(&slcN_per,&global_slcN,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

   if(!myid) printf("%d slice is reconstructed by BPT!\n",global_slcN);*/

   return TRUE;

}


//SIRT
/**********************************************************************************/
int SIRT_Slice_reproj_nearest (Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_calc_num_all,float *slc_data, int Z_start, int Z_end)
{

   int n=2;//order
   int j;
   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));


   int pixel_index,prj_index;


   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
           pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
           for(j=0;j<prj.AngN;j++)
             {
              val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //val_coef_recur(n, pixel, vol, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);
          
             if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
             if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;


             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                {
                
                  prj_index=comp_prj->x_min+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj index
                  prj_calc_all[prj_index]+=slc_data[pixel_index];
                  prj_calc_num_all[prj_index]++;

                 }//end if
              }//end for j
         }//end for pixel

    //  printf("slice %d is finished to reproj.\n",pixel.Z);
     }//end for pixel.Z
   
   return TRUE;
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

     printf("slice %d is finished to reproj.\n",pixel.Z);
     }//end for pixel.Z
   
   return TRUE;

}




/**********************************************************************************/
int SIRT_one_update_nearest(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_real,float *slc_data,float ASART_ITR_STEP, int Z_start, int Z_end)
{

  int n=2;//order
  int j;
   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));

   int pixel_index,prj_cal_index,prj_real_index;
   float del;

   double s,c;

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
             pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
             s=0;
             c=0;

             for(j=0;j<prj.AngN;j++)
             {
             val_coef(pixel, j, x_coef, y_coef, comp_prj); 
            // val_coef_recur(n, pixel, vol, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);
            
             if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
             if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;

             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                 {
                    prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X+j*prj.X*prj.Y; //prj_calc index
                 //   prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                  //  del=prj_real[prj_real_index]-prj_calc[prj_cal_index];
                  //   del=prj_calc_all[prj_cal_index];
                     s+=prj_calc_all[prj_cal_index];
                     c++;

          if(c!=0.0f)
                slc_data[pixel_index] += (float)(s/c)*ASART_ITR_STEP;
                 } //end if
             }//end for prj.AngN
           
         }// end pixel
    //   printf("slice %d is completed to update.\n",pixel.Z);       
                
     }// end pixel.Z

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
              for(j=0;j<prj.AngN;j++)
            	//  for(j=0;j<1;j++)
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
float NCC(float * test, float * data, int width, int height) {

	float *tmp1, *tmp2, dataavg, testavg, t1, t2, t3;

	tmp1 = new float[width * height];
	tmp2 = new float[width * height];
	memset(tmp1, 0, sizeof(float) * width * height);
	memset(tmp2, 0, sizeof(float) * width * height);

	dataavg = testavg = 0;
	int count = 0;

	float blankv = 0;
	for (int i = 0; i < width * height; i++) {
		if (data[i] == 0 || test[i] == blankv) {
			continue;
		}
		dataavg += data[i];
		testavg += test[i];
		count++;
	}

	dataavg /= count;
	testavg /= count;

	for (int i = 0; i < width * height; i++) {
		if (data[i] == 0 || test[i] == blankv) {
			continue;
		}
		tmp1[i] = data[i] - dataavg;
		tmp2[i] = test[i] - testavg;
	}

	t1 = t2 = t3 = 0;
	for (int i = 0; i < width * height; i++) {
		t1 += tmp1[i] * tmp2[i];
		t2 += tmp1[i] * tmp1[i];
		t3 += tmp2[i] * tmp2[i];
	}

	delete[] tmp1;
	delete[] tmp2;

	return t1 / sqrt(t2 * t3);
}
/**********************************************************************************/
int SIRT_update_slice(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_calc_num_all, float *prj_real,float *slc_data, float *global_prj_calc_all, float *global_prj_calc_num_all, float ATOM_ITR_STEP,int Z_start, int Z_end,int myid)
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
int SIRT_update_slice_NCC(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc_all, float *prj_calc_num_all, float *prj_real,float *slc_data, float *global_prj_calc_all, float *global_prj_calc_num_all, float ATOM_ITR_STEP,int Z_start, int Z_end,int myid)
{

   int n;
   int line_allnum=prj.X*prj.Y*prj.AngN;

   int slcN_per=0;
   int global_slcN=0;



   SIRT_Slice_reproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all, slc_data, Z_start, Z_end);


    MPI_Allreduce(prj_calc_all,global_prj_calc_all,line_allnum,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(prj_calc_num_all,global_prj_calc_num_all,line_allnum,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

    float errordata1=0;
    float errordata2=0;
    for(n=0;n<line_allnum;n++)
    {
       if(global_prj_calc_num_all[n] != 0)
    	   global_prj_calc_all[n]/=global_prj_calc_num_all[n];
    }


    if(!myid)
    {
    	double average=0;
    	double ncc_arr[200];
    	printf("The Error : %f %f \n",errordata1,errordata2);
		for (int i = 0; i < prj.AngN; i++) {
			ncc_arr[i] = NCC(prj_real+(prj.X*prj.Y*i), global_prj_calc_all+(prj.X*prj.Y*i), prj.X, prj.Y);
			average += ncc_arr[i];
		    printf("%f\t",ncc_arr[i]);
		}
    	printf("\n");
		average=average/prj.AngN;
		printf("aver : %f \n",average);
    }

    for(n=0;n<line_allnum;n++)
    {

        if(global_prj_calc_all[n]!=0)
     	   global_prj_calc_all[n]=prj_real[n]-global_prj_calc_all[n];  //prj error
        errordata1+=global_prj_calc_all[n];
        errordata2+=abs(global_prj_calc_all[n]);
    }
     SIRT_one_update_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc_all, prj_real, slc_data, ATOM_ITR_STEP,Z_start,Z_end);


   MPI_Barrier(MPI_COMM_WORLD);




   return TRUE;

}
//SART
/**********************************************************************************/
int Sart_Slice_reproj_nearest(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *prj_calc_num,float *slc_data,int j, int Z_start, int Z_end, double *prj_calc_abs)
{
   
   int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));


   int pixel_index,prj_index;
   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
           pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index

           val_coef(pixel, j, x_coef, y_coef, comp_prj); 
          //val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);
          
           if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
           if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;


           if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                 {
                
                prj_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj index
                prj_calc[prj_index]+=slc_data[pixel_index];
                prj_calc_abs[prj_index] += abs(slc_data[pixel_index]);
                prj_calc_num[prj_index]++;

              //  if (pixel_index ==325 && pixel.Z==-62) printf("prj_index is %d in angle %d\n",prj_index, j);
              //  if (pixel_index==325 && pixel.Z==-62) printf("prj_calc[%d] is %f, prj_calc_num[%d] is %f in angle %d\n",prj_index, prj_calc[prj_index],prj_index,prj_calc_num[prj_index],j);
                 }//end if
         
         }//end for pixel

    //  printf("slice %d is finished to reproj.\n",pixel.Z);
     }//end for pixel.Z
   
   return TRUE;

}



/**********************************************************************************/
int Sart_Slice_reproj_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *prj_calc_num,float *slc_data,int j, int Z_start, int Z_end, double *prj_calc_abs)
{
   
   int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));


   int pixel_index,prj_index;
   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
           pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index

           val_coef(pixel, j, x_coef, y_coef, comp_prj); 
         // val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);


           if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                 {
                
                   prj_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj index
                   prj_calc[prj_index]+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*slc_data[pixel_index];
                   prj_calc_abs[prj_index] += (1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*abs(slc_data[pixel_index]);
                   prj_calc_num[prj_index] += (1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);
    
                 }//end if
           if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min+1, y_min)
                 {
                
                   prj_index=comp_prj->x_min+1+comp_prj->y_min*prj.X; //prj index
                   prj_calc[prj_index]+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*slc_data[pixel_index];
                   prj_calc_abs[prj_index] += comp_prj->x_min_del*(1-comp_prj->y_min_del)*abs(slc_data[pixel_index]);
                   prj_calc_num[prj_index] += comp_prj->x_min_del*(1-comp_prj->y_min_del);
    
                 }//end if
           if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                 {
                    prj_index=comp_prj->x_min+(comp_prj->y_min+1)*prj.X; //prj index
                    prj_calc[prj_index]+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*slc_data[pixel_index];
                    prj_calc_abs[prj_index] += (1-comp_prj->x_min_del)*comp_prj->y_min_del*abs(slc_data[pixel_index]);
                    prj_calc_num[prj_index] += (1-comp_prj->x_min_del)*comp_prj->y_min_del;           
                 }
           if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min+1, y_min+1)
                 {
                    prj_index=(comp_prj->x_min+1)+(comp_prj->y_min+1)*prj.X; //prj index
                    prj_calc[prj_index]+=comp_prj->x_min_del*comp_prj->y_min_del*slc_data[pixel_index];
                    prj_calc_abs[prj_index] += comp_prj->x_min_del*comp_prj->y_min_del*abs(slc_data[pixel_index]);
                    prj_calc_num[prj_index] += comp_prj->x_min_del*comp_prj->y_min_del;   
                 }                
         }//end for pixel

    //  printf("slice %d is finished to reproj.\n",pixel.Z);
     }//end for pixel.Z
   
   return TRUE;

}


/**********************************************************************************/
int proj_calc_sart_div  (Projection prj, float *prj_calc, float *prj_calc_num)
{
   int line_num=prj.X*prj.Y;
   int i;
   for(i=0;i<line_num;i++)
    if(prj_calc_num[i]!=0)
     prj_calc[i]=prj_calc[i]/prj_calc_num[i];

   return TRUE;
}

/**********************************************************************************/
int SART_one_update_nearest(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *slc_data,int j,float ASART_ITR_STEP, int Z_start, int Z_end)
{

  int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));

   int pixel_index,prj_cal_index,prj_real_index;
   float del;

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
             pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
             
             val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);
            
             if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
             if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;

             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                 {
                    prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj_calc index
                 //   prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                  //  del=prj_real[prj_real_index]-prj_calc[prj_cal_index];
                     del=prj_calc[prj_cal_index];

                     slc_data[pixel_index]+=del*ASART_ITR_STEP;
                 } //end if
             
           
         }// end pixel
    //   printf("slice %d is completed to update.\n",pixel.Z);       
                
     }// end pixel.Z

   return TRUE;

}


/**********************************************************************************/
int SART_one_update_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *slc_data,int j,float ASART_ITR_STEP, int Z_start, int Z_end)
{

  int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));

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

             val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);

            
            if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                {
                   prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);                    
                 } //end if
                
              if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)//(x_min+1, y_min)
                {
                   prj_cal_index=comp_prj->x_min+1+comp_prj->y_min*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index]; 
                   c+=comp_prj->x_min_del*(1-comp_prj->y_min_del);                    
                 } //end if
              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+(1+comp_prj->y_min)*prj.X; //prj_calc index
               //    prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*comp_prj->y_min_del;                    
                 } //end if
               if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y)   //(x_min+1, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+1+(comp_prj->y_min+1)*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   c+=comp_prj->x_min_del*comp_prj->y_min_del;                    
                 } //end if

           if(c!=0.0f)
                slc_data[pixel_index] += (float)(s/c)*ASART_ITR_STEP;

                
         }// end pixel
    //   printf("slice %d is completed to update.\n",pixel.Z);       
                
     }// end pixel.Z

   return TRUE;

}
/**********************************************************************************/
int SART_update_slice(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *prj_calc_num, float *prj_real,float *slc_data, float *global_prj_calc, float *global_prj_calc_num, float ATOM_ITR_STEP,int Z_start, int Z_end, double *prj_calc_abs,double *global_prj_calc_abs)
{
   /*int *a;
   a=(int *)malloc(sizeof(int)*AngN);
   MLS(a,AngN);*/
   int i,j;
   int n,m;
   int line_num=prj.X*prj.Y;

   int slcN_per=0;
   int global_slcN=0;

   for(j=0;j<prj.AngN;j++)
   { 
      //  j=a[i];

    //Sart_Slice_reproj_nearest(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc, prj_calc_num,slc_data,j,Z_start,Z_end,prj_calc_abs); 
    Sart_Slice_reproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc, prj_calc_num,slc_data,j,Z_start,Z_end,prj_calc_abs);  
   //proj_calc_sart_div(prj, prj_calc, prj_calc_num);

   
   
    MPI_Allreduce(prj_calc,global_prj_calc,line_num,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(prj_calc_num,global_prj_calc_num,line_num,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(prj_calc_abs,global_prj_calc_abs,line_num,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     
    for(n=0;n<line_num;n++)
      {
        if(global_prj_calc_num[n]!=0)
           {
             global_prj_calc[n]/=global_prj_calc_num[n];
            // global_prj_calc_abs[n]/=global_prj_calc_num[n];
           }             
        m=n+j*prj.X*prj.Y;
        global_prj_calc[n]=prj_real[m]-global_prj_calc[n]; //prj error
       }
 
  /* //filter to prj error      
   //prj3 
   int d=5;
   char *filter="ram-lak";
 
   //char *filter="shepp-logan";
   filter_prj(global_prj_calc, filter, d, prj); //filter-prj3 */

  // SART_one_update_nearest(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc, prj_real, slc_data, j, ATOM_ITR_STEP,Z_start,Z_end);
   SART_one_update_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc,slc_data, j, ATOM_ITR_STEP,  Z_start,  Z_end);
  
  // printf("Angle %d has been completed in SART_one_update\n", j);

   MPI_Barrier(MPI_COMM_WORLD);

   memset(prj_calc, 0 , sizeof(float)*line_num);
   memset(prj_calc_num, 0 , sizeof(float)*line_num);

   memset(global_prj_calc, 0 , sizeof(float)*line_num);
   memset(global_prj_calc_num, 0 , sizeof(float)*line_num);

   memset(prj_calc_abs, 0 , sizeof(double)*line_num);
   memset(global_prj_calc_abs, 0 , sizeof(double)*line_num);

    } //end for j

   /*MPI_Reduce(&slcN_per,&global_slcN,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   if(!myid) printf("%d slice is reconstruct in ATOM !\n",global_slcN);*/

   return TRUE;

}

/**********************************************************************************/
int ASART_one_update_nearest(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *slc_data,int j,float ASART_ITR_STEP, int Z_start, int Z_end, float *prj_calc_num, double *prj_calc_abs)
{

  int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));

   int pixel_index,prj_cal_index;
   float del=0.0f;
   long double t=0.0f;

   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
             pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index

             val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);

             if(comp_prj->x_min_del>0.5)
                  comp_prj->x_min += 1;
             if(comp_prj->y_min_del>0.5)
                  comp_prj->y_min += 1;

             if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)
                {
                   prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj_calc index
                   /*prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   del=prj_real[prj_real_index]-prj_calc[prj_cal_index];*/
                   del=prj_calc[prj_cal_index];
        
                /*if (prj_calc[prj_cal_index]!=0 && prj_calc_num[prj_cal_index]!=0)  
                   {
                    t=prj_calc[prj_cal_index]*prj_calc_num[prj_cal_index];                
                    slc_data[pixel_index]+=slc_data[pixel_index]*del/t*ASART_ITR_STEP;
                   }
                 */
                
                 t=prj_calc_abs[prj_cal_index];
                 if (t != 0)
                    slc_data[pixel_index]+=abs(slc_data[pixel_index])*del/t*ASART_ITR_STEP;
                 
               // if (pixel_index ==325 && pixel.Z==-62) printf("prj_calc[%d] is %f, prj_calc_num[%d] is %f, t is %Lf, del is %f, slc_data[%d] is %f in angle %d\n",prj_cal_index, prj_calc[prj_cal_index], prj_cal_index, prj_calc_num[prj_cal_index], t, del, pixel_index, slc_data[pixel_index], j);
                 } //end if
                
         }// end pixel

     //if (pixel_index ==21482 && pixel.Z==-62) printf("prj_calc[%d] is %f, prj_calc_num[%d] is %f, t is %Lf, del is %f, slc_data[%d] is %f in angle %d\n",prj_cal_index, prj_calc[prj_cal_index], prj_cal_index, prj_calc_num[prj_cal_index], t, del, pixel_index, slc_data[pixel_index], j);
      // printf("slice %d is completed to update.\n",pixel.Z);       
     // exit(0);     
       //     if (pixel.Z==-62 && j == (prj.AngN-1)) exit (0);
     }// end pixel.Z
   
   return TRUE;


}


/**********************************************************************************/
int ASART_one_update_bilinear(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *slc_data,int j,float ASART_ITR_STEP, int Z_start, int Z_end, float *prj_calc_num, double *prj_calc_abs)
{

  int n=2;//order

   double *x_val;
   x_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_pre_val;
   qx_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qx_val;
   qx_val=(double *)malloc(sizeof(double)*(n+1));

   double *y_val;
   y_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_pre_val;
   qy_pre_val=(double *)malloc(sizeof(double)*(n+1));

   double *qy_val;
   qy_val=(double *)malloc(sizeof(double)*(n+1));

   int pixel_index,prj_cal_index;
   double s,c,t;


   for(pixel.Z=Z_start;pixel.Z<Z_end;pixel.Z++)
    {
           
      for (pixel.Y=vol.Ystart;pixel.Y<vol.Yend;pixel.Y++)
       for (pixel.X=vol.Xstart;pixel.X<vol.Xend;pixel.X++)
         {
             pixel_index=(pixel.X-vol.Xstart)+(pixel.Y-vol.Ystart)*vol.X+(pixel.Z-Z_start)*vol.X*vol.Y;   //pixel index
             s=0;
             t=0;
             c=0;

             val_coef(pixel, j, x_coef, y_coef, comp_prj); 
             //val_coef_recur(n, pixel, j, x_coef, y_coef, comp_prj,x_val,qx_pre_val,qx_val,y_val,qy_pre_val,qy_val);

            
            if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y) //(x_min, y_min)
                {
                   prj_cal_index=comp_prj->x_min+comp_prj->y_min*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index];
                   t+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del)*prj_calc_abs[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*(1-comp_prj->y_min_del);                    
                 } //end if
                
              if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && comp_prj->y_min>=0 && comp_prj->y_min<prj.Y)//(x_min+1, y_min)
                {
                   prj_cal_index=comp_prj->x_min+1+comp_prj->y_min*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*prj_calc[prj_cal_index]; 
                   t+=comp_prj->x_min_del*(1-comp_prj->y_min_del)*prj_calc_abs[prj_cal_index]; 
                   c+=comp_prj->x_min_del*(1-comp_prj->y_min_del);                    
                 } //end if
              if(comp_prj->x_min>=0 && comp_prj->x_min<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y) //(x_min, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+(1+comp_prj->y_min)*prj.X; //prj_calc index
               //    prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   t+=(1-comp_prj->x_min_del)*comp_prj->y_min_del*prj_calc_abs[prj_cal_index];
                   c+=(1-comp_prj->x_min_del)*comp_prj->y_min_del;                    
                 } //end if
               if((comp_prj->x_min+1)>=0 && (comp_prj->x_min+1)<prj.X && (comp_prj->y_min+1)>=0 && (comp_prj->y_min+1)<prj.Y)   //(x_min+1, y_min+1)
                {
                   prj_cal_index=comp_prj->x_min+1+(comp_prj->y_min+1)*prj.X; //prj_calc index
                 //  prj_real_index=prj_cal_index+j*prj.X*prj.Y;  //prj_real index
                  
                   s+=comp_prj->x_min_del*comp_prj->y_min_del*prj_calc[prj_cal_index];
                   t+=comp_prj->x_min_del*comp_prj->y_min_del*prj_calc_abs[prj_cal_index];
                   c+=comp_prj->x_min_del*comp_prj->y_min_del;                    
                 } //end if


           if (t != 0)
                    slc_data[pixel_index]+=abs(slc_data[pixel_index])*s/t*ASART_ITR_STEP;

                
         }// end pixel
    //   printf("slice %d is completed to update.\n",pixel.Z);       
                
     }// end pixel.Z

   return TRUE;

}


/**********************************************************************************/
int ASART_update_slice(Pixel pixel, Projection prj, Volume vol, Weight *comp_prj, double *x_coef, double *y_coef, float *prj_calc, float *prj_calc_num, float *prj_real,float *slc_data,float *global_prj_calc, float *global_prj_calc_num, float ATOM_ITR_STEP,int Z_start, int Z_end, double *prj_calc_abs, double *global_prj_calc_abs)
{
   /*int *a;
   a=(int *)malloc(sizeof(int)*AngN);
   MLS(a,AngN);*/
   int i,j;
   int n,m;
   int line_num=prj.X*prj.Y;
   for(j=0;j<prj.AngN;j++)
   { 
      //  j=a[i];

  // Sart_Slice_reproj_nearest(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc, prj_calc_num,slc_data,j,Z_start,Z_end,prj_calc_abs);
   Sart_Slice_reproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc, prj_calc_num,slc_data,j,Z_start,Z_end,prj_calc_abs);
   //proj_calc_sart_div(prj, prj_calc, prj_calc_num);

   MPI_Allreduce(prj_calc,global_prj_calc,line_num,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(prj_calc_num,global_prj_calc_num,line_num,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(prj_calc_abs,global_prj_calc_abs,line_num,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   
   for(n=0;n<line_num;n++)
    {
      if(global_prj_calc_num[n]!=0)
         global_prj_calc[n]/=global_prj_calc_num[n];
         
      m=n+j*prj.X*prj.Y;
      global_prj_calc[n]=prj_real[m]-global_prj_calc[n]; //prj error
    }
 
   /* if (j==0)
        for (n=0;n<line_num;n++)
            printf("global_prj_calc[%d] is %f\n",n,global_prj_calc[n]);*/

 
   //filter to prj error      
 /*  //prj3 
   int d=5;
  // char *filter="ram-lak";
 
   char *filter="shepp-logan";
   //filter_prj(global_prj_calc, filter, d, prj); //filter-prj3 
   filter_prj_sym(global_prj_calc, filter, d, prj);*/

/*   //prj4   
   //int filtlength=prj.Y;
   int filtlength=64;
   //char *filter="RamLak";
   char *filter="SheppLogan";
   filter_prj(global_prj_calc, filter, filtlength, prj);//filter-prj4*/
   
  // ASART_one_update_nearest(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc, prj_real, slc_data, j, ATOM_ITR_STEP,Z_start,Z_end, global_prj_calc_num, global_prj_calc_abs);
   
   ASART_one_update_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, global_prj_calc, slc_data, j, ATOM_ITR_STEP,Z_start,Z_end, global_prj_calc_num, global_prj_calc_abs);

  /*if (Z_start==0)
  {
   printf("Angle %d has been completed in ASART_one_update\n", j);
   printf("prj_calc_num[7965] is %f\n", global_prj_calc_num[7965]);
  }*/

   MPI_Barrier(MPI_COMM_WORLD);  

   memset(prj_calc, 0 , sizeof(float)*line_num);
   memset(prj_calc_num, 0 , sizeof(float)*line_num);

   memset(global_prj_calc, 0 , sizeof(float)*line_num);
   memset(global_prj_calc_num, 0 , sizeof(float)*line_num);

   memset(prj_calc_abs, 0 , sizeof(double)*line_num);
   memset(global_prj_calc_abs, 0 , sizeof(double)*line_num);

    
    } //end for j
   return TRUE;

}
/**********************************************************************************/
int  ATOM(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag) //Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
{
  
   Volume vol;
   Projection prj;

/************************File and file head processing  ****************************************/
   MPI_File fout, f_fbp;


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

   vol.Xstart=-23;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=0;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-32;
   vol.Zend=vol.Zstart+vol.Z;
  



   printf("%d %d %d \n",vol.Xstart,vol.Ystart,vol.Zstart);
   printf("%d %d %d \n",vol.X,vol.Y,vol.Z);
   mrc_init_head(outhead);

   outhead->nxstart=vol.Xstart;
   outhead->nystart=vol.Ystart;
   outhead->nzstart=vol.Zstart;

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
     // printf("%f %f \n",x_coef[0],y_coef[0]);

      fclose(f_coef);
   }

    MPI_Bcast(x_coef, prj.AngN*10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_coef, prj.AngN*10, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);  

  //  printf("x_coef[35] is %lf, y_coef[35] is %lf in process %d\n", x_coef[35], y_coef[35], myid);
   int line_allnum=prj.X*prj.Y*prj.AngN;
   printf("line_allnum is %d in process %d\n",line_allnum,myid);

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
    if(initflag==0)
    {
    	for(int i=0;i<pixel_num_per;i++)
    		slc_data[i]=0;
    }
    else
    {
    	mrc_read_block(f_fbp, fbphead, start, end, 'Z', slc_data);
    }



    MPI_Barrier(MPI_COMM_WORLD);

    int i;

    if(!strcmp(Method,"SIRT")) 
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

              SIRT_update_slice(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all,
            		  prj_real,slc_data, global_prj_calc_all, global_prj_calc_num_all,
            		  ATOM_ITR_STEP, Z_start,Z_end,myid);

              MPI_Barrier(MPI_COMM_WORLD);


              memset(prj_calc_all, 0 , sizeof(float)*line_allnum);
              memset(prj_calc_num_all, 0 , sizeof(float)*line_allnum);

              memset(global_prj_calc_all, 0 , sizeof(float)*line_allnum);  
              memset(global_prj_calc_num_all, 0 , sizeof(float)*line_allnum); 
              if(myid==0)
            	  printf("iter %d finsihed \n",i);
               
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



int  ATOM_NCC(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag) //Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
{

   Volume vol;
   Projection prj;

/************************File and file head processing  ****************************************/
   MPI_File fout, f_fbp;


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

   vol.Xstart=-5;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=14;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-33;
   vol.Zend=vol.Zstart+vol.Z;




   printf("%d %d %d \n",vol.Xstart,vol.Ystart,vol.Zstart);
   printf("%d %d %d \n",vol.X,vol.Y,vol.Z);
   mrc_init_head(outhead);

   outhead->nxstart=vol.Xstart;
   outhead->nystart=vol.Ystart;
   outhead->nzstart=vol.Zstart;

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
    if(initflag==0)
    {
    	for(int i=0;i<pixel_num_per;i++)
    		slc_data[i]=0;
    }
    else
    {
    	mrc_read_block(f_fbp, fbphead, start, end, 'Z', slc_data);
    }



    MPI_Barrier(MPI_COMM_WORLD);

    int i;




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
        	SIRT_update_slice_NCC(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all,
            		  prj_real,slc_data, global_prj_calc_all, global_prj_calc_num_all,
            		  ATOM_ITR_STEP, Z_start,Z_end,myid);

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

int  ATOM_SIMPLE(char *inf,char *outf,  char *fbpf, char *coef, float ATOM_ITR_STEP, int ATOM_ITR_NUM,char *Method,int myid,int mypro,int initflag) //Method ="SIRT" means sirt, Method = "ART" means ART;Method ="SART" means sart
{

   Volume vol;
   Projection prj;

/************************File and file head processing  ****************************************/
   MPI_File fout, f_fbp;


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

   vol.Xstart=-5;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=14;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-33;
   vol.Zend=vol.Zstart+vol.Z;




   printf("%d %d %d \n",vol.Xstart,vol.Ystart,vol.Zstart);
   printf("%d %d %d \n",vol.X,vol.Y,vol.Z);
   mrc_init_head(outhead);

   outhead->nxstart=vol.Xstart;
   outhead->nystart=vol.Ystart;
   outhead->nzstart=vol.Zstart;

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

   int line_allnum=prj.X*prj.Y*prj.AngN;

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
     volZ_per=vol.Z/mypro+1;
   }
   else
   {
     Z_start=vol.Zstart+(vol.Z/mypro)*myid+volZ_add;
     Z_end = vol.Zstart+(vol.Z/mypro)*(myid+1)+volZ_add;
     volZ_per=vol.Z/mypro;
   }
   printf("Z_start is %d, Z_end is %d in process %d\n",Z_start,Z_end,myid);


   int pixel_num_per=vol.X*vol.Y*volZ_per;

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
    if(initflag==0)
    {
    	for(int i=0;i<pixel_num_per;i++)
    		slc_data[i]=0;
    }
    else
    {
    	mrc_read_block(f_fbp, fbphead, start, end, 'Z', slc_data);
    }



    MPI_Barrier(MPI_COMM_WORLD);

    int i;

    if(!strcmp(Method,"SIRT"))
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
              SIRT_update_slice(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_calc_all, prj_calc_num_all,
            		  prj_real,slc_data, global_prj_calc_all, global_prj_calc_num_all,
            		  ATOM_ITR_STEP, Z_start,Z_end,myid);

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



