#include "mrcfile_atom.h"
#include "mpi.h"

/*******************************************************************************************/
long get_file_size(FILE *fin)
{
	fseek(fin,0,SEEK_END);

	return ftell(fin);

}

/*******************************************************************************************/
//head

/*******************************************************************************************/
int mrc_init_head(MrcHeader *head)
{
  head->nx=0;
  head->ny=0;
  head->nz=0;

  head->mode=MRC_MODE_FLOAT;

  head->nxstart=0;
  head->nystart=0;
  head->nzstart=0;

  head->mx=1;
  head->my=1;
  head->mz=1;

  head->xlen=1;
  head->ylen=1;
  head->zlen=1;

  head->alpha=90;
  head->beta=90;
  head->gamma=90;

  head->mapc=1;
  head->mapr=2;
  head->maps=3;

  head->amin=0;
  head->amax=255;
  head->amean=128;

  head->ispg=0;
  head->nsymbt=0;

  head->next=0;

  head->creatid=1000;
  head->cmap[0]='M';
  head->cmap[1]='A';
  head->cmap[2]='P';

  head->stamp[0]='D';
  return 0;
}


/*******************************************************************************************/

int mrc_read_head_MPI (MPI_File fin,  MrcHeader *head)
{
 // if(ftello64(fin)!=0)rewind(fin);

  MPI_File_seek( fin, 0, MPI_SEEK_SET ); 
  MPI_File_read( fin, head, 1024, MPI_CHAR, MPI_STATUS_IGNORE);
//  fread(head,1024,1,fin);

/*  if(!head||!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
  {
    printf("Error with function 'mrc_read_head()'! Can not Read the MrcHeader!");
    return -1;
  }*/
  return 0;
}

/*******************************************************************************************/
int mrc_read_head (FILE *fin,  MrcHeader *head)
{
  if(ftello64(fin)!=0)rewind(fin);
  fread(head,1024,1,fin);

  /*if(!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
  {
     printf("Error with function 'mrc_read_head()'! Warning: Not MRC format! \n");
     return -1;
  }*/

  return 0;
}


/*******************************************************************************************/
int mrc_write_head(MPI_File fout, MrcHeader *head)
{
 /* if(ftello64(fout)!=0)rewind(fout);
  MPI_offset *offset;
  MPI_File_get_position(fout, *offset);
  if*/ 

 /* if(!head||!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
  {
    printf("Error with function 'mrc_write_head()'! Can not write the MrcHeader!");
    return -1;
  }*/
  MPI_File_seek( fout, 0, MPI_SEEK_SET ); 
  MPI_File_write(fout, head, 1024, MPI_CHAR, MPI_STATUS_IGNORE);

  
 // fwrite(head,1024,1,fout);
  
  return 0;
}




/*******************************************************************************************/
int mrc_update_head(MPI_File fout)
{


  MrcHeader  *head;
  head=(MrcHeader *)malloc(sizeof(MrcHeader));

 /* FILE *finout;

  if((finout=fopen(inoutf,"r+"))==NULL)
  {
    printf("Cannot open file strike any key exit!\n");
  }*/

  mrc_read_head_MPI(fout,head);
  
  if(!head||!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
  {
    printf("Fatal error! The file is not a valid mrc file! in mrc_update_head\n");
    return -1;
  }

  long double sum;
  long double amin,amax,amean;
  int k,pNum;
  unsigned long site;
  unsigned char *p_uchar;
  short *p_short;
  float *p_float;

  //fseek(finout,1024,SEEK_SET);
  MPI_File_seek( fout, 1024, MPI_SEEK_SET ); 

  switch(head->mode)
  {      //switch start


    /**********case MRC_MODE_BYTE ***********/
    case MRC_MODE_BYTE :
    
    pNum=head->nx*head->ny;
   
    if((p_uchar=(unsigned char *)malloc(sizeof(unsigned char)*pNum))==NULL)
    {
       printf("Function 'malloc' erro, while updating head!\n");
       return -1;
    }

    printf("updating head!\n");
    MPI_File_read( fout, p_uchar, pNum, MPI_CHAR, MPI_STATUS_IGNORE);
    //fread(p_uchar,sizeof(unsigned char),pNum,finout);

    amin=amax=amean=p_uchar[0];
    sum=0;
    for(site=1;site<pNum;site++)
    {
       if(p_uchar[site]>amax)amax=p_uchar[site];
       if(p_uchar[site]<amin)amin=p_uchar[site];
       amean=amean+p_uchar[site];
    }
    amean/=pNum;
    sum=amean;
    for(k=1;k<head->nz;k++)
    {
       amean=0;
      // fread(p_uchar,sizeof(unsigned char),pNum,finout);
       MPI_File_read( fout, p_uchar, pNum, MPI_CHAR, MPI_STATUS_IGNORE);
       for(site=0;site<pNum;site++)
       {
          if(p_uchar[site]>amax)amax=p_uchar[site];
          if(p_uchar[site]<amin)amin=p_uchar[site];
          amean=amean+p_uchar[site];
       }
       amean/=pNum;
       sum+=amean;
     }
     amean=sum/head->nz;

     free(p_uchar);

    break;



    /**********case MRC_MODE_SHORT ***********/
    case MRC_MODE_SHORT :

    pNum=head->nx*head->ny;

    if((p_short=(short *)malloc(sizeof(short)*pNum))==NULL)
    {
       printf("Function 'malloc' error, while updating head!\n");

       return -1;
    }

    printf("updating head!\n");
    MPI_File_read( fout, p_short, pNum, MPI_SHORT, MPI_STATUS_IGNORE);
   // fread(p_short,sizeof(short),pNum,finout);
    amin=amax=amean=p_short[0];
    sum=0;
    for(site=1;site<pNum;site++)
    {
       if(p_short[site]>amax)amax=p_short[site];
       if(p_short[site]<amin)amin=p_short[site];
       amean=amean+p_short[site];
    }
    amean/=pNum;
    sum=amean;
    for(k=1;k<head->nz;k++)
    {
       amean=0;
       //fread(p_short,sizeof(short),pNum,finout);
       MPI_File_read( fout, p_short, pNum, MPI_SHORT, MPI_STATUS_IGNORE);
       for(site=0;site<pNum;site++)
       {
          if(p_short[site]>amax)amax=p_short[site];
          if(p_short[site]<amin)amin=p_short[site];
          amean=amean+p_short[site];
       }
       amean/=pNum;
       sum+=amean;
     }
     amean=sum/head->nz;

     free(p_short);

 
    break;


    /**********case MRC_MODE_FLOAT ***********/
    case MRC_MODE_FLOAT :

    pNum=head->nx*head->ny;
    if((p_float=(float *)malloc(sizeof(float)*pNum))==NULL)
    {

       printf("Function malloc erro!\n");
       return -1;
    }

    printf("updating head!\n");
    MPI_File_read( fout, p_float, pNum, MPI_FLOAT, MPI_STATUS_IGNORE);
    //fread(p_float,sizeof(float),pNum,finout);

    amin=amax=amean=p_float[0];
    sum=0;
    for(site=1;site<pNum;site++)
    {
       if(p_float[site]>amax)amax=p_float[site];
       if(p_float[site]<amin)amin=p_float[site];
       amean=amean+p_float[site];
    }
    amean/=pNum;
    sum=amean;
    for(k=1;k<head->nz;k++)
    {
       amean=0;
       //fread(p_float,sizeof(float),pNum,finout);
       MPI_File_read( fout, p_float, pNum, MPI_FLOAT, MPI_STATUS_IGNORE);
       for(site=0;site<pNum;site++)
       {
          if(p_float[site]>amax)amax=p_float[site];
          if(p_float[site]<amin)amin=p_float[site];
          amean=amean+p_float[site];
       }
       amean/=pNum;
       sum+=amean;
     }
     amean=sum/head->nz;

     free(p_float);

     break;
  }     //switch end
  
  


  head->amin=amin;
  head->amax=amax;
  head->amean=amean;

 // fclose(finout);

  printf("head->amin is %f, head->amax is %f, head->amean is %f\n",head->amin, head->amax, head->amean);
 // mrc_replace_head(inoutf,head);

  mrc_write_head(fout, head);
  free(head);
  printf("updating finished!\n");
  return 0;
}




/*******************************************************************************************/
int mrc_update_head_MPI(MPI_File fout, int id, int p)
{


  MrcHeader  *head;
  head=(MrcHeader *)malloc(sizeof(MrcHeader));


  mrc_read_head_MPI(fout,head);
  
  if(!head||!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
  {
    printf("Fatal error! The file is not a valid mrc file! in mrc_update_head\n");
    return -1;
  }


   int Z_start,Z_end,volZ_per;
   int volZ_add=head->nz%p;
   if(id<volZ_add)
   {
       volZ_per=head->nz/p+1;	   
	   Z_start=volZ_per*id;
	   Z_end=volZ_per*(id+1);
   }
   else
   {
       volZ_per=head->nz/p;	   
	   Z_start=volZ_per*id+volZ_add;	   
	   Z_end=volZ_per*(id+1)+volZ_add;
   }

   printf("Z_start is %d, volZ_per is %d in process %d\n",Z_start, volZ_per,id);
   
  int i,j,line_num;
  line_num=head->nx*head->ny;

  float *mrc_data_slice;
  if((mrc_data_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
   {
       printf("Can't malloc memery for 'mrc_data_slice'!");
       return FALSE;
   }
   memset(mrc_data_slice, 0 , sizeof(float)*line_num);

   
   if(id==0)
	   printf("updating head!\n");

   float amin,amin_all;
   float amax,amax_all;
   long double amean_slice, amean, amean_all;

//read the first slice
   mrc_read_slice(fout, head, Z_start, 'Z', mrc_data_slice);
   
   amin=amax=mrc_data_slice[0];
   amean_slice=mrc_data_slice[0];

   amin_all=0;
   amax_all=0;
   amean=0;
   amean_all=0;
   
   for(i=1;i<line_num;i++)
   {
	   if(mrc_data_slice[i]>amax) amax=mrc_data_slice[i];

	   if(mrc_data_slice[i]<amin) amin=mrc_data_slice[i];

       amean_slice=amean_slice+mrc_data_slice[i];
	}
   
   amean_slice/=line_num;
 //  printf("amean_slice is %Lf in slice %d \n",amean_slice, Z_start);

   amean+=amean_slice;

//read the rest slices
  for(j=1;j<volZ_per;j++)
  {
	  amean_slice=0;

	  mrc_read_slice(fout, head, (Z_start+j), 'Z', mrc_data_slice);

	 
	  for(i=0;i<line_num;i++)
	  {
		  if(mrc_data_slice[i]>amax) amax=mrc_data_slice[i];
	      if(mrc_data_slice[i]<amin) amin=mrc_data_slice[i];
          amean_slice=amean_slice+mrc_data_slice[i];
	  }
	  amean_slice/=line_num;
  // printf("amean_slice is %Lf in slice %d\n",amean_slice, Z_start+j);
	  

	  amean+=amean_slice;
  }


  amean/=volZ_per;


  printf("amin is %f, amax is %f, amean is %Lf in process %d\n",amin, amax, amean,id);


  MPI_Barrier(MPI_COMM_WORLD);
	

  MPI_Allreduce(&amin,&amin_all,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&amax,&amax_all,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&amean,&amean_all,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  amean_all/=p;

  head->amin=amin_all;
  head->amax=amax_all;
  head->amean=(float)amean_all;

  if(id==0)
  {
	  printf("head->amin is %f, head->amax is %f, head->amean is %f\n",head->amin, head->amax, head->amean);  
	  mrc_write_head(fout, head);
  }

  free(head);
  free(mrc_data_slice);
  if(id==0)
	  printf("updating finished!\n");

  return true;
}



/*******************************************************************************************/
//slice

/*******************************************************************************************/
int mrc_read_slice(MPI_File fin, MrcHeader  *head, int slcN, char axis,float *mrc_data_slice)
{

        int i,j,k,psize;

        unsigned char *data_byte;
        short *data_short;
        unsigned short *data_unshort;
	    float *data_float;		
		
		int Zslice_length;
		Zslice_length=head->nx*head->ny;

        MPI_Offset offset;
		offset=1024;
	    MPI_File_seek( fin, offset, MPI_SEEK_SET ); 


        switch(head->mode)
        {
          case MRC_MODE_BYTE :
          psize=sizeof(unsigned char);

          break;

          case MRC_MODE_SHORT :
          case MRC_MODE_USHORT:
          psize=sizeof(short);

          break;

          case MRC_MODE_FLOAT :
          psize=sizeof(float);

          break;
         }

         switch(axis)
         {

/***********************************X************************************/
          case 'x':
          case 'X': 
  
         switch(head->mode)
		 {
			 case MRC_MODE_BYTE:
				 if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_byte'!");
					 return FALSE;
				 }
				 memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 mrc_data_slice[j+k*head->ny]=(float)data_byte[slcN+j*head->nx];
				 }

				 free(data_byte);
				 
			 break;
			 
			 case MRC_MODE_SHORT:
				 if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_short'!");
					 return FALSE;
				 }
				 memset(data_short, 0 , sizeof(short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 mrc_data_slice[j+k*head->ny]=(float)data_short[slcN+j*head->nx];
				 }
				 free(data_short);

			 break;

			 case MRC_MODE_USHORT:
				 if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_unshort'!");
					 return FALSE;
				 }
				 memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 mrc_data_slice[j+k*head->ny]=(float)data_unshort[slcN+j*head->nx];
				 }
				 free(data_unshort);

			 break;

			 			 

			 
			 case MRC_MODE_FLOAT:
				 if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_float'!");
					 return FALSE;
				 }
				 memset(data_float, 0 , sizeof(float)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_float,Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 mrc_data_slice[j+k*head->ny]=data_float[slcN+j*head->nx];
				 }
				 free(data_float);

			 break;

		 }
		 
		 break;

/***********************************Y************************************/
        case 'y':
        case 'Y':
  	   
	   		
		switch(head->mode)
        {
			case MRC_MODE_BYTE:
				if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
				{
					printf("Error with Function 'mrc_read_slice() in axis-Y'!Can't malloc memery for 'data_byte'!");
					return FALSE;
				}
				memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);
				
				for(k=0;k<head->nz;k++)
				{
					MPI_File_read( fin, data_byte, Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
			 
					for(i=0;i<head->nx;i++)
					  mrc_data_slice[i+k*head->nx]=(float)data_byte[i+slcN*head->nx];
				}
				free(data_byte);

			break;

			case MRC_MODE_SHORT:
				if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
				{
					printf("Error with Function 'mrc_read_slice() in axis-Y'!Can't malloc memery for 'data_short'!");
					return FALSE;
				}
				memset(data_short, 0 , sizeof(short)*Zslice_length);
				
				for(k=0;k<head->nz;k++)
				{
					MPI_File_read( fin, data_short, Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
			 
					for(i=0;i<head->nx;i++)
					  mrc_data_slice[i+k*head->nx]=(float)data_short[i+slcN*head->nx];
				}
				free(data_short);

			break;

			case MRC_MODE_USHORT:
				if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
				{
					printf("Error with Function 'mrc_read_slice() in axis-Y'!Can't malloc memery for 'data_unshort'!");
					return FALSE;
				}
				memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);
				
				for(k=0;k<head->nz;k++)
				{
					MPI_File_read( fin, data_unshort, Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
			 
					for(i=0;i<head->nx;i++)
					  mrc_data_slice[i+k*head->nx]=(float)data_unshort[i+slcN*head->nx];
				}
				free(data_unshort);

			break;

			
			case MRC_MODE_FLOAT:
				if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
				{
					printf("Error with Function 'mrc_read_slice() in axis-Y'!Can't malloc memery for 'data_float'!");
					return FALSE;
				}
				memset(data_float, 0 , sizeof(float)*Zslice_length);
				
				for(k=0;k<head->nz;k++)
				{
					MPI_File_read( fin, data_float, Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
			 
					for(i=0;i<head->nx;i++)
					  mrc_data_slice[i+k*head->nx]=(float)data_float[i+slcN*head->nx];
				}
				free(data_float);

			break;

        }//end switch
		
		break;

/***********************************Z************************************/
        case 'z':
        case 'Z': 
  
        offset = head->nx*head->ny*psize;

        for (k=0;k<slcN;k++)
          MPI_File_seek( fin, offset, MPI_SEEK_CUR );


        switch(head->mode)
        {
          case MRC_MODE_BYTE:     

             if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
              {
                printf("Error with Function 'mrc_read_slice()'!Can't malloc memery for 'data_byte'!");
                return FALSE;
              }
             memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);
             MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
             
             for(i=0;i<Zslice_length;i++)
                  mrc_data_slice[i]=(float)data_byte[i];
             free(data_byte);
         
          break;

          case MRC_MODE_SHORT:

             if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
              {
                 printf("Error with Function 'mrc_read_slice()'!Can't malloc memery for 'data_short'!");
                 return FALSE;
              }
             memset(data_short, 0 , sizeof(short)*Zslice_length);
            
			 MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
             
             for(i=0;i<Zslice_length;i++)
                 mrc_data_slice[i]=(float)data_short[i];
             free(data_short);

		  break;

          case MRC_MODE_USHORT:  
             
             if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
              {
                 printf("Error with Function 'mrc_read_slice()'!Can't malloc memery for 'data_unshort'!");
                 return FALSE;
              }
             memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);
			 
			 MPI_File_read( fin, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
             
			  for(i=0;i<Zslice_length;i++)
                 mrc_data_slice[i]=(float)data_unshort[i];
             
			 free(data_unshort);
          break;

          case MRC_MODE_FLOAT:
               MPI_File_read( fin, mrc_data_slice,Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
          break;
         }

    break;

    }


	return true;
}



/*******************************************************************************************/
int mrc_write_slice(MPI_File fout, MrcHeader  *head, int slcN,char axis,float *slcdata) //only Z
{
  int psize;
  if (head->mode==MRC_MODE_FLOAT)
  psize=sizeof(float);
  else {
        printf ("outfile headmode is error!\n");
        return false;
        }

  int i,j,k;
  MPI_Offset offset;
  offset = 1024;
  MPI_File_seek( fout, offset, MPI_SEEK_SET ); 
	  

  int Zslice_length;
  Zslice_length=head->nx*head->ny;

 /* float *Zslice_data;
  if((Zslice_data=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
   {
       printf("Can't malloc memery for 'Zslice_data'!");
       return FALSE;
   }
   memset(Zslice_data, 0 , sizeof(float)*length);*/

  switch(axis)
  {

/***********************************X************************************/
    case 'x':
    case 'X':

	  offset = slcN*psize;
	  MPI_File_seek( fout, offset, MPI_SEEK_CUR);

	  for(k=0;k<head->nz;k++)
	  {
		  for(j=0;j<head->ny;j++)
		  {
		  	MPI_File_write(fout, slcdata+k*head->ny+j, 1, MPI_FLOAT, MPI_STATUS_IGNORE);

			offset=(head->nx-1)*psize;
			MPI_File_seek(fout,offset,MPI_SEEK_CUR);
		  }

	  }
	break;

	case 'y':
	case 'Y':

	  offset = head->nx*slcN*psize;
	  MPI_File_seek( fout, offset, MPI_SEEK_CUR);
	  
	  for(k=0;k<head->nz;k++)
	  {
		 MPI_File_write(fout, slcdata+k*head->nx, head->nx, MPI_FLOAT, MPI_STATUS_IGNORE);
		 
		 offset=head->nx*(head->ny-1)*psize;
		 MPI_File_seek(fout,offset,MPI_SEEK_CUR);
	  }
	break;

	case 'z':
	case 'Z':
	  
      offset=head->nx*head->ny*psize;
	  
	  for(k=0;k<slcN;k++)
		  MPI_File_seek( fout, offset, MPI_SEEK_CUR);

      MPI_File_write(fout, slcdata, Zslice_length, MPI_FLOAT, MPI_STATUS_IGNORE);
	  
    break;

    }

 // free(Zslice_data);
  return true;
}



/*******************************************************************************************/
int mrc_add_sliceN(FILE *fout, MrcHeader  *headout, float *slcdata, int slcN)
{
  off_t length;
  //fseeko(fout,1024+sizeof(float)*headout->nx*headout->ny*slcN,SEEK_SET);
  length=headout->nx*headout->ny;
  length*=(sizeof(float)*slcN);
  length+=1024;
  fseeko(fout,length,SEEK_SET);

  fwrite(slcdata,sizeof(float),headout->nx*headout->ny,fout);
  return 0;
}



  
/*******************************************************************************************/
int mrc_add_slice(FILE *fout , MrcHeader  *headout, float *slcdata)
{

  fseeko(fout,0,SEEK_END);
  fwrite(slcdata,sizeof(float),headout->nx*headout->ny,fout);
  return 0;
}


/*****************************************************************************************************/
//pixel

/*****************************************************************************************************/

int mrc_read_pixel( MPI_File fin, MrcHeader  *head, int x, int y, int z ,float *pixel_gray)

{
	unsigned char *data_byte;
    short *data_short;
    unsigned short *data_unshort;
	float *data_float;

	int Zslice_length;
	Zslice_length=head->nx*head->ny;

    MPI_Offset  offset;
	offset=1024;
	MPI_File_seek( fin, offset, MPI_SEEK_SET);

	int psize;
	switch(head->mode)
    {
		case MRC_MODE_BYTE :
        psize=sizeof(unsigned char);

        break;

        case MRC_MODE_SHORT :
        case MRC_MODE_USHORT:
        psize=sizeof(short);

        break;

        case MRC_MODE_FLOAT :
        psize=sizeof(float);

        break;
    }

	offset=head->nx*head->ny*psize;

	int k;
	
	for(k=0;k<z;k++)
		MPI_File_seek( fin, offset, MPI_SEEK_CUR);


    switch(head->mode)
    
    {
      case MRC_MODE_BYTE:
		  if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_read_pixel()'!Can't malloc memory for 'data_byte'!");
			  return FALSE;
		  }
		  memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);

          MPI_File_read( fin, data_byte, Zslice_length, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

		  pixel_gray[0]=(float)data_byte[x+y*head->nx];

		  free(data_byte);

	  break;

      case MRC_MODE_SHORT:
		  if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_read_pixel()'!Can't malloc memory for 'data_short'!");
			  return FALSE;
		  }
		  memset(data_short, 0 , sizeof(short)*Zslice_length);

          MPI_File_read( fin, data_short, Zslice_length, MPI_SHORT, MPI_STATUS_IGNORE);

		  pixel_gray[0]=(float)data_short[x+y*head->nx];

		  free(data_short);

	  break;

	  case MRC_MODE_USHORT:
		  if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_read_pixel()'!Can't malloc memory for 'data_unshort'!");
			  return FALSE;
		  }
		  memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

          MPI_File_read( fin, data_unshort, Zslice_length, MPI_UNSIGNED_SHORT, MPI_STATUS_IGNORE);

		  pixel_gray[0]=(float)data_unshort[x+y*head->nx];

		  free(data_unshort);

	  break;
	  
	  case MRC_MODE_FLOAT:
		  if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_read_pixel()'!Can't malloc memory for 'data_float'!");
			  return FALSE;
		  }
		  memset(data_float, 0 , sizeof(float)*Zslice_length);

          MPI_File_read( fin, data_float, Zslice_length, MPI_FLOAT, MPI_STATUS_IGNORE);

		  pixel_gray[0]=data_float[x+y*head->nx];

		  free(data_float);

	  break;


    }
return true;
}



/*****************************************************************************************************/

int mrc_write_pixel( MPI_File fout, MrcHeader  *head, int x, int y, int z ,float *pixel_gray) //serial

{
	char *data_byte;
    short *data_short;
	unsigned short *data_unshort;
	float *data_float;
	
	int Zslice_length;
	Zslice_length=head->nx*head->ny;

    MPI_Offset  offset;
	offset=1024;
	MPI_File_seek( fout, offset, MPI_SEEK_SET);

	int psize;
	switch(head->mode)
    {
		case MRC_MODE_BYTE :
        psize=sizeof(unsigned char);

        break;

        case MRC_MODE_SHORT :
        case MRC_MODE_USHORT:
        psize=sizeof(short);

        break;
        case MRC_MODE_FLOAT :
        psize=sizeof(float);

        break;
    }

	offset=head->nx*head->ny*psize;

	int k;
	
	for(k=0;k<z;k++)
		MPI_File_seek( fout, offset, MPI_SEEK_CUR);


    switch(head->mode)
    
    {
      case MRC_MODE_BYTE:

		  if((data_byte=( char *)malloc(sizeof( char)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_write_pixel()'!Can't malloc memory for 'data_byte'!");
			  return FALSE;
		  }
		  memset(data_byte, 0 , sizeof( char)*Zslice_length);

          MPI_File_read( fout, data_byte, Zslice_length, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
		  
		  data_byte[x+y*head->nx]=( char)pixel_gray[0];

		  MPI_File_seek(fout, -1*offset, MPI_SEEK_CUR);

		  MPI_File_write( fout, data_byte,Zslice_length, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

		  free(data_byte);

      break;

      case MRC_MODE_SHORT:
	      if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_write_pixel()'!Can't malloc memory for 'data_short'!");
			  return FALSE;
		  }
		  memset(data_short, 0 , sizeof(short)*Zslice_length);

          MPI_File_read( fout, data_short, Zslice_length, MPI_SHORT, MPI_STATUS_IGNORE);
		  
		  data_short[x+y*head->nx]=(short)pixel_gray[0];

		  MPI_File_seek(fout, -1*offset, MPI_SEEK_CUR);

		  MPI_File_write( fout, data_short,Zslice_length, MPI_SHORT, MPI_STATUS_IGNORE);

		  free(data_short);


	  break;

	  case MRC_MODE_USHORT:
	      if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_write_pixel()'!Can't malloc memory for 'data_unshort'!");
			  return FALSE;
		  }
		  memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

          MPI_File_read( fout, data_unshort, Zslice_length, MPI_UNSIGNED_SHORT, MPI_STATUS_IGNORE);
		  
		  data_unshort[x+y*head->nx]=(short)pixel_gray[0];

		  MPI_File_seek(fout, -1*offset, MPI_SEEK_CUR);		  

		  MPI_File_write( fout, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT, MPI_STATUS_IGNORE);

		  free(data_unshort);

	  break;


      case MRC_MODE_FLOAT:
	      if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
		  {
			  printf("Error with Function 'mrc_write_pixel()'!Can't malloc memory for 'data_float'!");
			  return FALSE;
		  }
		  memset(data_float, 0 , sizeof(float)*Zslice_length);

          MPI_File_read( fout, data_float, Zslice_length, MPI_FLOAT, MPI_STATUS_IGNORE);
		  
		  data_float[x+y*head->nx]=pixel_gray[0];

		  MPI_File_seek(fout, -1*offset, MPI_SEEK_CUR);		  

		  MPI_File_write( fout, data_float,Zslice_length, MPI_FLOAT, MPI_STATUS_IGNORE);

		  free(data_float);

	          
	  break;

    }
return 0;
}



/*****************************************************************************************************/

int mrc_write_pixel_mod( MPI_File fout, MrcHeader  *head, int x, int y, int z ,float *pixel_gray) //parallel

{
    MPI_Offset  offset;
	int k;

	char *pixel_byte;
    pixel_byte=(char *)malloc(sizeof(char));
	
    short *pixel_short;
	pixel_short=(short *)malloc(sizeof(short));
	

	unsigned short *pixel_unshort;
	pixel_unshort=(unsigned short *)malloc(sizeof(unsigned short));
	

	int psize;
	switch(head->mode)
    {
		case MRC_MODE_BYTE :
        psize=sizeof(unsigned char);

        break;

        case MRC_MODE_SHORT :
        case MRC_MODE_USHORT:
        psize=sizeof(short);

        break;

        case MRC_MODE_FLOAT :
        psize=sizeof(float);

        break;
    }

	offset=1024+(x+y*head->nx)*psize;
	MPI_File_seek( fout, offset, MPI_SEEK_SET );

	offset=head->nx*head->ny*psize;
	for (k=0;k<z;k++)
		MPI_File_seek( fout, offset, MPI_SEEK_CUR);

    switch(head->mode)
    
    {
      case MRC_MODE_BYTE:
		  
		  pixel_byte[0]=(char)pixel_gray[0];
		   
		  MPI_File_write( fout, pixel_byte,1, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);

      break;

      case MRC_MODE_SHORT:
	  
		  pixel_short[0]=(short)pixel_gray[0];

	      MPI_File_write( fout, pixel_short,1, MPI_SHORT,MPI_STATUS_IGNORE);
           
	  break;

	  case MRC_MODE_USHORT:
	  
		  pixel_unshort[0]=(unsigned short)pixel_gray[0];

	      MPI_File_write( fout, pixel_unshort,1, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
           
	  break;


      case MRC_MODE_FLOAT:

		  MPI_File_write( fout, pixel_gray,1, MPI_FLOAT,MPI_STATUS_IGNORE);

     
	  break;

    }

	free(pixel_byte);
	free(pixel_short);
	free(pixel_unshort);

return 0;
}


/*****************************************************************************************************/
//block

/*****************************************************************************************************/

int mrc_read_block(MPI_File fin, MrcHeader  *head, int start, int end, char axis,float *mrc_data_block)
{

    
        int i,j,k,psize;

        unsigned char *data_byte;
        short *data_short;
        unsigned short *data_unshort;
        int *data_int;
        float *data_float;

        MPI_Offset offset;
		
		offset=1024;
	    MPI_File_seek( fin, offset, MPI_SEEK_SET); 

        int Zslice_length;
		Zslice_length=head->nx*head->ny;


        int length;
        length=end-start;


        switch(head->mode)
        {
          case MRC_MODE_BYTE :
          psize=sizeof(unsigned char);

          break;

          case MRC_MODE_SHORT :
          case MRC_MODE_USHORT:
          psize=sizeof(short);

          break;

          case MRC_MODE_FLOAT :
          psize=sizeof(float);

          break;
         }

         switch(axis)
         {

/***********************************X************************************/
          case 'x':
          case 'X':
 
		  switch(head->mode)
		  {
			 
			 case MRC_MODE_BYTE:
				 if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_byte'!");
					 return FALSE;
				 }
				 memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 for(i=start;i<end;i++)
							 mrc_data_block[i-start+j*length+k*length*head->ny]=(float)data_byte[i+j*head->nx];
				 }

				 free(data_byte);
				 
			 break;
			 
			 case MRC_MODE_SHORT:
				 if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_short'!");
					 return FALSE;
				 }
				 memset(data_short, 0 , sizeof(short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 for(i=start;i<end;i++)
							 mrc_data_block[i-start+j*length+k*length*head->ny]=(float)data_short[i+j*head->nx];
				 }
				 free(data_short);

			 break;

			 case MRC_MODE_USHORT:
				 if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_unshort'!");
					 return FALSE;
				 }
				 memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 for(i=start;i<end;i++)
							 mrc_data_block[i-start+j*length+k*length*head->ny]=(float)data_unshort[i+j*head->nx];
				 }
				 free(data_unshort);

			 break;

			 
			 case MRC_MODE_FLOAT:
				 if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_float'!");
					 return FALSE;
				 }
				 memset(data_float, 0 , sizeof(float)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_float,Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
					 for(j=0;j<head->ny;j++)
						 for(i=start;i<end;i++)
							 mrc_data_block[i-start+j*length+k*length*head->ny]=data_float[i+j*head->nx];
				 }
				 free(data_float);

			 break;

		 }
		 
		 break;

/***********************************Y************************************/
        case 'y':
        case 'Y':
		switch(head->mode)
		  {
			 
			 case MRC_MODE_BYTE:
				 if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_byte'!");
					 return FALSE;
				 }
				 memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
					 for(j=start;j<end;j++)
						 for(i=0;i<head->nx;i++)
							 mrc_data_block[i+(j-start)*head->nx+k*length*head->nx]=(float)data_byte[i+j*head->nx];
				 }

				 free(data_byte);
				 
			 break;
			 
			 case MRC_MODE_SHORT:
				 if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_short'!");
					 return FALSE;
				 }
				 memset(data_short, 0 , sizeof(short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
					 for(j=start;j<end;j++)
						 for(i=0;i<head->nx;i++)
							 mrc_data_block[i+(j-start)*head->nx+k*length*head->nx]=(float)data_short[i+j*head->nx];
				 }
				 free(data_short);

			 break;

			 case MRC_MODE_USHORT:
				 if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_unshort'!");
					 return FALSE;
				 }
				 memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
					 for(j=start;j<end;j++)
						 for(i=0;i<head->nx;i++)
							 mrc_data_block[i+(j-start)*head->nx+k*length*head->nx]=(float)data_unshort[i+j*head->nx];
				 }
				 free(data_unshort);

			 break;

			 
			 case MRC_MODE_FLOAT:
				 if((data_float=(float *)malloc(sizeof(float)*Zslice_length))==NULL)
                 {
					 printf("Error with Function 'mrc_read_slice() in axis-X'!Can't malloc memery for 'data_float'!");
					 return FALSE;
				 }
				 memset(data_float, 0 , sizeof(float)*Zslice_length);

				 for(k=0;k<head->nz;k++)
				 {
					 MPI_File_read( fin, data_float,Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
					 for(j=start;j<end;j++)
						 for(i=0;i<head->nx;i++)
							 mrc_data_block[i+(j-start)*head->nx+k*length*head->nx]=data_float[i+j*head->nx];
				 }
				 free(data_float);

			 break;

		 }

          break;

/***********************************Z************************************/
        case 'z':
        case 'Z': 
		  
		offset = head->nx*head->ny*psize;

        for (k=0;k<start;k++)
			MPI_File_seek( fin, offset, MPI_SEEK_CUR );

        switch(head->mode)
        {
          case MRC_MODE_BYTE:     

             if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
              {
                printf("Error with Function 'mrc_read_block()'!Can't malloc memery for 'data_byte'!");
                return FALSE;
              }
             memset(data_byte, 0 , sizeof(unsigned char)*(head->nx*head->ny));

			 for(k=start;k<end;k++)
			 {
				 MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
				 for(j=0;j<head->ny;j++)
				   for(i=0;i<head->nx;i++)
					 mrc_data_block[i+j*head->nx+(k-start)*Zslice_length]=(float)data_byte[i+j*head->nx];
				}

             free(data_byte);
         
          break;

          case MRC_MODE_SHORT:

             if((data_short=(short *)malloc(sizeof(short)*(head->nx*head->ny)))==NULL)
              {
                 printf("Error with Function 'mrc_read_block()'!Can't malloc memery for 'data_short'!");
                 return FALSE;
              }
             memset(data_short, 0 , sizeof(short)*(head->nx*head->ny));
          
			 for(k=start;k<end;k++)
             {
                 MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
				 for(j=0;j<head->ny;j++)
				   for(i=0;i<head->nx;i++)
					 mrc_data_block[i+j*head->nx+(k-start)*Zslice_length]=(float)data_short[i+j*head->nx];
             
             }
             free(data_short);
		  break;

          case MRC_MODE_USHORT:  
             
             if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*(head->nx*head->ny)))==NULL)
              {
                 printf("Error with Function 'mrc_read_block()'!Can't malloc memery for 'data_unshort'!");
                 return FALSE;
              }
             memset(data_unshort, 0 , sizeof(unsigned short)*(head->nx*head->ny));
			 int k;
			 for(k=start;k<end;k++)
             {
                 MPI_File_read( fin, data_unshort,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
                 for(j=0;j<head->ny;j++)
				   for(i=0;i<head->nx;i++)
					 mrc_data_block[i+j*head->nx+(k-start)*Zslice_length]=(float)data_unshort[i+j*head->nx];

                
			 }
			 free(data_unshort);
           break;

           case MRC_MODE_FLOAT:
               MPI_File_read( fin, mrc_data_block,length*Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
           break;
         }

    break;

    }


	return 0;
}
	

/*****************************************************************************************************/


int mrc_write_block(MPI_File fout, MrcHeader  *head, int start, int end, char axis,float *mrc_data_block)
{

    
        int i,k,j,psize;

		psize=sizeof(float);

        int length;
		int del=end-start;
       
        MPI_Offset offset;
 
		offset=1024;
		MPI_File_seek( fout, offset, MPI_SEEK_SET);

       
         switch(axis)
         {

/***********************************X************************************/
          case 'x':
          case 'X':
			  
			  length=del*head->ny*head->nz;

			  offset=start*psize;
			  MPI_File_seek( fout, offset, MPI_SEEK_CUR);
			  
			  for(k=0;k<head->nz;k++)
				for(j=0;j<head->ny;j++)
				{
					  MPI_File_write( fout, mrc_data_block+j*del+k*del*head->ny, del, MPI_FLOAT, MPI_STATUS_IGNORE);
							
					  offset=(head->nx-del)*psize;
					  MPI_File_seek( fout, offset, MPI_SEEK_CUR);
				}
				
		break;

/***********************************Y************************************/
        case 'y':
        case 'Y':
		
		      length=del*head->nx*head->nz;

			  offset=head->nx*start*psize;
			  MPI_File_seek( fout, offset, MPI_SEEK_CUR);
			 
			  for(k=0;k<head->nz;k++)
			  {
					  MPI_File_write( fout, mrc_data_block+k*del*head->nx, del*head->nx, MPI_FLOAT, MPI_STATUS_IGNORE);
							  
					  offset=head->nx*(head->ny-del)*psize;
					  MPI_File_seek( fout, offset, MPI_SEEK_CUR);
			  }

        break;

/***********************************Z************************************/
        case 'z':
        case 'Z': 
       
		length=head->nx*head->ny*(end-start);

        offset = head->nx*head->ny*psize;

        for (i=0;i<start;i++)
          MPI_File_seek( fout, offset, MPI_SEEK_CUR );


		for (k=start;k<end;k++)
		   {
               MPI_File_write( fout, mrc_data_block+(k-start)*head->nx*head->ny,(head->nx*head->ny), MPI_FLOAT,MPI_STATUS_IGNORE);
		   }

    break;

    }


	return 0;
}
	



/*****************************************************************************************************/


int mrc_read_all_MPI(MPI_File fin, MrcHeader  *head, float *mrc_data_all)
{

        MPI_Offset filesize;
        MPI_File_get_size (fin, &filesize);
        //check the mrc file to make sure the size is exact in register with the head
	switch(head->mode)
	{
		case MRC_MODE_BYTE:
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(char))
			{
				printf("Error with Function 'mrc_read_all()'! File size error!");
			}
			break;

		case MRC_MODE_SHORT:
		case MRC_MODE_USHORT:
                        printf("hello4\n");
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(short))
			{
				printf("Error with Function 'mrc_read_all()'! File size error!");
			}
			break;

		case MRC_MODE_FLOAT:
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(float))
			{
				printf("Error with Function 'mrc_read_all()'! File size error!");
			}
			break;
			
		default:
			printf("Error with Function 'mrc_read_all()'! File type unknown!");

			break;
	}
 


	int i,k;
	unsigned char *data_byte;
	short *data_short;
	unsigned short *data_unshort;

    int Zslice_length;
    Zslice_length=head->nx*head->ny;
        
	MPI_File_seek( fin, 1024, MPI_SEEK_SET ); 
	
	switch(head->mode)
	{
		case MRC_MODE_BYTE:
			if((data_byte=(unsigned char *)malloc(sizeof(unsigned char)*Zslice_length))==NULL)
              {
                printf("Error with Function 'mrc_read_all_MPI()'!Can't malloc memery for 'data_byte'!");
                return FALSE;
              }
             memset(data_byte, 0 , sizeof(unsigned char)*Zslice_length);


			for(k=0;k<head->nz;k++)
            {          
               MPI_File_read( fin, data_byte,Zslice_length, MPI_UNSIGNED_CHAR,MPI_STATUS_IGNORE);
			   for(i=0;i<Zslice_length;i++)
				   mrc_data_all[i+k*Zslice_length]=(float)data_byte[i];
            }
        
			free(data_byte);

	        break;
		
		case MRC_MODE_SHORT:
			if((data_short=(short *)malloc(sizeof(short)*Zslice_length))==NULL)
              {
                printf("Error with Function 'mrc_read_all_MPI()'!Can't malloc memery for 'data_short'!");
                return FALSE;
              }
             memset(data_short, 0 , sizeof(short)*Zslice_length);

			for(k=0;k<head->nz;k++)
            {          
               MPI_File_read( fin, data_short,Zslice_length, MPI_SHORT,MPI_STATUS_IGNORE);
			   for(i=0;i<Zslice_length;i++)
				   mrc_data_all[i+k*Zslice_length]=(float)data_short[i];
            }

			free(data_short);
			break;

		case MRC_MODE_USHORT:
			if((data_unshort=(unsigned short *)malloc(sizeof(unsigned short)*Zslice_length))==NULL)
              {
                printf("Error with Function 'mrc_read_all_MPI()'!Can't malloc memery for 'data_unshort'!");
                return FALSE;
              }
             memset(data_unshort, 0 , sizeof(unsigned short)*Zslice_length);

			for(k=0;k<head->nz;k++)
            {          
               MPI_File_read( fin, data_unshort,Zslice_length, MPI_UNSIGNED_SHORT,MPI_STATUS_IGNORE);
			   for(i=0;i<Zslice_length;i++)
				   mrc_data_all[i+k*Zslice_length]=(float)data_unshort[i];
            }

			free(data_unshort);
			break;

		case MRC_MODE_FLOAT:

			for(k=0;k<head->nz;k++)         
               MPI_File_read( fin, mrc_data_all+k*Zslice_length,Zslice_length, MPI_FLOAT,MPI_STATUS_IGNORE);
			  
	
			break;

	}

	return 0;
}

/*****************************************************************************************************/


int mrc_read_all(FILE *fin, MrcHeader  *head, float *mrc_data_all)
{


	int headsize = get_file_size(fin) - (long long) head->nx * head->ny * head->nz * sizeof(short);
//check the mrc file to make sure the size is exact in register with the head


	long i;
	unsigned char buf_byte;
	short buf_short;
	short buf_ushort;

	fseek(fin,(headsize),SEEK_SET);

	switch(head->mode)
	{
		case MRC_MODE_BYTE:


			for(i=0;i<head->nx*head->ny*head->nz;i++)
			{
				fread(&buf_byte,sizeof(char),1,fin);
				mrc_data_all[i] = (float)buf_byte;
			}
			break;
		
		case MRC_MODE_SHORT:

			for(i=0;i<head->nx*head->ny*head->nz;i++)
			{
				fread(&buf_short,sizeof(short),1,fin);
				mrc_data_all[i] = (float)buf_short;
			}
			break;

		case MRC_MODE_USHORT:

			for(i=0;i<head->nx*head->ny*head->nz;i++)
			{
				fread(&buf_ushort,sizeof(short),1,fin);
				mrc_data_all[i] = (float)buf_ushort;
			}
			break;

		case MRC_MODE_FLOAT:

			if((fread(mrc_data_all,head->nx*head->ny*head->nz*sizeof(float),1,fin)==0))
			{
				printf("Error with Function 'mrc_read_all()'! Reading file failed!");
				return -1;
			}
			break;

		default:
			printf("Error with Function 'mrc_read_all()'! File type unknown!");
			break;

	}

	return 0;
}
	


/*****************************************************************************************************/
int mrc_write_all(MPI_File fout, MrcHeader  *head,  int Z_start,int Z_end, float *mrc_data_all)
{
    MPI_Offset offset;
    offset=head->nx*head->ny*Z_start*sizeof(float);
    offset+=1024;
    MPI_File_seek( fout, offset, MPI_SEEK_SET ); 
    MPI_File_write(fout, mrc_data_all, head->nx*head->ny*(Z_end-Z_start), MPI_FLOAT, MPI_STATUS_IGNORE);

return true;
}


/*****************************************************************************************************/
int mrc_flipYZ_block(MPI_File fin, MPI_File fout, MrcHeader *inhead, MrcHeader *outhead, int start, int end)
{
  printf("\nBegin flipping:");
  
 /* outhead->nx=inhead->nx;
  outhead->ny=inhead->nz;
  outhead->nz=inhead->ny;
  outhead->mode=MRC_MODE_FLOAT;

  mrc_write_head(fout,outhead);*/

  float *buf;
  buf=(float *)malloc(sizeof(float)*inhead->nx*inhead->nz);
  int j;
  for(j=start;j<end;j++)
  {
    mrc_read_slice(fin,inhead,j,'y',buf);
	mrc_write_slice(fout, outhead, j,'Z',buf);
  }
  
  free (buf);
  printf("\nflipping finished!\n");
  return 0;
}


