#include "mrcfile_atom.h"
#include "mpi.h"

/*******************************************************************************************/
long get_file_size(FILE *fin)
{
	fseek(fin,0,SEEK_END);

	return ftell(fin);

}


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


/*****************************************************************************************************/


int mrc_read_block(MPI_File fin, MrcHeader  *head, int start, int end, char axis,float *mrc_data_block)
{

        MPI_Offset filesize;
        MPI_File_get_size (fin, &filesize);
//check the mrc file to make sure the size is exact in register with the head
	switch(head->mode)
	{
		case MRC_MODE_BYTE:
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(char))
			{
				printf("Error with Function 'mrc_read_block()'! File size error!");
			}
			break;

		case MRC_MODE_SHORT:
		case MRC_MODE_USHORT:
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(short))
			{
				printf("Error with Function 'mrc_read_block()'! File size error!");
			}
			break;

		case MRC_MODE_FLOAT:
			if(filesize - 1024 != head->nx*head->ny*head->nz*sizeof(float))
			{
				printf("Error with Function 'mrc_read_block()'! File size error!");
			}
			break;
			
		default:
			printf("Error with Function 'mrc_read_block()'! File type unknown!");

			break;
	}


	int i;
	unsigned char buf_byte;
	short buf_short;
	short buf_ushort;

       int psize;
  int length;
  length=head->nx*head->ny*(end-start);

  MPI_Offset offset;
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

      /*fseeko(fin,1024+slcN*psize,SEEK_SET);

      switch(head->mode)
      {
        case MRC_MODE_BYTE:
          for(i=0;i<head->ny*head->nz;i++)
              {
                fread(&buf_byte,psize,1,fin);
                slcdata[i]=(float)buf_byte;
                fseeko(fin,(head->nx-1)*psize,SEEK_CUR);
              }

          break;

        case MRC_MODE_SHORT:
          for(i=0;i<head->ny*head->nz;i++)
             {
                fread(&buf_short,psize,1,fin);
                slcdata[i]=(float)(buf_short);
                fseeko(fin,(head->nx-1)*psize,SEEK_CUR);
              }

          break;

        case MRC_MODE_FLOAT:
          for(i=0;i<head->ny*head->nz;i++)
            {
            fread(&buf_float,psize,1,fin);
            slcdata[i]=buf_float;
            fseeko(fin,(head->nx-1)*psize,SEEK_CUR);
            }
          break;

       }*/

    break;

/***********************************Y************************************/
    case 'y':
    case 'Y':
   //   fseeko(fin,1024+(k*head->nx*head->ny+head->nx*slcN)*psize,SEEK_SET);
     /* fseeko(fin,1024+slcN*head->nx*psize,SEEK_SET);
      for(k=0;k<head->nz;k++)
      {

      switch(head->mode)
      {
        case MRC_MODE_BYTE:
        for(i=0;i<head->nx;i++)
              {
                fread(&buf_byte,psize,1,fin);
                slcdata[k*head->nx+i]=(float)buf_byte;
              }

          break;

        case MRC_MODE_SHORT:
        for(i=0;i<head->nx;i++)
             {
                fread(&buf_short,psize,1,fin);
                slcdata[k*head->nx+i]=(float)(buf_short);
              }

          break;

        case MRC_MODE_FLOAT:
        fread(slcdata+k*head->nx,psize,head->nx,fin);


          break;

      }//end switch
      fseeko(fin,head->nx*(head->ny-1)*psize,SEEK_CUR);
      }//end for*/
    break;

/***********************************Z************************************/
    case 'z':
    case 'Z': 
      offset=head->nx*head->ny;
      offset*=(start*psize);
      offset+=(1024);

      MPI_File_seek( fin, offset, MPI_SEEK_SET ); 

     // fseeko(fin, offset, SEEK_SET);
      //fseeko(fin,1024+slcN*head->nx*head->ny*psize,SEEK_SET);

     switch(head->mode)
      {
        case MRC_MODE_BYTE:     
        for(i=0;i<length;i++)
        {
          // MPI_File_read( fin, head, 1024, MPI_CHAR, MPI_STATUS_IGNORE);
          
          MPI_File_read( fin, &buf_byte,1, MPI_CHAR,MPI_STATUS_IGNORE);
        //  fread(&buf_byte,psize,1,fin);
          mrc_data_block[i]=(float)buf_byte;
         }
        break;

        case MRC_MODE_SHORT:
        case MRC_MODE_USHORT:  
        for(i=0;i<length;i++)
        {
          //fread(&buf_short,psize,1,fin);
          MPI_File_read( fin, &buf_short,1, MPI_SHORT,MPI_STATUS_IGNORE);
          mrc_data_block[i]=(float)buf_short;
        }
         // MPI_File_read( fin, mrc_data_block,length, MPI_SHORT, MPI_STATUS_IGNORE);
        break;

        case MRC_MODE_FLOAT:
       // fread(mrc_data_block,psize,head->nx*head->ny*(end-start),fin); 
        MPI_File_read(fin, mrc_data_block,length, MPI_FLOAT,MPI_STATUS_IGNORE);
        break;
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
 


	long i;
	unsigned char buf_byte;
	short buf_short;
	short buf_ushort;

        int length;
        length=head->nx*head->ny*head->nz;
	//fseek(fin,(1024),SEEK_SET);
        MPI_File_seek( fin, 1024, MPI_SEEK_SET ); 
        printf("length is %d\n",length);
	switch(head->mode)
	{
		case MRC_MODE_BYTE:



                for(i=0;i<length;i++)
                {
          
                   MPI_File_read( fin, &buf_byte,1, MPI_CHAR,MPI_STATUS_IGNORE);
                   mrc_data_all[i]=(float)buf_byte;
                 }
        

	        break;
		
		case MRC_MODE_SHORT:

			/*for(i=0;i<head->nx*head->ny*head->nz;i++)
			{
				fread(&buf_short,sizeof(short),1,fin);
				mrc_data_all[i] = (float)buf_short;
			}*/
                       for(i=0;i<length;i++)
                       {
          
                           MPI_File_read( fin, &buf_short,1, MPI_SHORT,MPI_STATUS_IGNORE);
                           mrc_data_all[i]=(float)buf_short;
                           if (i==length-1)
                              printf("mrc_data_all[%d] is %f\n",i,mrc_data_all[i]);
                        }
			break;

		case MRC_MODE_USHORT:

			for(i=0;i<length;i++)
			{
                           MPI_File_read( fin, &buf_short,1, MPI_SHORT,MPI_STATUS_IGNORE);
                           mrc_data_all[i]=(float)buf_ushort;
			}
			break;

		case MRC_MODE_FLOAT:

                        MPI_File_read( fin, mrc_data_all,length, MPI_FLOAT,MPI_STATUS_IGNORE);                           
			break;

		default:
			printf("Error with Function 'mrc_read_all()'! File type unknown!");
			break;

	}

	return 0;
}

/*****************************************************************************************************/


int mrc_read_all(FILE *fin, MrcHeader  *head, float *mrc_data_all)
{

//check the mrc file to make sure the size is exact in register with the head
	int headsize ;
	switch(head->mode)
	{
		case MRC_MODE_BYTE:
			headsize=get_file_size(fin) - (long long) head->nx * head->ny * head->nz * sizeof(char);
			break;

		case MRC_MODE_SHORT:
		case MRC_MODE_USHORT:
			headsize=get_file_size(fin) - (long long) head->nx * head->ny * head->nz * sizeof(short);
			break;

		case MRC_MODE_FLOAT:
			headsize=get_file_size(fin) - (long long) head->nx * head->ny * head->nz * sizeof(float);

			break;
			
		default:
			printf("Error with Function 'mrc_read_all()'! File type unknown!");

			break;
	}


	long i;
	unsigned char buf_byte;
	short buf_short;
	short buf_ushort;
	float buf_float;

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

			for(i=0;i<head->nx*head->ny*head->nz;i++)
			{
				fread(&buf_float,sizeof(float),1,fin);
				mrc_data_all[i] = (float)buf_float;
			}
			break;
/*			if((fread(mrc_data_all,head->nx*head->ny*head->nz*sizeof(float),1,fin)==0))
			{
				printf("Error with Function 'mrc_read_all()'! Reading file failed!");
				return -1;
			}*/

		default:
			printf("Error with Function 'mrc_read_all()'! File type unknown!");
			break;

	}

	return 0;
}
	


/*****************************************************************************************************/
int mrc_write_all(MPI_File fout, MrcHeader  *head,  int Z_start,int Z_end, float *mrc_data_all)
{
   /* int psize;

    switch(head->mode)
    
    {
      case MRC_MODE_BYTE:
       psize=sizeof(char);
       
      break;

      case MRC_MODE_SHORT:
       psize=sizeof(short);

      break;

      case MRC_MODE_FLOAT:
       psize=sizeof(float);

      break;
    }*/
    MPI_Offset offset;
    offset=head->nx*head->ny*Z_start*sizeof(float);
    offset+=1024;
    MPI_File_seek( fout, offset, MPI_SEEK_SET ); 
    MPI_File_write(fout, mrc_data_all, head->nx*head->ny*(Z_end-Z_start), MPI_FLOAT, MPI_STATUS_IGNORE);
    

return true;
}
