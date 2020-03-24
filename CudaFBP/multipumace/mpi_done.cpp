

/**********************************************************************************main.c*/



#include "atom.h"



void datadeepcopy(float *dst,float *src,int datasize)
{
	for(int i=0;i<datasize;i++)
		dst[i]=src[i];
}
void writefilefromzero(char *initfile,float *data,int datasize,Volume vol)
{
	   MPI_File  fout;
	   MPI_File_open(MPI_COMM_WORLD, initfile, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fout);
	   MrcHeader *outhead;

	   outhead=(MrcHeader *)malloc(sizeof(MrcHeader));

	   mrc_init_head(outhead);
	   outhead->nx=vol.X;
	   outhead->ny=vol.Y;
	   outhead->nz=vol.Z;
	   mrc_write_head(fout,outhead);
	   int initoffset=outhead->nx*outhead->ny;
	   int offset=0;
	   for(int i=0;i<outhead->nz;i++)
	   {
		   mrc_write_slice(fout, outhead, i, 'z', data+offset);
		   offset=initoffset+offset;
	   }
		mrc_update_head(fout);
		free(outhead);
		MPI_File_close(&fout);

}

void readdatafromfile(char *filename,float *data)
{
	   MPI_File  fin;
	   MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fin);
	   MrcHeader *inhead;

	   inhead=(MrcHeader *)malloc(sizeof(MrcHeader));
	   mrc_read_head_MPI(fin, inhead);


	   int initoffset=inhead->nx*inhead->ny;
	   int offset=0;
	   for(int i=0;i<inhead->nz;i++)
	   {
		   mrc_read_slice(fin, inhead, i, 'z', data+offset);
		   offset=initoffset+offset;
	   }
		free(inhead);
		MPI_File_close(&fin);
}


int main(int argc, char *argv[])
{
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
   int ATOM_Out_NUM=1;
   int ATOM_Iner_NUM=1;

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
		   ATOM_Iner_NUM=atoi(argv[i+1]);
		   break;
		 case 't':
		   ATOM_Out_NUM=atof(argv[i+1]);
		   break;
	  }
   }



   char *fbpf=argv[3];

   MPI_File  f_fbp;
   MPI_File_open(MPI_COMM_WORLD, fbpf, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &f_fbp);
   MrcHeader  *fbphead;
   fbphead=(MrcHeader *)malloc(sizeof(MrcHeader));
   mrc_read_head_MPI(f_fbp,fbphead);

   Volume vol;
   vol.X=fbphead->nx;
   vol.Y=fbphead->ny;
   vol.Z=fbphead->nz;
   vol.Xstart=-37;
   vol.Xend=vol.Xstart+vol.X;
   vol.Ystart=-61;
   vol.Yend=vol.Ystart+vol.Y;
   vol.Zstart=-33;
   vol.Zend=vol.Zstart+vol.Z;

   int num_subset=2;



	int wholelength=vol.X*vol.Y*vol.Z;
	float *zhatdata=(float* )malloc(sizeof(float)*wholelength*2);
	float *Gw=(float* )malloc(sizeof(float)*wholelength*2);
    //memset(GW,0,sizeof(float)*wholelength*2);
	float *v=(float* )malloc(sizeof(float)*wholelength*2);
	float *wdata=(float* )malloc(sizeof(float)*wholelength*2);
	float *FV=(float* )malloc(sizeof(float)*wholelength*2);
	int max_iter=ATOM_Out_NUM;
	float gamma=0.5;
	//init read in
    int start=vol.Zstart-vol.Zstart;
    int end=vol.Zend-vol.Zstart;
    mrc_read_block(f_fbp, fbphead, start, end, 'Z', wdata);

    MPI_File_close(&f_fbp);
    datadeepcopy(wdata+wholelength,wdata,wholelength);
    datadeepcopy(zhatdata,wdata,wholelength*2);
    if(!id)
    	printf("After read init data\n");

	char subset1init[20];
	char subset2init[20];
	char subset1recon[20];
	char subset2recon[20];

	char sumrecon[20];


	//eel-tomo4.st  eel-tomo4.txbr
	for(i=0;i<max_iter;i++)
	{
	    //%==== Update v ====%
		int k=0;
		datadeepcopy(Gw,zhatdata,wholelength*2); //    Gw = zhat;
		for(k=0;k<wholelength*2;k++)             //	 v  = 2*Gw - wdata;
		{
			Gw[k]=Gw[k]*2;
			v[k]=Gw[k]-wdata[k];
		}
		sprintf(subset1init,"%03dsubset1%03d.mrc",max_iter,i);
		sprintf(subset2init,"%03dsubset2%03d.mrc",max_iter,i);
		sprintf(subset1recon,"%03dsubset1recon%03d.mrc",max_iter,i);
		sprintf(subset2recon,"%03dsubset2recon%03d.mrc",max_iter,i);
		sprintf(sumrecon,"%03dsumrecon%03d.mrc",max_iter,i);
		writefilefromzero(subset1init,v,wholelength,vol);
		writefilefromzero(subset2init,v+wholelength,wholelength,vol);
		for(k=0;k<num_subset;k++)
		{
			if(k==0)
			{
				int flagpart=k;
				ATOM_MPI_SIRT(subset1recon, subset1init, ATOM_ITR_STEP, ATOM_Iner_NUM,id,p,flagpart);
				MPI_Barrier(MPI_COMM_WORLD);
				readdatafromfile(subset1recon,FV);
			    if(!id)
			    	printf("After first data\n");
			}

			if(k==1)
			{
				int flagpart=k;
				ATOM_MPI_SIRT(subset2recon, subset2init, ATOM_ITR_STEP, ATOM_Iner_NUM,id,p,flagpart);
				MPI_Barrier(MPI_COMM_WORLD);
				readdatafromfile(subset2recon,FV+wholelength);
			    if(!id)
			    	printf("After second data\n");
			}
		}
		float subgamma=1-gamma;
		for(k=0;k<wholelength*2;k++)
		{
			wdata[k]=subgamma*wdata[k]+ gamma*(2*FV[k]-v[k]);

		}//对整体的w进行更新
		for(k=0;k<wholelength;k++)
		{
			zhatdata[k]=0.5*wdata[k]+0.5*(wdata+wholelength)[k];
			(zhatdata+wholelength)[k]=zhatdata[k];
		}//这是取w平均值的操作
		if(id==0)
			writefilefromzero(sumrecon,zhatdata,wholelength,vol);

		printf("The Iteration finish");

	}

	MPI_Barrier(MPI_COMM_WORLD );

	MPI::Finalize();
	printf("Sucess \n");

	return 0;




}
/*
2233.158936 2257.893799 2280.263916 2318.859863 2347.979004 2354.653564 2363.364990 2405.657959 2412.657715 2424.730225
After read init data
updating head!
head->amin is 0.000000, head->amax is 3190.257568, head->amean is 2720.039307
updating finished!
updating head!
head->amin is 0.000000, head->amax is 3190.257568, head->amean is 2720.039307
updating finished!
Z_start is -33, Z_end is 33 in process 0
slice -33 is finished to reproj.
slice -32 is finished to reproj.
slice -31 is finished to reproj.
slice -30 is finished to reproj.
slice -29 is finished to reproj.
slice -28 is finished to reproj.
slice -27 is finished to reproj.
slice -26 is finished to reproj.
slice -25 is finished to reproj.
slice -24 is finished to reproj.
slice -23 is finished to reproj.
slice -22 is finished to reproj.
slice -21 is finished to reproj.
slice -20 is finished to reproj.
slice -19 is finished to reproj.
slice -18 is finished to reproj.
slice -17 is finished to reproj.
slice -16 is finished to reproj.
slice -15 is finished to reproj.
slice -14 is finished to reproj.
slice -13 is finished to reproj.
slice -12 is finished to reproj.
slice -11 is finished to reproj.
slice -10 is finished to reproj.
slice -9 is finished to reproj.
slice -8 is finished to reproj.
slice -7 is finished to reproj.
slice -6 is finished to reproj.
slice -5 is finished to reproj.
slice -4 is finished to reproj.
slice -3 is finished to reproj.
slice -2 is finished to reproj.
slice -1 is finished to reproj.
slice 0 is finished to reproj.
slice 1 is finished to reproj.
slice 2 is finished to reproj.
slice 3 is finished to reproj.
slice 4 is finished to reproj.
slice 5 is finished to reproj.
slice 6 is finished to reproj.
slice 7 is finished to reproj.
slice 8 is finished to reproj.
slice 9 is finished to reproj.
slice 10 is finished to reproj.
slice 11 is finished to reproj.
slice 12 is finished to reproj.
slice 13 is finished to reproj.
slice 14 is finished to reproj.
slice 15 is finished to reproj.
slice 16 is finished to reproj.
slice 17 is finished to reproj.
slice 18 is finished to reproj.
slice 19 is finished to reproj.
slice 20 is finished to reproj.
slice 21 is finished to reproj.
slice 22 is finished to reproj.
slice 23 is finished to reproj.
slice 24 is finished to reproj.
slice 25 is finished to reproj.
slice 26 is finished to reproj.
slice 27 is finished to reproj.
slice 28 is finished to reproj.
slice 29 is finished to reproj.
slice 30 is finished to reproj.
slice 31 is finished to reproj.
slice 32 is finished to reproj.
slice -33 is completed to update.
slice -32 is completed to update.
slice -31 is completed to update.
slice -30 is completed to update.
slice -29 is completed to update.
slice -28 is completed to update.
slice -27 is completed to update.
slice -26 is completed to update.
slice -25 is completed to update.
slice -24 is completed to update.
slice -23 is completed to update.
slice -22 is completed to update.
slice -21 is completed to update.
slice -20 is completed to update.
slice -19 is completed to update.
slice -18 is completed to update.
slice -17 is completed to update.
slice -16 is completed to update.
slice -15 is completed to update.
slice -14 is completed to update.
slice -13 is completed to update.
slice -12 is completed to update.
slice -11 is completed to update.
slice -10 is completed to update.
slice -9 is completed to update.
slice -8 is completed to update.
slice -7 is completed to update.
slice -6 is completed to update.
slice -5 is completed to update.
slice -4 is completed to update.
slice -3 is completed to update.
slice -2 is completed to update.
slice -1 is completed to update.
slice 0 is completed to update.
slice 1 is completed to update.
slice 2 is completed to update.
slice 3 is completed to update.
slice 4 is completed to update.
slice 5 is completed to update.
slice 6 is completed to update.
slice 7 is completed to update.
slice 8 is completed to update.
slice 9 is completed to update.
slice 10 is completed to update.
slice 11 is completed to update.
slice 12 is completed to update.
slice 13 is completed to update.
slice 14 is completed to update.
slice 15 is completed to update.
slice 16 is completed to update.
slice 17 is completed to update.
slice 18 is completed to update.
slice 19 is completed to update.
slice 20 is completed to update.
slice 21 is completed to update.
slice 22 is completed to update.
slice 23 is completed to update.
slice 24 is completed to update.
slice 25 is completed to update.
slice 26 is completed to update.
slice 27 is completed to update.
slice 28 is completed to update.
slice 29 is completed to update.
slice 30 is completed to update.
slice 31 is completed to update.
slice 32 is completed to update.
updating head!
head->amin is 0.000000, head->amax is 3336.772217, head->amean is 2716.203613
updating finished!
After first data
Z_start is -33, Z_end is 33 in process 0
slice -33 is finished to reproj.
slice -32 is finished to reproj.
slice -31 is finished to reproj.
slice -30 is finished to reproj.
slice -29 is finished to reproj.
slice -28 is finished to reproj.
slice -27 is finished to reproj.
slice -26 is finished to reproj.
slice -25 is finished to reproj.
slice -24 is finished to reproj.
slice -23 is finished to reproj.
slice -22 is finished to reproj.
slice -21 is finished to reproj.
slice -20 is finished to reproj.
slice -19 is finished to reproj.
slice -18 is finished to reproj.
slice -17 is finished to reproj.
slice -16 is finished to reproj.
slice -15 is finished to reproj.
slice -14 is finished to reproj.
slice -13 is finished to reproj.
slice -12 is finished to reproj.
slice -11 is finished to reproj.
slice -10 is finished to reproj.
slice -9 is finished to reproj.
slice -8 is finished to reproj.
slice -7 is finished to reproj.
slice -6 is finished to reproj.
slice -5 is finished to reproj.
slice -4 is finished to reproj.
slice -3 is finished to reproj.
slice -2 is finished to reproj.
slice -1 is finished to reproj.
slice 0 is finished to reproj.
slice 1 is finished to reproj.
slice 2 is finished to reproj.
slice 3 is finished to reproj.
slice 4 is finished to reproj.
slice 5 is finished to reproj.
slice 6 is finished to reproj.
slice 7 is finished to reproj.
slice 8 is finished to reproj.
slice 9 is finished to reproj.
slice 10 is finished to reproj.
slice 11 is finished to reproj.
slice 12 is finished to reproj.
slice 13 is finished to reproj.
slice 14 is finished to reproj.
slice 15 is finished to reproj.
slice 16 is finished to reproj.
slice 17 is finished to reproj.
slice 18 is finished to reproj.
slice 19 is finished to reproj.
slice 20 is finished to reproj.
slice 21 is finished to reproj.
slice 22 is finished to reproj.
slice 23 is finished to reproj.
slice 24 is finished to reproj.
slice 25 is finished to reproj.
slice 26 is finished to reproj.
slice 27 is finished to reproj.
slice 28 is finished to reproj.
slice 29 is finished to reproj.
slice 30 is finished to reproj.
slice 31 is finished to reproj.
slice 32 is finished to reproj.
slice -33 is completed to update.
slice -32 is completed to update.
slice -31 is completed to update.
slice -30 is completed to update.
slice -29 is completed to update.
slice -28 is completed to update.
slice -27 is completed to update.
slice -26 is completed to update.
slice -25 is completed to update.
slice -24 is completed to update.
slice -23 is completed to update.
slice -22 is completed to update.
slice -21 is completed to update.
slice -20 is completed to update.
slice -19 is completed to update.
slice -18 is completed to update.
slice -17 is completed to update.
slice -16 is completed to update.
slice -15 is completed to update.
slice -14 is completed to update.
slice -13 is completed to update.
slice -12 is completed to update.
slice -11 is completed to update.
slice -10 is completed to update.
slice -9 is completed to update.
slice -8 is completed to update.
slice -7 is completed to update.
slice -6 is completed to update.
slice -5 is completed to update.
slice -4 is completed to update.
slice -3 is completed to update.
slice -2 is completed to update.
slice -1 is completed to update.
slice 0 is completed to update.
slice 1 is completed to update.
slice 2 is completed to update.
slice 3 is completed to update.
slice 4 is completed to update.
slice 5 is completed to update.
slice 6 is completed to update.
slice 7 is completed to update.
slice 8 is completed to update.
slice 9 is completed to update.
slice 10 is completed to update.
slice 11 is completed to update.
slice 12 is completed to update.
slice 13 is completed to update.
slice 14 is completed to update.
slice 15 is completed to update.
slice 16 is completed to update.
slice 17 is completed to update.
slice 18 is completed to update.
slice 19 is completed to update.
slice 20 is completed to update.
slice 21 is completed to update.
slice 22 is completed to update.
slice 23 is completed to update.
slice 24 is completed to update.
slice 25 is completed to update.
slice 26 is completed to update.
slice 27 is completed to update.
slice 28 is completed to update.
slice 29 is completed to update.
slice 30 is completed to update.
slice 31 is completed to update.
slice 32 is completed to update.
updating head!
head->amin is 0.000000, head->amax is 3267.509277, head->amean is 2719.735107
updating finished!
After second data
updating head!
head->amin is 0.000000, head->amax is 3302.140625, head->amean is 2717.969238
updating finished!
The Iteration finishupdating head!
head->amin is 0.000000, head->amax is 3267.509033, head->amean is 2719.735107
updating finished!
updating head!
head->amin is 0.000000, head->amax is 3336.771973, head->amean is 2716.203613
updating finished!
Z_start is -33, Z_end is 33 in process 0
slice -33 is finished to reproj.
slice -32 is finished to reproj.
slice -31 is finished to reproj.
slice -30 is finished to reproj.
slice -29 is finished to reproj.
slice -28 is finished to reproj.
slice -27 is finished to reproj.
slice -26 is finished to reproj.
slice -25 is finished to reproj.
slice -24 is finished to reproj.
slice -23 is finished to reproj.
slice -22 is finished to reproj.
slice -21 is finished to reproj.
slice -20 is finished to reproj.
slice -19 is finished to reproj.
slice -18 is finished to reproj.
slice -17 is finished to reproj.
slice -16 is finished to reproj.
slice -15 is finished to reproj.
slice -14 is finished to reproj.
slice -13 is finished to reproj.
slice -12 is finished to reproj.
slice -11 is finished to reproj.
slice -10 is finished to reproj.
slice -9 is finished to reproj.
slice -8 is finished to reproj.
slice -7 is finished to reproj.
slice -6 is finished to reproj.
slice -5 is finished to reproj.
slice -4 is finished to reproj.
slice -3 is finished to reproj.
slice -2 is finished to reproj.
slice -1 is finished to reproj.
slice 0 is finished to reproj.
slice 1 is finished to reproj.
slice 2 is finished to reproj.
slice 3 is finished to reproj.
slice 4 is finished to reproj.
slice 5 is finished to reproj.
slice 6 is finished to reproj.
slice 7 is finished to reproj.
slice 8 is finished to reproj.
slice 9 is finished to reproj.
slice 10 is finished to reproj.
slice 11 is finished to reproj.
slice 12 is finished to reproj.
slice 13 is finished to reproj.
slice 14 is finished to reproj.
slice 15 is finished to reproj.
slice 16 is finished to reproj.
slice 17 is finished to reproj.
slice 18 is finished to reproj.
slice 19 is finished to reproj.
slice 20 is finished to reproj.
slice 21 is finished to reproj.
slice 22 is finished to reproj.
slice 23 is finished to reproj.
slice 24 is finished to reproj.
slice 25 is finished to reproj.
slice 26 is finished to reproj.
slice 27 is finished to reproj.
slice 28 is finished to reproj.
slice 29 is finished to reproj.
slice 30 is finished to reproj.
slice 31 is finished to reproj.
slice 32 is finished to reproj.
slice -33 is completed to update.
slice -32 is completed to update.
slice -31 is completed to update.
slice -30 is completed to update.
slice -29 is completed to update.
slice -28 is completed to update.
slice -27 is completed to update.
slice -26 is completed to update.
slice -25 is completed to update.
slice -24 is completed to update.
slice -23 is completed to update.
slice -22 is completed to update.
slice -21 is completed to update.
slice -20 is completed to update.
slice -19 is completed to update.
slice -18 is completed to update.
slice -17 is completed to update.
slice -16 is completed to update.
slice -15 is completed to update.
slice -14 is completed to update.
slice -13 is completed to update.
slice -12 is completed to update.
slice -11 is completed to update.
slice -10 is completed to update.
slice -9 is completed to update.
slice -8 is completed to update.
slice -7 is completed to update.
slice -6 is completed to update.
slice -5 is completed to update.
slice -4 is completed to update.
slice -3 is completed to update.
slice -2 is completed to update.
slice -1 is completed to update.
slice 0 is completed to update.
slice 1 is completed to update.
slice 2 is completed to update.
slice 3 is completed to update.
slice 4 is completed to update.
slice 5 is completed to update.
slice 6 is completed to update.
slice 7 is completed to update.
slice 8 is completed to update.
slice 9 is completed to update.
slice 10 is completed to update.
slice 11 is completed to update.
slice 12 is completed to update.
slice 13 is completed to update.
slice 14 is completed to update.
slice 15 is completed to update.
slice 16 is completed to update.
slice 17 is completed to update.
slice 18 is completed to update.
slice 19 is completed to update.
slice 20 is completed to update.
slice 21 is completed to update.
slice 22 is completed to update.
slice 23 is completed to update.
slice 24 is completed to update.
slice 25 is completed to update.
slice 26 is completed to update.
slice 27 is completed to update.
slice 28 is completed to update.
slice 29 is completed to update.
slice 30 is completed to update.
slice 31 is completed to update.
slice 32 is completed to update.
updating head!
head->amin is 0.000000, head->amax is 3487.013184, head->amean is 2715.938232
updating finished!
After first data
Z_start is -33, Z_end is 33 in process 0
slice -33 is finished to reproj.
slice -32 is finished to reproj.
slice -31 is finished to reproj.
slice -30 is finished to reproj.
slice -29 is finished to reproj.
slice -28 is finished to reproj.
slice -27 is finished to reproj.
slice -26 is finished to reproj.
slice -25 is finished to reproj.
slice -24 is finished to reproj.
slice -23 is finished to reproj.
slice -22 is finished to reproj.
slice -21 is finished to reproj.
slice -20 is finished to reproj.
slice -19 is finished to reproj.
slice -18 is finished to reproj.
slice -17 is finished to reproj.
slice -16 is finished to reproj.
slice -15 is finished to reproj.
slice -14 is finished to reproj.
slice -13 is finished to reproj.
slice -12 is finished to reproj.
slice -11 is finished to reproj.
slice -10 is finished to reproj.
slice -9 is finished to reproj.
slice -8 is finished to reproj.
slice -7 is finished to reproj.
slice -6 is finished to reproj.
slice -5 is finished to reproj.
slice -4 is finished to reproj.
slice -3 is finished to reproj.
slice -2 is finished to reproj.
slice -1 is finished to reproj.
slice 0 is finished to reproj.
slice 1 is finished to reproj.
slice 2 is finished to reproj.
slice 3 is finished to reproj.
slice 4 is finished to reproj.
slice 5 is finished to reproj.
slice 6 is finished to reproj.
slice 7 is finished to reproj.
slice 8 is finished to reproj.
slice 9 is finished to reproj.
slice 10 is finished to reproj.
slice 11 is finished to reproj.
slice 12 is finished to reproj.
slice 13 is finished to reproj.
slice 14 is finished to reproj.
slice 15 is finished to reproj.
slice 16 is finished to reproj.
slice 17 is finished to reproj.
slice 18 is finished to reproj.
slice 19 is finished to reproj.
slice 20 is finished to reproj.
slice 21 is finished to reproj.
slice 22 is finished to reproj.
slice 23 is finished to reproj.
slice 24 is finished to reproj.
slice 25 is finished to reproj.
slice 26 is finished to reproj.
slice 27 is finished to reproj.
slice 28 is finished to reproj.
slice 29 is finished to reproj.
slice 30 is finished to reproj.
slice 31 is finished to reproj.
slice 32 is finished to reproj.
slice -33 is completed to update.
slice -32 is completed to update.
slice -31 is completed to update.
slice -30 is completed to update.
slice -29 is completed to update.
slice -28 is completed to update.
slice -27 is completed to update.
slice -26 is completed to update.
slice -25 is completed to update.
slice -24 is completed to update.
slice -23 is completed to update.
slice -22 is completed to update.
slice -21 is completed to update.
slice -20 is completed to update.
slice -19 is completed to update.
slice -18 is completed to update.
slice -17 is completed to update.
slice -16 is completed to update.
slice -15 is completed to update.
slice -14 is completed to update.
slice -13 is completed to update.
slice -12 is completed to update.
slice -11 is completed to update.
slice -10 is completed to update.
slice -9 is completed to update.
slice -8 is completed to update.
slice -7 is completed to update.
slice -6 is completed to update.
slice -5 is completed to update.
slice -4 is completed to update.
slice -3 is completed to update.
slice -2 is completed to update.
slice -1 is completed to update.
slice 0 is completed to update.
slice 1 is completed to update.
slice 2 is completed to update.
slice 3 is completed to update.
slice 4 is completed to update.
slice 5 is completed to update.
slice 6 is completed to update.
slice 7 is completed to update.
slice 8 is completed to update.
slice 9 is completed to update.
slice 10 is completed to update.
slice 11 is completed to update.
slice 12 is completed to update.
slice 13 is completed to update.
slice 14 is completed to update.
slice 15 is completed to update.
slice 16 is completed to update.
slice 17 is completed to update.
slice 18 is completed to update.
slice 19 is completed to update.
slice 20 is completed to update.
slice 21 is completed to update.
slice 22 is completed to update.
slice 23 is completed to update.
slice 24 is completed to update.
slice 25 is completed to update.
slice 26 is completed to update.
slice 27 is completed to update.
slice 28 is completed to update.
slice 29 is completed to update.
slice 30 is completed to update.
slice 31 is completed to update.
slice 32 is completed to update.
updating head!
head->amin is 0.000000, head->amax is 3472.899902, head->amean is 2716.652588
updating finished!
After second data
updating head!
head->amin is 0.000000, head->amax is 3479.956543, head->amean is 2716.295410
updating finished!
The Iteration finishSucess

inf outf doublebpt.mrc coef SIRT -n 1 -t 0.2
 */
