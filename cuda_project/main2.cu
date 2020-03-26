#include "read_write_mrc.h"
#include "atom.h"
#include <iostream>
#include <cuda_runtime.h>
#include <sys/time.h>

#define FALSE 0
#define TRUE 1
using namespace std;

int PrjXYAngN;
int vol_pixel_num;
double iStart,iElaps;

int iLen = 256;

double cpuSecond(){
	struct timeval tp;
	gettimeofday(&tp,NULL);
	//sce + msec
	return (double)tp.tv_sec +(double )tp.tv_usec*1e-6;
}

__host__ int read_data(Volume *vol,Projection *prj,MrcHeader *in_head,MrcHeader *out_head,double *x_coef,double *y_coef,float *prj_real,float *vol_real,char *in_addr,char *angle_addr)
{
/***********************Malloc space**********************************************/
	//cudaMallocManaged((void **)&vol,sizeof(Volume));
	//cudaMallocManaged((void **)&prj,sizeof(Projection));
	//cudaMallocManaged((void **)&in_head,sizeof(MrcHeader));
	//cudaMallocManaged((void **)&out_head,sizeof(MrcHeader));
	vol = (Volume *)malloc(sizeof(Volume));
	prj = (Projection *)malloc(sizeof(Projection));
	in_head = (MrcHeader *)malloc(sizeof(MrcHeader));
	out_head = (MrcHeader *)malloc(sizeof(MrcHeader));

/*****************read head-file ande angle-file*************************************************************/
	FILE *in_file;
	in_file = fopen(in_addr,"r");
	if(!in_file){
		printf("Can not open in_file");
		return FALSE;	
	}
	mrc_read_head(in_file,in_head);
	fclose(in_file);
	printf("%d %d %d\n",in_head->nx,in_head->ny,in_head->nz);
	prj->X = in_head->nx;
	prj->Y = in_head->ny;
	prj->AngN = in_head->nz;
	
	vol->X=1467;
	vol->Y=1521;
	vol->Z=58;

	vol->Xstart=-259;
	vol->Xend=vol->Xstart+vol->X;
	vol->Ystart=-209;
	vol->Yend=vol->Ystart+vol->Y;
	vol->Zstart=-32;
	vol->Zend=vol->Zstart+vol->Z;

	mrc_init_head(out_head);
	out_head->nx=vol->X;
	out_head->ny=vol->Y;
	out_head->nz=vol->Z;

	out_head->nxstart=vol->Xstart;
	out_head->nystart=vol->Ystart;
	out_head->nzstart=vol->Zstart;

	out_head->mx=vol->X;
	out_head->my=vol->Y;
	out_head->mz=vol->Z;

	//cudaMallocManaged((void **)&x_coef,sizeof(double)*prj->AngN*10);
	//cudaMallocManaged((void **)&y_coef,sizeof(double)*prj->AngN*10);
	x_coef = (double *)malloc(sizeof(double)*prj->AngN*10);
	y_coef = (double *)malloc(sizeof(double)*prj->AngN*10);
	memset(x_coef, 0 , sizeof(double)*prj->AngN*10);
	memset(y_coef, 0 , sizeof(double)*prj->AngN*10);

	printf("%d\n",sizeof(double)*prj->AngN*10);

	FILE *angle_file;
	angle_file = fopen(angle_addr,"r");
	if(!angle_file){
		printf("Can not open angle_file");
		return FALSE;	
	}
	read_coef(x_coef, y_coef, angle_file);
	fclose(angle_file);
/***********************************************************************************************************/

/*****************************************read all data*****************************************************/
	PrjXYAngN = prj->X*prj->Y*prj->AngN;
	vol_pixel_num = vol->X*vol->Y*vol->Z;

	//cudaMallocManaged((void **)&prj_real,sizeof(float)*PrjXYAngN);
	//cudaMallocManaged((void **)&vol_real,sizeof(float)*vol_pixel_num);
	prj_real=(float *)malloc(sizeof(float)*PrjXYAngN);
	vol_real = (float *)malloc(sizeof(float)*vol_pixel_num);
	memset(prj_real, 0 , sizeof(float)*PrjXYAngN);
	memset(vol_real, 0 , sizeof(float)*vol_pixel_num);
	
	FILE *in_context=fopen(in_addr,"r");
	if(!in_context){
		printf("Can not open in_context");
		return FALSE;	
	}
	mrc_read_all(in_context,in_head,prj_real);
	fclose(in_context);
/*******************************************************************************************************/

	return TRUE;
}

/**
	Pixel pixel;
	Weight *comp_prj;
	comp_prj=(Weight *)malloc(sizeof(Weight));
	if((comp_prj=(Weight *)malloc(sizeof(Weight)))==NULL)
	{
		printf("Error with Function 'read_data()'!Can't malloc memery for 'comp_prj'!");
		return FALSE;
	}
        Slice_backproj_bilinear(pixel, prj, vol, comp_prj, x_coef, y_coef, prj_real,slc_data,Z_start,Z_end,slcN_per); //four-weight 
**/

/*
	int sizeofZ_per_block = vol->Z/process_num+1;
	int* Z_start = (int *)malloc(sizeof(int)*process_num);
	int* Z_end = (int *)malloc(sizeof(int)*process_num);//the start or end slice of reproject per process
	int* Z_per = (int *)malloc(sizeof(int)*process_num);
	printf("%d / %d = %d \n",vol->Z,process_num,sizeofZ_per_block);
	for(int i=0;i<process_num;i++)
	{
		Z_start[i] = vol->Zstart+i*sizeofZ_per_block;
		Z_end[i] = min(Z_start[i]+sizeofZ_per_block,vol->Zstart+vol->Z);
		Z_per[i] = Z_end[i]-Z_start[i];
		printf("For the %dth block, Z_start is %d, Z_end is %d,Z_per is %d\n",i,Z_start[i],Z_end[i],Z_per[i]);
	}
*/

__device__ void computeWeight(int x,int y,int z,int angle,double *x_coef,double *y_coef,Weight *weight)
{
	double res_x,res_y;
	int index = 4*angle;
	res_x = x_coef[index]+x_coef[index+1]*x+x_coef[index+2]*y+x_coef[index+3]*z;
	res_y = y_coef[index]+y_coef[index+1]*x+y_coef[index+2]*y+y_coef[index+3]*z;	

	weight->x_min = floor(res_x);
	weight->y_min = floor(res_y);

	weight->x_min_del = res_x - weight->x_min;
	weight->y_min_del = res_y - weight->y_min;
}

__global__ void backProjOnGPU(Projection *prj,Volume *vol,double *x_coef,double *y_coef,float *prj_real,float *vol_real,Weight *weight,int *test)
{
	double divisor;//分子
	double dividend;//分母
	int z = threadIdx.x+blockIdx.x*blockDim.x;
	if(z>=vol->Z) return;
	int x,y,index,angle,n;
	//test[z] = 1;
	//printf("Slice %d is beginning.\n",z);

	for(y=vol->Ystart;y<vol->Yend;y++)
	{
		for(x=vol->Xstart;x<vol->Xend;x++)
		{
		/*	divisor = 0;
			dividend = 0;
			index = (x-vol->Xstart)+(y-vol->Ystart)*vol->X;//+(z-vol->Zstart)*vol->X*vol->Y;
			//要计算的像素在三维结构中的坐标
			//if(x-vol->Xstart<2&&y-vol->Ystart<2&&index<vol->X)
			//printf("%d %d %d %d %d\n",index,x,vol->Xstart,y,vol->Ystart);
			for(angle=0;angle<prj->AngN;angle++)
			{
				computeWeight(x,y,z,angle,x_coef,y_coef,weight);
				//printf("%lf %lf \n",weight->x_min,weight->y_min);
				if(weight->x_min>=0 && weight->x_min<prj->X && weight->y_min>=0 && weight->y_min<prj->Y)//(x_min,y_min)
				{
					n = weight->x_min + weight->y_min*prj->X + angle*prj->X*prj->Y;
					divisor += (1-weight->x_min_del)*(1-weight->y_min_del)*prj_real[n];
					dividend += (1-weight->x_min_del)*(1-weight->y_min_del);
				}
				if(weight->x_min>=0 && weight->x_min<prj->X && weight->y_min>=0 && weight->y_min<prj->Y)//(x_min+1,y_min)
				{
					n = (weight->x_min+1) + weight->y_min*prj->X + angle*prj->X*prj->Y;
					divisor += weight->x_min_del*(1-weight->y_min_del)*prj_real[n];
					dividend += weight->x_min_del*(1-weight->y_min_del);
				}
				if(weight->x_min>=0 && weight->x_min<prj->X && weight->y_min>=0 && weight->y_min<prj->Y)//(x_min,y_min+1)
				{
					n = weight->x_min + (weight->y_min+1)*prj->X + angle*prj->X*prj->Y;
					divisor += (1-weight->x_min_del)*weight->y_min_del*prj_real[n];
					dividend += (1-weight->x_min_del)*weight->y_min_del;
				}
				if(weight->x_min>=0 && weight->x_min<prj->X && weight->y_min>=0 && weight->y_min<prj->Y)//(x_min+1,y_min+1)
				{
					n = (weight->x_min+1)+ (weight->y_min+1)*prj->X + angle*prj->X*prj->Y;
					divisor += weight->x_min_del*weight->y_min_del*prj_real[n];
					dividend += weight->x_min_del*weight->y_min_del;
				}
			}
			if(dividend!=0.0f)
			{
				vol_real[index] = (float)(divisor/dividend);
				//printf("vol_read[%d]:%f\n",index,vol_real[index]);
			}
			
			//if(19000<index&&index<19010)
			//if(index<10)
			//printf("%d-------------%lf\n",index,vol_real[index]);
		*/
		}
		
	}
	//test[z] = 2;
	//printf("Slice %d has been finished.\n",z);
}

void write_data(char *out_addr,MrcHeader *out_head,int Z_end,float *vol_real)
{
	FILE *out_file;
	out_file = fopen(out_addr,"w");
	if(!out_file){
		printf("Can not open in_file");
		return;	
	}
	mrc_write_head(out_file,out_head);
	printf("siezof out_head %ld \n",sizeof(out_head));
	mrc_write_all(out_file,out_head,Z_end,vol_real);
	mrc_update_head(out_file);
	fclose(out_file);
	return;
}

int main(int argc,char *argv[])
{
	iStart = cpuSecond();

	iLen = atoi(argv[1]);
	char* in_addr = argv[2];
	char* out_addr = argv[3];
	char* angle_addr = argv[4];
	//cout<<"arg1:"<<process_num<<" arg2:"<<in_addr<<"  arg3:"<<out_addr<<"  arg4:"<<angle_addr<<endl;
	
	/*******************Read data_in_host********************/
	Volume *h_vol;
	Projection *h_prj;
	MrcHeader *in_head,*out_head;
	double *h_x_coef;
	double *h_y_coef;
	float *h_prj_real,*h_vol_real;//pri_real is inputted data ande vol_real is calculated data

	read_data(h_vol,h_prj,in_head,out_head,h_x_coef,h_y_coef,h_prj_real,h_vol_real,in_addr,angle_addr);
	/******************************************************/

	/********************copy data from host to device************************/
	Volume *vol;
	Projection *prj;
	double *x_coef;
	double *y_coef;
	float *prj_real,*vol_real;
	Weight *weight;
	cudaMalloc((void **)&vol,sizeof(Volume));
	cudaMalloc((void **)&prj,sizeof(Projection));
	cudaMalloc((void **)&x_coef,sizeof(double)*h_prj->AngN*10);
	cudaMalloc((void **)&y_coef,sizeof(double)*h_prj->AngN*10);
	cudaMalloc((void **)&prj_real,sizeof(float)*PrjXYAngN);
	cudaMalloc((void **)&vol_real,sizeof(float)*vol_pixel_num);
	cudaMalloc((void **)&weight,sizeof(Weight));
	
	cudaMemcpy((void *)vol,(void *)h_vol,sizeof(Volume),cudaMemcpyHostToDevice);
	cudaMemcpy((void *)prj,(void *)h_prj,sizeof(Projection),cudaMemcpyHostToDevice);
	cudaMemcpy((void *)x_coef,(void *)h_x_coef,sizeof(double)*h_prj->AngN*10,cudaMemcpyHostToDevice);
	cudaMemcpy((void *)y_coef,(void *)h_y_coef,sizeof(double)*h_prj->AngN*10,cudaMemcpyHostToDevice);
	cudaMemcpy((void *)prj_real,(void *)h_prj_real,sizeof(float)*PrjXYAngN,cudaMemcpyHostToDevice);
	/********************************************************************/
	
	/*********************Debug**********************/
	printf("%d\n",h_vol->Z);
	int *test = (int *)malloc(sizeof(int)*h_vol->Z);
	//cudaMallocManaged((void **)&test,vol->Z*sizeof(int)*10);
	for(int i=0;i<h_vol->Z;i++) test[i] = 0;
	int *d_test;
	cudaMalloc((void **)&d_test,sizeof(int)*h_vol->Z);
	cudaMemcpy(d_test,test,sizeof(int)*h_vol->Z,cudaMemcpyHostToDevice);
	/************************************************/


	/************run on GPU*********************/
	dim3 block(iLen);
	dim3 grid((h_vol->Z+block.x-1)/block.x);
	backProjOnGPU<<<grid,block>>>(prj,vol,x_coef,y_coef,prj_real,vol_real,weight,test);	
	//cudaDeviceSynchronize();
	/******************************************/

	/************copy result from Device to host***********/
	cudaMemcpy((void *)h_vol_real,(void *)vol_real,sizeof(float)*vol_pixel_num,cudaMemcpyDeviceToHost);	
	/*******/

	//write_data(out_addr,out_head,h_vol->Z,h_vol_real);

	/************Debug******************/
	cudaMemcpy((void *)test,(void *)d_test,sizeof(int)*h_vol->Z,cudaMemcpyDeviceToHost);
	for(int i=0;i<vol->Z;i++) printf("%d ",test[i]);
	/******************************/

	/*************free space****************************/
	cudaFree(vol);
	cudaFree(prj);
	cudaFree(x_coef);
	cudaFree(y_coef);
	cudaFree(prj_real);
	cudaFree(vol_real);
	free(h_vol);
	free(h_prj);
	free(h_x_coef);
	free(h_y_coef);
	free(prj_real);
	free(vol_real);
	free(in_head);
	free(out_head);
	/****************************************************/

	iElaps = cpuSecond()-iStart;
	printf("Host time elapsed:%lfsec\n",iElaps);
	return 0;
}

