#include "read_write_mrc.h"
#include "atom.h"
#include <iostream>
#include <cuda_runtime.h>
#include <sys/time.h>

#define FALSE 0
#define TRUE 1
#define checkCudaErrors( a ) do { \
	if (cudaSuccess != (a)) { \
	fprintf(stderr, "Cuda runtime error in line %d of file %s \
	: %s \n", __LINE__, __FILE__, cudaGetErrorString(cudaGetLastError()) ); \
	exit(EXIT_FAILURE); \
	} \
	} while(0);
using namespace std;

int PrjXYAngN;
int vol_pixel_num;
double iStart,iElaps;

int iLen = 32;
int SIRT_ITER_NUM = 4;
float ITER_STEP_LENGTH = 0.2;

double cpuSecond(){
	struct timeval tp;
	gettimeofday(&tp,NULL);
	//sce + msec
	return (double)tp.tv_sec +(double )tp.tv_usec*1e-6;
}

int read_head_data(Volume *vol,Projection *prj,MrcHeader *in_head,MrcHeader *out_head,char *in_addr)
{
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

	printf("%d %d %d 0\n",out_head->nx,out_head->ny,out_head->nz);
	return TRUE;
}


int read_txbr_data(double *x_coef,double *y_coef,char *angle_addr)
{
	FILE *angle_file;
	angle_file = fopen(angle_addr,"r");
	if(!angle_file){
		printf("Can not open angle_file");
		return FALSE;	
	}
	read_coef(x_coef, y_coef, angle_file);
	fclose(angle_file);
	return true;
}

int read_all_data(MrcHeader *in_head,float *prj_real,char *in_addr)
{	
	FILE *in_context=fopen(in_addr,"r");
	if(!in_context){
		printf("Can not open in_context");
		return FALSE;	
	}
	mrc_read_all(in_context,in_head,prj_real);
	fclose(in_context);
	return true;
}

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


__global__ void backProjOnGPU(Projection *prj,Volume *vol,double *x_coef,double *y_coef,float *prj_real,float *vol_real,float iter_step_length)
{
	double divisor;//分子
	double dividend;//分母
	int x = threadIdx.x+blockIdx.x*blockDim.x +vol->Xstart;
	int y = threadIdx.y+blockIdx.y*blockDim.y +vol->Ystart;
	//printf("%d %d\n ",y,z);
	if(x>=vol->Xend || y>=vol->Yend) return;
	int z,index,angle,n;	
	for(z=vol->Zstart;z<vol->Zend;z++)
	{
		divisor = 0;
		dividend = 0;
		for(angle=0;angle<prj->AngN;angle++)
		{
			double res_x,res_y,x_min_del,y_min_del;
			int id = 4*angle,x_min,y_min;
			res_x = x_coef[id]+x_coef[id+1]*x+x_coef[id+2]*y+x_coef[id+3]*z;
			res_y = y_coef[id]+y_coef[id+1]*x+y_coef[id+2]*y+y_coef[id+3]*z;	
			x_min = floor(res_x);
			y_min = floor(res_y);
			x_min_del = res_x - x_min;
			y_min_del = res_y - y_min;
			
			if(x_min>=0 && x_min<prj->X && y_min>=0 && y_min<prj->Y)//(x_min,y_min)
			{
				n = x_min + y_min*prj->X + angle*prj->X*prj->Y;
				divisor += (1-x_min_del)*(1-y_min_del)*prj_real[n];
				dividend += (1-x_min_del)*(1-y_min_del);
			}
			if(x_min+1>=0 && x_min+1<prj->X && y_min>=0 && y_min<prj->Y)//(x_min+1,y_min)
			{
				n = (x_min+1) + y_min*prj->X + angle*prj->X*prj->Y;
				divisor += x_min_del*(1-y_min_del)*prj_real[n];
				dividend += x_min_del*(1-y_min_del);
			}
			if(x_min>=0 && x_min<prj->X && y_min+1>=0 && y_min+1<prj->Y)//(x_min,y_min+1)
			{
				n = x_min + (y_min+1)*prj->X + angle*prj->X*prj->Y;
				divisor += (1-x_min_del)*y_min_del*prj_real[n];
				dividend += (1-x_min_del)*y_min_del;
			}
			if(x_min+1>=0 && x_min+1<prj->X && y_min+1>=0 && y_min+1<prj->Y)//(x_min+1,y_min+1)
			{
				n = (x_min+1)+ (y_min+1)*prj->X + angle*prj->X*prj->Y;
				divisor += x_min_del*y_min_del*prj_real[n];
				dividend += x_min_del*y_min_del;
			}
		}
		if(dividend!=0.0f)
		{
			index = (x-vol->Xstart)+(y-vol->Ystart)*vol->X+(z-vol->Zstart)*vol->X*vol->Y;
			atomicAdd(&vol_real[index], (float)(divisor/dividend)*iter_step_length);
			//vol_real[index] = (float)(divisor/dividend);
			//if(index>vol->X*vol->Y+90000&&index<vol->X*vol->Y+90500)
			//	printf("vol_real[%d]:%f\n",index,vol_real[index]);
		}
	}	
}

__global__ void reProjOnGPU(Projection *prj,Volume *vol,double *x_coef,double *y_coef,float *vol_real,float *iter_prj_divisor,float *iter_prj_dividend)
{
	
	int x = threadIdx.x+blockIdx.x*blockDim.x +vol->Xstart;
	int y = threadIdx.y+blockIdx.y*blockDim.y +vol->Ystart;
	if(x>=vol->Xend || y>=vol->Yend) return;
	int z,index,angle,n;	
	for(z=vol->Zstart;z<vol->Zend;z++)
	{
		index = (x-vol->Xstart)+(y-vol->Ystart)*vol->X+(z-vol->Zstart)*vol->X*vol->Y;
		for(angle=0;angle<prj->AngN;angle++)
		{
			double res_x,res_y,x_min_del,y_min_del;
			int id = 4*angle,x_min,y_min;
			res_x = x_coef[id]+x_coef[id+1]*x+x_coef[id+2]*y+x_coef[id+3]*z;
			res_y = y_coef[id]+y_coef[id+1]*x+y_coef[id+2]*y+y_coef[id+3]*z;	
			x_min = floor(res_x);
			y_min = floor(res_y);
			x_min_del = res_x - x_min;
			y_min_del = res_y - y_min;
			
			if(x_min>=0 && x_min<prj->X && y_min>=0 && y_min<prj->Y)//(x_min,y_min)
			{
				n = x_min + y_min*prj->X + angle*prj->X*prj->Y;
				atomicAdd(&iter_prj_divisor[n], (1-x_min_del)*(1-y_min_del)*vol_real[index]);
				atomicAdd(&iter_prj_dividend[n], (1-x_min_del)*(1-y_min_del));
			}
			if(x_min+1>=0 && x_min+1<prj->X && y_min>=0 && y_min<prj->Y)//(x_min+1,y_min)
			{
				n = (x_min+1) + y_min*prj->X + angle*prj->X*prj->Y;
				atomicAdd(&iter_prj_divisor[n], x_min_del*(1-y_min_del)*vol_real[index]);
				atomicAdd(&iter_prj_dividend[n], x_min_del*(1-y_min_del));
			}
			if(x_min>=0 && x_min<prj->X && y_min+1>=0 && y_min+1<prj->Y)//(x_min,y_min+1)
			{
				n = x_min + (y_min+1)*prj->X + angle*prj->X*prj->Y;
				atomicAdd(&iter_prj_divisor[n], (1-x_min_del)*y_min_del*vol_real[index]);
				atomicAdd(&iter_prj_dividend[n], (1-x_min_del)*y_min_del);
			}
			if(x_min+1>=0 && x_min+1<prj->X && y_min+1>=0 && y_min+1<prj->Y)//(x_min+1,y_min+1)
			{
				n = (x_min+1)+ (y_min+1)*prj->X + angle*prj->X*prj->Y;
				atomicAdd(&iter_prj_divisor[n], x_min_del*y_min_del*vol_real[index]);
				atomicAdd(&iter_prj_dividend[n], x_min_del*y_min_del);
			}
		}
	}
}

int update_head(float *vol_real,MrcHeader *head)
{
	long double sum=0,amin,amax,amean;
	int prj_size=head->nx*head->ny,i,j;
	printf("updating head(FLOAT)...\n");
	amax = amin = vol_real[0];
	for(j = 0;j<head->nz;j++)
	{
		amean = 0;
		//printf("%d :%f\n",j,vol_real[90499]);
		for(i = 0;i<prj_size;i++)
		{
			int tmp_index = i+j*prj_size;
			if(vol_real[tmp_index]>amax) amax = vol_real[tmp_index];
			if(vol_real[tmp_index]<amin) amin = vol_real[tmp_index];
			amean+=vol_real[tmp_index];
		}
		amean/=prj_size;
		sum += amean;
	}
	amean = sum/head->nz;
	head->amin=amin;
	head->amax=amax;
	head->amean=amean;
	printf("head->amin is %f, head->amax is %f, head->amean is %f\n",head->amin, head->amax, head->amean);
	return true;
}

__global__ void computePrjError(Projection *prj,float *prj_real,float *iter_prj_divisor,float *iter_prj_dividend)
{
	int y = threadIdx.x+blockIdx.x*blockDim.x;
	int z = threadIdx.y+blockIdx.y*blockDim.y;
	if(y>=prj->Y || z>=prj->AngN) return;
	int x,index;	
	for(x=0;x<prj->X;x++)
	{
		index = x+y*prj->X+z*prj->X*prj->Y;
		if(iter_prj_dividend[index]!=0)
			iter_prj_divisor[index] /= iter_prj_dividend[index];
		iter_prj_divisor[index] = prj_real[index]-iter_prj_divisor[index];
	}
}

void write_data(char *out_addr,MrcHeader *out_head,float *vol_real)
{
	clean_file(out_addr);
	FILE *out_file;
	out_file = fopen(out_addr,"r+");
	if(!out_file){
		printf("Can not open out_file!\n");
		return;	
	}
	mrc_write_head(out_file,out_head);

	//printf("siezof out_head %ld \n",sizeof(MrcHeader));
	mrc_write_all(out_file,out_head,vol_real);
	printf("%d %d %d 1\n",out_head->nx,out_head->ny,out_head->nz);
	
	//mrc_update_head(out_file);
	fclose(out_file);
	return;
}

int main(int argc,char *argv[])
{
	iStart = cpuSecond();

	iLen = atoi(argv[1]);
	SIRT_ITER_NUM = atoi(argv[2]);
	ITER_STEP_LENGTH = atof(argv[3]);
	char* in_addr = argv[4];
	char* out_addr = argv[5];
	char* angle_addr = argv[6];
	//cout<<"arg1:"<<process_num<<" arg2:"<<in_addr<<"  arg3:"<<out_addr<<"  arg4:"<<angle_addr<<endl;
	
	Volume *vol;
	Projection *prj;
	MrcHeader *in_head,*out_head;
	double *x_coef;
	double *y_coef;
	float *prj_real,*vol_real;//pri_real is inputted data ande vol_real is calculated data

/*************Head file read ande malloc space*******************/
	cudaMallocManaged((void **)&vol,sizeof(Volume));
	cudaMallocManaged((void **)&prj,sizeof(Projection));
	cudaMallocManaged((void **)&in_head,sizeof(MrcHeader));
	cudaMallocManaged((void **)&out_head,sizeof(MrcHeader));

	read_head_data(vol,prj,in_head,out_head,in_addr);
/********************************************************************/


/*************TXBR file read ande malloc space*******************/
	cudaMallocManaged((void **)&x_coef,sizeof(double)*prj->AngN*10);
	memset(x_coef, 0 , sizeof(double)*prj->AngN*10);
	//printf("%d",sizeof(double)*prj->AngN*10);
	cudaMallocManaged((void **)&y_coef,sizeof(double)*prj->AngN*10);
	memset(y_coef, 0 , sizeof(double)*prj->AngN*10);

	read_txbr_data(x_coef,y_coef,angle_addr);
/********************************************************************/


/*************Reminding data read ande malloc space*******************/
	PrjXYAngN = prj->X*prj->Y*prj->AngN;
	vol_pixel_num = vol->X*vol->Y*vol->Z;
	printf("vol_pixel_num:%d\n",vol_pixel_num);
	printf("PriXYAngN:%d\n",PrjXYAngN);
	/*for input file*/
	cudaMallocManaged((void **)&prj_real,sizeof(float)*PrjXYAngN);
	memset(prj_real, 0 , sizeof(float)*PrjXYAngN);
	/*for output file*/
	cudaMallocManaged((void **)&vol_real,sizeof(float)*vol_pixel_num);
	memset(vol_real, 0 , sizeof(float)*vol_pixel_num);	
	read_all_data(in_head,prj_real, in_addr);
/********************************************************************/

/*****************back projection (initial modle)*****************************/
	printf("%d 1\n",iLen);
	dim3 block(iLen,iLen);
	dim3 grid((vol->X+block.x-1)/block.x,(vol->Y+block.y-1)/block.y);
	dim3 grid_prj((prj->X+block.x-1)/block.x,(prj->Y+block.y-1)/block.y);
	//cudaDeviceSynchronize();
	backProjOnGPU<<<grid,block>>>(prj,vol,x_coef,y_coef,prj_real,vol_real,1);	
	cudaDeviceSynchronize();
/*****************************************************************************/

/*********************DATA space needed by SIRT***********************/
	float *iter_prj_divisor,*iter_prj_dividend;
	checkCudaErrors(cudaMallocManaged((void **)&iter_prj_divisor,sizeof(float)*PrjXYAngN));
	checkCudaErrors(cudaMallocManaged((void **)&iter_prj_dividend,sizeof(float)*PrjXYAngN));
/*********************************************************************/

/*******************SIRT************************/
	for(int i=0;i<SIRT_ITER_NUM;i++)
	{
		memset(iter_prj_divisor,0,sizeof(float)*PrjXYAngN);
		memset(iter_prj_dividend,0,sizeof(float)*PrjXYAngN);
		reProjOnGPU<<<grid,block>>>(prj,vol,x_coef,y_coef,vol_real,iter_prj_divisor,iter_prj_dividend);
		cudaDeviceSynchronize();
	/*	for(int index=0;index<PrjXYAngN;index++)
		{
			if(iter_prj_dividend[index]!=0)
				iter_prj_divisor[index] /= iter_prj_dividend[index];
			iter_prj_divisor[index] = prj_real[index]-iter_prj_divisor[index];
		}
	 */
		computePrjError<<<grid_prj,block>>>(prj,prj_real,iter_prj_divisor,iter_prj_dividend);
		cudaDeviceSynchronize();
		backProjOnGPU<<<grid,block>>>(prj,vol,x_coef,y_coef,iter_prj_divisor,vol_real,ITER_STEP_LENGTH);
		cudaDeviceSynchronize();
		printf("Iteration %d finished..\n",i);
	}
/**********************************************/

	update_head(vol_real,out_head);
	write_data(out_addr,out_head,vol_real);
	
	cudaDeviceReset();//重置CUDA设备释放程序占用的资源

	iElaps = cpuSecond()-iStart;
	printf("Host time elapsed:%lfsec\n",iElaps);
	return 0;
}

