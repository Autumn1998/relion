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

	printf("%d 0\n",vol->Z);
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


__global__ void backProjOnGPU(Projection *prj,Volume *vol,double *x_coef,double *y_coef,float *prj_real,float *vol_real,int *test)
{
	double divisor;//分子
	double dividend;//分母
	int z = threadIdx.x+blockIdx.x*blockDim.x + vol->Zstart;
	if(z!=vol->Zstart) return;
	if(z>=vol->Zend) return;
	int x,y,index,angle,n;	
	test[z] = 1;
	for(y=vol->Ystart;y<vol->Ystart+vol->Y;y++)
	{
		for(x=vol->Xstart;x<0;x++)
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
				if(x_min>=0 && x_min<prj->X && y_min>=0 && y_min<prj->Y)//(x_min+1,y_min)
				{
					n = (x_min+1) + y_min*prj->X + angle*prj->X*prj->Y;
					divisor += x_min_del*(1-y_min_del)*prj_real[n];
					dividend += x_min_del*(1-y_min_del);
				}
				if(x_min>=0 && x_min<prj->X && y_min>=0 && y_min<prj->Y)//(x_min,y_min+1)
				{
					n = x_min + (y_min+1)*prj->X + angle*prj->X*prj->Y;
					divisor += (1-x_min_del)*y_min_del*prj_real[n];
					dividend += (1-x_min_del)*y_min_del;
				}
				if(x_min>=0 && x_min<prj->X && y_min>=0 && y_min<prj->Y)//(x_min+1,y_min+1)
				{
					n = (x_min+1)+ (y_min+1)*prj->X + angle*prj->X*prj->Y;
					divisor += x_min_del*y_min_del*prj_real[n];
					dividend += x_min_del*y_min_del;
				}
			}
			index = (x-vol->Xstart)+(y-vol->Ystart)*vol->X+(z-vol->Zstart)*vol->X*vol->Y;
			printf("%d\n",index);
			if(dividend!=0.0f)
			{
				vol_real[index] = (float)(divisor/dividend);
				//printf("vol_read[%d]:%f\n",index,vol_real[index]);
			}
		}
		
	}
	test[z] = 2;	
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

__global__ void testOnGPU(int *test)
{
	int z = threadIdx.x+blockIdx.x*blockDim.x;
	if(z>=58) return;
	test[z] = 1;
}

int main(int argc,char *argv[])
{
	iStart = cpuSecond();

	iLen = atoi(argv[1]);
	char* in_addr = argv[2];
	char* out_addr = argv[3];
	char* angle_addr = argv[4];
	//cout<<"arg1:"<<process_num<<" arg2:"<<in_addr<<"  arg3:"<<out_addr<<"  arg4:"<<angle_addr<<endl;
	
	Volume *vol;
	Projection *prj;
	MrcHeader *in_head,*out_head;
	double *x_coef;
	double *y_coef;
	float *prj_real,*vol_real;//pri_real is inputted data ande vol_real is calculated data
	int *test;

/*************Head file read ande malloc space*******************/
	cudaMallocManaged((void **)&vol,sizeof(Volume));
	cudaMallocManaged((void **)&prj,sizeof(Projection));
	cudaMallocManaged((void **)&in_head,sizeof(MrcHeader));
	cudaMallocManaged((void **)&out_head,sizeof(MrcHeader));

	read_head_data(vol,prj,in_head,out_head,in_addr);

/********************************************************************/


/*************TXBR file read ande malloc space*******************/
	cudaMallocManaged((void **)&test,sizeof(int)*vol->Z);
	memset(test,0,sizeof(int)*vol->Z);
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
	/*for input file*/
	cudaMallocManaged((void **)&prj_real,sizeof(float)*PrjXYAngN);
	memset(prj_real, 0 , sizeof(float)*PrjXYAngN);
	/*for output file*/
	cudaMallocManaged((void **)&vol_real,sizeof(float)*vol_pixel_num);
	memset(vol_real, 0 , sizeof(float)*vol_pixel_num);	
	read_all_data(in_head,prj_real, in_addr);
/********************************************************************/

	printf("%d 1\n",vol->Z);
	//for(int i=0;i<vol->Z;i++ ) printf("%d",test[i]);
	dim3 block(iLen);
	dim3 grid((vol->Z+block.x-1)/block.x);
	//cudaDeviceSynchronize();
	backProjOnGPU<<<grid,block>>>(prj,vol,x_coef,y_coef,prj_real,vol_real,test);	
	//testOnGPU<<<grid,block>>>(test);
	cudaDeviceSynchronize();
	
	//write_data(out_addr,out_head,vol->Zend-vol->Zstart,vol_real);
	
	printf("%d 2\n",vol->Z);
	for(int qwe = 0;qwe<vol->Z;qwe++) printf("%d",test[qwe]);
	//cudaDeviceReset();//重置CUDA设备释放程序占用的资源

	iElaps = cpuSecond()-iStart;
	printf("Host time elapsed:%lfsec\n",iElaps);
	return 0;
}

