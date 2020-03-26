#ifndef MRCFILES_H
#define MRCFILES_H

#define TEX_LINE_MAX 500
#define HEAD_SIZE 1024
#define INF 99999999999999

#define PI_180 0.01745329252f
#ifndef PI
#define     PI  3.14159265358979323846
#endif
#define D2R(__ANGLE__) ((__ANGLE__) * PI_180) // PI / 180

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif



#define MRC_MODE_BYTE          0
#define MRC_MODE_SHORT         1
#define MRC_MODE_FLOAT         2
#define MRC_MODE_COMPLEX_SHORT 3
#define MRC_MODE_COMPLEX_FLOAT 4
#define MRC_MODE_USHORT        6
#define MRC_MODE_RGB           16


#define MRC_LABEL_SIZE         80
#define MRC_NEXTRA             16
#define MRC_NLABELS            10
#define MRC_HEADER_SIZE        1024   /* Length of Header is 1024 Bytes. */
#define MRC_MAXCSIZE           3

#pragma pack(1)
typedef struct MRCheader
{
  int   nx;         /*  # of Columns                  */
  int   ny;         /*  # of Rows                     */
  int   nz;         /*  # of Sections.                */
  int   mode;       /*  given by #define MRC_MODE...  */

  int   nxstart;    /*  Starting point of sub image.  */
  int   nystart;
  int   nzstart;

  int   mx;         /* Grid size in x, y, and z       */
  int   my;
  int   mz;

  float   xlen;       /* length of x element in um.     */
  float   ylen;       /* get scale = xlen/nx ...        */
  float   zlen;

  float   alpha;      /* cell angles, ignore */
  float   beta;
  float   gamma;

  int   mapc;       /* map coloumn 1=x,2=y,3=z.       */
  int   mapr;       /* map row     1=x,2=y,3=z.       */
  int   maps;       /* map section 1=x,2=y,3=z.       */

  float   amin;
  float   amax;
  float   amean;
  
  short   ispg;       /* image type */
  short   nsymbt;     /* space group number */


  /* 64 bytes */

  int   next;
  short   creatid;  /* Creator id, hvem = 1000, DeltaVision = -16224 */

  
  char    blank[30];
  
  short   nint;
  short   nreal;

  short   sub;
  short   zfac;

  float   min2;
  float   max2;
  float   min3;
  float   max3;
  float   min4;
  float   max4;


  short   idtype;
  short   lens;
  short   nd1;     /* Devide by 100 to get float value. */
  short   nd2;
  short   vd1;
  short   vd2;
  float   tiltangles[6];  /* 0,1,2 = original:  3,4,5 = current */


  float   xorg;
  float   yorg;
  float   zorg;
  char    cmap[4];
  char    stamp[4];
  float   rms;

  int nlabl;
  char  labels[10][80];


} MrcHeader;
#pragma pack()
/* END_CODE */



/******************************** Header functions or useful io functions **************************/
int mrc_init_head(MrcHeader *head);

int mrc_read_head (FILE *fin,  MrcHeader *head);

int read_coef(double *x_coef, double *y_coef, FILE *f_coef);

int mrc_read_all(FILE *fin,MrcHeader *head,float *mrc_data_all);

int mrc_write_all(FILE *fout,MrcHeader *head,int data_length, float *mrc_data_all);

int mrc_update_head(FILE *fout);

int mrc_write_head(FILE *fout, MrcHeader *head);
#endif

