#ifndef MRCFILES_H
#define MRCFILES_H

#ifndef PI
#define     PI  3.14159265358979323846
#endif

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>


//#include <dirent.h>
//#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#ifndef FALSE
#define FALSE       1           /*false for boolean*/
#endif
#ifndef TRUE
#define TRUE        0           /*true for boolean*/
#endif

#ifndef false
#define false       1
#endif
#ifndef true
#define true        0
#endif

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



/*typedef struct  //complex floating number
{
  float a;
  float b;

} ComplexFloat;

typedef struct  //complex short number
{
  short a;
  short b;

} ComplexShort;
*/

/* The header structure for MRC files */
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

int mrc_read_head_MPI(MPI_File fin,  MrcHeader *head);

int mrc_read_head (FILE *fin,  MrcHeader *head);

int mrc_write_head(MPI_File fout, MrcHeader *head);

int mrc_update_head(MPI_File fout);

int mrc_read_block(MPI_File fin, MrcHeader  *head, int start, int end, char axis,float *mrc_data_block);

int mrc_read_all_MPI(MPI_File fin, MrcHeader  *head, float *mrc_data_all);

int mrc_read_all(FILE *fin, MrcHeader  *head, float *mrc_data_all);

int mrc_write_all(MPI_File fout, MrcHeader  *head,  int Z_start,int Z_end, float *mrc_data_all);

/////new for PUMACE
int mrc_write_slice(MPI_File fout, MrcHeader  *head, int slcN,char axis,float *slcdata);

int mrc_read_slice(MPI_File fin, MrcHeader  *head, int slcN, char axis,float *mrc_data_slice);

#endif

