#ifndef ATOM_H
#define ATOM_H


#include "mrcfile_atom.h"

#ifndef TRUE
#define TRUE  0
#define FALSE 1
#endif

#ifndef TEXT_LINE_MAX
#define TEXT_LINE_MAX 500
#endif



/*Volume stores the parameter(not the data) of a volum data*/
typedef struct
{
  int Xstart;
  int Xend;
  int Ystart;
  int Yend;
  int Zstart;
  int Zend;

  int X;                       //X,Y,Z are the pixel size of the three dimesions in a volum data;
  int Y;
  int Z;                               //it also equals the thickness of a volum data;

} Volume;

/*Proj stores the parameter(not the data) of a projection data*/
typedef struct
{
  int X;                       //X,Y,Z are the pixel size of the three dimesions in a volum data;
  int Y;
  int AngN;                               //it also equals the thickness of a volum data;

} Projection;

/*Pixel is the coordinate number of a 3d map*/
typedef struct
{
  int    X;
  int    Y;
  int    Z;

} Pixel;

/*computing proj by the coordinate of a 3D pixel*/
typedef struct
{
  int    x_min;//x coordinate of the proj
  int    y_min;//y coordinate of the proj
  
  double x_min_del;
  double y_min_del; //weight of the proj

} Weight;


int  ATOM_MPI_SIRT(char *outf,char *fbpf, float ATOM_ITR_STEP, int ATOM_ITR_NUM,int myid,int mypro,int flagpart);

#endif
