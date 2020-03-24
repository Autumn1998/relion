
#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>
//#include <opencv.hpp>
#include "mrcfile_atom.h"
#include "atom.h"
//#include "mrcfiles.h"
//#include "mrcslice.h"

//#define PI 3.1415927f
CvMat *create_CvMat_32F_from_float_data(int rows, int cols, float *data);
int RamLak(int width, float *ramlak);
int SheppLogan(int length,  float *shepp);


int max(int i, int j);
void Show_Mat2D(CvMat *mat,int row,int col);
void symmetrize_IplImage_values_CvRect(IplImage *image, CvRect rect);

int filter_prj_sym(float *prj_real, char *filter, int filtlength, Projection prj, int row_pad, int col_pad, int symmetrize_2D_flag, int angle, int id);


int filter_prj_X_sym(float *prj_real, char *filter, int filtlength, Projection prj, int row_pad, int col_pad, int symmetrize_2D_flag, int angle, int id);

