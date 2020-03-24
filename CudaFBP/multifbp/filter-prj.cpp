//gcc -I /home/x5wan/software/OpenCV-2.2.0/include -o opencv_test3 opencv_test3.c -L /home/x5wan/software/OpenCV-2.2.0/lib -lopencv_core -lopencv_highgui -lopencv_imgproc

#include "filter-prj.h"

///////////////////////////////// OpenCV utilities: CvMat.
CvMat *create_CvMat_32F_from_float_data(int rows, int cols, float *data)
{
    CvMat *matrix = cvCreateMatHeader(rows, cols, CV_32FC1);

    int n_bytes_row = cols * sizeof(float);
    cvSetData(matrix, data, n_bytes_row);
    
    return matrix;
}


/////////////////////////////
void write_MRC_image_from_CvMat(CvMat *Mat, char *filepath, int angle)
{
    int i,j;
    
    MPI_File fout;
    MPI_File_open(MPI_COMM_WORLD, filepath, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fout);

    MrcHeader  *outhead;
   
    outhead=(MrcHeader *)malloc(sizeof(MrcHeader));
    if (angle==0)
        printf("Mat->cols is %d, Mat->rows is %d\n",Mat->cols, Mat->rows);


    mrc_init_head(outhead);
    outhead->nx=Mat->cols;
    outhead->ny=Mat->rows;
    outhead->nz=61;

    if(angle==0)    
        mrc_write_head(fout,outhead);

    int line_num=Mat->cols*Mat->rows;

    float *prj_slice;
    if((prj_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
    {
       printf("Error with Function 'filter_prj()'!Can't malloc memery for 'prj_slice'!");
       exit(1);
    }
    memset(prj_slice, 0 , sizeof(float)*line_num);



  
    for(j=0;j<Mat->rows;j++)
    
        for(i=0;i<Mat->cols;i++)
            prj_slice[i+j*Mat->cols]=cvGetReal2D(Mat,j,i);
           //prj_slice[i+j*Mat->cols]=0;

     

    mrc_write_slice(fout, outhead, angle, 'Z', prj_slice);

    free(outhead); 
    outhead=NULL;

    MPI_File_close(&fout);
    free(prj_slice);
    

}


/////////////////////////////
void write_MRC_image_from_IplImage(IplImage *image, char *filepath, int angle)
{
    int i,j;
    
    MPI_File fout;
    MPI_File_open(MPI_COMM_WORLD, filepath, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fout);

    MrcHeader  *outhead;
   
    outhead=(MrcHeader *)malloc(sizeof(MrcHeader));
    mrc_init_head(outhead);
    outhead->nx=image->width;
    outhead->ny=image->height;
    outhead->nz=61;
    outhead->xorg=0.0;
    outhead->yorg=0.0;
    outhead->zorg=0.0;

    if(angle==0)
    {
        mrc_write_head(fout,outhead);
        printf("outhead->nx is %d, outhead->ny is %d, outhead->nz is %d\n",outhead->nx, outhead->ny,outhead->nz);
        
    }
 

    unsigned int line_num=image->width*image->height;



    float *prj_slice;
    if((prj_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
    {
       printf("Error with Function 'filter_prj()'!Can't malloc memery for 'prj_slice'!");
       exit(1);
    }
    memset(prj_slice, 0 , sizeof(float)*line_num);


  /*  float *p_image;
  
    for(j=0;j<image->height;j++)
    {
        p_image = (float *) (image->imageData+j*image->widthStep);
        for(i=0;i<image->width;i++)
        {
            prj_slice[i+j*image->width]=*p_image++;
        }
    }*/

    for(j=0;j<image->height;j++)
    {
               
        for(i=0;i<image->width;i++)
            prj_slice[i+j*image->width]=cvGetReal2D(image,j,i);

    }


    mrc_write_slice(fout, outhead, angle, 'Z', prj_slice);


     free(outhead); 
    outhead=NULL;

    MPI_File_close(&fout);

    free(prj_slice);



}



/////////////////////////////RamLak
int RamLak(int width, float *ramlak)
{
  int n;
  n=-(width-1)/2;
  
  int i;
  for(i=0;i<width;i++)
    {
       if (n==0)
          ramlak[i]=M_PI/4;
       else if (abs(n%2)==1)
                ramlak[i]=-1.0/(M_PI*n*n);
            else
                ramlak[i]=0;
       n++;
     }
   return TRUE;
}

//////////////////////////////SheppLogan
int SheppLogan(int length,  float *shepp)
{
    if (length % 2 != 1)
    {
        printf("length %d is not Odd\n", length);
        exit(1);
    }

    double n;
    n= -((double) length - 1.0) / 2.0;

    int i;
    for(i=0;i<length;i++)
    {
       shepp[i]=(float)(-2.0/(M_PI*(4.0*n*n-1.0)));
       n = n+1.0;

     //  printf("shepp[%d] is %f\n",i,shepp[i]);
    }

    return TRUE;
}

//////////////////////////////////
int max(int i, int j)
{ 
  if (i>=j) return i;
  else return j;
}

double lg2(double n)
{
   return log(n)/log(2);
}

/*****************************************************************************************************/      
void Show_Mat2D(CvMat *mat,int row,int col)
{
 int i,j;
 for(i=0;i<row;i++)
 {
  for(j=0;j<col;j++)
  {
   printf("%lf  ",cvGet2D(mat,i,j).val[0]);
  }
  printf("\n");
 }
}

/*****************************************************************************************************/      
void symmetrize_IplImage_values_CvRect(IplImage *image, CvRect rect)
{
    if ((image->depth != IPL_DEPTH_32F) && (image->depth != IPL_DEPTH_64F))
    {
        printf("Expected image->depth = %i (IPL_DEPTH_32F) or %i (IPL_DEPTH_64F), but image->depth = %i.\n",
            IPL_DEPTH_32F, IPL_DEPTH_64F, image->depth);
        exit(1);
    }

    //write_MRC_image_from_IplImage(image, "/home/akulo/filter_1D_images/image.mrc");

    int n_x = image->width;
    int n_y = image->height;

    CvMat matrix_src;
    CvMat matrix_dst;

    // BEGIN: Top center.
    // NOTE: A flip about x.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, (rect.y + rect.height) - (n_y - (rect.y + rect.height)), rect.width, (n_y - (rect.y + rect.height))));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, rect.y + rect.height, rect.width, (n_y - (rect.y + rect.height))));
    cvFlip(&matrix_src, &matrix_dst, 0);
    // END: Top center.

    // BEGIN: Bottom center.
    // NOTE: A flip about x.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, rect.y, rect.width, rect.y));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x, 0, rect.width, rect.y));
    cvFlip(&matrix_src, &matrix_dst, 0);
    // END: Bottom center.

    // BEGIN: Left top, center, and bottom.
    // NOTE: A flip about y.
    cvGetSubRect(image, &matrix_src, cvRect(rect.x, 0, rect.x, n_y));
    cvGetSubRect(image, &matrix_dst, cvRect(0, 0, rect.x, n_y));
    cvFlip(&matrix_src, &matrix_dst, 1);
    // END: Left top, center, and bottom.

    // BEGIN: Right top, center, and bottom.
    // NOTE: A flip about y.
    cvGetSubRect(image, &matrix_src, cvRect((rect.x + rect.width) - (n_x - (rect.x + rect.width)), 0, n_x - (rect.x + rect.width), n_y));
    cvGetSubRect(image, &matrix_dst, cvRect(rect.x + rect.width, 0, n_x - (rect.x + rect.width), n_y));
    cvFlip(&matrix_src, &matrix_dst, 1);
    // END: Right top, center, and bottom.
}


/*****************************************************************************************************/
int filter_prj_sym(float *prj_real, char *filter, int filtlength, Projection prj, int row_pad, int col_pad, int symmetrize_2D_flag, int angle, int id)
{
 unsigned int i,j;

 int n_x_pad=prj.X+row_pad; int n_y_pad=prj.Y+col_pad;

 int n_x_pad_conv=filtlength+n_x_pad-1;
 int n_x_pad_conv_opt=cvGetOptimalDFTSize(n_x_pad_conv);

 int n_y_pad_opt=cvGetOptimalDFTSize(n_y_pad);

 if(id==0 && angle==0)
 printf("n_x_pad is %d, n_y_pad is %d, n_x_pad_conv is %d, n_x_pad_conv_opt is %d, n_y_pad_opt is %d\n", n_x_pad, n_y_pad, n_x_pad_conv, n_x_pad_conv_opt, n_y_pad_opt);


 int line_num=prj.X*prj.Y;
 float *prj_slice;
 if((prj_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
   {
       printf("Error with Function 'filter_prj()'!Can't malloc memery for 'prj_slice'!");
       return FALSE;
   }
 memset(prj_slice, 0 , sizeof(float)*line_num);

 for(i=0;i<line_num;i++)
     prj_slice[i]=prj_real[i];

 // original projection (prj.Y*prj.X)
 IplImage *projection_org = cvCreateImageHeader(cvSize(prj.X, prj.Y), IPL_DEPTH_32F,1); 

 int n_bytes_row = prj.X*sizeof(float);
 cvSetData(projection_org, prj_slice, n_bytes_row);

// DEBUG:	write_MRC_image_from_IplImage(projection_org, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_org.mrc");
// printf("projection_org is done\n");

 // pad projection (n_y_pad*n_x_pad)
 IplImage *projection_pad = cvCreateImage(cvSize(n_x_pad, n_y_pad), IPL_DEPTH_32F,1); 
 cvSetZero(projection_pad);

 CvRect projection_pad_center_ROI = cvRect(row_pad/2, col_pad/2, prj.X, prj.Y);

 cvSetImageROI(projection_org, cvRect(0, 0, prj.X, prj.Y));
 cvSetImageROI(projection_pad, projection_pad_center_ROI); // NOTE: The extents of this ROI are the same as image0's ROI.
 cvCopy(projection_org, projection_pad, NULL);
 cvResetImageROI(projection_org);
 cvResetImageROI(projection_pad);


 
// DEBUG:	
/* if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(projection_pad, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }*/

/* 
  //DEBUG:
 CvMat *Matrix1;
 Matrix1=cvCreateMatHeader(n_y_pad,n_x_pad,CV_32FC1);
 cvGetMat(projection_pad,Matrix1, NULL,0);

if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(Matrix1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }*/




 
 //padding projection_pad using symmetrize
 if (symmetrize_2D_flag)
     symmetrize_IplImage_values_CvRect(projection_pad, projection_pad_center_ROI);

/* // DEBUG:	
 if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(projection_pad, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }

 
  //DEBUG:
 CvMat *Matrix1;
 Matrix1=cvCreateMatHeader(n_y_pad,n_x_pad,CV_32FC1);
 cvGetMat(projection_pad,Matrix1, NULL,0);

if (id ==0) 
 {
     printf("Begin to write projection_pad_mat.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(Matrix1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad_mat.mrc",angle);
 }*/




 //filter
 float *filter_rs = (float *) malloc(filtlength * sizeof(float));  //filter_rs (1*251)
    if (filter_rs == NULL)
    {
        printf("Memory request for filter_rs failed .\n");
        exit(1);
    }
    
 if (!strcmp(filter,"RamLak"))
      RamLak(filtlength, filter_rs);
 else if (!strcmp(filter,"SheppLogan"))
      SheppLogan(filtlength,filter_rs); 
      else
          printf("Invalid filter selected.\n");

// filter_rs_cv (1*251)
 CvMat *filter_rs_cv = create_CvMat_32F_from_float_data(1, filtlength, filter_rs);  

//Debug:
/* Show_Mat2D(filter_rs_cv,1,filtlength);

 exit(0);*/

// filter_fs_cv (1*n_x_pad_cov_opt)
 CvMat *filter_fs_cv = cvCreateMat(1, n_x_pad_conv_opt, CV_32FC1);
 cvZero(filter_fs_cv);

 CvMat tmp;
 cvGetSubRect(filter_fs_cv, &tmp, cvRect(0, 0, filtlength, 1));
 cvCopy(filter_rs_cv, &tmp, NULL);

 cvDFT(filter_fs_cv, filter_fs_cv, CV_DXT_FORWARD | CV_DXT_ROWS, 1);

// filter_fs_cv (n_y_pad_opt*n_x_pad_cov_opt)
 CvMat *filter_2D_fs_cv= cvCreateMat(n_y_pad_opt, n_x_pad_conv_opt,CV_32FC1);
 cvZero(filter_2D_fs_cv);

 cvRepeat(filter_fs_cv,filter_2D_fs_cv);


//DEBUG:
/*IplImage *filter_2D_fs_cv_img=cvCreateImage(cvSize(n_x_pad_conv_opt, n_y_pad_opt), IPL_DEPTH_32F,1); 
cvSetZero(filter_2D_fs_cv_img);*/
/* IplImage filter_2D_fs_cv_img;

 if (id==0)
 {
     printf("Begin to write filter_2D_fs_cv for %d angle\n",angle);     
     cvGetImage (filter_2D_fs_cv, &filter_2D_fs_cv_img);
    // write_MRC_image_from_CvMat(filter_2D_fs_cv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/filter_2D_fs_cv.mrc",angle);
     write_MRC_image_from_IplImage(&filter_2D_fs_cv_img, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/filter_2D_fs_cv_img.mrc",angle); 
 }*/

 //convolution projection_pad and filter_2D_fs_cv into projection_conv
 CvMat *projection_conv=cvCreateMat(n_y_pad_opt, n_x_pad_conv_opt, CV_32FC1);
 cvZero(projection_conv);

 

 CvMat matrix_src;
 CvMat matrix_dst;

 cvGetSubRect(projection_conv, &matrix_dst, cvRect(0, 0, n_x_pad, n_y_pad));
 cvCopy(projection_pad, &matrix_dst, NULL);
// cvGetSubRect(projection_conv, &matrix_dst, cvRect(n_x_pad, 0, n_x_pad_conv_opt-n_x_pad, n_y_pad));

/* // DEBUG:	
 IplImage *Image1;
 Image1=cvCreateImageHeader(cvSize(n_x_pad_conv_opt,n_y_pad_opt), IPL_DEPTH_32F,1);
 cvGetImage(projection_conv,Image1);
 
 if (id ==0) 
 {
     printf("Begin to write projection_conv_img.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(Image1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv_img.mrc",angle);
 }

 
  //DEBUG:

 if (id ==0) 
 {
     printf("Begin to write projection_conv.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(projection_conv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv.mrc",angle);
 }
*/

 if (symmetrize_2D_flag)
 {
   /*  cvGetSubRect(projection_pad, &matrix_src, cvRect(n_x_pad - (n_x_pad_conv_opt - n_x_pad), 0, n_x_pad_conv_opt - n_x_pad, n_y_pad));
     cvFlip(&matrix_src, &matrix_dst, 1);// END: Symmetrize right.  Only a single mirrored copy is necessary.*/

          
       
     int i,j;
     if(n_x_pad>=n_x_pad_conv_opt-n_x_pad)
     {
         for(j=0;j<n_y_pad_opt;j++)
             for(i=n_x_pad;i<n_x_pad_conv_opt;i++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,j,2*n_x_pad-1-i));
     }
     else
     {
         
         for(j=0;j<n_y_pad_opt;j++)
             for(i=n_x_pad;i<2*n_x_pad;i++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,j,2*n_x_pad-1-i));
     }

     
 }

/*// DEBUG:	
 IplImage *Image1;
 Image1=cvCreateImageHeader(cvSize(n_x_pad_conv_opt,n_y_pad_opt), IPL_DEPTH_32F,1);
 cvGetImage(projection_conv,Image1);
 
 if (id ==0) 
 {
     printf("Begin to write projection_conv_img.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(Image1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv_img.mrc",angle);
 }

 
  //DEBUG:

 if (id ==0) 
 {
     printf("Begin to write projection_conv.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(projection_conv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv.mrc",angle);
 }

 */
 


 cvDFT(projection_conv, projection_conv, CV_DXT_FORWARD | CV_DXT_ROWS, n_y_pad);

 cvMulSpectrums(projection_conv, filter_2D_fs_cv, projection_conv, CV_DXT_ROWS);

 cvDFT(projection_conv, projection_conv, CV_DXT_INV_SCALE | CV_DXT_ROWS, n_y_pad);
 
 

 cvGetSubRect(projection_conv, &matrix_src, cvRect(filtlength/2, 0, projection_pad->width, projection_pad->height));
 cvCopy(&matrix_src, projection_pad, NULL);



 float *p_image;
  
 for(j=0;j<projection_pad_center_ROI.height;j++)
 {
     p_image = (float *) (projection_pad->imageData+projection_pad_center_ROI.x*sizeof(float)+(projection_pad_center_ROI.y+j)*projection_pad->widthStep);
     
     for(i=0;i<projection_pad_center_ROI.width;i++)
     {
         prj_slice[i+j*prj.X]=*p_image++;
     }
 }

 
          
 for(i=0;i<line_num;i++)
     prj_real[i]=prj_slice[i];


 free(prj_slice);

 
 cvReleaseImageHeader(&projection_org);
 cvReleaseImage(&projection_pad);

 
 free(filter_rs);

 cvReleaseMat(&filter_rs_cv);
 cvReleaseMat(&filter_fs_cv);
 cvReleaseMat(&filter_2D_fs_cv);
 cvReleaseMat(&projection_conv);


return TRUE;
}






/*****************************************************************************************************/
int filter_prj_X_sym(float *prj_real, char *filter, int filtlength, Projection prj, int row_pad, int col_pad, int symmetrize_2D_flag, int angle, int id)
{
 unsigned int i,j;

 int n_x_pad=prj.X+row_pad; int n_y_pad=prj.Y+col_pad;

 int n_x_pad_conv=filtlength+n_x_pad-1;
 int n_x_pad_conv_opt=cvGetOptimalDFTSize(n_x_pad_conv);

 int n_y_pad_opt=cvGetOptimalDFTSize(n_y_pad);

 if(id==0 && angle==0)
 printf("n_x_pad is %d, n_y_pad is %d, n_x_pad_conv is %d, n_x_pad_conv_opt is %d, n_y_pad_opt is %d\n", n_x_pad, n_y_pad, n_x_pad_conv, n_x_pad_conv_opt, n_y_pad_opt);


 int line_num=prj.X*prj.Y;
 float *prj_slice;
 if((prj_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
   {
       printf("Error with Function 'filter_prj()'!Can't malloc memery for 'prj_slice'!");
       return FALSE;
   }
 memset(prj_slice, 0 , sizeof(float)*line_num);

/* for(i=0;i<line_num;i++)
     prj_slice[i]=prj_real[i+angle*line_num];*/

 for(i=0;i<prj.X;i++)
     for(j=0;j<prj.Y;j++)
         prj_slice[j+i*prj.Y]=prj_real[i+j*prj.X];

 // original projection (prj.Y*prj.X)
 IplImage *projection_org = cvCreateImageHeader(cvSize(prj.X, prj.Y), IPL_DEPTH_32F,1); 

 int n_bytes_row = prj.X*sizeof(float);
 cvSetData(projection_org, prj_slice, n_bytes_row);

// DEBUG:	write_MRC_image_from_IplImage(projection_org, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_org.mrc");
// printf("projection_org is done\n");

 // pad projection (n_y_pad*n_x_pad)
 IplImage *projection_pad = cvCreateImage(cvSize(n_x_pad, n_y_pad), IPL_DEPTH_32F,1); 
 cvSetZero(projection_pad);

 CvRect projection_pad_center_ROI = cvRect(row_pad/2, col_pad/2, prj.X, prj.Y);

 cvSetImageROI(projection_org, cvRect(0, 0, prj.X, prj.Y));
 cvSetImageROI(projection_pad, projection_pad_center_ROI); // NOTE: The extents of this ROI are the same as image0's ROI.
 cvCopy(projection_org, projection_pad, NULL);
 cvResetImageROI(projection_org);
 cvResetImageROI(projection_pad);


 
// DEBUG:	
/* if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(projection_pad, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }*/

/* 
  //DEBUG:
 CvMat *Matrix1;
 Matrix1=cvCreateMatHeader(n_y_pad,n_x_pad,CV_32FC1);
 cvGetMat(projection_pad,Matrix1, NULL,0);

if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(Matrix1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }*/




 
 //padding projection_pad using symmetrize
 if (symmetrize_2D_flag)
     symmetrize_IplImage_values_CvRect(projection_pad, projection_pad_center_ROI);

/* // DEBUG:	
 if (id ==0) 
 {
     printf("Begin to write projection_pad.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(projection_pad, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad.mrc",angle);
 }

 
  //DEBUG:
 CvMat *Matrix1;
 Matrix1=cvCreateMatHeader(n_y_pad,n_x_pad,CV_32FC1);
 cvGetMat(projection_pad,Matrix1, NULL,0);

if (id ==0) 
 {
     printf("Begin to write projection_pad_mat.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(Matrix1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_pad_mat.mrc",angle);
 }*/




 //filter
 float *filter_rs = (float *) malloc(filtlength * sizeof(float));  //filter_rs (1*251)
    if (filter_rs == NULL)
    {
        printf("Memory request for filter_rs failed .\n");
        exit(1);
    }
    
 if (!strcmp(filter,"RamLak"))
      RamLak(filtlength, filter_rs);
 else if (!strcmp(filter,"SheppLogan"))
      SheppLogan(filtlength,filter_rs); 
      else
          printf("Invalid filter selected.\n");

// filter_rs_cv (1*251)
 CvMat *filter_rs_cv = create_CvMat_32F_from_float_data(1, filtlength, filter_rs);  

//Debug:
/* Show_Mat2D(filter_rs_cv,1,filtlength);

 exit(0);*/

// filter_fs_cv (1*n_x_pad_cov_opt)
 CvMat *filter_fs_cv = cvCreateMat(1, n_x_pad_conv_opt, CV_32FC1);
 cvZero(filter_fs_cv);

 CvMat tmp;
 cvGetSubRect(filter_fs_cv, &tmp, cvRect(0, 0, filtlength, 1));
 cvCopy(filter_rs_cv, &tmp, NULL);

 cvDFT(filter_fs_cv, filter_fs_cv, CV_DXT_FORWARD | CV_DXT_ROWS, 1);

// filter_fs_cv (n_y_pad_opt*n_x_pad_cov_opt)
 CvMat *filter_2D_fs_cv= cvCreateMat(n_y_pad_opt, n_x_pad_conv_opt,CV_32FC1);
 cvZero(filter_2D_fs_cv);

 cvRepeat(filter_fs_cv,filter_2D_fs_cv);


//DEBUG:
/*IplImage *filter_2D_fs_cv_img=cvCreateImage(cvSize(n_x_pad_conv_opt, n_y_pad_opt), IPL_DEPTH_32F,1); 
cvSetZero(filter_2D_fs_cv_img);*/
/* IplImage filter_2D_fs_cv_img;

 if (id==0)
 {
     printf("Begin to write filter_2D_fs_cv for %d angle\n",angle);     
     cvGetImage (filter_2D_fs_cv, &filter_2D_fs_cv_img);
    // write_MRC_image_from_CvMat(filter_2D_fs_cv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/filter_2D_fs_cv.mrc",angle);
     write_MRC_image_from_IplImage(&filter_2D_fs_cv_img, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/filter_2D_fs_cv_img.mrc",angle); 
 }*/

 //convolution projection_pad and filter_2D_fs_cv into projection_conv
 CvMat *projection_conv=cvCreateMat(n_y_pad_opt, n_x_pad_conv_opt, CV_32FC1);
 cvZero(projection_conv);

 

 CvMat matrix_src;
 CvMat matrix_dst;

 cvGetSubRect(projection_conv, &matrix_dst, cvRect(0, 0, n_x_pad, n_y_pad));
 cvCopy(projection_pad, &matrix_dst, NULL);
// cvGetSubRect(projection_conv, &matrix_dst, cvRect(n_x_pad, 0, n_x_pad_conv_opt-n_x_pad, n_y_pad));

/* // DEBUG:	
 IplImage *Image1;
 Image1=cvCreateImageHeader(cvSize(n_x_pad_conv_opt,n_y_pad_opt), IPL_DEPTH_32F,1);
 cvGetImage(projection_conv,Image1);
 
 if (id ==0) 
 {
     printf("Begin to write projection_conv_img.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(Image1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv_img.mrc",angle);
 }

 
  //DEBUG:

 if (id ==0) 
 {
     printf("Begin to write projection_conv.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(projection_conv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv.mrc",angle);
 }
*/

 if (symmetrize_2D_flag)
 {
   /*  cvGetSubRect(projection_pad, &matrix_src, cvRect(n_x_pad - (n_x_pad_conv_opt - n_x_pad), 0, n_x_pad_conv_opt - n_x_pad, n_y_pad));
     cvFlip(&matrix_src, &matrix_dst, 1);// END: Symmetrize right.  Only a single mirrored copy is necessary.*/

          
       
     int i,j;
     if(n_x_pad>=n_x_pad_conv_opt-n_x_pad)
     {
         for(j=0;j<n_y_pad_opt;j++)
             for(i=n_x_pad;i<n_x_pad_conv_opt;i++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,j,2*n_x_pad-1-i));
     }
     else
     {
         
         for(j=0;j<n_y_pad_opt;j++)
             for(i=n_x_pad;i<2*n_x_pad;i++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,j,2*n_x_pad-1-i));
     }

     
 }

/*// DEBUG:	
 IplImage *Image1;
 Image1=cvCreateImageHeader(cvSize(n_x_pad_conv_opt,n_y_pad_opt), IPL_DEPTH_32F,1);
 cvGetImage(projection_conv,Image1);
 
 if (id ==0) 
 {
     printf("Begin to write projection_conv_img.mrc for %d angle\n",angle);
     write_MRC_image_from_IplImage(Image1, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv_img.mrc",angle);
 }

 
  //DEBUG:

 if (id ==0) 
 {
     printf("Begin to write projection_conv.mrc for %d angle\n",angle);
     write_MRC_image_from_CvMat(projection_conv, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_conv.mrc",angle);
 }

 */
 


 cvDFT(projection_conv, projection_conv, CV_DXT_FORWARD | CV_DXT_ROWS, n_y_pad);

 cvMulSpectrums(projection_conv, filter_2D_fs_cv, projection_conv, CV_DXT_ROWS);

 cvDFT(projection_conv, projection_conv, CV_DXT_INV_SCALE | CV_DXT_ROWS, n_y_pad);
 
 

 cvGetSubRect(projection_conv, &matrix_src, cvRect(filtlength/2, 0, projection_pad->width, projection_pad->height));
 cvCopy(&matrix_src, projection_pad, NULL);



 float *p_image;
  
 for(j=0;j<projection_pad_center_ROI.height;j++)
 {
     p_image = (float *) (projection_pad->imageData+projection_pad_center_ROI.x*sizeof(float)+(projection_pad_center_ROI.y+j)*projection_pad->widthStep);
     
     for(i=0;i<projection_pad_center_ROI.width;i++)
     {
         prj_slice[i+j*prj.X]=*p_image++;
     }
 }

 
          
/* for(i=0;i<line_num;i++)
     prj_real[i+angle*line_num]=prj_slice[i];*/

 for(i=0;i<prj.X;i++)
     for(j=0;j<prj.Y;j++)
         prj_real[i+j*prj.X]=prj_slice[j+i*prj.Y];


 free(prj_slice);

 
 cvReleaseImageHeader(&projection_org);
 cvReleaseImage(&projection_pad);

 
 free(filter_rs);

 cvReleaseMat(&filter_rs_cv);
 cvReleaseMat(&filter_fs_cv);
 cvReleaseMat(&filter_2D_fs_cv);
 cvReleaseMat(&projection_conv);


return TRUE;
}


/*****************************************************************************************************/
/*int filter_prj_X_sym(float *prj_real, char *filter, int filtlength, Projection prj, int row_pad, int col_pad, int symmetrize_2D_flag, int angle, int id)
{
 unsigned int i,j;

 int n_x_pad=prj.X+row_pad; 
 int n_y_pad=prj.Y+col_pad;

 int n_y_pad_conv=filtlength+n_y_pad-1;
 int n_y_pad_conv_opt=cvGetOptimalDFTSize(n_y_pad_conv);

 int n_x_pad_opt=cvGetOptimalDFTSize(n_x_pad);

 if(id==0 && angle==0)
 printf("n_x_pad is %d, n_y_pad is %d, n_y_pad_conv is %d, n_y_pad_conv_opt is %d, n_x_pad_opt is %d\n", n_x_pad, n_y_pad, n_y_pad_conv, n_y_pad_conv_opt, n_x_pad_opt);


 int line_num=prj.X*prj.Y;
 float *prj_slice;
 if((prj_slice=(float *)malloc(sizeof(float)*line_num))==NULL)
   {
       printf("Error with Function 'filter_prj_X_sym()'!Can't malloc memery for 'prj_slice'!");
       return FALSE;
   }
 memset(prj_slice, 0 , sizeof(float)*line_num);

 for(i=0;i<line_num;i++)
     prj_slice[i]=prj_real[i+angle*line_num];

 // original projection (prj.Y*prj.X)
 IplImage *projection_org = cvCreateImageHeader(cvSize(prj.X, prj.Y), IPL_DEPTH_32F,1); 

 int n_bytes_row = prj.X*sizeof(float);
 cvSetData(projection_org, prj_slice, n_bytes_row);

// DEBUG:	write_MRC_image_from_IplImage(projection_org, "/home/wanxiaohua/exp/3D_phantom/reconstruct/txbr-test/projection_org.mrc");
// printf("projection_org is done\n");

 // pad projection (n_y_pad*n_x_pad)
 IplImage *projection_pad = cvCreateImage(cvSize(n_x_pad, n_y_pad), IPL_DEPTH_32F,1); 
 cvSetZero(projection_pad);

 CvRect projection_pad_center_ROI = cvRect(row_pad/2, col_pad/2, prj.X, prj.Y);

 cvSetImageROI(projection_org, cvRect(0, 0, prj.X, prj.Y));
 cvSetImageROI(projection_pad, projection_pad_center_ROI); // NOTE: The extents of this ROI are the same as image0's ROI.
 cvCopy(projection_org, projection_pad, NULL);
 cvResetImageROI(projection_org);
 cvResetImageROI(projection_pad);


 



 
 //padding projection_pad using symmetrize
 if (symmetrize_2D_flag)
     symmetrize_IplImage_values_CvRect(projection_pad, projection_pad_center_ROI);




 //filter
 float *filter_rs = (float *) malloc(filtlength * sizeof(float));  //filter_rs (1*251)
    if (filter_rs == NULL)
    {
        printf("Memory request for filter_rs failed .\n");
        exit(1);
    }
    
 if (!strcmp(filter,"RamLak"))
      RamLak(filtlength, filter_rs);
 else if (!strcmp(filter,"SheppLogan"))
      SheppLogan(filtlength,filter_rs); 
      else
          printf("Invalid filter selected.\n");

// filter_rs_cv (251*1)
 CvMat *filter_rs_cv = create_CvMat_32F_from_float_data(filtlength, 1, filter_rs);  


 // filter_fs_cv (n_y_pad_cov_opt*1)
 CvMat *filter_fs_cv = cvCreateMat(n_y_pad_conv_opt, 1, CV_32FC1);
 cvZero(filter_fs_cv);

 CvMat tmp;
 cvGetSubRect(filter_fs_cv, &tmp, cvRect(0, 0, 1, filtlength));
 cvCopy(filter_rs_cv, &tmp, NULL);

 cvDFT(filter_fs_cv, filter_fs_cv, CV_DXT_FORWARD, n_y_pad_conv_opt);

// filter_fs_cv (n_y_pad_conv_opt*n_x_pad_opt)
 CvMat *filter_2D_fs_cv= cvCreateMat(n_y_pad_conv_opt, n_x_pad_opt,CV_32FC1);
 cvZero(filter_2D_fs_cv);

 cvRepeat(filter_fs_cv,filter_2D_fs_cv);


 //convolution projection_pad and filter_2D_fs_cv into projection_conv
 CvMat *projection_conv=cvCreateMat(n_y_pad_conv_opt, n_x_pad_opt, CV_32FC1);
 cvZero(projection_conv);

 

 CvMat matrix_src;
 CvMat matrix_dst;

 cvGetSubRect(projection_conv, &matrix_dst, cvRect(0, 0, n_x_pad, n_y_pad));
 cvCopy(projection_pad, &matrix_dst, NULL);
// cvGetSubRect(projection_conv, &matrix_dst, cvRect(n_x_pad, 0, n_x_pad_conv_opt-n_x_pad, n_y_pad));


 if (symmetrize_2D_flag)
 {
            
       
     int i,j;
     if(n_y_pad>=n_y_pad_conv_opt-n_y_pad)
     {
         for(i=0;i<n_x_pad_opt;i++)
             for(j=n_y_pad;j<n_y_pad_conv_opt;j++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,2*n_y_pad-1-j,i));
     }
     else
     {
         for(i=0;i<n_x_pad_opt;i++)
             for(j=n_y_pad;j<2*n_y_pad;j++)
                 cvmSet(projection_conv,j,i,cvmGet(projection_conv,2*n_x_pad-1-j,i));
     }

     
 }



 cvDFT(projection_conv, projection_conv, CV_DXT_FORWARD,n_y_pad_conv_opt);

 cvMulSpectrums(projection_conv, filter_2D_fs_cv, projection_conv, CV_DXT_ROWS);

 cvDFT(projection_conv, projection_conv, CV_DXT_INV_SCALE,n_y_pad_conv_opt);
 
 

 cvGetSubRect(projection_conv, &matrix_src, cvRect(0,filtlength/2, projection_pad->width, projection_pad->height));
 cvCopy(&matrix_src, projection_pad, NULL);



 float *p_image;
  
 for(j=0;j<projection_pad_center_ROI.height;j++)
 {
     p_image = (float *) (projection_pad->imageData+projection_pad_center_ROI.x*sizeof(float)+(projection_pad_center_ROI.y+j)*projection_pad->widthStep);
     
     for(i=0;i<projection_pad_center_ROI.width;i++)
     {
         prj_slice[i+j*prj.X]=*p_image++;
     }
 }

 
          
 for(i=0;i<line_num;i++)
     prj_real[i+angle*line_num]=prj_slice[i];


 free(prj_slice);

 
 cvReleaseImageHeader(&projection_org);
 cvReleaseImage(&projection_pad);

 
 free(filter_rs);

 cvReleaseMat(&filter_rs_cv);
 cvReleaseMat(&filter_fs_cv);
 cvReleaseMat(&filter_2D_fs_cv);
 cvReleaseMat(&projection_conv);


return TRUE;
}

*/
