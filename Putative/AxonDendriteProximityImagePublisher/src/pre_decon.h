

#ifndef PREDECON
#define PREDECON

#include "typedefs.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>


#define MAX_DECON_SIZE 	3000


    //------------------------------------------------------------------------------------
    //Type Definitions for 2D Inputfiles
    //------------------------------------------------------------------------------------

    
    typedef unsigned char                                         dim2_PixelType;
  
    typedef itk::Image< dim2_PixelType, 2 >                       dim2_ImageType;
    
    typedef itk::ImageFileReader< dim2_ImageType >                dim2_ReaderType;
    typedef itk::ImageFileWriter< dim2_ImageType >                dim2_WriterType;

    typedef itk::ImageRegionConstIterator< dim2_ImageType >       DeconConstIteratorType;
    typedef itk::ImageRegionIterator< dim2_ImageType >            DeconIteratorType;

    //------------------------------------------------------------------------------------
    //Type Definitions for 3D Outputfiles
    //------------------------------------------------------------------------------------

    typedef itk::Image< dim2_PixelType, 3 >			StackImageType;
    typedef itk::ImageSeriesReader< StackImageType >		StackReaderType;
    typedef itk::ImageFileWriter< StackImageType >		StackWriterType;

    typedef itk::NumericSeriesFileNames				NameGeneratorType;
    
    typedef itk::InvertIntensityImageFilter< dim2_ImageType >	InvertFilterType;
    
    

class PreDecon 
{  
public:
    char norm_buffer[1024];
    char inv_buffer[1024];
    int bricks_x;
    int bricks_y;
    int dimension_matrix [32][32][4];
    
    PreDecon(char* inputfilename, char* outputfilename, int dim_z);
    static void readImageStats(unsigned int* bricks_x, unsigned int* bricks_y, unsigned int* dim_x, unsigned int* dim_y);
    static void writeImageStats(int bricks_x, int bricks_y, int dim_x, int dim_y);
    
private:
    
    
  
};


#endif


