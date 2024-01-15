#include "typedefs.h"
#include "utility.h"
#include "basics.h"
#include "amiraReader.h"


//#define ROI_SIZE                25
#ifdef CONFOCAL
  #define BRICK_DIMENSIONS      3000  //confocal
  #define BUFFER_DIMENSIONS     50
#else
  #define BRICK_DIMENSIONS      1500    //brightfield
  #define BUFFER_DIMENSIONS     25 
#endif

class BrickImageReader{
public:
    BrickImageReader(const char * input, unsigned int start, unsigned int end, const char * loadhx_file);
    void RoiGenerator( int section_num, const char * output_path, const char * spatial_graph, const char * loadhx_file);
    ImageType::Pointer GetImage(){return input_image;};
    void writeProximityImage(ImageType::Pointer image, const char * output_file, int section_num, int index_num, ImageType::RegionType source_region);
    void writeImagePlanes(ImageType::Pointer input_image, const char* output_file, ImageType::RegionType source_region);
    void writeInfoFile(const char* output_file, const unsigned int first_slice, const unsigned int last_slice);
    void writeHXFile( const char* outputFilename, ImageType::RegionType source_region);
    //bool getTransformationandInverseTransformation(const char* inputfilename);
    bool readAmiraTransformations(const char * inputfilename);
    void getSpatialGraphWithinRoi(ImageType::RegionType source_region);
    
private:
  ImageType::Pointer input_image;
  const char * inputfilename;
  Reader * amiraInputGraphReader;
  AmiraSpatialGraph * roiSpatialGraph;
  unsigned int bricks_x, bricks_y, start_index, end_index, dim_x, dim_y;
  //ImageType::Pointer *Roi;
  double   sectionTranslation[4][4];// = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double   sectionRotation[4][4];// =  {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double   transformation[4][4];//{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double   inverse_transformation[4][4];//= {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  
};

