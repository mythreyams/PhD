#include "typedefs.h"
#include "utility.h"
#include "pre_decon.h"

#ifdef CONFOCAL
  #define BRICK_DIMENSIONS      3000  //confocal
  #define BUFFER_DIMENSIONS     50
#else
  #define BRICK_DIMENSIONS      1500    //brightfield
  #define BUFFER_DIMENSIONS     25 
#endif

class BrickImageReader{
public:
    BrickImageReader(const char * input, unsigned int start, unsigned int end);
    // To only load a brick if necessary
    BrickImageReader(const char * input, unsigned int start, unsigned int end, std::vector<double*>*prox_points);
    ImageType::Pointer GetImage(){return input_image;};
    
private:
  ImageType::Pointer input_image;
  const char * inputfilename;
  unsigned int bricks_x, bricks_y, start_index, end_index, dim_x, dim_y;
  
  bool isAnyProxInThisBrick(std::vector<double*> *prox_points, ImageType::IndexType start_index, ImageType::SizeType brick_size );
  
};

