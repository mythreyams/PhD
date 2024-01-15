#include "axon_dendrite_proximity_finder.h"

ImageType::Pointer readOriginalImage(char* file_name, int start_index, int end_index);

int main( int argc , char * argv[])
{
  Reader * amiraReader;
  AmiraSpatialGraph * input_graph;
  ImageType::Pointer input_image;
  Image2DType::Pointer image_projection;
  
  
  //std::cout << "Begin reading" << std::endl;
  amiraReader = new Reader(argv[1]);  
  amiraReader->readSpatialGraphFile(false);
  input_graph = amiraReader->getSpatialGraph();
  
  
  AxonDendriteProximityFinder * prox_finder;
  
  
  if(argc == 9)
    {   BrickImageReader *BrickImage = new BrickImageReader(argv[2], atoi(argv[4]), atoi(argv[5]) );               
        input_image = BrickImage->GetImage();
        input_image->Update();
        

                //image_projection = reader->GetZProjection(1);
                ////std::cout <<std::endl << "brick image successfully read!" <<std::endl << std::endl;
                
         prox_finder = new AxonDendriteProximityFinder(input_graph, input_image, argv[1], atoi(argv[6]), atoi(argv[7]), atof(argv[8]));
         
  
  