/****************************************************************************/
/*                                                                          */
/* File:      proximity_finder_client.cpp                                   */
/*                                                                          */
/* Purpose:                                                                 */
/*                                                                          */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   16 January 2014                                               */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


#include "proximity_finder.h"




//*******************************************************

ImageType::Pointer readOriginalImage(char* file_name, int start_index, int end_index);



//*******************************************************

int main( int argc , char * argv[])
{
  Reader * amiraReader;
  AmiraSpatialGraph * input_graph;
  ImageType::Pointer input_image;
  Image2DType::Pointer image_projection;
  

  
  
  //Validate command line arguements
  if(argc != 8 && argc != 9)
  {
      ////std::cout << std::endl << "Proper format:  ./ProximityFinder 'amira_file.am' 'dendrite_image%03d.png' 'output_name' startIndex endIndex [maxDist (default 5)]" << std::endl << std::endl;
      //////std::cout << std::endl << "For bricks: ./ProximityFinder 'amira_file.am' 'dendrite_image%03d.png' 'output_name' startIndex endIndex" << std::endl << std::endl;
      return -1;
      // " For bricks: ./ProximityFinder 'amira_file.am' 'input_filename' 'dendrite_image%03d.png' 'output_name' startIndex endIndex'
  }
  
  
  //Read in amira file and decon image
  ////std::cout << "Begin reading" << std::endl;
  amiraReader = new Reader(argv[1]);  
  amiraReader->readSpatialGraphFile(false);
  input_graph = amiraReader->getSpatialGraph();
  

//   input_image = readOriginalImage(argv[2], atoi(argv[4]), atoi(argv[5]));
// 
//     
//     
//     IteratorType2::SizeType size;
//     size[X_COORD] = input_image->GetLargestPossibleRegion().GetSize()[X_COORD];
//     size[Y_COORD] = input_image->GetLargestPossibleRegion().GetSize()[Y_COORD];
//     size[Z_COORD] = input_image->GetLargestPossibleRegion().GetSize()[Z_COORD];
//     ////std::cout<< "input size "<<size[0]<<" "<<size[1]<<" "<<size[2] << std::endl;
//     IteratorType2::IndexType index;
//     index[X_COORD] = 0;
//     index[Y_COORD] = 0;
//     index[Z_COORD] = 0;
//     
//     IteratorType2::RegionType region;
//     region.SetIndex(index);
//     region.SetSize(size);

    ProximityFinder * prox_finder;
    if(argc == 9)
    {   BrickImageReader *BrickImage = new BrickImageReader(argv[2], atoi(argv[4]), atoi(argv[5]) );               
        input_image = BrickImage->GetImage();
        input_image->Update();
        

                //image_projection = reader->GetZProjection(1);
                //////std::cout <<std::endl << "brick image successfully read!" <<std::endl << std::endl;
                
         prox_finder = new ProximityFinder(input_graph, input_image, argv[1], atoi(argv[6]), atoi(argv[7]), atof(argv[8]));
         
//         prox_finder->readAmiraSectionTransformations(inputfilename);
         prox_finder->findProximities();
         prox_finder->readAmiraTransformations();   
         prox_finder->writeProximityImages(argv[3]);
         prox_finder->writeProximityLandmarks(argv[3]);

         delete BrickImage;
         
         
    }
    
    if(argc == 8)
    {
                BrickImageReader *BrickImage = new BrickImageReader(argv[2], atoi(argv[4]), atoi(argv[5]));               
                input_image = BrickImage->GetImage();
                input_image->Update();

                //image_projection = reader->GetZProjection(1);
                //////std::cout <<std::endl << "brick image successfully read!" <<std::endl << std::endl;
                
                prox_finder = new ProximityFinder(input_graph, input_image, argv[1], atoi(argv[6]), atoi(argv[7]));
//                const char* inputfilename = argv[2];
//                prox_finder->readAmiraSectionTransformations();
                prox_finder->findProximities();
                prox_finder->readAmiraTransformations();
                prox_finder->writeProximityImages(argv[3]);
                prox_finder->writeProximityLandmarks(argv[3]);

                delete BrickImage;
    }
    
    //else
    //{
    //     prox_finder = new ProximityFinder(input_graph, input_image);
    //}
  
  return 0;
};



//********************************************************



ImageType::Pointer readOriginalImage(char* file_name, int start_index, int end_index)
{

        char input_file[1024];
        strcpy(input_file, file_name);


        NameGeneratorType::Pointer name_gen = NameGeneratorType::New();         //DECLARE AND INITIALIZE NAME_GENERATOR
        name_gen->SetSeriesFormat( input_file );
        name_gen->SetStartIndex( start_index );
        name_gen->SetEndIndex( end_index );
        name_gen->SetIncrementIndex( 1 );
                        
        SeriesReaderType::Pointer input_reader = SeriesReaderType::New();       //DECLARE AND INITIALIZE INPUT_READER
        input_reader->SetImageIO( itk::PNGImageIO::New() );
        input_reader->SetFileNames( name_gen->GetFileNames() );

        try
        {
                input_reader->Update();
        }
        catch( itk::ExceptionObject & err )
        {
                std::cerr << "ImageReaderExceptionObject caught !" << std::endl;
                std::cerr << err << std::endl;
        }


        ImageType::Pointer original_image;

        original_image = input_reader->GetOutput();
        original_image->Update();
        
        return original_image;
};


        
    







