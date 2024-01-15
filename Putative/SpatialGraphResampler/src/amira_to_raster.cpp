

#include "bouton_finder.h"
#include "typedefs.h"
#include "typedefs_two.h"
#include "amiraReader.h"

#define DEBUG

void WriteImagePlanes(ImageType::Pointer input_image, char output_file[1024]);

int main( int argc , char * argv[])
{
  if(argc != 5)
  {
      ////////////////std:://cout << "FORMAT:  ./program amira_file.am  pixels_x  pixels_y  dim_z" << std::endl;
      return -1;
  }
  
  int dim_z = atoi(argv[4]);
  
  Reader * amiraReader = new Reader(argv[1]);  
  
  ////////////////std:://cout << "Begin reading" << std::endl;

  amiraReader->readSpatialGraphFile(false);
  
  
  AmiraSpatialGraph* graph = amiraReader->getSpatialGraph();
  std::vector< Vertex * > * vertices = graph->verticesPointer();
  int numOfVerts = vertices->size();
  
  ImageType::Pointer output_image = ImageType::New(); 
  
  ImageType::IndexType start;
  start[0] =   0;  // first index on X
  start[1] =   0;  // first index on Y
  start[2] =   0;  // first index on Z
  
  ImageType::SizeType  size;
  size[0]  = std::atoi(argv[2]);  // size along X
  size[1]  = std::atoi(argv[3]);  // size along Y
  size[2]  = dim_z;  // size along Z
  
  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  
  output_image->SetRegions( region );
  output_image->Allocate();
  
  
   ////////////////std:://cout << "Begin processing" << std::endl;
  
  for(int i=0; i<numOfVerts; i++)
  {
    
    Vertex * currentVertex = vertices->at(i);
    
    ImageType::IndexType pixelIndex;
    int x_pos = rint( currentVertex->coordinates[0] / .094 );
    int y_pos = rint( currentVertex->coordinates[1] / .094 );
    int z_pos = rint( currentVertex->coordinates[2] / .5 );
    
    pixelIndex[0] = x_pos;
    pixelIndex[1] = y_pos;
    pixelIndex[2] = z_pos;
    
    ////////////////std:://cout << "Vert:" << i << " x=" << x_pos << " y=" << y_pos << " z=" << z_pos << std::endl;
    
    output_image->SetPixel(pixelIndex, 255);
    
  }
  
  std::vector< Edge * > * edges = graph->edgesPointer();
  int numOfEdges = edges->size();
  
  for(int i=0; i<numOfEdges; i++)  //for each edge
  {
    
    Edge * currentEdge = edges->at(i);
    std::list< double * >::iterator it;
    
    for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
    {
      double * coords = *it;
           
      ImageType::IndexType pixelIndex;
      int x_pos = rint( coords[0] / .094 );
      int y_pos = rint( coords[1] / .094 );
      int z_pos = rint( coords[2] / .5);
      
      pixelIndex[0] = x_pos;
      pixelIndex[1] = y_pos;
      pixelIndex[2] = z_pos;
      
      ////////////////std:://cout << "Edge:" << i << " x=" << x_pos << " y=" << y_pos << " z=" << z_pos << std::endl;
     
      output_image->SetPixel(pixelIndex, 255);
    }
    
  }

  ////////////////std:://cout << "Begin writing" << std::endl;

  
  WriteImagePlanes(output_image, "raster_amira_slice");
  
  
  return 0;
}



/****************************************************************************/
/*WriteImagePlanes(output filename) evokes writing of 2D png planes         */
/****************************************************************************/

void WriteImagePlanes(ImageType::Pointer input_image, char output_file[1024])
{
#ifdef DEBUG
  ////////////////std:://cout<< "Generate Writer !!" << std::endl;
#endif

  Writer2DType::Pointer writer = Writer2DType::New();
  writer->SetInput( input_image );

  NameGeneratorType::Pointer writer_name_gen = NameGeneratorType::New();

  std::string format = output_file;
  //format += "_seg";
  format += "%03d.";
  format += "png";

  writer_name_gen->SetSeriesFormat( format.c_str() );

  ImageType::RegionType region = input_image->GetLargestPossibleRegion();
  ImageType::IndexType start = region.GetIndex();
  ImageType::SizeType size = region.GetSize();

  const unsigned int first_slice = start[2];
  const unsigned int last_slice  = start[2] + size[2] -1;

  writer_name_gen->SetStartIndex( first_slice );
  writer_name_gen->SetEndIndex( last_slice );
  writer_name_gen->SetIncrementIndex( 1 );

#ifdef DEBUG
  ////////////////std:://cout << "-------------------------------------------------------------------" << std::endl;
  ////////////////std:://cout<< "Start Writing !!" << std::endl;
  ////////////////std:://cout << "-------------------------------------------------------------------" << std::endl;
#endif

  writer->SetFileNames( writer_name_gen->GetFileNames() );

  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & err )
    {
      std::cerr << "SegWriterExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
    }

#ifdef DEBUG
  ////////////////std:://cout<< "Writing Done !!" << std::endl;
#endif
}; 

