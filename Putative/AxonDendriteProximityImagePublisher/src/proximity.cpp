/****************************************************************************/
/*                                                                          */
/* File:      proximity.cpp                                                 */
/*                                                                          */
/* Purpose:   Defines the methods for the Proximity class                   */
/*                                                                          */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   31.03.2014                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include "proximity.h"
#include "typedefs.h"





/****************************************************************************/
/*                                                                          */
/* writeProximityImage( char * output_file)                                 */
/*                                                                          */
/*      writes the region in the vicinity of the bounding box to an         */
/*      imagestack in a separate folder                                     */
/*                                                                          */
/****************************************************************************/
void Proximity::writeProximityImage(const char * output_folder, const char * output_file, double ** amira_transformation, int prox_indx)
{
    ////std::cout<< " In Writing Proximity Image  "<< std::endl;
    
    //get the cropped roi from the original image
    this->getProximityImage();
    
    std::string folder_name = output_folder;

    
        
//         if(Axon_ID == AXON_ID_1 && Dendrite_ID == DENDRITE_ID_2)
//         folder_name += "_CELL#1-->CELL#2";
//         if(Axon_ID == AXON_ID_2 && Dendrite_ID == DENDRITE_ID_1)
//         folder_name += "_CELL#2-->CELL#1";
//         if(Axon_ID == AXON_ID_1 && Dendrite_ID == DENDRITE_ID_3)
//         folder_name += "_CELL#1-->CELL#3";
//         if(Axon_ID == AXON_ID_3 && Dendrite_ID == DENDRITE_ID_1)
//         folder_name += "_CELL#3-->CELL#1";
//         if(Axon_ID == AXON_ID_2 && Dendrite_ID == DENDRITE_ID_3)
//         folder_name += "_CELL#2-->CELL#3";
//         if(Axon_ID == AXON_ID_3 && Dendrite_ID == DENDRITE_ID_2)
//         folder_name += "_CELL#3-->CELL#2";
        
//     folder_name += "_roi";
//     folder_name += Utility::makeString(index_num);

    std::cout<<"Proximity image folder name: "<<folder_name<<std::endl;   
    
//     std::string mkdir_command = "mkdir " + folder_name;
//     system(mkdir_command.c_str());
    if(mkdir(folder_name.c_str(),
            S_IRWXU| //full access for owner
            S_IRWXG| //full access for grp
            S_IRWXO  //full access for other
            ))
        {
	   std::cout<<"Creating image folders "<<std::endl;
//         if(errno == EEXIST){
//             std::cerr << folder_name << " already exists" << std::endl;
//         }
//         else {
//             perror("Error creating directory");
//             exit(-1);
//         }
        }
    
//     std::string cd_command = "cd " + folder_name; 
//     system(cd_command.c_str());
    
    
    
        
    writeImagePlanes( folder_name.c_str(), output_file);
    
    writeLandmarkFile(folder_name.c_str(), output_file, amira_transformation, prox_indx);
    
    writeInfoFile(folder_name.c_str(), output_file, prox_indx);
    
    writeHXFile( folder_name.c_str(), output_file,  amira_transformation, prox_indx);
  
//     cd_command = "cd ..";
//     system(cd_command.c_str());
};


/*****************************************************************************************
Method:  getProximityImage()

Copies all pixel values from the IteratorRegion in the original image into a second
OutputRegion in the output image.

Return: A pointer to the filled output image
******************************************************************************************/
void Proximity::getProximityImage()
{
        ////std::cout<< " In getProximityImage  "<< std::endl;
    this->iterator_region = getIteratorRegion();
    this->output_region = getOutputRegion();
        
    ImageType::Pointer output_image = ImageType::New(); 
    output_image->SetRegions( output_region );
    output_image->Allocate(); 
    
    IteratorType2 og_iterator(this->original_image, iterator_region);
    IteratorType2 output_iterator(output_image, output_image->GetLargestPossibleRegion());
    
    for (og_iterator.GoToBegin(), output_iterator.GoToBegin();!og_iterator.IsAtEnd(); ++og_iterator, ++output_iterator ) 
    {
        unsigned char value = og_iterator.Get();      
        output_iterator.Set((unsigned char)value);
    }
    
    this->proximity_image = output_image;
};


/*******************************************************************************************************
Method: getIteratorRegion()

Defines a Region in the original image that contains all pixels in the bounding_box plus a buffer 
around the edges of the box. The size of this buffer is defined by the global constant XY_BUFFER
********************************************************************************************************/
ImageType::RegionType Proximity::getIteratorRegion()
{
    
    ////std::cout<< " In getIteratorRegion  "<< std::endl;
    int max_x  = original_image ->GetLargestPossibleRegion().GetSize(0);  // size along X
    int max_y  = original_image->GetLargestPossibleRegion().GetSize(1);  // size along Y
    int max_z  = original_image->GetLargestPossibleRegion().GetSize(2);  // size along Z
    
   ////std::cout<< max_x<<" "<<max_y<<" "<<max_z<< std::endl;

    ImageType::SizeType  size;   
    ImageType::IndexType start;
    
    if((this->getRealCentroidInImageCoords()[0] - MAX_PROX_BOX_SIZE_IMAGE_COORD_XY) > 0)
        start[0] = this->getRealCentroidInImageCoords()[0] - MAX_PROX_BOX_SIZE_IMAGE_COORD_XY;
    else
        start[0] = 0;
    
   if((this->getRealCentroidInImageCoords()[1] - MAX_PROX_BOX_SIZE_IMAGE_COORD_XY) > 0)
        start[1] = this->getRealCentroidInImageCoords()[1] - MAX_PROX_BOX_SIZE_IMAGE_COORD_XY;
    else
        start[1] = 0;
    
    if((this->getRealCentroidInImageCoords()[2] - MAX_PROX_BOX_SIZE_IMAGE_COORD_Z) > 0)
        start[2] = this->getRealCentroidInImageCoords()[2] - MAX_PROX_BOX_SIZE_IMAGE_COORD_Z;
    else
        start[2] = 0;
    
    
    if( (start[0] + MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2) >= max_x)
        size[0] = max_x - start[0] - 1;
    else 
        size[0] = MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2;
    
    if( (start[1] + MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2) >= max_y)
        size[1] = max_y - start[1] - 1;
    else 
        size[1] = MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2;
    
    if( (start[2] + MAX_PROX_BOX_SIZE_IMAGE_COORD_Z*2) >= max_z)
        size[2] = max_z - start[2] - 1;
    else 
        size[2] = MAX_PROX_BOX_SIZE_IMAGE_COORD_Z*2;
       
    
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    //region.Print(//////std::cout);

    return region;
};


/*******************************************************************************************************
Method: getOutputRegion()

Defines a Region that will define the size of the output image. This should be the same size as the 
IteratorRegion
********************************************************************************************************/
ImageType::RegionType Proximity::getOutputRegion()
{
    
    ////std::cout<< " In getOutputRegion  "<< std::endl;
    int max_x  = original_image ->GetLargestPossibleRegion().GetSize(0);  // size along X
    int max_y  = original_image->GetLargestPossibleRegion().GetSize(1);  // size along Y
    int max_z  = original_image->GetLargestPossibleRegion().GetSize(2);  // size along Z
    
   

    ImageType::SizeType  size;   
    ImageType::IndexType start;
    
    start.Fill(0);
    
    
    if( (start[0] + MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2) >= max_x)
        size[0] = max_x - start[0] - 1;
    else 
        size[0] = MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2;
    
    if( (start[1] + MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2) >= max_y)
        size[1] = max_y - start[1] - 1;
    else 
        size[1] = MAX_PROX_BOX_SIZE_IMAGE_COORD_XY*2;
    
    if( (start[2] + MAX_PROX_BOX_SIZE_IMAGE_COORD_Z*2) >= max_z)
        size[2] = max_z - start[2] - 1;
    else 
        size[2] = MAX_PROX_BOX_SIZE_IMAGE_COORD_Z*2;
       
    
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    //region.Print(//////std::cout);

    return region;
};


    
/****************************************************************************/
/*                                                                          */
/* writeImagePlanes(ImageType::Pointer input_image, char * output_file)     */
/*                                                                          */
/*Writes the 3D ITK image into an image stack of individual PNG slices      */
/*                                                                          */
/****************************************************************************/
void Proximity::writeImagePlanes( const char* folder_name, const char* file_name)
{
     ////std::cout<< " In writeImagePlanes  "<< std::endl;
        #ifdef DEBUG
          //////std::cout<< "Generate Writer !!" << std::endl;
        #endif

          Writer2DType::Pointer writer = Writer2DType::New();
          writer->SetInput( this->proximity_image );

          NameGeneratorType::Pointer writer_name_gen = NameGeneratorType::New();

          std::string format = folder_name;
          format += "/"; //+ output_file;
	  format += file_name;
          //std::string format = output_file;
          //format += "_seg";
          format += "%03d.";
          format += "tif";

          writer_name_gen->SetSeriesFormat( format.c_str() );

          ImageType::RegionType region = this->proximity_image->GetLargestPossibleRegion();
          ImageType::IndexType start = region.GetIndex();
          ImageType::SizeType size = region.GetSize();

          const unsigned int first_slice = start[2];
          const unsigned int last_slice  = start[2] + size[2] -1;

          writer_name_gen->SetStartIndex( first_slice );
          writer_name_gen->SetEndIndex( last_slice );
          writer_name_gen->SetIncrementIndex( 1 );

        #ifdef DEBUG
          //////std::cout << "-------------------------------------------------------------------" << std::endl;
          //////std::cout<< "Start Writing !!" << std::endl;
          //////std::cout << "-------------------------------------------------------------------" << std::endl;
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
          //////std::cout<< "Writing Done !!" << std::endl;
        #endif

}; 


void Proximity::writeLandmarkFile( const char* outputFilename, const char* file_name, double ** amira_transformation, int prox_indx)
{
    ////std::cout<< " In writeLandmarkFile  "<< std::endl;
  //////std::cout<< "in writeLandmarkFile" << std::endl;
  
// transformToWorldCoordinates(list);
  
  std::string format = outputFilename;
  format += "/";
  format += file_name;
  //format += Utility::inttochar(prox_indx+1);
  format += "_landmark_locations.landmarkAscii";
  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << 1                << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "Parameters {"                        << std::endl;
  LandMarkData << "    NumSets 1,"                      << std::endl;
  LandMarkData << "    ContentType \"LandmarkSet\""     << std::endl;
  LandMarkData << "}"                                   << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "Markers { float[3] Coordinates } @1" << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "# Data section follows"              << std::endl;
  LandMarkData << "@1"                                  << std::endl;
  
  //for(int i=0; i < list->size(); i++)
  {
    Utility::applyTransformationToPoint(getRealCentroid(), amira_transformation);
    LandMarkData << this->getRealCentroid()[X_COORD] << " " << this->getRealCentroid()[Y_COORD] << " " << this->getRealCentroid()[Z_COORD] << std::endl;
  } 
  
  LandMarkData.close();
};


void Proximity::writeInfoFile(const char* folder_name, const char* file_name, int prox_indx)
{
    ////std::cout<< " In writeInfoFile  "<< std::endl;
    //writes an InfoFile to ensure right positioning of the proximity images in Amira
    std::string format = folder_name;
    format += "/";
    format += file_name;
    //format += Utility::inttochar(prox_indx+1);
    format += "_InfoFile.info";
    
    ////std::cout<<"info Filename: "<<format<<std::endl;
    
    std::ofstream ListData( format.c_str() );
    
    ImageType::SizeType region_size;
    
    ImageType::RegionType region = this->proximity_image->GetLargestPossibleRegion();
    double pixel_size = XYSAMPLING;
    
    region_size = region.GetSize();
    ImageType::IndexType start = region.GetIndex();
    

    
    const unsigned int first_slice = 0;
    const unsigned int last_slice  = region_size[2] -1;
    
    ListData << "# Amira Stacked Slices"<< std::endl;
    ListData << " "<< std::endl;
    ListData << "pixelsize "<<pixel_size<<"  "<<pixel_size<< std::endl;
    ListData << "       "<<std::endl;
    
    for(int i = 0 ; i <= last_slice; i++)
    {
        
        std::string imagename = file_name;
           
           
        
        ListData <<" "<<imagename<<std::setfill('0')<< std::setw(3)<<i<<".tif"<<" "<<first_slice + 0.5*i <<std::endl;
    }
    
    ListData.close();  
    

    
};

void Proximity::writeHXFile( const char* outputFoldername, const char* file_name, double** transformation, int prox_indx)
{
    ////std::cout<< " In writeHXFile  "<< std::endl;
    ImageType::Pointer prox_image = this->proximity_image;//this->getProximityImage();        
    ImageType::RegionType output_region = this->output_region;//this->getOutputRegion();
    ImageType::RegionType iterator_region = this->iterator_region;//this->getIteratorRegion();
    ImageType::SizeType region_size;
        
     
        region_size = output_region.GetSize();
        ImageType::IndexType start = iterator_region.GetIndex();
        
        std::string format2 = outputFoldername;
        format2 += "/Load_Proximity_Zone_";
        format2 += Utility::inttochar(prox_indx+1);    
        format2 += ".hx";
        //////std::cout<<"HX Filename: "<<format2<<std::endl;
        std::ofstream ListData( format2.c_str() );
        
//         std::string proximity_label = "_Axon_";
//         proximity_label += Utility::makeString(Axon_ID);
//         proximity_label += "_Dendrite_";
//         proximity_label += Utility::makeString(Dendrite_ID);
//         proximity_label += "_";
        
	
        
        ListData << "# Amira Script"<< std::endl;
        ListData << " "<< std::endl;
        //ListData << "[ load ${SCRIPTDIR}/Proximity_Zone_"<< prox_indx+1 <<"_InfoFile.info ] setLabel {Proximity_Zone_"<< prox_indx+1 <<"_Images.info}"<< std::endl;
	ListData << "[ load ${SCRIPTDIR}/"<< file_name <<"_InfoFile.info ] setLabel {"<< file_name <<"_Images.info}"<< std::endl;
        ListData << "       "<<std::endl;
        
        //ListData << "Proximity_Zone_"<< prox_indx+1 <<"_Images.info setTranslation " <<start[X_COORD]*XYSAMPLING<<" "<<start[Y_COORD]*XYSAMPLING<<" "<<start[Z_COORD]*ZSAMPLING<<" "<<std::endl;
        ListData << file_name <<"_Images.info setTranslation " <<start[X_COORD]*XYSAMPLING<<" "<<start[Y_COORD]*XYSAMPLING<<" "<<start[Z_COORD]*ZSAMPLING<<" "<<std::endl;
	//ListData << "proximity"<<proximity_label<<i<<".info setTranslation " <<list->at(i)[X_COORD] - region_size[0]/2*XYSAMPLING<<" "<<list->at(i)[Y_COORD] - region_size[1]/2*XYSAMPLING<<" "<<list->at(i)[Z_COORD] - region_size[2]/2*ZSAMPLING<<std::endl;
        ListData << file_name <<"_Images.info applyTransform "<<  std::endl;
        ListData << file_name <<"_Images.info setTransform " <<transformation[0][0]<< " "<<transformation[1][0]<< " "<<transformation[2][0]<< " "<<transformation[3][0]<< " "<<transformation[0][1]<< " "<<transformation[1][1]<< " "<<transformation[2][1]<< " "<<transformation[3][1]<< " "<<transformation[0][2]<< " "<<transformation[1][2]<< " "<<transformation[2][2]<< " "<<transformation[3][2]<< " "<<transformation[0][3]<< " "<<transformation[1][3]<< " "<<transformation[2][3]<< " "<<transformation[3][3]<< " "<< std::endl;
        ListData << "       "<<std::endl;
	
        ListData << "[ load ${SCRIPTDIR}/"<<file_name<<"_landmark_locations.landmarkAscii ]"<< std::endl;
        
            
        ListData.close(); 
   
}; 

// void Proximity::calculateConfidenceValue(std::vector<double*> * confidenceList, std::vector<double>  distances){
//     
//     double confidence_val;
//     
//     double proximity_avg_bright;
//     double proximity_avg_std_dev;
//     double proximity_avg_radius;
//     double proximity_avg_distance;
//     double min_distance = 5;
//     int counter = 0;
//     int counter2 = 0;
//     double bouton_confidence_val;
//     int count = 0;
//     int boutons_detected = confidenceList.size();
//     double confidence;
//     
//    
//     
//     for(int i = 0; i < distances.size(); i++){
//         
//             proximity_avg_distance += distances.at(i);
//             counter2++;
//             if(distances.at(i) < min_distance)
//             min_distance = distances.at(i);
//             
//             
//             
//             
//         }
//         
//         
//         
//         
//         
//         proximity_avg_distance = proximity_avg_distance/counter2;
// //        //////std::cout<<"average_distance:  "<<proximity_avg_distance<<std::endl;
// //        //////std::cout<<"min_distance:  "<<min_distance<<std::endl;
// //         
//         
//    for(int i = 0; i< confidenceList->size(); i++){ 
//        
//        
//         double * coords = confidenceList->at(i);
// //     if(coords[IS_BOUTON]){
//      float brightness = coords[LOCAL_BRIGHTNESS];
//      float std_dev = coords[LOCAL_SIGMA];
//      float threshold = coords[THRESHOLD];
//      float radius = coords[SURFACE];
//      
// //     if(coords[IS_BOUTON])
// //        //////std::cout<<"has bouton!"<<std::endl;
// //     bouton_confidence_val = brightness/threshold - std_dev + radius;
// //     }
// //     if( (brightness > (avg_bright + 10)) && (radius > avg_rad) && (radius >= 3.5) ) 
// //                coords[IS_BOUTON] = 1;
// //              else
// //                coords[IS_BOUTON] = 0;
//     
//     
// //     proximity_avg_bright = proximity_avg_bright/counter;
// //     proximity_avg_std_dev = proximity_avg_std_dev/counter;
// //     proximity_avg_radius = proximity_avg_radius/counter;
// //     counter++;
//    }
//     
// //    confidence_val = proximity_avg_bright + proximity_avg_std_dev + proximity_avg_radius /*+ 2* min_distance*/;
// //      //////std::cout<<"confidence list size: "<<confidenceList->size()<<std::endl;
//     for(int i= 0;  i< boutons_detected; i++){
// //       //////std::cout<<"confidence : "<<confidenceList->at(i)<<std::endl;
//       confidence = confidenceList.at(i);
//     bouton_confidence_val += confidence;
// //     //////std::cout<<"bouton_confidence_value: "<<bouton_confidence_val<<std::endl;
//     ++count;
// //     //////std::cout<<"counter: "<<count<<std::endl;
//     }
//     bouton_confidence_val = bouton_confidence_val/count;
//     if(boutons_detected == 0)
//       bouton_confidence_val = 0;
//     
// //     //////std::cout<<"bouton_confidence_value: "<<bouton_confidence_val<<std::endl;
// //     //////std::cout<<"min_distance:  "<<min_distance<<std::endl;
//     //this->confidence_value = (bouton_confidence_val + 1/(0.1 + min_distance))/2;
//     this->confidence_value = bouton_confidence_val;
//     //////std::cout<<"confidence value: "<<this->confidence_value<<std::endl;
//     
//     
//     
//     
// 
//     
//     
// }








