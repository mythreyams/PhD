#include "brick_image_reader.h"
#include "amiraReader.h"


BrickImageReader::BrickImageReader(const char * input, unsigned int start, unsigned int end, const char * loadhx_file)
{
  
  //TODO: take offset region into consideration
        //readAmiraTransformations(loadhx_file);
  
  
        input_image = ImageType::New();
        inputfilename  = "deconvolved_z%03d.png";
        start_index = start;
        end_index = end;
        
        //PreDecon::readImageStats(&bricks_x, &bricks_y, &dim_x, &dim_y);

        std::cout << "Inputfilename:" <<  inputfilename << std::endl;
        
        if(chdir(inputfilename) != 0)
                  perror("Couldn't open section drirectory!");
        
        char input_file[1024];
        strcpy(input_file, inputfilename);      

        unsigned int bricks, brick_dim_x, brick_dim_y, offset_x, offset_y; 
        
        
        std::ifstream fp;
        
        //fp=fopen("image_stats.txt", "r");
        
         fp.open("image_stats.txt", std::fstream::in);
    
        
        if(!fp.good())
        { 
            perror("Error opening 'image_stats.txt'\nFile formatting: bricks_x bricks_y pixels_x pixels_y\n");
            exit(-1);
        }
        
        //fscanf(fp, "%d %d %d %d ", bricks_x, bricks_y, dim_x, dim_y);
        fp>>bricks_x>> bricks_y>> dim_x>> dim_y;
        
        std::cout<<bricks_x<< bricks_y<< dim_x<< dim_y<<std::endl;
        
        bricks = bricks_x * bricks_y;
        
        
        
        brick_dim_x = dim_x / bricks_x;
        brick_dim_y = dim_y / bricks_y;
        
        offset_x = dim_x % bricks_x;
        offset_y = dim_y % bricks_y;
        
        
        ImageType::IndexType input_index;               //INITIALIZE VECTOR VALUES FOR REGION
        ImageType::SizeType input_size;

        input_index[0] = 0;
        input_index[1] = 0;
        input_index[2] = 0;

        input_size[0] = dim_x;
        input_size[1] = dim_y;
        input_size[2] = end_index - start_index +1;
        
        
        std::cout << "input_size: " << input_size[0] << " " << input_size[1] << " " << input_size[2] << std::endl;
        
        ImageType::RegionType input_region;             //DECLARE REGION
        input_region.SetSize( input_size );
        input_region.SetIndex( input_index );

        input_image->SetRegions( input_region );        //ALLOCATE REGION
        input_image->Allocate();
        
        
        const char * directoryname = "brick";           //DIRECTORY NAMES FOR ACCESSING DIFFERENT BRICKS
        char buf[1024];
        strcpy(buf, directoryname);

        const char * dots = "..";
        char dotbuf[1024];
        strcpy(dotbuf, dots);
        
        int count = 1;
        
        for(int i = 0; i < bricks_y; i++)               //FOR ALL BRICKS IN BRAIN SECTION
        for(int j = 0; j < bricks_x; j++)
        {
                char buffer[1024];
                sprintf(buffer,"%s%02d", buf, count);                   //OPEN DIRECTORY OF CURRENT BRICK
                if(chdir(buffer) != 0)
                  perror("Couldn't open brick drirectory!");

                std::cout<< "Reading brick y=" << i << " x=" << j << std::endl;
                
                ImageType::IndexType paste_index;                       //INITIALIZE VECTOR VALUES FOR REGION
                ImageType::SizeType paste_size;

                //paste_index[0] = j * BRICK_DIMENSIONS;
                //paste_index[1] = i * BRICK_DIMENSIONS;
                paste_index[0] = j * brick_dim_x;
                paste_index[1] = i * brick_dim_y;
                paste_index[2] = 0;

                

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
                        std::cerr << "BrickReaderExceptionObject caught !" << std::endl;
                        std::cerr << err << std::endl;
                }

                ImageType::IndexType brick_index;
                ImageType::SizeType brick_size;
                
                ImageType::RegionType brick_region;

                ImageType::Pointer brick_image;
                brick_image = input_reader->GetOutput();
                brick_image->Update();
                
		ImageType::RegionType region = brick_image->GetLargestPossibleRegion();
		ImageType::IndexType start = region.GetIndex();
		ImageType::SizeType size = region.GetSize();

                std::cout<<"read brick"<<std::endl;
                std::cout<<"read brick index"<<start[0]<<" "<<start[1]<<" "<<start[2]<<std::endl;
                std::cout<<"read brick size"<<size[0]<<" "<<size[1]<<" "<<size[2]<<std::endl;
               
                
                
                if(i != bricks_y - 1 && j != bricks_x - 1 )
                {
                  brick_size[0] = brick_dim_x;
                  brick_size[1] = brick_dim_y;
                }
                else if(i != bricks_y - 1 && j == bricks_x - 1 )
                {
                  brick_size[0] = brick_dim_x + offset_x;
                  brick_size[1] = brick_dim_y;
                }
                else if(i == bricks_y - 1 && j != bricks_x - 1 )
                {
                  brick_size[0] = brick_dim_x;
                  brick_size[1] = brick_dim_y + offset_y;
                }
                else if(i == bricks_y - 1 && j == bricks_x - 1 )
                {
                  brick_size[0] = brick_dim_x + offset_x;
                  brick_size[1] = brick_dim_y + offset_y;       
                }       
                brick_size[2] = end_index - start_index +1;
                
                
                
                if(bricks == 1)
                {
                  brick_index[0] = 0;
                  brick_index[1] = 0;
                  brick_index[2] = 0;
                }
                else if(bricks_x == 1)
                {
                  brick_index[0] = 0;
                  brick_index[1] = BUFFER_DIMENSIONS;
                  brick_index[2] = 0;
                }
                else if(bricks_y == 1)
                {
                  brick_index[0] = BUFFER_DIMENSIONS;
                  brick_index[1] = 0;
                  brick_index[2] = 0;
                }
                else
                {
                  brick_index[0] = BUFFER_DIMENSIONS;
                  brick_index[1] = BUFFER_DIMENSIONS;
                  brick_index[2] = 0;   
                  
                }
                

                brick_region.SetSize( brick_size );
                brick_region.SetIndex( brick_index );
                
                


                PasteFilterType::Pointer paster = PasteFilterType::New();

                paster->SetSourceImage( brick_image );
                paster->SetSourceRegion( brick_region );
                paster->SetDestinationImage( input_image );
                paster->SetDestinationIndex( paste_index );
                paster->Update();

                input_image = paster->GetOutput();

                input_image->Update();
                
                
                std::cout<<"pasted brick"<<std::endl;
                
                
                
                char dotbuffer[1024];
                sprintf(dotbuffer,"%s",dotbuf);
                if(chdir(dotbuffer) != 0)
                  perror("Couldn't switch to Previous Directory!");
                                
                count++;
               
        }
        
        
};
 










void BrickImageReader::RoiGenerator( int section_num, const char * output_path, const char * spatial_graph, const char * loadhx_file)
{

        //input_image = ImageType::New();
        //inputfilename  = "deconvolved_z%03d.png";
        //start_index = start;
        //end_index = end;
        
        roiSpatialGraph = new AmiraSpatialGraph();
        
        std::cout<<"before"<<std::endl;
        // Read the whole cell graph
        amiraInputGraphReader = new Reader(spatial_graph, "RoiSpatialgraph");//"amira_file_withRad_centerline");  
        amiraInputGraphReader->readSpatialGraphFile(false);
        
        // Get the  transform from the file
        //getTransformationandInverseTransformation(loadhx_file);
        readAmiraTransformations(loadhx_file);
        
        // Apply inverse transform to match the image
        (amiraInputGraphReader->getSpatialGraph())->applyInverseTransformation(inverse_transformation);
         //std::vector< Edge * > * edges = (amiraInputGraphReader->getSpatialGraph())->edgesPointer();
    
        std::cout<<"after"<<std::endl;
        
        
#if 0
        
        //PreDecon::readImageStats(&bricks_x, &bricks_y, &dim_x, &dim_y);

        std::cout << "Inputfilename:" <<  inputfilename << std::endl;

        char input_file[1024];
        strcpy(input_file, inputfilename);      

        //unsigned int bricks, brick_dim_x, brick_dim_y, offset_x, offset_y; 
//         FILE *fp;
//     
        if(chdir(input) != 0)
                  perror("Couldn't open brick drirectory!");

        std::cout<<"Opened Success"<<input<<std::endl;
        

        
        const char * directoryname = "brick";           //DIRECTORY NAMES FOR ACCESSING DIFFERENT BRICKS
        char buf[1024];
        strcpy(buf, directoryname);

        const char * dots = "..";
        char dotbuf[1024];
        strcpy(dotbuf, dots);
        
        int count = 1;
        
        /*for(int i = 0; i < bricks_y; i++)               //FOR ALL BRICKS IN BRAIN SECTION
        for(int j = 0; j < bricks_x; j++)*/
        for(unsigned int i = 0; i < num_bricks; i++)    
        {
                char buffer[1024];
                sprintf(buffer,"%s%02d", buf, i+1);                   //OPEN DIRECTORY OF CURRENT BRICK
                std::cout<<buf<<std::endl;
                if(chdir(buffer) != 0)
                  perror("Couldn't open brick drirectory!");

                std::cout<< "Reading brick =" << i <<  std::endl;
                
                NameGeneratorType::Pointer name_gen = NameGeneratorType::New();         //DECLARE AND INITIALIZE NAME_GENERATOR
                name_gen->SetSeriesFormat( input_file );
                name_gen->SetStartIndex( 0 );
                name_gen->SetEndIndex( num_zslices-1 );
                name_gen->SetIncrementIndex( 1 );
                
                
                                
                SeriesReaderType::Pointer input_reader = SeriesReaderType::New();       //DECLARE AND INITIALIZE INPUT_READER
                input_reader->SetImageIO( itk::PNGImageIO::New() );
                input_reader->SetFileNames( name_gen->GetFileNames() );

                std::cout<<"ok"<<std::endl;
                try
                {
                        input_reader->Update();
                }
                catch( itk::ExceptionObject & err )
                {
                        std::cerr << "BrickReaderExceptionObject caught !" << std::endl;
                        std::cerr << err << std::endl;
                }
                
                

                ImageType::IndexType brick_index;
                ImageType::SizeType brick_size;
                
                ImageType::RegionType brick_region;

                ImageType::Pointer brick_image;
                brick_image = input_reader->GetOutput();
                brick_image->Update();
                
                
#endif
                
                std::cout<<"size"<<std::endl;
                
                std::cout<< input_image->GetLargestPossibleRegion().GetSize(0) << " " << input_image->GetLargestPossibleRegion().GetSize(1) << " " << input_image->GetLargestPossibleRegion().GetSize(2) << std::endl;
                // Determine the number of rois to be split out of this brick
                // if we want split it into 8 rois = 2*2*2
                unsigned int xdelta = (input_image->GetLargestPossibleRegion().GetSize(0)) / (ROI_SIZE_XY / XYSAMPLING);
                unsigned int ydelta = (input_image->GetLargestPossibleRegion().GetSize(1)) / (ROI_SIZE_XY / XYSAMPLING);
                unsigned int zdelta = (input_image->GetLargestPossibleRegion().GetSize(2)) / (ROI_SIZE_Z / ZSAMPLING);
                
                std::cout<< xdelta << " " << ydelta << " " << zdelta << std::endl;
                
                int roi_num = 0;
                
                for(unsigned int z = 0 ; z < zdelta; z++)
                    for(unsigned int y = 0 ; y < ydelta; y++)
                        for(unsigned int x = 0 ; x < xdelta; x++)
                            {
                                std::cout<<"in loop"<<std::endl;
                                roi_num++;
                                // Split the brick image_stats
                                ImageType::IndexType paste_index;                       //INITIALIZE VECTOR VALUES FOR REGION
                                ImageType::SizeType paste_size;

                                //paste_index[0] = j * BRICK_DIMENSIONS;
                                //paste_index[1] = i * BRICK_DIMENSIONS;
                                paste_index[0] = x* (ROI_SIZE_XY / XYSAMPLING + 1);
                                std::cout<<paste_index[0]<<std::endl;
                                paste_index[1] = y* (ROI_SIZE_XY / XYSAMPLING + 1);
                                std::cout<<paste_index[1]<<std::endl;
                                paste_index[2] = z* (ROI_SIZE_Z / ZSAMPLING + 1);
                                std::cout<<paste_index[2]<<std::endl;
                                
                                paste_size[0] = (ROI_SIZE_XY / XYSAMPLING );
                                std::cout<<paste_size[0]<<std::endl;
                                paste_size[1] = (ROI_SIZE_XY / XYSAMPLING );
                                std::cout<<paste_size[1]<<std::endl;
                                paste_size[2] = (ROI_SIZE_Z / ZSAMPLING );
                                std::cout<<paste_size[2]<<std::endl;
                                
                                
                                
                                
                                // Copy the this roi image
                                 
                                // Set the target region start and size
                                ImageType::Pointer roi_image = ImageType::New(); 
                                
                                ImageType::RegionType roi_region;  
                                
                                ImageType::IndexType roi_index;                       //INITIALIZE VECTOR VALUES FOR REGION
                                //ImageType::SizeType roi_size;
                                
                                roi_index[0] = 0;
                                roi_index[1] = 0;
                                roi_index[2] = 0;
                                
                                roi_region.SetSize( paste_size );
                                roi_region.SetIndex( roi_index );
                                
                                roi_image->SetRegions( roi_region );
                                roi_image->Allocate(); 
                                
                                
                                // Set the source size and index
                                
                                ImageType::RegionType source_region;
                                //ImageType::IndexType source_index;                       //INITIALIZE VECTOR VALUES FOR REGION
                                //ImageType::SizeType roi_size;
                                source_region.SetSize( paste_size );
                                source_region.SetIndex( paste_index );
                                
                                std::cout<<"ok allocating"<<std::endl;
                                
                                
                                IteratorType2 og_iterator(input_image, source_region);
                                IteratorType2 output_iterator(roi_image, roi_region);
                                
                                for (og_iterator.GoToBegin(), output_iterator.GoToBegin();!og_iterator.IsAtEnd(); ++og_iterator, ++output_iterator ) 
                                {
                                    unsigned char value = og_iterator.Get();      
                                    output_iterator.Set((unsigned char)value);
                                }
                                
                                std::cout<<"ok copying"<<std::endl;
                                
                                
                                // write this image output_iterator
                                writeProximityImage( roi_image, output_path,  section_num, roi_num, source_region);
                                
                                std::cout<<"end loop"<<std::endl;
                                
                                // Cut pieces of axon tht fall within this roi region 
                                //getSpatialGraphWithinRoi(source_region);
                               
                            }
                            
                            
//                 char dotbuffer[1024];
//                 sprintf(dotbuffer,"%s",dotbuf);
//                 if(chdir(dotbuffer) != 0)
//                   perror("Couldn't switch to Previous Directory!");
//                                 
//                 count++;
        
               
};
        
        


void BrickImageReader::writeProximityImage(ImageType::Pointer image, const char * output_file, int section_num, int index_num, ImageType::RegionType source_region)
{
    std::cout<< "Writing roi image number "<< index_num <<std::endl;
    
    //get the cropped roi from the original image
    
    std::string folder_name = output_file;

    
        
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
        
//     folder_name += "_roi_s";
//     folder_name += Utility::makeString(section_num);
// //     folder_name += "_b";
// //     folder_name += Utility::makeString(brick_num);
     folder_name += "_";
    folder_name += Utility::makeString(index_num);

    std::cout<<"Proximity image folder name: "<<folder_name<<std::endl;   
    
//     std::string mkdir_command = "mkdir " + folder_name;
//     system(mkdir_command.c_str());

    if(mkdir(folder_name.c_str(),
            S_IRWXU| //full access for owner
            S_IRWXG| //full access for grp
            S_IRWXO  //full access for other
            ))
        {
            
           
 
//         if(errno == EEXIST){
//             std::cerr << folder_name << " already exists" << std::endl;
//         }
//         else {
//             perror("Error creating directory");
//             exit(-1);
//         }
        }
        
//         if(chdir(folder_name.c_str()) != 0)
//                   perror("Couldn't open brick drirectory!");
    
//     std::string cd_command = "cd " + folder_name; 
//     system(cd_command.c_str());
    
    std::string output_fileName = folder_name + "/roi"; //+ output_file;
    
        
    writeImagePlanes(image, output_fileName.c_str(),  source_region);
    
    

//     const char * dots = "..";
//     char dotbuf[1024];
//     strcpy(dotbuf, dots);
// 
//     if(chdir(dotbuf) != 0)
//                 perror("Couldn't switch to Previous Directory!");
                
    
//     cd_command = "cd ..";
//     system(cd_command.c_str());
};



void BrickImageReader::writeImagePlanes(ImageType::Pointer input_image, const char* output_file, ImageType::RegionType source_region)
{
        #ifdef DEBUG
          std::cout<< "Generate Writer !!" << std::endl;
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
          std::cout << "-------------------------------------------------------------------" << std::endl;
          std::cout<< "Start Writing !!" << std::endl;
          std::cout << "-------------------------------------------------------------------" << std::endl;
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
          std::cout<< "Writing Done !!" << std::endl;
        #endif
          
        writeInfoFile(output_file, first_slice, last_slice);
	
	writeHXFile(  output_file,  source_region);
          
          

}; 


void BrickImageReader::getSpatialGraphWithinRoi(ImageType::RegionType source_region)
{
    
    // scale and clip it to image boundary
    std::cout<<"ok here"<<std::endl;
    std::vector< Edge * > * edges = (amiraInputGraphReader->getSpatialGraph())->edgesPointer();
    //std::vector< Vertex * > * vertices = (amiraInputGraphReader->getSpatialGraph())->verticesPointer();
    std::cout<<"ok here"<<std::endl;
    std::vector< Edge * >  destedges;// = roiSpatialGraph->edgesPointer();
    std::vector< Vertex * >  destvertices;// = roiSpatialGraph->verticesPointer();
    
    int numOfEdges = edges->size();
    //int numofVertices = vertices->size();
        
    //int pixel_value = 0;
    /*for(int i=0; i<numofVertices; i++)  //for each edge
    {
        
        Vertex * currentVertex = vertices->at(i);
        
        double * coords = (currentVertex)->coordinates;
        coords[X_COORD] = coords[X_COORD] - imageTranslation[X_COORD];
        coords[Y_COORD] = coords[Y_COORD] - imageTranslation[Y_COORD];
        coords[Z_COORD] = coords[Z_COORD] - imageTranslation[Z_COORD];
        
    }*/
    
    for(int i=0; i<numOfEdges; i++)  //for each edge
    {
        
        Edge * currentEdge = edges->at(i);
        std::list< double * >::iterator it;
        std::list< double * > tempEdgePoints;
                
        for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
        {
            double * coords = *it;
            coords[0]  = coords[0]/XYSAMPLING;
            coords[1]  = coords[1]/XYSAMPLING;
            coords[2]  = coords[2]/ZSAMPLING;
            
            //std::cout<<coords[0]<<" "<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
            //std::cout<<source_region.GetIndex(0)<<" "<<" "<<source_region.GetIndex(1)<<" "<<source_region.GetIndex(2)<<std::endl;
            //std::cout<<source_region.GetSize(0)<<" "<<" "<<source_region.GetSize(1)<<" "<<source_region.GetSize(2)<<std::endl;
            
            // if this point is within the roi bounds
            if( ((coords[X_COORD] >= source_region.GetIndex(0)  ) && (coords[X_COORD] <= source_region.GetSize(0))) && 
                ((coords[Y_COORD] >= source_region.GetIndex(1)  ) && (coords[Y_COORD] <= source_region.GetSize(1))) && 
                ((coords[Z_COORD] >= source_region.GetIndex(2)  ) && (coords[Z_COORD] <= source_region.GetSize(2))) )
            {
                
                // Apply image translation before scaling
                coords[X_COORD] = coords[X_COORD] - source_region.GetIndex(0);
                coords[Y_COORD] = coords[Y_COORD] - source_region.GetIndex(1);
                coords[Z_COORD] = coords[Z_COORD] - source_region.GetIndex(2);
            
                std::cout<<"after translation"<<std::endl;
                std::cout<<coords[0]<<" "<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
            
                // Add this one t
                tempEdgePoints.push_back(coords);
                
                std::cout<<"after push"<<std::endl;
                
                /*
                int x_pos = rint( coords[X_COORD] / XYSAMPLING );
                int y_pos = rint( coords[Y_COORD] / XYSAMPLING );
                int z_pos = rint( coords[Z_COORD] / ZSAMPLING);
                
                coords[X_COORD] = x_pos;
                coords[Y_COORD] = y_pos;
                coords[Z_COORD] = z_pos;*/
                
            }
        }
        
        if(tempEdgePoints.size())
        {
            //std::list< double * >::iterator tempEdgePointIterator;
            
            Edge * tmpEdge = new Edge(currentEdge->edgeConnectivity, currentEdge->numEdgePoints, currentEdge->label, tempEdgePoints);
            destedges.push_back(tmpEdge);
        } 
        
        std::cout<<"after new edge push_back"<<std::endl;
                        
        
        // All the relevant points have been added to temp list
        // Now populate the new edges with these points and the other edge attributes
        
    }
    
    std::cout<<"after edge loop"<<std::endl;
                
    if(destedges.size())
    {
        std::vector< Edge * >::iterator edgeIter;
        for(edgeIter = destedges.begin(); edgeIter != destedges.end(); ++edgeIter)
        roiSpatialGraph->addEdge(*edgeIter);
        
        
        std::cout<<"New Edges"<<std::endl;
        // print only
        std::vector< Edge * > * newedges = roiSpatialGraph->edgesPointer();
        
        numOfEdges = newedges->size();
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
            
            Edge * currentEdge = newedges->at(i);
            std::list< double * >::iterator it;
            std::list< double * > tempEdgePoints;
                    
            for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
            {
                std::cout<<(*it)[0]<<" "<<" "<<(*it)[1]<<" "<<(*it)[2]<<std::endl;
                
            }
            
        }
                
        // print only
        std::cout<<"before write "<<std::endl;
        amiraInputGraphReader->writeSpatialGraphFile2(roiSpatialGraph);
        std::cout<<"before exit "<<std::endl;
        
    }
    
    
};
                               

void BrickImageReader::writeInfoFile(const char* output_file, const unsigned int first_slice, const unsigned int last_slice)
{
    std::string format = output_file;
    format += ".info";
    std::ofstream ListData( format.c_str() );

    ListData << "# Amira Stacked Slices"<< std::endl;
    ListData << " "<< std::endl;
    ListData << "pixelsize "<<XYSAMPLING<<"  "<<XYSAMPLING<< std::endl;
    ListData << "       "<<std::endl;
    
    for(unsigned int i = 0 ; i <= last_slice; i++)
    {
        
        std::string imagename = "roi";
         
        ListData <<" "<<imagename<<std::setfill('0')<< std::setw(3)<<i<<".png"<<" "<<first_slice + 0.5*i <<std::endl;
    }
    
    ListData.close();  

};

void BrickImageReader::writeHXFile( const char* outputFilename, ImageType::RegionType source_region)
{
  //writes an HXFile to ensure right positioning of the proximity images in Amira  
    
    
        std::string format2 = outputFilename;
        format2 += "_load.hx";
        
        std::ofstream ListData( format2.c_str() );
        
        

        
        ListData << "# Amira Script"<< std::endl;
        ListData << " "<< std::endl;
        ListData << "[ load ${SCRIPTDIR}/roi.info ] setLabel {roi.info}"<< std::endl;
        ListData << "       "<<std::endl;
        
        ListData << "roi.info"<<" setTranslation " <<source_region.GetIndex(0)*XYSAMPLING<<" "<<source_region.GetIndex(1)*XYSAMPLING<<" "<<source_region.GetIndex(2)*ZSAMPLING<<" "<<std::endl;
        //ListData << "proximity"<<proximity_label<<i<<".info setTranslation " <<list->at(i)[X_COORD] - region_size[0]/2*XYSAMPLING<<" "<<list->at(i)[Y_COORD] - region_size[1]/2*XYSAMPLING<<" "<<list->at(i)[Z_COORD] - region_size[2]/2*ZSAMPLING<<std::endl;
        ListData << "roi.info applyTransform "<<  std::endl;
        ListData << "roi.info setTransform " <<transformation[0][0]<< " "<<transformation[1][0]<< " "<<transformation[2][0]<< " "<<transformation[3][0]<< " "<<transformation[0][1]<< " "<<transformation[1][1]<< " "<<transformation[2][1]<< " "<<transformation[3][1]<< " "<<transformation[0][2]<< " "<<transformation[1][2]<< " "<<transformation[2][2]<< " "<<transformation[3][2]<< " "<<transformation[0][3]<< " "<<transformation[1][3]<< " "<<transformation[2][3]<< " "<<transformation[3][3]<< " "<< std::endl;
        
        
            
        ListData.close(); 
   
        
        
    
    
};

#if 0
bool BrickImageReader::getTransformationandInverseTransformation(const char* inputfilename)
{
    std::cout<<"In getTranformfromfile"<<std::endl;
    // intialize transformation matrices
    
    //this->sectionTranslation = new double *[4];
    //this->sectionRotation = new double *[4];

    
            for(int ii = 0; ii < 4; ++ii)
            {
                    //this->sectionTranslation[ii] = new double[4];
                    //this->sectionRotation[ii] = new double[4];
                    for(int jj = 0; jj < 4; ++jj)
                    {
                            this->sectionTranslation[ii][jj] = 0;
                            this->sectionRotation[ii][jj] = 0;
                    }
                    this->sectionTranslation[ii][ii] = 1;
                    this->sectionRotation[ii][ii] = 1;
            }
        
            //this->transformation = new double *[4];
                    for(int ii = 0; ii < 4; ++ii)
                    {
                            //this->transformation[ii] = new double[4];
                            for(int jj = 0; jj < 4; ++jj)
                            {
                                    if(ii != jj)
                                            this->transformation[ii][jj] = 0;
                                    else
                                            this->transformation[ii][jj] = 1;
                            }
                    }
                
            //this->inverse_transformation = new double *[4];
                    for(int ii = 0; ii < 4; ++ii)
                    {
                            //this->inverse_transformation[ii] = new double[4];
                            for(int jj = 0; jj < 4; ++jj)
                            {
                                    if(ii != jj)
                                            this->inverse_transformation[ii][jj] = 0;
                                    else
                                            this->inverse_transformation[ii][jj] = 1;
                            }
                    }
                    
    //std::cout<<"inputfilename: " << inputfilename <<std::endl;
        std::ifstream inputStream(inputfilename);
     
        if(!inputStream.fail())
        {
                //const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                //const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t\n ";
                
                std::string currentLine;
                //unsigned int line = 0;
                
               
                
                //bool parameters = 1;
                //bool transform = 0;
 //               bool correctSection = 1;
 //               bool correctPrevSection = 0;
 //               int sectionID = 0;
                //unsigned int brackets = 0, transformBrackets = 0;
                //unsigned int currentIndex = 0;
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
                        
                        if(currentLine.size())
                            
                            
                                if(currentLine.find("setTransform ", 0) != std::string::npos)
                                {
                                    
                                              //std::cout << "found correct section transform parameters!" << std::endl;
                                        unsigned int count = 0;
                                        std::string::size_type loc1, loc2, loc3;
                                        loc1 = currentLine.find("setTransform ", 0);
                                        loc1 += 13;// Lenth of setTransform
                                        loc2 = currentLine.find_first_of(signs, 0);
                                        if(loc2 != std::string::npos)
                                                if(loc2 < loc1)
                                                        loc1 = loc2;
                                        loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                        while(loc2 != std::string::npos && count < 16)
                                        {
//                                             std::cout << loc1 << '\t'<<loc2<<std::endl;
                                            
                                            //std::cout<<"May crash here if the compiler is different.. seee comments in code"<<std::endl;
                                            
//##################WARNING#############################################################
                                            // eventhough we are creating the buffer at run time
                                            // the variable still has 4-5 bytes allocated and when we do atof the 
                                            // junk value (if a number) will get converted and hence the read values will be different
                                            // Therefore for this version of the compiler, I am initializing 5 bytes eventhough the variable is supposed to have 
                                            // as many number of bytes as required at runtime.. may cause issues if another compiler does not work like this
                                            
                                                char * tmp = new char[20];
                                                
                                                for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                {
                                                    tmp[i] = 'f';
                                                }
                                                
//                                                 std::cout << tmp<<std::endl;
                                                
                                                currentLine.copy(tmp, loc2 - loc1, loc1);
                                                
//                                                 std::cout << tmp<<std::endl;
                                                
                                                
                                                double ftmp1 = atof(tmp);
                                                
                                                //std::cout << "before!" << std::endl;
                                                //sectionRotation[0][0]= 1;        // amira files are columns after each other
                                                this->sectionRotation[count%4][count/4]= ftmp1;        // amira files are columns after each other
//                                                 std::cout << ftmp1 << '\t'<<tmp<<std::endl;
                                                loc3 = loc2;
                                                loc1 = currentLine.find_first_of(numbers, loc3);
                                                loc2 = currentLine.find_first_of(signs, loc3);
                                                if(loc2 != std::string::npos)
                                                        if(loc2 < loc1)
                                                                loc1 = loc2;
                                                loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                                                ++count;
                                                delete [] tmp;
                                        }
                                        
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        
                                                               
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
                                        //sectionRotation = transformation;
                                        //remove numeric artifacts from z-axis:
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         sectionRotation[2][ii] = 0;
//                                                         sectionRotation[ii][2] = 0;
//                                                 }
//                                                 sectionRotation[2][2] = 1;

                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                this->sectionTranslation[ii][3] = this->sectionRotation[ii][3];
                                                this->sectionRotation[ii][3] = 0;
                                        }
                                        this->sectionRotation[3][3] = 1;
                                        
//                                                 sectionTranslation[2][3] += manual_z_scale[1]/**manual_z_scale[0]*/;
//                                                 std::cout<<"shift : "<<sectionTranslation[2][3]<<std::endl;
                                        
                                        std::cout << "translation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionTranslation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionTranslation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
                                        std::cout << "rotation matrix with scaling:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionRotation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        double square_scaling[3];//myth increased to 3 from 2
                                        double scaling[4];
                                        scaling[3] = 1;
                                        
                                        
                                        
                                        for(int jj = 0; jj< 3; ++jj){
                                        
                                        square_scaling[jj] = this->sectionRotation[0][jj]*this->sectionRotation[0][jj] + this->sectionRotation[1][jj]*this->sectionRotation[1][jj] + this->sectionRotation[2][jj]*this->sectionRotation[2][jj];
                                        scaling[jj] = sqrt(square_scaling[jj]);
                                        std::cout<<"scaling : "<<scaling[jj]<<std::endl;
                                        
                                        }
                                        //scaling[2] = manual_z_scale[0];
                                        //std::cout<<"scaling : "<<scaling[2]<<std::endl;
                                        
                                        for(int ii = 0; ii < 3; ++ii)
                                        {                                                        
                                                for(int jj = 0; jj < 3; ++jj)                                                        
                                                    this->sectionRotation[ii][jj] = this->sectionRotation[ii][jj]/scaling[jj];                                                       
                                        } 
                                        
                                        
                                        std::cout << "rotation matrix without scaling:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionRotation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
//                                                 double ** mInverse = new double *[4];
//                                                 for(int ii = 0; ii < 4; ++ii)
//                                                 {
//                                                         mInverse[ii] = new double[4];
//                                                         for(int jj = 0; jj < 4; ++jj)
//                                                                 mInverse[ii][jj] = 0;
//                                                 }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
//                                                 mInverse[0][3] = -1*(sectionRotation[0][0]*sectionTranslation[0][3] + sectionRotation[1][0]*sectionTranslation[1][3]);
//                                                 mInverse[1][3] = -1*(sectionRotation[0][1]*sectionTranslation[0][3] + sectionRotation[1][1]*sectionTranslation[1][3]);
//                                                 mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                 mInverse[3][3] = 1;
//                                                 
//                                                 this->inverse_transformation = mInverse;
                                        
                                        
                                        
                                        
                                        double ** mProduct = new double *[4];
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                mProduct[ii] = new double[4];
                                                for(int jj = 0; jj < 4; ++jj)
                                                        mProduct[ii][jj] = 0;
                                        }
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                                for(int jj = 0; jj < 4; ++jj)
                                                        for(int kk = 0; kk < 4; ++kk){
                                                            
                                                            //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                                mProduct[ii][jj] += this->sectionTranslation[ii][kk]*this->sectionRotation[kk][jj]*scaling[jj];
                                                        }
                                               
                                        /*for(int ii = 0; ii < 4; ++ii)
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                    
                                                    //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                        this->transformation[ii][jj] = mProduct[ii][jj];
                                        }*/
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                //std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        this->transformation[ii][jj] = mProduct[ii][jj];
                                                }
                                                
                                                
                                        }
                                        //this->transformation = mProduct;
                                        std::cout << "transformation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->transformation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->transformation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                                
                                        }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         transformation[2][ii] = 0;
//                                                         transformation[ii][2] = 0;
//                                                 }
//                                                 transformation[2][2] = 1;

                                        double ** mInverse = new double *[4];
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                mInverse[ii] = new double[4];
                                                for(int jj = 0; jj < 4; ++jj)
                                                        mInverse[ii][jj] = 0;
                                        }
                                        
                                        double a1 = this->sectionRotation[0][0]*this->sectionRotation[1][1]*this->sectionRotation[2][2]; 
                                        double a2 = this->sectionRotation[0][1]*this->sectionRotation[1][2]*this->sectionRotation[2][0];  
                                        double a3 = this->sectionRotation[0][2]*this->sectionRotation[1][0]*this->sectionRotation[2][1];
                                        double a4 = this->sectionRotation[2][0]*this->sectionRotation[1][1]*this->sectionRotation[0][2];
                                        double a5 = this->sectionRotation[2][1]*this->sectionRotation[1][2]*this->sectionRotation[0][0];
                                        double a6 = this->sectionRotation[2][2]*this->sectionRotation[1][0]*this->sectionRotation[0][1];
                                        
                                        double det = (a1 + a2 + a3 -a4 -a5 -a6)*scaling[0]*scaling[1]*scaling[2];
                                        
                                        mInverse[0][0] = (this->sectionRotation[1][1]*this->sectionRotation[2][2] - this->sectionRotation[1][2]*this->sectionRotation[2][1])*scaling[1]*scaling[2]/det;
                                        mInverse[0][1] = (this->sectionRotation[0][2]*this->sectionRotation[2][1] - this->sectionRotation[0][1]*this->sectionRotation[2][2])*scaling[1]*scaling[2]/det;
                                        mInverse[0][2] = (this->sectionRotation[0][1]*this->sectionRotation[1][2] - this->sectionRotation[0][2]*this->sectionRotation[1][1])*scaling[1]*scaling[2]/det;
                                        mInverse[1][0] = (this->sectionRotation[1][2]*this->sectionRotation[2][0] - this->sectionRotation[1][0]*this->sectionRotation[2][2])*scaling[0]*scaling[2]/det;
                                        mInverse[1][1] = (this->sectionRotation[0][0]*this->sectionRotation[2][2] - this->sectionRotation[0][2]*this->sectionRotation[2][0])*scaling[0]*scaling[2]/det;
                                        mInverse[1][2] = (this->sectionRotation[0][2]*this->sectionRotation[1][0] - this->sectionRotation[0][0]*this->sectionRotation[1][2])*scaling[0]*scaling[2]/det;
                                        mInverse[2][0] = (this->sectionRotation[1][0]*this->sectionRotation[2][1] - this->sectionRotation[1][1]*this->sectionRotation[2][0])*scaling[1]*scaling[0]/det;
                                        mInverse[2][1] = (this->sectionRotation[0][1]*this->sectionRotation[2][0] - this->sectionRotation[0][0]*this->sectionRotation[2][1])*scaling[1]*scaling[0]/det;
                                        mInverse[2][2] = (this->sectionRotation[0][0]*this->sectionRotation[1][1] - this->sectionRotation[0][1]*this->sectionRotation[1][0])*scaling[1]*scaling[0]/det;
                                        
                                        
                                        
//                                                mInverse[0][0] = sectionRotation[1][1]/scaling[1]*sectionRotation[2][2]/scaling[2] - sectionRotation[1][2]/scaling[2]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[0][1] = sectionRotation[0][2]/scaling[2]*sectionRotation[2][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[0][2] = sectionRotation[0][1]/scaling[1]*sectionRotation[1][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[1][1]/scaling[1];
//                                                mInverse[1][0] = sectionRotation[1][2]/scaling[2]*sectionRotation[2][0]/scaling[0] - sectionRotation[1][0]/scaling[0]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[1][1] = sectionRotation[0][0]/scaling[0]*sectionRotation[2][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[1][2] = sectionRotation[0][2]/scaling[2]*sectionRotation[1][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[1][2]/scaling[2];
//                                                mInverse[2][0] = sectionRotation[1][0]/scaling[0]*sectionRotation[2][1]/scaling[1] - sectionRotation[1][1]/scaling[1]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[2][1] = sectionRotation[0][1]/scaling[1]*sectionRotation[2][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[2][2] = sectionRotation[0][0]/scaling[0]*sectionRotation[1][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[1][0]/scaling[0];
//                                                
                                        
//                                                for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
                                        //scaling[2] = 1;
                                        mInverse[0][3] = -1*(this->sectionRotation[0][0]/scaling[0]*this->sectionTranslation[0][3] + this->sectionRotation[1][0]/scaling[0]*this->sectionTranslation[1][3] + this->sectionRotation[2][0]/scaling[0]*this->sectionTranslation[2][3]);
                                        mInverse[1][3] = -1*(this->sectionRotation[0][1]/scaling[1]*this->sectionTranslation[0][3] + this->sectionRotation[1][1]/scaling[1]*this->sectionTranslation[1][3] + this->sectionRotation[2][1]/scaling[1]*this->sectionTranslation[2][3]);
                                        mInverse[2][3] = -1*(this->sectionRotation[0][2]/scaling[2]*this->sectionTranslation[0][3] + this->sectionRotation[1][2]/scaling[2]*this->sectionTranslation[1][3] + this->sectionRotation[2][2]/scaling[2]*this->sectionTranslation[2][3]);
                                        mInverse[3][3] = 1;
                                        /*
                                        for(int ii = 0; ii < 4; ++ii)
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                    
                                                    //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                        this->inverse_transformation[ii][jj] = mInverse[ii][jj];
                                        }*/
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                //std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        this->inverse_transformation[ii][jj] = mInverse[ii][jj];
                                                }
                                                
                                                
                                        }
                                        //this->inverse_transformation = mInverse;
                                        
//                                                mInverse[0][3] = -1*(sectionTranslation[0][3]);
//                                                mInverse[1][3] = -1*(sectionTranslation[1][3]);
//                                                mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                mInverse[3][3] = 1;
                                        
                                        //this->inverse_transformation = mInverse;
                                        
                                        std::cout << "inverse transformation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->inverse_transformation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->inverse_transformation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                                
                                        } 
//                                                for(int ii = 0; ii < 3; ++ii)
//                                                 {
//                                                         
//                                                         for(int jj = 0; jj < 3; ++jj){
//                                                             
//                                                             mIn 
//                                                             
//                                                         }
//                                                 }
                                }
                       
                }
        }
        std::cout << "out of transformReading" << std::endl;
        
        inputStream.close();
        return 1;
    
}
#endif

bool BrickImageReader::readAmiraTransformations(const char * inputfilename)
{
    //reads Amira transformations from spatial graph file in order to transform the merged graph into it's original position in to the current section's image planes
  
    

#ifdef DEBUG    
std::cout<<"in readAmiraTransformations"<<std::endl;
#endif

	    for(int ii = 0; ii < 4; ++ii)
            {
                    //this->sectionTranslation[ii] = new double[4];
                    //this->sectionRotation[ii] = new double[4];
                    for(int jj = 0; jj < 4; ++jj)
                    {
                            this->sectionTranslation[ii][jj] = 0;
                            this->sectionRotation[ii][jj] = 0;
                    }
                    this->sectionTranslation[ii][ii] = 1;
                    this->sectionRotation[ii][ii] = 1;
            }
        
            //this->transformation = new double *[4];
	    for(int ii = 0; ii < 4; ++ii)
	    {
		    //this->transformation[ii] = new double[4];
		    for(int jj = 0; jj < 4; ++jj)
		    {
			    if(ii != jj)
				    this->transformation[ii][jj] = 0;
			    else
				    this->transformation[ii][jj] = 1;
		    }
	    }
	
	    //this->inverse_transformation = new double *[4];
	    for(int ii = 0; ii < 4; ++ii)
	    {
		    //this->inverse_transformation[ii] = new double[4];
		    for(int jj = 0; jj < 4; ++jj)
		    {
			    if(ii != jj)
				    this->inverse_transformation[ii][jj] = 0;
			    else
				    this->inverse_transformation[ii][jj] = 1;
		    }
	    }


        std::cout<<"inputfilename: " << inputfilename <<std::endl;
        std::ifstream inputStream(inputfilename);
     
        if(!inputStream.fail())
        {
                const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t ";
                
                std::string currentLine;
                unsigned int line = 0;
                
               
                
                bool parameters = 1;
                bool transform = 0;
 //               bool correctSection = 1;
 //               bool correctPrevSection = 0;
 //               int sectionID = 0;
                unsigned int brackets = 0, transformBrackets = 0;
                unsigned int currentIndex = 0;
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
                        
                        if(currentLine.size())
                                if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
                                        {
                                              std::cout << "found correct section transform parameters!" << std::endl;
                                                unsigned int count = 0;
                                                std::string::size_type loc1, loc2, loc3;
                                                loc1 = currentLine.find_first_of(numbers, 0);
                                                loc2 = currentLine.find_first_of(signs, 0);
                                                if(loc2 != std::string::npos)
                                                        if(loc2 < loc1)
                                                                loc1 = loc2;
                                                loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                                while(loc2 != std::string::npos && count < 16)
                                                {
                                                        char * tmp1 = new char[loc2 - loc1];
                                                        currentLine.copy(tmp1, loc2 - loc1, loc1);
                                                        double ftmp1 = atof(tmp1);
                                                        sectionRotation[count%4][count/4]= ftmp1;        // amira files are columns after each other
                                                        loc3 = loc2;
                                                        loc1 = currentLine.find_first_of(numbers, loc3);
                                                        loc2 = currentLine.find_first_of(signs, loc3);
                                                        if(loc2 != std::string::npos)
                                                                if(loc2 < loc1)
                                                                        loc1 = loc2;
                                                        loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                                                        ++count;
                                                        delete [] tmp1;
                                                }
                                              
                                                //sectionRotation = transformation;
                                                //remove numeric artifacts from z-axis:
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         sectionRotation[2][ii] = 0;
//                                                         sectionRotation[ii][2] = 0;
//                                                 }
//                                                 sectionRotation[2][2] = 1;

                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        sectionTranslation[ii][3] = sectionRotation[ii][3];
                                                        sectionRotation[ii][3] = 0;
                                                }
                                                sectionRotation[3][3] = 1;
                                                
//                                                 sectionTranslation[2][3] += manual_z_scale[1]/**manual_z_scale[0]*/;
//                                                 std::cout<<"shift : "<<sectionTranslation[2][3]<<std::endl;
                                                
                                                std::cout << "translation matrix:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                if(jj < 3)
                                                                        std::cout << sectionTranslation[ii][jj] << ",\t";
                                                                else
                                                                        std::cout << sectionTranslation[ii][jj];
                                                        }
                                                        std::cout << "]" << std::endl;
                                                }
                                                
                                                std::cout << "rotation matrix with scaling:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                if(jj < 3)
                                                                        std::cout << sectionRotation[ii][jj] << ",\t";
                                                                else
                                                                        std::cout << sectionRotation[ii][jj];
                                                        }
                                                        std::cout << "]" << std::endl;
                                                }
                                                double square_scaling[2];
                                                double scaling[4];
                                                scaling[3] = 1;
                                                
                                                
                                                
                                                for(int jj = 0; jj< 3; ++jj){
                                                
                                                square_scaling[jj] = sectionRotation[0][jj]*sectionRotation[0][jj] + sectionRotation[1][jj]*sectionRotation[1][jj] + sectionRotation[2][jj]*sectionRotation[2][jj];
                                                scaling[jj] = sqrt(square_scaling[jj]);
                                                std::cout<<"scaling : "<<scaling[jj]<<std::endl;
                                                
                                                }
                                                //scaling[2] = manual_z_scale[0];
                                                //std::cout<<"scaling : "<<scaling[2]<<std::endl;
                                                
                                                for(int ii = 0; ii < 3; ++ii)
                                                {                                                        
                                                        for(int jj = 0; jj < 3; ++jj)                                                        
                                                            sectionRotation[ii][jj] = sectionRotation[ii][jj]/scaling[jj];                                                       
                                                } 
                                                
                                                
                                                std::cout << "rotation matrix without scaling:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                if(jj < 3)
                                                                        std::cout << sectionRotation[ii][jj] << ",\t";
                                                                else
                                                                        std::cout << sectionRotation[ii][jj];
                                                        }
                                                        std::cout << "]" << std::endl;
                                                }
                                                
//                                                 double ** mInverse = new double *[4];
//                                                 for(int ii = 0; ii < 4; ++ii)
//                                                 {
//                                                         mInverse[ii] = new double[4];
//                                                         for(int jj = 0; jj < 4; ++jj)
//                                                                 mInverse[ii][jj] = 0;
//                                                 }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
//                                                 mInverse[0][3] = -1*(sectionRotation[0][0]*sectionTranslation[0][3] + sectionRotation[1][0]*sectionTranslation[1][3]);
//                                                 mInverse[1][3] = -1*(sectionRotation[0][1]*sectionTranslation[0][3] + sectionRotation[1][1]*sectionTranslation[1][3]);
//                                                 mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                 mInverse[3][3] = 1;
//                                                 
//                                                 this->inverse_transformation = mInverse;
                                                
                                                
                                                
                                                
                                                double ** mProduct = new double *[4];
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        mProduct[ii] = new double[4];
                                                        for(int jj = 0; jj < 4; ++jj)
                                                                mProduct[ii][jj] = 0;
                                                }
                                                
                                                for(int ii = 0; ii < 4; ++ii)
                                                        for(int jj = 0; jj < 4; ++jj)
                                                                for(int kk = 0; kk < 4; ++kk){
                                                                    
                                                                    //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                                        mProduct[ii][jj] += sectionTranslation[ii][kk]*sectionRotation[kk][jj]*scaling[jj];
                                                                }
                                                                
                                                //this->transformation = mProduct;
                                                for(int ii = 0; ii < 4; ++ii)
						{
							//std::cout << "[";
							for(int jj = 0; jj < 4; ++jj)
							{
								this->transformation[ii][jj] = mProduct[ii][jj];
							}
							
							
						}
						
						
                                                std::cout << "transformation matrix:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                if(jj < 3)
                                                                        std::cout << transformation[ii][jj] << ",\t";
                                                                else
                                                                        std::cout << transformation[ii][jj];
                                                        }
                                                        std::cout << "]" << std::endl;
                                                        
                                                }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         transformation[2][ii] = 0;
//                                                         transformation[ii][2] = 0;
//                                                 }
//                                                 transformation[2][2] = 1;

                                               double ** mInverse = new double *[4];
                                               
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        mInverse[ii] = new double[4];
                                                        for(int jj = 0; jj < 4; ++jj)
                                                                mInverse[ii][jj] = 0;
                                                }
                                                
                                               double a1 = sectionRotation[0][0]*sectionRotation[1][1]*sectionRotation[2][2]; 
                                               double a2 = sectionRotation[0][1]*sectionRotation[1][2]*sectionRotation[2][0];  
                                               double a3 = sectionRotation[0][2]*sectionRotation[1][0]*sectionRotation[2][1];
                                               double a4 = sectionRotation[2][0]*sectionRotation[1][1]*sectionRotation[0][2];
                                               double a5 = sectionRotation[2][1]*sectionRotation[1][2]*sectionRotation[0][0];
                                               double a6 = sectionRotation[2][2]*sectionRotation[1][0]*sectionRotation[0][1];
                                               
                                               double det = (a1 + a2 + a3 -a4 -a5 -a6)*scaling[0]*scaling[1]*scaling[2];
                                               
                                               mInverse[0][0] = (sectionRotation[1][1]*sectionRotation[2][2] - sectionRotation[1][2]*sectionRotation[2][1])*scaling[1]*scaling[2]/det;
                                               mInverse[0][1] = (sectionRotation[0][2]*sectionRotation[2][1] - sectionRotation[0][1]*sectionRotation[2][2])*scaling[1]*scaling[2]/det;
                                               mInverse[0][2] = (sectionRotation[0][1]*sectionRotation[1][2] - sectionRotation[0][2]*sectionRotation[1][1])*scaling[1]*scaling[2]/det;
                                               mInverse[1][0] = (sectionRotation[1][2]*sectionRotation[2][0] - sectionRotation[1][0]*sectionRotation[2][2])*scaling[0]*scaling[2]/det;
                                               mInverse[1][1] = (sectionRotation[0][0]*sectionRotation[2][2] - sectionRotation[0][2]*sectionRotation[2][0])*scaling[0]*scaling[2]/det;
                                               mInverse[1][2] = (sectionRotation[0][2]*sectionRotation[1][0] - sectionRotation[0][0]*sectionRotation[1][2])*scaling[0]*scaling[2]/det;
                                               mInverse[2][0] = (sectionRotation[1][0]*sectionRotation[2][1] - sectionRotation[1][1]*sectionRotation[2][0])*scaling[1]*scaling[0]/det;
                                               mInverse[2][1] = (sectionRotation[0][1]*sectionRotation[2][0] - sectionRotation[0][0]*sectionRotation[2][1])*scaling[1]*scaling[0]/det;
                                               mInverse[2][2] = (sectionRotation[0][0]*sectionRotation[1][1] - sectionRotation[0][1]*sectionRotation[1][0])*scaling[1]*scaling[0]/det;
                                               
                                               
                                                
//                                                mInverse[0][0] = sectionRotation[1][1]/scaling[1]*sectionRotation[2][2]/scaling[2] - sectionRotation[1][2]/scaling[2]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[0][1] = sectionRotation[0][2]/scaling[2]*sectionRotation[2][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[0][2] = sectionRotation[0][1]/scaling[1]*sectionRotation[1][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[1][1]/scaling[1];
//                                                mInverse[1][0] = sectionRotation[1][2]/scaling[2]*sectionRotation[2][0]/scaling[0] - sectionRotation[1][0]/scaling[0]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[1][1] = sectionRotation[0][0]/scaling[0]*sectionRotation[2][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[1][2] = sectionRotation[0][2]/scaling[2]*sectionRotation[1][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[1][2]/scaling[2];
//                                                mInverse[2][0] = sectionRotation[1][0]/scaling[0]*sectionRotation[2][1]/scaling[1] - sectionRotation[1][1]/scaling[1]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[2][1] = sectionRotation[0][1]/scaling[1]*sectionRotation[2][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[2][2] = sectionRotation[0][0]/scaling[0]*sectionRotation[1][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[1][0]/scaling[0];
//                                                
                                               
//                                                for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
                                               //scaling[2] = 1;
                                                mInverse[0][3] = -1*(sectionRotation[0][0]/scaling[0]*sectionTranslation[0][3] + sectionRotation[1][0]/scaling[0]*sectionTranslation[1][3] + sectionRotation[2][0]/scaling[0]*sectionTranslation[2][3]);
                                                mInverse[1][3] = -1*(sectionRotation[0][1]/scaling[1]*sectionTranslation[0][3] + sectionRotation[1][1]/scaling[1]*sectionTranslation[1][3] + sectionRotation[2][1]/scaling[1]*sectionTranslation[2][3]);
                                                mInverse[2][3] = -1*(sectionRotation[0][2]/scaling[2]*sectionTranslation[0][3] + sectionRotation[1][2]/scaling[2]*sectionTranslation[1][3] + sectionRotation[2][2]/scaling[2]*sectionTranslation[2][3]);
                                                mInverse[3][3] = 1;
                                                
                                                
                                                //this->inverse_transformation = mInverse;
                                                for(int ii = 0; ii < 4; ++ii)
						{
							//std::cout << "[";
							for(int jj = 0; jj < 4; ++jj)
							{
								this->inverse_transformation[ii][jj] = mInverse[ii][jj];
							}
							
							
						}
//                                                mInverse[0][3] = -1*(sectionTranslation[0][3]);
//                                                mInverse[1][3] = -1*(sectionTranslation[1][3]);
//                                                mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                mInverse[3][3] = 1;
                                               
                                               //this->inverse_transformation = mInverse;
                                               
                                               std::cout << "inverse transformation matrix:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                if(jj < 3)
                                                                        std::cout << inverse_transformation[ii][jj] << ",\t";
                                                                else
                                                                        std::cout << inverse_transformation[ii][jj];
                                                        }
                                                        std::cout << "]" << std::endl;
                                                        
                                                } 
//                                                for(int ii = 0; ii < 3; ++ii)
//                                                 {
//                                                         
//                                                         for(int jj = 0; jj < 3; ++jj){
//                                                             
//                                                             mIn 
//                                                             
//                                                         }
//                                                 }
                                        }
                       
                }
        }
        
        inputStream.close();
        return 0;
}



        