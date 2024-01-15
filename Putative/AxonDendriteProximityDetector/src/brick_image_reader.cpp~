#include "brick_image_reader.h"
#include "pre_decon.h"

BrickImageReader::BrickImageReader(const char * input, unsigned int start, unsigned int end)
{
    input_image = ImageType::New();
        inputfilename  = input;
        start_index = start;
        end_index = end;
        
        PreDecon::readImageStats(&bricks_x, &bricks_y, &dim_x, &dim_y);

        //std::cout << "Inputfilename:" <<  inputfilename << std::endl;

        char input_file[1024];
        strcpy(input_file, inputfilename);      

        unsigned int bricks, brick_dim_x, brick_dim_y, offset_x, offset_y; 
//         FILE *fp;
//     
//     fp=fopen("image_stats.txt", "r");
//     
//     if(fp == NULL)
//     { 
//       perror("Error opening 'image_stats.txt'\nFile formatting: bricks_x bricks_y pixels_x pixels_y\n");
//       exit(-1);
//     }
//     
//     fscanf(fp, "%d %d %d %d ", bricks_x, bricks_y, dim_x, dim_y);
        
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
        
        
        //std::cout << "input_size: " << input_size[0] << " " << input_size[1] << " " << input_size[2] << std::endl;
        
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

                //std::cout<< "Reading brick y=" << i << " x=" << j << std::endl;
                
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
                char dotbuffer[1024];
                sprintf(dotbuffer,"%s",dotbuf);
                if(chdir(dotbuffer) != 0)
                  perror("Couldn't switch to Previous Directory!");
                                
                count++;
               
        };
        
        
}
        