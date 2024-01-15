/*
After entering the dimension in x y and z and the stackname, 
which will be carried through the whole working pipeline, the 
2D Planes are read, cut into the specified size (dim x/y), its 
grayvalues are inverted and all planes are saved as a 3D Tif stack!

Date Modified: 14.2.2013
*/

#include "pre_decon.h"
#include "typedefs.h"


  PreDecon::PreDecon(char* inputfilename, char* outputfilename, int max_z)
  { 
    bool info_written = false;

    unsigned int x = 0, y = 0, Dimension = 2;

    char input_buffer[1024], output_buffer[1024], buf[1024];

    const char * directoryname = "brick";

    strcpy(buf, directoryname);	
	
   
    //------------------------------------------------------------------------------------
//#pragma omp parallel for
      for(int i = 0; i <= max_z; i++)							//for each image plane
      { 
	sprintf(input_buffer, inputfilename,i);
	sprintf(output_buffer, outputfilename,i);
	
	const char * inputfilename = input_buffer;
	const char * outputfilename = output_buffer;

	dim2_ReaderType::Pointer dim2_reader = dim2_ReaderType::New();
	dim2_WriterType::Pointer dim2_writer = dim2_WriterType::New();
  
	dim2_reader->SetFileName( inputfilename );
  #ifdef DEBUG
	//std::cout << "Konverting in tif and Inverting of file: " << inputfilename << std::endl;
  #endif
	dim2_reader->Update();
    
	dim2_ImageType::RegionType input_region;

	dim2_ImageType::RegionType::IndexType input_start;
	dim2_ImageType::RegionType::SizeType output_size, input_size;

	unsigned int dim_x = 0, dim_y = 0, brick_dim_x = 0, brick_dim_y = 0; 

	dim_x = dim2_reader->GetOutput()->GetRequestedRegion().GetSize(0);
	dim_y = dim2_reader->GetOutput()->GetRequestedRegion().GetSize(1);

	unsigned int bricks = 0;
	unsigned int offset_x, offset_y;
	
	if(dim_x % MAX_DECON_SIZE == 0)
	    bricks_x = dim_x / MAX_DECON_SIZE;
	else
	    bricks_x = (dim_x / MAX_DECON_SIZE) + 1;
	
	
	if(dim_y % MAX_DECON_SIZE == 0)
	    bricks_y = dim_y / MAX_DECON_SIZE;
	else
	    bricks_y = (dim_y / MAX_DECON_SIZE) + 1;
	
	if( !info_written )
	{
	  writeImageStats(bricks_x, bricks_y, dim_x, dim_y);
	  info_written = true;
	}

	x = bricks_x;
	y = bricks_y;

	brick_dim_x = dim_x / bricks_x;
	brick_dim_y = dim_y / bricks_y;
	offset_x = dim_x % bricks_x;
	offset_y = dim_y % bricks_y;
	
	//std::cout << brick_dim_x << " " << brick_dim_y << " " << dim_x << " " << dim_y << " " << offset_x << " " << offset_y << std::endl;

	bricks = bricks_x * bricks_y;

	if(i == 0)
	  {
	    for(int j = 1; j <= bricks; j++)
	      {
		char buffer[1024];
		sprintf(buffer,"%s%02d", buf, j);
		if(mkdir(buffer,
			S_IRWXU| //full access for owner
			S_IRWXG| //full access for grp
			S_IROTH|S_IXOTH  //read/chdir access for other
			))
		  {
		    if(errno == EEXIST){
		      std::cerr << buffer << " already exists" << std::endl;
		    }
		    else {
		      perror("Error creating directory");
		      exit(-1);
		    }
		  }
	      }
	  }
		    
	input_start[0] = 0;
	input_start[1] = 0;

	int count = 1;

	for(int k = 0; k < bricks_y; k++)
	  for(int l = 0; l < bricks_x; l++)
	  {	
	    
	    dim2_ImageType::RegionType output_region;
	    dim2_ImageType::RegionType::IndexType output_start;
	
	    
	      //calculate dimensions of brick
	      if(k == 0 && l == 0)				//top left corner
	      {
		input_start[0] = 0;
	        input_start[1] = 0;
		
		output_start[0] = 50;
		output_start[1] = 50;
		
		input_size[0] = brick_dim_x + 50;
		input_size[1] = brick_dim_y + 50;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + 100;
	      }
	      else if(k == 0 && l != bricks_x - 1) 			 //top edge
	      {
		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = 0;
		
		output_start[0] = 0;
		output_start[1] = 50;
		
		input_size[0] = brick_dim_x + 100;
		input_size[1] = brick_dim_y + 50;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + 100;
	      }
	      else if(k != bricks_y - 1 && l == 0) 			 //left edge
	      {
		input_start[0] = 0;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 50;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + 50;
		input_size[1] = brick_dim_y + 100;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + 100;
	      }
	      else if(k == 0 && l == bricks_x - 1) 			 //top right corner
	      {
		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = 0;
		
		output_start[0] = 0;
		output_start[1] = 50;
		
		input_size[0] = brick_dim_x + offset_x + 50;
		input_size[1] = brick_dim_y + 50;
		
		output_size[0] = brick_dim_x + offset_x + 50;
		output_size[1] = brick_dim_y + 100;
	      }
	      else if(k == bricks_y - 1 && l == 0) 			 //bottom left corner
	      {
		input_start[0] = 0;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 50;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + 50;
		input_size[1] = brick_dim_y + offset_y + 50;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + offset_y + 50;
	      }
	      else if(k != bricks_y - 1 && l == bricks_x - 1 )   //right edge
	      {

		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 0;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + offset_x + 50;
		input_size[1] = brick_dim_y + 100;
		
		output_size[0] = brick_dim_x + offset_x + 50;
		output_size[1] = brick_dim_y + 100;
	      }
	      else if(k == bricks_y - 1 && l != bricks_x - 1 )   //bottom edge
	      {
		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 0;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + 100;
		input_size[1] = brick_dim_y + offset_y + 50;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + offset_y + 50;
	      }
	      else if(k == bricks_y - 1 && l == bricks_x - 1 )    //bottom right corner
	      {
		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 0;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + offset_x + 50;
		input_size[1] = brick_dim_y + offset_y + 50;
		
		output_size[0] = brick_dim_x + offset_x + 50;
		output_size[1] = brick_dim_y + offset_y + 50;
	      }
	      else						  //center bricks
	      {
		input_start[0] = l*brick_dim_x - 50;
	        input_start[1] = k*brick_dim_y - 50;
		
		output_start[0] = 0;
		output_start[1] = 0;
		
		input_size[0] = brick_dim_x + 100;
		input_size[1] = brick_dim_y + 100;
		
		output_size[0] = brick_dim_x + 100;
		output_size[1] = brick_dim_y + 100;
		
	      }
	      
	    
	      input_region.SetIndex( input_start );
	      input_region.SetSize( input_size );	    
	    
	      output_region.SetSize( input_size );
	      output_region.SetIndex( output_start );
	      
	      
	      

	      dim2_ImageType::Pointer output_image = dim2_ImageType::New();
	      dim2_ImageType::RegionType output_image_region;
	      dim2_ImageType::IndexType output_image_index;
	      
	      output_image_index[0] = 0;
	      output_image_index[1] = 0;
	      output_image_region.SetIndex( output_image_index );
	      output_image_region.SetSize( output_size );
		    
	      output_image->SetRegions( output_image_region );
	      const dim2_ImageType::SpacingType& spacing = dim2_reader->GetOutput()->GetSpacing();
	      const dim2_ImageType::PointType& input_origin = dim2_reader->GetOutput()->GetOrigin();
	      double output_origin[ Dimension ];
		    
	      for( unsigned int i = 0; i < Dimension; i++)
		{
		  output_origin[i] = input_origin[i] + spacing[i] * input_start[i];
		}
		

	      output_image->SetSpacing( spacing );
	      output_image->SetOrigin( output_origin );
	      output_image->Allocate();
	      output_image->FillBuffer(0);

	      DeconConstIteratorType input_it( dim2_reader->GetOutput(), input_region );
	      DeconIteratorType output_it( output_image, output_region );

	      for( input_it.GoToBegin(), output_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it,++output_it)
		{
		  output_it.Set( input_it.Get() );
		}

	      
	      sprintf(norm_buffer,"%s%02d/norm_%s", buf, count, output_buffer);

	      dim2_writer->SetFileName( norm_buffer );
	      dim2_writer->SetInput( output_image );

	      try
	      {
			    dim2_writer->Update();
	      }
	      catch( itk::ExceptionObject & err )
	      {
		  //std::cout << "ExceptionObject caught !" << std::endl;
		  std::cerr << err << std::endl;
		  //return EXIT_FAILURE;
	      }
	      
	      
	      InvertFilterType::Pointer invert_filter = InvertFilterType::New();
	  
	      invert_filter->SetMaximum( 255 );
	      invert_filter->SetInput( output_image );
	      
	      output_image = invert_filter->GetOutput();
	      output_image->Update();
		
	      sprintf(inv_buffer,"%s%02d/inv_%s", buf, count, output_buffer);

	      dim2_writer->SetFileName( inv_buffer );
	      dim2_writer->SetInput( output_image );

	      try
	      {
			    dim2_writer->Update();
	      }
	      catch( itk::ExceptionObject & err )
	      {
		  //std::cout << "ExceptionObject caught !" << std::endl;
		  std::cerr << err << std::endl;
		  //return EXIT_FAILURE;
	      }
		
	      count++;
	    }
      }
      
    sprintf(norm_buffer,"norm_%s", outputfilename);
    sprintf(inv_buffer,"inv_%s", outputfilename);
      
    //std::cout<< "x " << x << std::endl;
    //std::cout<< "y " << y << std::endl;
    
  };
  
  void PreDecon::writeImageStats(int bricks_x, int bricks_y, int dim_x, int dim_y)
  {
    std::string format = "image_stats.txt";
    std::ofstream ImageStats( format.c_str() );
    
    ImageStats << bricks_x << " " << bricks_y << " " << dim_x << " " << dim_y << std::endl;
    
  };
  
  void PreDecon::readImageStats(unsigned int* bricks_x, unsigned int* bricks_y, unsigned int* dim_x, unsigned int* dim_y)
  {
    FILE *fp;
    
    fp=fopen("image_stats.txt", "r");
    
    if(fp == NULL)
    { 
      perror("Error opening 'image_stats.txt'\nFile formatting: bricks_x bricks_y pixels_x pixels_y\n");
      exit(-1);
    }
    
    fscanf(fp, "%d %d %d %d ", bricks_x, bricks_y, dim_x, dim_y);
    
  };
  
  
  
  

