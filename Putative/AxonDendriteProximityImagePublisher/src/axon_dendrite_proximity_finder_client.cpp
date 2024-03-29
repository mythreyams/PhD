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


#include "axon_dendrite_proximity_finder.h"

using namespace std;



//*******************************************************

ImageType::Pointer readOriginalImage(char* file_name, int start_index, int end_index);


#if 0
void writeStatsXls(char * ouptputStatsXlsPath, int section_number, int roi_number, int axon_num, int dend_num, int num_boutons, int num_putatives)
{
    

    bool file_exsists = 0;
    std::string s1("PutativeContacts.xls");
    std::string s2("PutativeContacts.csv");
    
    ////////////////std:://cout<<"In writexls"<<std::endl;
    
    if(chdir(ouptputStatsXlsPath) != 0)
            perror("Couldn't open image drirectory!");
    
    // Iterate over the params list and write each param of each bouton one by one
    std::list< bouton_params * >::iterator paramIterator;
    //std::ofstream xlsParamWriter("BoutonParams.xls",+a);
    
    std::ifstream f1(s1.c_str());
    if (f1.good()) {
        f1.close();
        file_exsists = true;
    } else {
        f1.close();
        file_exsists = false;
    } 
    
    std::ifstream f2(s2.c_str());
    if (f2.good()) {
        f2.close();
        file_exsists = true;
    } else {
        f2.close();
        file_exsists = false;
    } 
    
    
    std::fstream xlsParamWriter;
    std::fstream csvParamWriter;
    
    
    xlsParamWriter.open ("PutativeContacts.xls", std::fstream::in | std::fstream::out | std::fstream::app);
    csvParamWriter.open ("PutativeContacts.csv", std::fstream::in | std::fstream::out | std::fstream::app);
    
    // write header if it does not exists
    if(!file_exsists)
    {
        xlsParamWriter<<"Section Number"<<'\t'<<"ROI Number"<<'\t'
        <<"From Axon of Cell #"<<'\t'<<"To Dendrite of Cell #"<<'\t'
        <<"Number of Boutons"<<'\t'<<"Number of Putative Contacts"<<'\n';
        
         csvParamWriter<<"Section Number"<<","<<"ROI Number"<<","
        <<"From Axon of Cell #"<<","<<"To Dendrite of Cell #"<<","
        <<"Number of Boutons"<<","<<"Number of Putative Contacts"<<'\n';
    }
    
    
    //for(int i = 0; i<NumOfLandmarks; i++)
    {
        //////////////////std:://cout<< (*paramIterator)->landmark[0]<<std::endl;
        
        xlsParamWriter<<section_number<<'\t'<<roi_number<<'\t'
        <<axon_num<<'\t'<<dend_num<<'\t'<<num_boutons<<'\t'<<num_putatives<<'\n';
        
        
       
        
        std::fstream xlsParamWriter;
        csvParamWriter<<section_number<<","<<roi_number<<","
        <<axon_num<<","<<dend_num<<","<<num_boutons<<","<<num_putatives<<'\n';
        
    }
    
    
    xlsParamWriter.close();
    csvParamWriter.close();
    
    
}
#endif

//*******************************************************

int main( int argc , char * argv[])
{
  Reader * amiraReader;
  Reader * mergedAmiraReader;
  AmiraSpatialGraph * merged_graph;
  AmiraSpatialGraph * input_graph;
  ImageType::Pointer input_image;
//   ImageType::Pointer nondecon_image;
  Image2DType::Pointer image_projection;
  

  std::cout<<"argc "<<argc<<std::endl;
  std::cout<<"Prox Image Size = "<<MAX_PROX_BOX_SIZE*2<<std::endl;
  
  
  //Validate command line arguements
  if(argc != 13 /*&& argc != 9  && argc != 10 && argc != 7 */)
  {
      //////std::cout << std::endl << "Proper format for whole cell:  ./ProximityFinder 'amira_file_with_transformations.am' 'deconvolved_image%03d.png' 'output_name' startIndex endIndex includeUnknownSegments  'amira_file_whole_cell.am'  'nondecon_images.png'" << std::endl << std::endl;
      //////std::cout << std::endl << "Proper format for whole cell with manual z scale :  ./ProximityFinder 'amira_file_with_transformations.am' 'deconvolved_image%03d.png' 'output_name' startIndex endIndex includeUnknownSegments maxDist (default 5) 'amira_file_whole_cell.am' additional_scaling.txt 'nondecon_images.png'" << std::endl << std::endl;
      //////std::cout << std::endl << "Proper format for single section:  ./ProximityFinder 'amira_file.am' 'deconvolved_image%03d.png' 'output_name' startIndex endIndex includeUnknownSegments pointRange maxDist 'nondecon_images.png'" << std::endl << std::endl;
      ////////std::cout << std::endl << "For bricks: ./ProximityFinder 'amira_file.am' 'dendrite_image%03d.png' 'output_name' startIndex endIndex" << std::endl << std::endl;
      return -1;
      // " For bricks: ./ProximityFinder 'amira_file.am' 'input_filename' 'dendrite_image%03d.png' 'output_name' startIndex endIndex'
      

  } 
 
    AxonDendriteProximityFinder * prox_finder;
    
    if(argc == 13)
    {
      std::cout<<"in main "<<argc<<std::endl;
      //for whole merged cell files without scaling errors  
        
      //argv[1] : amira transformation
      //argv[2] : deconvolved image files
      //argv[3] : output_name
      //argv[4] : image index start_index
      //argv[5] : image index end
      //argv[6] : include unknown: y/n?     
      //argv[7] : merged graph file      
      //argv[8] : nondeconvolved images
        const char * amiraTransformationName = argv[1];
	cout<< "amiraTransformationName" << amiraTransformationName<< endl;
        const char * deconvolvedImageFileName = argv[2];
        const char * outputFilename = argv[3];
        int imagestartindex = atoi(argv[4]);
        int imageendindex = atoi(argv[5]);
        int includeunknown = atoi(argv[6]);
        const char * amiraMergedGraphName = argv[7];
//        const char * nondeconvolvedImageFileName = argv[8];
        bool realBild = atoi(argv[8]);
        const char * landmarkfilename = argv[9];
        const char * xlsinputfilename = argv[10];
	const char * xlsoututfilename = argv[11];
        int section_number = atoi(argv[12]);
        
                //Read in amira file and decon image
                //////std::cout << "Begin Reading" << std::endl;  
                //std::cout<< "amiraTransformationName " << amiraTransformationName<<std::endl;
//                  amiraReader = new Reader(amiraTransformationName);  
//                  amiraReader->readSpatialGraphFile(false);
//                  input_graph = amiraReader->getSpatialGraph();
		 
		 
		 
//                 
// 		
//                 mergedAmiraReader = new Reader(amiraMergedGraphName);  
//                 mergedAmiraReader->readSpatialGraphFile(false);
//                 merged_graph = mergedAmiraReader->getSpatialGraph();

                //image_projection = reader->GetZProjection(1);
                                
                prox_finder = new AxonDendriteProximityFinder(landmarkfilename, xlsinputfilename, xlsoututfilename, section_number, deconvolvedImageFileName, imagestartindex, imageendindex, amiraTransformationName, amiraMergedGraphName, includeunknown,outputFilename/*, nondecon_image*/);
//                const char* inputfilename = argv[2];
//                prox_finder->readAmiraSectionTransformations();
//                 std::string folder_name = argv[3];
//                 
//                 if(mkdir(folder_name.c_str(),
//                         S_IRWXU| //full access for owner
//                         S_IRWXG| //full access for grp
//                         S_IRWXO  //full access for other
//                         ))
//                     {
//             //         if(errno == EEXIST){
//             //             std::cerr << folder_name << " already exists" << std::endl;
//             //         }
//             //         else {
//             //             perror("Error creating directory");
//             //             exit(-1);
//             //         }
//                     }
//                 std::string filenames = folder_name;
//                 filenames += "/";
//                 filenames += folder_name;
//                 
//                 
//                 
//                 char * filename = new char[filenames.size() + 1];
//                 std::copy(filenames.begin(), filenames.end(), filename);
//                 filename[filenames.size()] = '\0';
                
//                prox_finder->readAmiraTransformations();
  //            prox_finder->readAmiraTransformations(argv[1]);
//                 prox_finder->cutOffGraph(merged_graph);
                
 
               //prox_finder->findProximities1(realBild);
//                prox_finder->findProximities(false);
                
                
          /*if(realBild){  
         //prox_finder->writeProximityImages(outputFilename);
         prox_finder->writeProximityLandmarks(outputFilename);
         
         std::string outputFilename12(outputFilename);
         outputFilename12 += "cell_1_cell_2";
         std::string outputFilename13(outputFilename);
         outputFilename13 += "cell_1_cell_3";
         std::string outputFilename23(outputFilename);
         outputFilename23 += "cell_2_cell_3";
         
         prox_finder->writeMasterHxFile(outputFilename12.c_str(), amiraMergedGraphName,12);
         prox_finder->writeMasterHxFile(outputFilename13.c_str(), amiraMergedGraphName, 13);
         prox_finder->writeMasterHxFile(outputFilename23.c_str(), amiraMergedGraphName, 23);
         
         //prox_finder->writeMasterHxFile(outputFilename, amiraMergedGraphName);
                  
         const char * dummyName = "bla";
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_1) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_1))
         {
            std::string cell_1_name(outputFilename);
            cell_1_name += "cell_1";
            AmiraSpatialGraph *  cell_1_Graph = new AmiraSpatialGraph;
            Reader*  cell_1_GraphWriter = new Reader(dummyName, cell_1_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_1);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_1);
            cell_1_GraphWriter->setSpatialGraph(cell_1_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_1_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_2) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_2))
         {
            std::string cell_2_name(outputFilename);
            cell_2_name += "cell_2";         
            AmiraSpatialGraph *  cell_2_Graph = new AmiraSpatialGraph;
            Reader*  cell_2_GraphWriter = new Reader(dummyName, cell_2_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_2);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_2);
            cell_2_GraphWriter->setSpatialGraph(cell_2_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_2_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_3) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_3))
         {
            std::string cell_3_name(outputFilename);
            cell_3_name += "cell_3";         
            AmiraSpatialGraph *  cell_3_Graph = new AmiraSpatialGraph;
            Reader*  cell_3_GraphWriter = new Reader(dummyName, cell_3_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_3);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_3);
            cell_3_GraphWriter->setSpatialGraph(cell_3_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_3_GraphWriter->writeSpatialGraphFile2();
         }

         
          }*/
          
          

                
    }
    
    if(argc == 12)    
    {   
        
        //for whole cell files with scaling errors  
        
      //argv[1] : amira transformation
      //argv[2] : deconvolved image files
      //argv[3] : output_name
      //argv[4] : image index start_index
      //argv[5] : image index end
      //argv[6] : include unknown: y/n?
      //argv[7] : maximum distance
      //argv[8] : merged graph file
      //argv[9] : additional transformations from text file_name
//       //argv[10]: nondeconvolved images
        
        const char * amiraTransformationName = argv[1];
        const char * deconvolvedImageFileName = argv[2];
        const char * outputFilename = argv[3];
        int imagestartindex = atoi(argv[4]);
        int imageendindex = atoi(argv[5]);
        int includeunknown = atoi(argv[6]);
        int maxDist = atoi(argv[7]);
        const char * amiraMergedGraphName = argv[8];
        const char * additionalScalingFileName = argv[9];
//         const char * nondeconvolvedImageFileName = argv[10];
        bool realBild = atoi(argv[10]);
        //Read in amira file and decon image
        //////std::cout << "Begin Reading" << std::endl;
        
        
        BrickImageReader *BrickImage;
//         BrickImageReader *BrickImageNonDecon;
        if(realBild){
        
        BrickImage = new BrickImageReader(deconvolvedImageFileName, imagestartindex, imageendindex );               
        input_image = BrickImage->GetImage();
        input_image->Update();
        
//         BrickImageNonDecon = new BrickImageReader(nondeconvolvedImageFileName , imagestartindex, imageendindex);               
//         nondecon_image = BrickImageNonDecon->GetImage();
//         nondecon_image->Update();
        }
        else
        {
            ImageType::SizeType fakeSize;
            ImageType::RegionType fakeRegion;
            ImageType::IndexType fakeIndex;
            fakeSize[0] = 100;
            fakeSize[1] = 100;
            fakeSize[2] = 2;
            fakeIndex[0] = 0;  
            fakeIndex[1] = 0;   
            fakeIndex[2] = 0;  
            fakeRegion.SetIndex(fakeIndex);
            fakeRegion.SetSize(fakeSize);
            input_image = ImageType::New();
            input_image->SetRegions(fakeRegion);
            input_image->Allocate();
            input_image->FillBuffer(0);
//                     nondecon_image = ImageType::New();
//                     nondecon_image->SetRegions(fakeRegion);
//                     nondecon_image->Allocate();
//                     nondecon_image->FillBuffer(0);
        }
        

        
        amiraReader = new Reader(amiraTransformationName);  
        amiraReader->readSpatialGraphFile(false);
        input_graph = amiraReader->getSpatialGraph();
        
        mergedAmiraReader = new Reader(amiraMergedGraphName);  
        mergedAmiraReader->readSpatialGraphFile(false);
        merged_graph = mergedAmiraReader->getSpatialGraph();
        

                //image_projection = reader->GetZProjection(1);
                ////////std::cout <<std::endl << "brick image successfully read!" <<std::endl << std::endl;
                
        prox_finder = new AxonDendriteProximityFinder(mergedAmiraReader, input_graph, merged_graph, input_image, amiraTransformationName, includeunknown, maxDist, additionalScalingFileName,outputFilename);
         
//         prox_finder->readAmiraSectionTransformations(inputfilename);
//          std::string folder_name = argv[3];
//          
//                     if(mkdir(folder_name.c_str(),
//                         S_IRWXU| //full access for owner
//                         S_IRWXG| //full access for grp
//                         S_IRWXO  //full access for other
//                         ))
//                     {
//             //         if(errno == EEXIST){
//             //             std::cerr << folder_name << " already exists" << std::endl;
//             //         }
//             //         else {
//             //             perror("Error creating directory");
//             //             exit(-1);
//             //         }
//                     }
//                     
//                     std::string filenames = folder_name;
//                     filenames += "/";
//                     filenames += folder_name;
//                     
//                     char * filename = new char[filenames.size() + 1];
//                     std::copy(filenames.begin(), filenames.end(), filename);
//                     filename[filenames.size()] = '\0';
                    
//          prox_finder->readAmiraTransformations();
//          prox_finder->cutOffGraph(merged_graph);
//          
        

         prox_finder->findProximities1(realBild);
         
//         prox_finder->findProximities(false);
         /*
          if(realBild){  
         //prox_finder->writeProximityImages(outputFilename);
         //prox_finder->writeProximityLandmarks(outputFilename);
         std::string outputFilename12(outputFilename);
         outputFilename12 += "cell_1_cell_2";
         std::string outputFilename13(outputFilename);
         outputFilename13 += "cell_1_cell_3";
         std::string outputFilename23(outputFilename);
         outputFilename23 += "cell_2_cell_3";
         
         prox_finder->writeMasterHxFile(outputFilename12.c_str(), amiraMergedGraphName, 12);
         prox_finder->writeMasterHxFile(outputFilename13.c_str(), amiraMergedGraphName, 13);
         prox_finder->writeMasterHxFile(outputFilename23.c_str(), amiraMergedGraphName, 23);         
         
         const char * dummyName = "bla";
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_1) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_1))
         {
            std::string cell_1_name(outputFilename);
            cell_1_name += "cell_1";
            AmiraSpatialGraph *  cell_1_Graph = new AmiraSpatialGraph;
            Reader*  cell_1_GraphWriter = new Reader(dummyName, cell_1_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_1);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_1);
            cell_1_GraphWriter->setSpatialGraph(cell_1_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_1_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_2) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_2))
         {
            std::string cell_2_name(outputFilename);
            cell_2_name += "cell_2";         
            AmiraSpatialGraph *  cell_2_Graph = new AmiraSpatialGraph;
            Reader*  cell_2_GraphWriter = new Reader(dummyName, cell_2_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_2);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_2);
            cell_2_GraphWriter->setSpatialGraph(cell_2_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_2_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(merged_graph->isLabelInSpatialGraph(AXON_ID_3) || merged_graph->isLabelInSpatialGraph(DENDRITE_ID_3))
         {
            std::string cell_3_name(outputFilename);
            cell_3_name += "cell_3";         
            AmiraSpatialGraph *  cell_3_Graph = new AmiraSpatialGraph;
            Reader*  cell_3_GraphWriter = new Reader(dummyName, cell_3_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(merged_graph, AXON_ID_3);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(merged_graph, DENDRITE_ID_3);
            cell_3_GraphWriter->setSpatialGraph(cell_3_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_3_GraphWriter->writeSpatialGraphFile2();
         }
          }*/
          

        if(realBild) delete BrickImage;
         
         
    }
    
//#if 0
    
    if(argc == 10){
        
      //for single section files  
        
      //argv[1] : amira spatial graph
      //argv[2] : deconvolved image files
      //argv[3] : output_name
      //argv[4] : image index start_index
      //argv[5] : image index end
      //argv[6] : include unknown: y/n?
      //argv[7] : point range
      //argv[8] : maximum distance      
      //argv[9] : nondeconvolved images
        
        const char * amiraGraphName = argv[1];
        const char * deconvolvedImageFileName = argv[2];
        const char * outputFilename = argv[3];
        int imagestartindex = atoi(argv[4]);
        int imageendindex = atoi(argv[5]);
        int includeunknown = atoi(argv[6]);
        int pointRange = atoi(argv[7]);
        int maxDist = atoi(argv[8]);    
//         const char * nondeconvolvedImageFileName = argv[9];
        bool realBild = atoi(argv[9]);
        
        //Read in amira file and decon image
        //////std::cout << "Begin Reading" << std::endl;
        
        amiraReader = new Reader(amiraGraphName);  
        amiraReader->readSpatialGraphFile(false);
        input_graph = amiraReader->getSpatialGraph();
        
      
        BrickImageReader *BrickImage;
//         BrickImageReader *BrickImageNonDecon;
        
        if(realBild){
        BrickImage = new BrickImageReader(deconvolvedImageFileName, imagestartindex, imageendindex);               
        input_image = BrickImage->GetImage();
        input_image->Update();
        
//         BrickImageNonDecon = new BrickImageReader(nondeconvolvedImageFileName , imagestartindex, imageendindex);               
//         nondecon_image = BrickImageNonDecon->GetImage();
//         nondecon_image->Update();
        }
        else
                {
                    ImageType::SizeType fakeSize;
                    ImageType::RegionType fakeRegion;
                    ImageType::IndexType fakeIndex;
                    fakeSize[0] = 100;
                    fakeSize[1] = 100;
                    fakeSize[2] = 2;
                    fakeIndex[0] = 0;  
                    fakeIndex[1] = 0;   
                    fakeIndex[2] = 0;  
                    fakeRegion.SetIndex(fakeIndex);
                    fakeRegion.SetSize(fakeSize);
                    input_image = ImageType::New();
                    input_image->SetRegions(fakeRegion);
                    input_image->Allocate();
                    input_image->FillBuffer(0);
//                    nondecon_image = ImageType::New();
//                     nondecon_image->SetRegions(fakeRegion);
//                     nondecon_image->Allocate();
//                     nondecon_image->FillBuffer(0);
                }
        

                //image_projection = reader->GetZProjection(1);
                ////////std::cout <<std::endl << "brick image successfully read!" <<std::endl << std::endl;
                
         prox_finder = new AxonDendriteProximityFinder(amiraReader, input_graph, input_image, amiraGraphName, includeunknown, pointRange, maxDist,outputFilename);
         prox_finder->findProximities1(realBild);
//                prox_finder->findProximities(false);
                
                
          /*if(realBild){  
         //prox_finder->writeProximityImages(outputFilename);
         prox_finder->writeProximityLandmarks(outputFilename);
         std::string outputFilename12(outputFilename);
         outputFilename12 += "cell_1_cell_2";
         std::string outputFilename13(outputFilename);
         outputFilename13 += "cell_1_cell_3";
         std::string outputFilename23(outputFilename);
         outputFilename23 += "cell_2_cell_3";
         
         prox_finder->writeMasterHxFile(outputFilename12.c_str(), amiraGraphName,12);
         prox_finder->writeMasterHxFile(outputFilename13.c_str(), amiraGraphName, 13);
         prox_finder->writeMasterHxFile(outputFilename23.c_str(), amiraGraphName, 23);
         
         const char * dummyName = "bla";
         if(input_graph->isLabelInSpatialGraph(AXON_ID_1) || input_graph->isLabelInSpatialGraph(DENDRITE_ID_1))
         {
            std::string cell_1_name(outputFilename);
            cell_1_name += "cell_1";
            AmiraSpatialGraph *  cell_1_Graph = new AmiraSpatialGraph;
            Reader*  cell_1_GraphWriter = new Reader(dummyName, cell_1_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(input_graph, AXON_ID_1);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_1"<<std::endl);
            cell_1_Graph->copyEdgesWithLabel(input_graph, DENDRITE_ID_1);
            cell_1_GraphWriter->setSpatialGraph(cell_1_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_1_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(input_graph->isLabelInSpatialGraph(AXON_ID_2) || input_graph->isLabelInSpatialGraph(DENDRITE_ID_2))
         {
            std::string cell_2_name(outputFilename);
            cell_2_name += "cell_2";         
            AmiraSpatialGraph *  cell_2_Graph = new AmiraSpatialGraph;
            Reader*  cell_2_GraphWriter = new Reader(dummyName, cell_2_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(input_graph, AXON_ID_2);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_2"<<std::endl);
            cell_2_Graph->copyEdgesWithLabel(input_graph, DENDRITE_ID_2);
            cell_2_GraphWriter->setSpatialGraph(cell_2_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_2_GraphWriter->writeSpatialGraphFile2();
         }
         
         if(input_graph->isLabelInSpatialGraph(AXON_ID_3) || input_graph->isLabelInSpatialGraph(DENDRITE_ID_3))
         {
            std::string cell_3_name(outputFilename);
            cell_3_name += "cell_3";         
            AmiraSpatialGraph *  cell_3_Graph = new AmiraSpatialGraph;
            Reader*  cell_3_GraphWriter = new Reader(dummyName, cell_3_name.c_str());
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : AXON_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(input_graph, AXON_ID_3);
//             std::flush(//////std::cout<<"in copyEdgeswithLabels : DENDRITE_ID_3"<<std::endl);
            cell_3_Graph->copyEdgesWithLabel(input_graph, DENDRITE_ID_3);
            cell_3_GraphWriter->setSpatialGraph(cell_3_Graph);
            
            //std::flush(//////std::cout<<"in writeSpatialGraph"<<std::endl);
            cell_3_GraphWriter->writeSpatialGraphFile2();
         }
         
//          std::string cell_2_name(outputFilename);
//          cell_2_name += "cell_2";            
//          Reader*  cell_2_Graph = new Reader(dummyName, cell_2_name.c_str());
//          cell_2_Graph->setSpatialGraph(merged_graph);
//          cell_2_Graph->writeSpatialGraphFile(AXON_ID_2, DENDRITE_ID_2);
//          
//          std::string cell_3_name(outputFilename);
//          cell_3_name += "cell_3";           
//          Reader*  cell_3_Graph = new Reader(dummyName, cell_3_name.c_str());
//          cell_3_Graph->setSpatialGraph(merged_graph);
//          cell_3_Graph->writeSpatialGraphFile(AXON_ID_3, DENDRITE_ID_3);
         
          }*/
          
          

        if(realBild) delete BrickImage;
        
        
        
        
        
    }
    
    
    
    if(argc == 7){
        
        //for graph checking
        
        
        
      //argv[1] : amira transformation
      //argv[2] : deconvolved image files
      //argv[3] : output_name
      //argv[4] : image index start_index
      //argv[5] : image index end
      //argv[6] : nondeconvolved images    
        
        const char * amiraTransformationName = argv[1];
        const char * deconvolvedImageFileName = argv[2];
        const char * outputFilename = argv[3];
        int imagestartindex = atoi(argv[4]);
        int imageendindex = atoi(argv[5]);
//        const char * nondeconvolvedImageFileName = argv[6];
        bool realBild = atoi(argv[6]);
        
                //Read in amira file and decon image
                //////std::cout << "Begin Reading" << std::endl;        
        
                
                BrickImageReader *BrickImage;
//                 BrickImageReader *BrickImageNonDecon;
                
                if(realBild){
                BrickImage = new BrickImageReader(deconvolvedImageFileName, imagestartindex, imageendindex);               
                input_image = BrickImage->GetImage();
                input_image->Update();
                
//                 BrickImageNonDecon = new BrickImageReader(nondeconvolvedImageFileName , imagestartindex, imageendindex);               
//                 nondecon_image = BrickImageNonDecon->GetImage();
//                 nondecon_image->Update();
                }
                else
                {
                    ImageType::SizeType fakeSize;
                    ImageType::RegionType fakeRegion;
                    ImageType::IndexType fakeIndex;
                    fakeSize[0] = 100;
                    fakeSize[1] = 100;
                    fakeSize[2] = 2;
                    fakeIndex[0] = 0;  
                    fakeIndex[1] = 0;   
                    fakeIndex[2] = 0;  
                    fakeRegion.SetIndex(fakeIndex);
                    fakeRegion.SetSize(fakeSize);
                    input_image = ImageType::New();
                    input_image->SetRegions(fakeRegion);
                    input_image->Allocate();
                    input_image->FillBuffer(0);
//                     nondecon_image = ImageType::New();
//                     nondecon_image->SetRegions(fakeRegion);
//                     nondecon_image->Allocate();
//                     nondecon_image->FillBuffer(0);
                }
                

                amiraReader = new Reader(amiraTransformationName);  
                amiraReader->readSpatialGraphFile(false);
                input_graph = amiraReader->getSpatialGraph();
                
                //image_projection = reader->GetZProjection(1);
                                
                prox_finder = new AxonDendriteProximityFinder(mergedAmiraReader, input_graph, input_image, amiraTransformationName,outputFilename /*, nondecon_image*/);
                
                               prox_finder->findProximities1(realBild);
//                prox_finder->findProximities(false);
                
                
          /*if(realBild){  
         //prox_finder->writeProximityImages(outputFilename);
         prox_finder->writeProximityLandmarks(outputFilename);
         
          }*/
          

                if(realBild) delete BrickImage;
        
        
        
    }
    
    //else
    //{
    //     prox_finder = new ProximityFinder(input_graph, input_image);
    //}
    
//#endif

//writeStatsXls();
  
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


        
    







