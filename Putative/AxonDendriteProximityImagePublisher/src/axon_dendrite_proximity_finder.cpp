/****************************************************************************/
/*                                                                          */
/* File:      proximity_finder.cpp                                          */
/*                                                                          */
/* Purpose:                                                                 */
/*                                                                          */
/*                                                                          */
/* Authors:   Christopher Tull                                              */
/*            Alison Smyth                                                  */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   26.03.2014                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



#include "axon_dendrite_proximity_finder.h"
#include "typedefs.h"
#include "basics.h"
#include "utility.h"

using namespace std;

/********************************************************************************************************************
Method: findProximities()

Performs the bulk of the work for the ProximityFinder class.
Cycles through every axon-dendrite pairing and creates local regions called Proximities where the axon 
and dendrite are close to each other. This 'closeness' is determined by the global constant MAX_DIST
*********************************************************************************************************************/

AxonDendriteProximityFinder::AxonDendriteProximityFinder(const char* landmarkfilename, const char* xlsinputfilename, const char* xlsoutputfilename, int section_number, const char*deconvolvedImageFileName, int imagestartindex, int imageendindex,  const char * sectiongraphname, const char * mergedgraphname, int include_unknown,const char *outputFilename/*, ImageType::Pointer nondecon_image*/)
     //5 arguments----for whole merged cell files without scaling errors
     {
       
#ifdef MEGE_OVERLAPPING_PROX_IMAGES
       double * merge_box_coords = new double[6];
       vector<double *> merged_bb_list;
#endif
       
       
         ////std::cout<< "sec no "<< section_number << std::endl;
         //graphReaderPointer = ameraReaderPointer;
         
//          original_image_nondecon = nondecon_image;
//          original_image_nondecon->Update();
         std::cout<<"In correct constructor "<<std::endl;
	 
         this->inputfilename = sectiongraphname;
         this->eliminateDoubledProximites = true;
         this->include_unknown = int(include_unknown);
         this->numberOfPointsAway = 24;        
         this->max_dist = MAX_DIST;
         this->all_bouton_list = new std::vector<double *>;
         this->outputpath = outputFilename;
         
         
        //initialize transformation matrices and set to default
        this->sectionboundingbox = new double [6];
        this->sectionTranslation = new double *[4];
        this->sectionRotation = new double *[4];

        for(int ii = 0; ii < 4; ++ii)
        {
            this->sectionTranslation[ii] = new double[4];
            this->sectionRotation[ii] = new double[4];
            for(int jj = 0; jj < 4; ++jj)
            {
                    this->sectionTranslation[ii][jj] = 0;
                    this->sectionRotation[ii][jj] = 0;
            }
            this->sectionTranslation[ii][ii] = 1;
            this->sectionRotation[ii][ii] = 1;
        }
    
        this->transformation = new double *[4];
        for(int ii = 0; ii < 4; ++ii)
        {
                this->transformation[ii] = new double[4];
                for(int jj = 0; jj < 4; ++jj)
                {
                        if(ii != jj)
                                this->transformation[ii][jj] = 0;
                        else
                                this->transformation[ii][jj] = 1;
                }
        }
            
        this->inverse_transformation = new double *[4];
        for(int ii = 0; ii < 4; ++ii)
        {
                this->inverse_transformation[ii] = new double[4];
                for(int jj = 0; jj < 4; ++jj)
                {
                        if(ii != jj)
                                this->inverse_transformation[ii][jj] = 0;
                        else
                                this->inverse_transformation[ii][jj] = 1;
                }
        }   
        
        cout<< "done initializing"<<endl;
	
	cout<< "section graph name"<<sectiongraphname <<endl;
 
        Reader * sectionGraphReader = new Reader(sectiongraphname , "/home/mythreya/temp/sectionGraphAsis");  
        sectionGraphReader->readSpatialGraphFile(false);
        untransformed_section_graph = sectionGraphReader->getSpatialGraph();
	
	cout<< "read section graph"<<endl;
        
        Reader * mergedAmiraReader = new Reader(mergedgraphname /*, "/home/mythreya/Documents/outputprox/bricktest/DS_20140122/mergedGraph"*/);  
        mergedAmiraReader->readSpatialGraphFile(false);
        amira_graph = mergedAmiraReader->getSpatialGraph();
	
	cout<< "read merged graph"<<endl;

        Reader * transformedsectionGraphReader = new Reader(sectiongraphname , "/home/mythreya/temp/sectionGraphtransformed");  
        transformedsectionGraphReader->readSpatialGraphFile(true);
        transformed_section_graph = transformedsectionGraphReader->getSpatialGraph();
	
	cout<< "read transformed graph"<<endl;
        
        transformedsectionGraphReader->writeSpatialGraphFile2();
	sectionGraphReader->writeSpatialGraphFile2();
	
	cout<< "wrote transformed graph"<<endl;
        
         // section bounds can be only read after transforming the section
         transformed_section_graph->getBoundingBox(sectionboundingbox);
         
         //std::cout<<" trans: "<<sectionboundingbox[0]<<"  Y :  "<<sectionboundingbox[2]<<"  Z : "<<sectionboundingbox[4]<<std::endl;
         //std::cout<<" trans: "<<sectionboundingbox[1]<<"  Y :  "<<sectionboundingbox[3]<<"  Z : "<<sectionboundingbox[5]<<std::endl;
         
         // read the prox and find out which ones belong to this section landmarks
         readLandmarkfile(landmarkfilename, &ProximityLandmarks);
         for(int i = 0; i < ProximityLandmarks.size(); i++)
         {
             if(isInSectionBoundingBox(ProximityLandmarks.at(i)))
             {
                // modifiy the xls sheet with the respective section numbers
                //std::cout<<" prox num = "<<i+1<<" "<<ProximityLandmarks.at(i)[0]<<" "<<ProximityLandmarks.at(i)[1]<<" "<<ProximityLandmarks.at(i)[2]<<std::endl;
                
                //this->sectionLineNumbers.push_back(i+1);
                
             }
            
         }
         
	 cout<< "added section number to the xls"<<endl;
         
         //read transformation info from transformed graph
         readAmiraTransformations();
         
         Utility::applyTransformationToPoint(&sectionboundingbox[0], inverse_transformation);
         Utility::applyTransformationToPoint(&sectionboundingbox[3], inverse_transformation);
	 
	 
	//getMacthingGraphPoints(transformed_section_graph, amira_graph);
	 
	 
         //std::cout<<"inv X trans: "<<sectionboundingbox[0]<<"  Y :  "<<sectionboundingbox[2]<<"  Z : "<<sectionboundingbox[4]<<std::endl;
         //std::cout<<"inv X trans: "<<sectionboundingbox[1]<<"  Y :  "<<sectionboundingbox[3]<<"  Z : "<<sectionboundingbox[5]<<std::endl;
         //inverse transform merged graph so that the section to be processed is in its original place
         amira_graph->applyInverseTransformation(inverse_transformation);
         
         // Now use the untransformed sections bounding box that corresponds to the image
         // the section graph as well as the image are at the global 0,0,0 coordinates
         untransformed_section_graph->getBoundingBox(sectionboundingbox);
         //std::cout<<" untrans X : "<<sectionboundingbox[0]<<"  Y :  "<<sectionboundingbox[2]<<"  Z : "<<sectionboundingbox[4]<<std::endl;
         //std::cout<<" untrans X : "<<sectionboundingbox[1]<<"  Y :  "<<sectionboundingbox[3]<<"  Z : "<<sectionboundingbox[5]<<std::endl;
         
         bool need_to_read_images = false;
         
        for(int i = 0; i < ProximityLandmarks.size(); i++)
        {
            //std::cout<<"before "<<ProximityLandmarks.at(i)[0]<<" "<<ProximityLandmarks.at(i)[1]<<" "<<ProximityLandmarks.at(i)[2]<<std::endl;
            Utility::applyTransformationToPoint(ProximityLandmarks.at(i), inverse_transformation);
            //std::cout<<" affter "<<ProximityLandmarks.at(i)[0]<<" "<<ProximityLandmarks.at(i)[1]<<" "<<ProximityLandmarks.at(i)[2]<<std::endl;
            ////std::cout<<" "<<ProximityLandmarks.at(i)[0]<<" "<<ProximityLandmarks.at(i)[1]<<" "<<ProximityLandmarks.at(i)[2]<<std::endl;
            if(isInSectionBoundingBox(ProximityLandmarks.at(i)))
            {
                //std::cout<<" prox num = "<<i+1<<" "<<ProximityLandmarks.at(i)[0]<<" "<<ProximityLandmarks.at(i)[1]<<" "<<ProximityLandmarks.at(i)[2]<<std::endl;
                ////std::cout<<"found points inside"<<std::endl;
                need_to_read_images = true;
                // line number would be i +2 including the header line
                // column number is 
                this->sectionLineNumbers.push_back(i+1);
                this->inSectionLandmarks.push_back(ProximityLandmarks.at(i));
		
#ifdef MEGE_OVERLAPPING_PROX_IMAGES
		if()
		{
		  // 
		   
		}
		  
#endif
                //break;
                
            }
            
        }
         
	 if(this->sectionLineNumbers.size() > 0)
	 {
	     cout<< "appending to xls" << endl;
	     Utility::appendToXls(xlsinputfilename, xlsoutputfilename, &(this->sectionLineNumbers), section_number);
	 
	 }
	 
	 cout<< "got the section number"<<endl;
        
         cout<< "found out whether to read images "<< need_to_read_images <<endl;
	
	
        
         
         if(need_to_read_images)
         {
            BrickImageReader *BrickImage;
            BrickImage = new BrickImageReader(deconvolvedImageFileName, imagestartindex, imageendindex, &(this->inSectionLandmarks) );               
            original_image = BrickImage->GetImage();
            original_image->Update();
	    
#ifdef NONDECON
	    
	    BrickImageReader *BrickImage1;
            BrickImage1 = new BrickImageReader("nondecon%03d.png", imagestartindex, imageendindex, &(this->inSectionLandmarks) );               
            original_image_non_decon = BrickImage1->GetImage();
            original_image_non_decon->Update();
	    
#endif
           
            cout<< "read the required images "<< need_to_read_images <<endl;
         }
         
         
        for(int i = 0; i < ProximityLandmarks.size(); i++)
        {
            if(isInSectionBoundingBox(ProximityLandmarks.at(i)))
            {
                // write the images out
                std::string format(this->outputpath);
		std::string filename("Proximity_Zone_");
                format += "/Proximity_Zone_";
                format += Utility::inttochar(i+1);
		filename += Utility::inttochar(i+1);
		cout<< "Proximity_Zone: "<<format<<endl;
                Proximity * tmpprox = new Proximity(ProximityLandmarks.at(i),this->original_image);
                ////std::cout<<"before writing prox images"<<std::endl;
                tmpprox->writeProximityImage(format.c_str(), filename.c_str(), this->transformation, i);
		
#ifdef NONDECON
		std::string format1(this->outputpath);
                format1 += "/Proximity_Zone_";
		
		std::string filename1("Proximity_Zone_");
		filename1 += Utility::inttochar(i+1);
		filename1 += "_non_decon";
		
                format1 += Utility::inttochar(i+1);
		format1 += "_non_decon";
		cout<< "Proximity_Zone: "<<format1<<endl;
		Proximity * tmpprox1 = new Proximity(ProximityLandmarks.at(i),this->original_image_non_decon);
                ////std::cout<<"before writing prox images"<<std::endl;
                tmpprox1->writeProximityImage(format1.c_str(), filename1.c_str(), this->transformation, i);
		
#endif
		
#ifdef MEGE_OVERLAPPING_PROX_IMAGES
		
		
		
#endif
                
                
            }
            
        }
         
         //set_edge_lists(amira_graph, true);         
         
};

// this function matches the part of section graph thats in line with the final merged graph
void AxonDendriteProximityFinder::getMacthingGraphPoints(AmiraSpatialGraph * sectiongraph_to_be_matched, AmiraSpatialGraph * main_graph_to_be_matched)
{
  vector<double * > matching_point_list;
  vector<Edge *> *sectiongraph_edges = transformed_section_graph->edgesPointer();
  vector<Edge *> *maingraph_edges = amira_graph->edgesPointer(); 
  for(int sectiongraphEdgeCount = 0; sectiongraphEdgeCount < sectiongraph_edges->size(); sectiongraphEdgeCount++)
  {
      Edge * current_section_edge = sectiongraph_edges->at(sectiongraphEdgeCount);
      
      std::list<double*>::iterator current_section_point_it;// = current_section_edge->edgePointCoordinates;
      for(current_section_point_it = current_section_edge->edgePointCoordinates.begin(); current_section_point_it != current_section_edge->edgePointCoordinates.end(); current_section_point_it++)
      {
	  for(int maingraphEdgeCount = 0; maingraphEdgeCount < maingraph_edges->size(); maingraphEdgeCount++)
	  {
	      Edge * current_main_edge = maingraph_edges->at(maingraphEdgeCount);
	      
	      std::list<double*>::iterator current_main_point_it;// = current_main_edge->edgePointCoordinates;
	      for(current_main_point_it = current_main_edge->edgePointCoordinates.begin(); current_main_point_it != current_main_edge->edgePointCoordinates.end(); current_main_point_it++)
	      {
		  //Edge * current_section_edge = sectiongraph_edges->at(0);
      
		  //std::list<double*>::iterator current_section_point_it = current_section_edge->edgePointCoordinates.begin();
		  cout<< "current_maingraph_points: "<< (*current_main_point_it)[0] << " "<<(*current_main_point_it)[1] <<" "<< (*current_main_point_it)[2] << endl;
		  //if (Utility::Match3DPoints(*current_section_point_it, *current_main_point_it, 1))
		  {
		      // found a matching point,
		      double * temp = new double[3];
		      temp[0] = (*current_main_point_it)[0];
		      temp[1] = (*current_main_point_it)[1];
		      temp[2] = (*current_main_point_it)[2];
		      matching_point_list.push_back(temp);
		      cout<< "matching section point "<< (temp)[0] << " "<<(temp)[1] << " "<< (temp)[2] << " "<< endl;
		      
		  }
	      }
	  }
      }
  }
#if 0
  
  for(int sectiongraphEdgeCount = 0; sectiongraphEdgeCount < sectiongraph_edges->size(); sectiongraphEdgeCount++)
  {
      Edge * current_section_edge = sectiongraph_edges->at(sectiongraphEdgeCount);
      
      std::list<double*>::iterator current_section_point_it;// = current_section_edge->edgePointCoordinates;
      for(current_section_point_it = current_section_edge->edgePointCoordinates.begin(); current_section_point_it != current_section_edge->edgePointCoordinates.end(); current_section_point_it++)
      {
	  cout<< "current_section_points: "<< (*current_section_point_it)[0] << " "<<(*current_section_point_it)[1] << " "<<(*current_section_point_it)[2] << endl;
	  // for each point of the section graph see which points match the final merged main graph
	  
	  for(int maingraphEdgeCount = 0; maingraphEdgeCount < maingraph_edges->size(); maingraphEdgeCount++)
	  {
	      Edge * current_main_edge = maingraph_edges->at(maingraphEdgeCount);
	      
	      std::list<double*>::iterator current_main_point_it;// = current_main_edge->edgePointCoordinates;
	      for(current_main_point_it = current_main_edge->edgePointCoordinates.begin(); current_main_point_it != current_main_edge->edgePointCoordinates.end(); current_main_point_it++)
	      {
		  cout<< "current_maingraph_points: "<< (*current_main_point_it)[0] << " "<<(*current_main_point_it)[1] <<" "<< (*current_main_point_it)[2] << endl;
		  if (Utility::Match3DPoints(*current_section_point_it, *current_main_point_it, 1))
		  {
		      // found a matching point,
		      double * temp = new double[3];
		      temp[0] = (*current_main_point_it)[0];
		      temp[1] = (*current_main_point_it)[1];
		      temp[2] = (*current_main_point_it)[2];
		      matching_point_list.push_back(temp);
		      cout<< "matching section point "<< (temp)[0] << " "<<(temp)[1] << " "<< (temp)[2] << " "<< endl;
		      
		  }
	      }
	  }
      }
  
  }
  
   writeLandmarks("/home/mythreya/tmp/mathingpoints", &matching_point_list);
   
#endif

};


     
void AxonDendriteProximityFinder::readLandmarkfile(const char* filename, std::vector<double *>*prox_landmarks)
{
    
  ////////////std:://cout<<"reading"<<std::endl;
    ////////////std:://cout<<filepath<<std::endl;
//   
//   std::string format = filepath;
//   format += "/Landmarks.landmarkAscii";
//   ////////////std:://cout<<format<<std::endl;
  std::ifstream inputStream(filename);
  std::string currentLine = "";
  int count = 0;
  
  //////////////std:://cout<<format<<std::endl;
  if(!inputStream.fail())
  {
    ////std::cout<<"valid path"<<std::endl;
    while(!std::getline(inputStream, currentLine).eof() && currentLine.find("@1") != 0 ){}
    while(!std::getline(inputStream, currentLine).eof() && currentLine.size() > 0 )
    { 
      double * tmp = new double[ARRAY_LENGTH];
      
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[X_COORD], &tmp[Y_COORD], &tmp[Z_COORD]);
      ////std::cout<<" "<<tmp[X_COORD]<<" "<<tmp[Y_COORD]<<" "<<tmp[Z_COORD]<<std::endl;
      //Proximity * tmpprox = new Proximity(tmp,this->original_image);
      prox_landmarks->push_back(tmp);
      
    }
    ////std::cout<<"done reading prox landmarks"<<std::endl;
  }

};

     

     
     
void AxonDendriteProximityFinder::sortProxDist( std::vector<Proximity *> *prox_list_unsorted, std::vector<Proximity *> *prox_list_sorted )
{
    std::vector<Proximity*>::iterator it;
    
    while(!prox_list_unsorted->empty())
    {
        double min_dist = 0xffff;
        int min_ind = 0xffff;
        //for (it = prox_list_unsorted->begin(); it != prox_list_unsorted->end(); it++)
        for(int i = 0; i < prox_list_unsorted->size(); i++)
        {
            // find the least and push it to sorted list
            if(prox_list_unsorted->at(i)->getDistance() < min_dist)
            {
                min_dist = prox_list_unsorted->at(i)->getDistance();
                min_ind = i;
                
            }
        }
        
        prox_list_sorted->push_back( new Proximity(prox_list_unsorted->at(min_ind), this->original_image));
        prox_list_unsorted->erase(prox_list_unsorted->begin()+min_ind);
        
    }
    
    
};
     
void AxonDendriteProximityFinder::findValleys(std::vector<Proximity*>*inputProxVector, std::vector<Proximity*>*outputProxVector, double threshold, double dist_to_merge)
{
    // iterate over the input list and find the min with
     //std::vector<Proximity*> outProxList;
    int num_of_neighbouring_points = 1;//rint (dist_to_merge * 0.1 / 2);
    int max_interdist = 3;
    int start_index = 0;
    int end_index = 0;
    for(int i = 0; i < (inputProxVector->size()-1); i++)
    {
        double dist = Utility::euclidean_distance(inputProxVector->at(i)->getAxonPoint(), inputProxVector->at(i+1)->getAxonPoint(),3,1);
        
        if((dist > max_interdist) || (i == (inputProxVector->size()-2)))
        {
            double local_min = 0xffff;
            int min_indx = 0;
            // let find minima for the points so far
            for(int k = start_index; k <= i ; k++)
            {
                if(inputProxVector->at(k)->getDistance() < local_min)
                {

                    local_min = inputProxVector->at(k)->getDistance();
                    min_indx = k;
                    
                }
                
            }
            
            Proximity * prox = new Proximity(inputProxVector->at(min_indx)->getAxonPoint(), inputProxVector->at(min_indx)->getDendritePoint(), inputProxVector->at(min_indx)->getDistance(), inputProxVector->at(min_indx)->Axon_ID, inputProxVector->at(min_indx)->Dendrite_ID, this->original_image);
            outputProxVector->push_back(prox); 
            
            start_index = i+1;
        }
        
    }
        
        
        
        
       
    //std::vector<Proximity*>::iterator it,it2;
    
    /*for(it = inputProxVector->begin()+num_of_neighbouring_points; it != inputProxVector->end() - num_of_neighbouring_points; it++)
    {
        double local_min = 0xffff;
        for(it2 = (it -num_of_neighbouring_points); it2 < (it + num_of_neighbouring_points); it2++)
        {
            if(it != it2)
            {
                if((*it2)->getDistance() < local_min)
                {
                    local_min = (*it2)->getDistance();
                }
                
            }
        }
        
        if((*it)->getDistance() == local_min)
        {
            // this prox zone is a local minima
            // lets add this to the final prox zone
            outputProxVector->push_back((*it));
        }
        
        
    }*/
    
    /*for(int i = num_of_neighbouring_points; i < (inputProxVector->size() - num_of_neighbouring_points); i++)
    {
        ////std::cout<<"   "<< i << " " << num_of_neighbouring_points<< " "<< inputProxVector->size()<<std::endl;
        double local_min = 0xffff;
        for(int j = (i -num_of_neighbouring_points ); j < (i + num_of_neighbouring_points); j++)
        {
            ////std::cout<< " "<< j << " "<< local_min<< inputProxVector->at(j)->getDistance()<<std::endl;
            if((i != j) && (inputProxVector->at(j)->getDistance() < local_min))
            {
                
                local_min = inputProxVector->at(j)->getDistance();
                
            }
        }
        ////std::cout<< "   " << i << " "<< local_min<< inputProxVector->at(i)->getDistance()<<std::endl;
        if(inputProxVector->at(i)->getDistance() == local_min)
        {
            Proximity * prox = new Proximity(inputProxVector->at(i)->getAxonPoint(), inputProxVector->at(i)->getDendritePoint(), inputProxVector->at(i)->getDistance(), inputProxVector->at(i)->Axon_ID, inputProxVector->at(i)->Dendrite_ID);
            ////std::cout<< "min "<<local_min << std::endl;
            // this prox zone is a local minima
            // lets add this to the final prox zone
            ////std::cout<<inputProxVector->at(i)->getRealCentroid()[X_COORD] << " " << inputProxVector->at(i)->getRealCentroid()[Y_COORD] << " " << inputProxVector->at(i)->getRealCentroid()[Z_COORD] << std::endl;
            outputProxVector->push_back(prox);
        }
        
        
    }*/
    
    
    
};

void AxonDendriteProximityFinder::writeLandmarks( const char* outputFilename, std::vector<double*> *landmark_vector)
{
    ////std::cout<< " In writeLandmarkFile  "<< std::endl;
  //////std::cout<< "in writeLandmarkFile" << std::endl;
  
// transformToWorldCoordinates(list);
  
//   std::string format = outputFilename;
//   format += "/landmark_locations.landmarkAscii";
  std::ofstream LandMarkData( outputFilename );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << landmark_vector->size()                << std::endl;
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
  
  for(int i=0; i < landmark_vector->size(); i++)
  {
    double * point = landmark_vector->at(i);
    LandMarkData << point[X_COORD] << " " << point[Y_COORD] << " " << point[Z_COORD] << std::endl;
     cout << "writing landmarks to file " << point[X_COORD] << " " << point[Y_COORD] << " " << point[Z_COORD] << std::endl;
    
  } 
  
  LandMarkData.close();
};

void AxonDendriteProximityFinder::writeProxLandmarks( const char* outputFilename, std::vector<Proximity*> *sorted_prox_vector)
{
    ////std::cout<< " In writeLandmarkFile  "<< std::endl;
  //////std::cout<< "in writeLandmarkFile" << std::endl;
  
// transformToWorldCoordinates(list);
  
//   std::string format = outputFilename;
//   format += "/landmark_locations.landmarkAscii";
  std::ofstream LandMarkData( outputFilename );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << sorted_prox_vector->size()                << std::endl;
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
  
  for(int i=0; i < sorted_prox_vector->size(); i++)
  {
    LandMarkData << sorted_prox_vector->at(i)->getRealCentroid()[X_COORD] << " " << sorted_prox_vector->at(i)->getRealCentroid()[Y_COORD] << " " << sorted_prox_vector->at(i)->getRealCentroid()[Z_COORD] << std::endl;
  } 
  
  LandMarkData.close();
};

void AxonDendriteProximityFinder::writeProxXls(const char* filename, std::vector<Proximity*> *sorted_prox_vector)
{
   
    bool file_exsists = 0;
    
    std::string s1(filename);
    //std::string s2("PutativeContacts.csv");
    
    ////////////////std:://cout<<"In writexls"<<std::endl;
    
    if(chdir(this->outputpath) != 0)
            perror("Couldn't open output drirectory!");
    
    
    
    std::ifstream f1(s1.c_str());
    if (f1.good()) {
        f1.close();
        file_exsists = true;
    } else {
        f1.close();
        file_exsists = false;
    } 
    
//     std::ifstream f2(s2.c_str());
//     if (f2.good()) {
//         f2.close();
//         file_exsists = true;
//     } else {
//         f2.close();
//         file_exsists = false;
//     } 
    
    
    std::fstream xlsParamWriter;
//     std::fstream csvParamWriter;
    
    
    xlsParamWriter.open (s1.c_str(),  std::fstream::out );
    //csvParamWriter.open ("PutativeContacts.csv", std::fstream::in | std::fstream::out | std::fstream::app);
    
    // write header if it does not exists
    if(!file_exsists)
    {
        xlsParamWriter<<"ID"<<'\t'<<"Proximity Location X"<<'\t'<<"Proximity Location Y"<<'\t'<<"Proximity Location Z"<<'\t'
        //<<"From Axon of Cell #"<<'\t'<<"To Dendrite of Cell #"<<'\t'
        <<"Distance"<<'\t'<<"From Axon of Cell #"<<'\t'<<"To Dendrite of Cell #"<<'\n';
        
//          csvParamWriter<<"Section Number"<<","<<"ROI Number"<<","
//         <<"From Axon of Cell #"<<","<<"To Dendrite of Cell #"<<","
//         <<"Number of Boutons"<<","<<"Number of Putative Contacts"<<'\n';
    }
    
    
    for(int i = 0; i<sorted_prox_vector->size(); i++)
    {
        //////////////////std:://cout<< (*paramIterator)->landmark[0]<<std::endl;
        
        xlsParamWriter<<i<<'\t'<<sorted_prox_vector->at(i)->getRealCentroid()[0]<<'\t'<<sorted_prox_vector->at(i)->getRealCentroid()[1]<<'\t'<<sorted_prox_vector->at(i)->getRealCentroid()[2]<<'\t'
        <<sorted_prox_vector->at(i)->getDistance()<<'\t'
        <<(amira_graph->getAxonCellIndex(sorted_prox_vector->at(i)->Axon_ID))+1<<'\t'<<(amira_graph->getDendriteCellIndex(sorted_prox_vector->at(i)->Dendrite_ID))+1<<'\n';
        
        
       
        
//         std::fstream xlsParamWriter;
//         csvParamWriter<<section_number<<","<<roi_number<<","
//         <<axon_num<<","<<dend_num<<","<<num_boutons<<","<<num_putatives<<'\n';
        
    }
    
    
    xlsParamWriter.close();
    //csvParamWriter.close();
    
    
    
};
     
     
void AxonDendriteProximityFinder::findProximities1(bool notestrun)
{
    #if 0
    ////std::cout<<"In findProximities1()"<<std::endl;
    int num_of_cells = amira_graph->getNumberofCells();
    ////std::cout<<"# of cells"<<num_of_cells<<std::endl;
    std::vector<Proximity *>* ProxMasterArray[num_of_cells][num_of_cells];
    std::vector<Proximity *>* ProxMasterArrayCleaned[num_of_cells][num_of_cells];
    std::vector<Proximity *>* ProxMasterArrayCleanedSorted[num_of_cells][num_of_cells];
    std::vector<Proximity *> ProxLinearArray;
    std::vector<Proximity *> ProxLinearArraySorted;
    std::vector<Proximity *> ProxLinearArrayUnSorted;
     
    std::vector<Proximity *> PrelimProxVector;
    bool isfirstone = true;
    int prox_counter = 0;
    
    // Init the Prox master array with valid vector pointers
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
              ProxMasterArray[i][j] = new std::vector<Proximity *>;
        }
    }
    
    // Init the Prox master array with valid vector pointers
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
              ProxMasterArrayCleaned[i][j] = new std::vector<Proximity *>;
        }
    }
    
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
              ProxMasterArrayCleanedSorted[i][j] = new std::vector<Proximity *>;
        }
    }
    
    
    ////std::cout<<"axon edge count"<<axon_list.size()<<std::endl;
    ////std::cout<<"dend edge count"<<dendrite_list.size()<<std::endl;
    
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    ////std::cout<<"edge count"<<edges->size()<<std::endl;
    
    
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
                        
            if(i != j)
            {
                // For each axon-dendrite pairing
                for(int axonCount = 0; axonCount < edges->size(); axonCount++)    
                {
                    //////std::cout<< edges->at(axonCount)->label<<std::endl;
                    if(isAxon(edges->at(axonCount)->label) && (amira_graph->getAxonCellIndex(edges->at(axonCount)->label) == i))
                    {
                        Edge * currentAxon = edges->at(axonCount);
                        int AxonID = currentAxon->label;
                        std::list<double*>::iterator axon_point_it ;//= currentAxon->edgePointCoordinates;
                        for(axon_point_it = currentAxon->edgePointCoordinates.begin(); axon_point_it != currentAxon->edgePointCoordinates.end(); axon_point_it++)
                        {
                                        
                            double * axonPoint = *axon_point_it;//currentAxon.getNext();
                            //////std::cout<<"axon edge index"<<axonCount<<std::endl;
                            
                            //if(isInSectionBoundingBox(axonPoint))
                            {
                                
                                double min_dist = 0xffff;
                                double *min_dend_point = new double[3];
                                int min_dend_id = 0xff;
                                            
                                for(int dendriteCount = 0; dendriteCount < edges->size(); dendriteCount++)
                                {
                                    if(isDendrite(edges->at(dendriteCount)->label) && (amira_graph->getDendriteCellIndex(edges->at(dendriteCount)->label) == j)) 
                                    {
                                        //////std::cout<<"dend edge index"<<dendriteCount<<std::endl;
                                        //SimpleDendrite currentDendrite = dendrite_list.at(dendriteCount);
                                        Edge * currentDendrite = edges->at(dendriteCount);
                                        
                                        int DendriteID = currentDendrite->label;
                                        
                                        
                                        //if(DendriteID == DENDRITE_ID_1 && AxonID == AXON_ID_2 || DendriteID == DENDRITE_ID_1 && AxonID == AXON_ID_3 || DendriteID == DENDRITE_ID_2 && AxonID == AXON_ID_1 || DendriteID == DENDRITE_ID_2 && AxonID == AXON_ID_3 || DendriteID == DENDRITE_ID_3 && AxonID == AXON_ID_1 || DendriteID == DENDRITE_ID_3 && AxonID == AXON_ID_2)
                                        //if(isTheSameCell(AxonID,DendriteID) == false)
                                        {
                                            
                                            // find the closest dendrite point for this one 
                                            //while(currentDendrite.hasNext())
                                            std::list<double*>::iterator dend_point_it ;//= currentAxon->edgePointCoordinates;
                                            //currentAxon.setToBegin();
                                            //while(currentAxon.hasNext())
                                            for(dend_point_it = currentDendrite->edgePointCoordinates.begin(); dend_point_it != currentDendrite->edgePointCoordinates.end(); dend_point_it++)
                                            {
                                            
                                                double * dendritePoint = *dend_point_it;//currentAxon.getNext();
                                                
                                                //if(isInSectionBoundingBox(dendritePoint))
                                                {
                                                    
                                                    double distance = Utility::euclidean_distance(axonPoint, dendritePoint, 3, 1);
                                                    if((distance < this->max_dist) && (distance < min_dist)) 
                                                    {   
                                                        
                                                        min_dist = distance;
                                                        min_dend_point[0] = dendritePoint[0];
                                                        min_dend_point[1] = dendritePoint[1];
                                                        min_dend_point[2] = dendritePoint[2];
                                                        min_dend_id = currentDendrite->label;
                                                    }
                                                    
                                                
                                                }// if end dendrite point in section
                                            }
                                            
                                            //////std::cout<<"done with current dend"<<std::endl;
                                            
                                            // if there is a close dendrite point lets add it to the prilim prox list
                                            
                                                
                                    
                                        }
                                    
                                    }
                                }
                                
                                if(min_dist < this->max_dist)
                                {
//                                     ////std::cout<<"inside prox"<<std::endl;
//                                     ////std::cout<< "axon " << i << "dend "<< j << std::endl;
//                                     ////std::cout<< axonPoint[0] <<" "<<axonPoint[1] <<" "<< axonPoint[2]<< std::endl;
                                    Proximity * prox = new Proximity(axonPoint, min_dend_point, min_dist, currentAxon->label, min_dend_id, this->original_image);
                                    //PrelimProxVector.push_back();
                                    //ProxMasterArray[amira_graph->getAxonCellIndex(AxonID)][amira_graph->getDendriteCellIndex(min_dend_id)]->push_back(prox);
                                    ProxMasterArray[i][j]->push_back(prox);
                                    ProxLinearArray.push_back(prox);
                                    
                                }
                            }
                            
                            
                        }
                    }
                }
            }
        }
    }
    
    //sortProxDist(&ProxLinearArray, &ProxLinearArraySorted);
    
    //////std::cout<< "out of loop"<<std::endl;
    
    std::string masterhxfilename(this->outputpath);
    
    ////std::cout<<masterhxfilename<<std::endl;
    
    masterhxfilename += "/Masterhxload.hx";
    
    ////std::cout<<masterhxfilename<<std::endl;
    std::ofstream MasterHxData(masterhxfilename.c_str());
    
    MasterHxData << "# Amira Script"<< std::endl;
    MasterHxData << " "<< std::endl;
    
    
    findValleys(&ProxLinearArray, &ProxLinearArrayUnSorted, this->max_dist, 1.0 );
    sortProxDist(&ProxLinearArrayUnSorted, &ProxLinearArraySorted);
    
    std::string s1 = this->outputpath;
    s1+="/Proximities.xls";
    writeProxXls(s1.c_str(),&ProxLinearArraySorted);
    
    std::string s2 = this->outputpath;
    s2+="/Prox_Lanmarks.landmarkAscii";
    writeProxLandmarks(s2.c_str(), &ProxLinearArraySorted);
    
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
            if((i!=j) && (!ProxMasterArray[i][j]->empty()))
            {
                
               findValleys(ProxMasterArray[i][j], ProxMasterArrayCleaned[i][j], this->max_dist, 1.0 );
               sortProxDist(ProxMasterArrayCleaned[i][j], ProxMasterArrayCleanedSorted[i][j]);
               
               // for each merged prox list write out the images and related amira files 
               
               std::string format(this->outputpath);
               format += "/CELL#";
               format += Utility::inttochar(i+1);
              
               format += "-CELL#";
               format += Utility::inttochar(j+1);
               
               std::string xlsfilename = format + "Proxlist.xls";
               writeProxXls(xlsfilename.c_str(), ProxMasterArrayCleanedSorted[i][j]);
              
               std::string landmarkfilename = format + "Proximity_Landmarks.landmarkAscii";
               writeProxLandmarks(landmarkfilename.c_str(), ProxMasterArrayCleanedSorted[i][j]);
              
              //format += "_roi_";
              //////std::cout<<"inttochar "<<j+1<<" "<<Utility::inttochar(j+1)<<std::endl;
              //////std::cout<<format<<std::endl;
              //format += "_CELL#1-->CELL#2";
//               for(int k=0; k < ProxMasterArrayCleaned[i][j]->size(); k++)
//               {
//                 //std::string format(this->outputpath);
//                 //format += "_CELL#1-->CELL#2";
//                 
//                 
//                 format += Utility::inttochar(k);
//                // (ProxMasterArrayCleaned[i][j]->at(k))->writeProximityImage(format.c_str(), this->transformation);
//                 MasterHxData<<"[ load "<<format<<"/loadROI.hx ]"<<std::endl;
//                
//             
//               }
               
            }
        }
    }
    
    
    MasterHxData.close();
    

   #endif  
};
  
void AxonDendriteProximityFinder::findProximities(bool notestrun)
{   
    ////std::cout<<"In findProximities()"<<std::endl;
    int num_of_cells = amira_graph->getNumberofCells();
    std::vector<Proximity *>* ProxMasterArray[num_of_cells][num_of_cells];
    bool isfirstone = true;
    int prox_counter = 0;
    
    // Init the Prox master array with valid vector pointers
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
              ProxMasterArray[i][j] = new std::vector<Proximity *>;
        }
    }
    
    // For each axon-dendrite pairing
//    for(int axonCount = 0; axonCount < axon_list.size(); axonCount++)
    for(int axonCount = 0; axonCount < axon_list.size(); axonCount++)    
    {
        SimpleAxon currentAxon = axon_list.at(axonCount);
        int AxonID = currentAxon.getID();
        bool hasBouton = false;
        
        
        for(int dendriteCount = 0; dendriteCount < dendrite_list.size(); dendriteCount++)
        {
            
            
            SimpleDendrite currentDendrite = dendrite_list.at(dendriteCount);
            std::vector<double *> * point_list = new std::vector<double *>;
            std::vector<double> distance_list /*= new std::vector<double >*/;
            std::vector<double* > * bouton_list = new std::vector<double * >;
            
            int DendriteID = currentDendrite.getID();
            
            // This variable ensures that edge segments whose distances hover around max_dist 
            // are not accidentally grouped into separate proximities when they should be one
            int numPointsAway = 0;
            
            //if(DendriteID == DENDRITE_ID_1 && AxonID == AXON_ID_2 || DendriteID == DENDRITE_ID_1 && AxonID == AXON_ID_3 || DendriteID == DENDRITE_ID_2 && AxonID == AXON_ID_1 || DendriteID == DENDRITE_ID_2 && AxonID == AXON_ID_3 || DendriteID == DENDRITE_ID_3 && AxonID == AXON_ID_1 || DendriteID == DENDRITE_ID_3 && AxonID == AXON_ID_2)
            if(isTheSameCell(AxonID,DendriteID) == false)
            {
                
                //////std::cout<<AxonID<<" "<<DendriteID<<" "<<std::endl;
            
            currentAxon.setToBegin();
            while(currentAxon.hasNext())
            {
                
                double * axonPoint = currentAxon.getNext();
                
                if(isInSectionBoundingBox(axonPoint))
                {
                    
                    bool isMatched = false;
                    bool firstLoop = true;
                    currentDendrite.setToBegin();
                    while(currentDendrite.hasNext())
                    {
                        double * dendritePoint = currentDendrite.getNext();
                        
                        if(isInSectionBoundingBox(dendritePoint))
                        {
                            
                            
                            
        //                     //////std::cout<<"axonID: "<<currentAxon.edgeID<< "  dendriteID: "<< currentDendrite.edgeID<<std::endl;
                            if(currentAxon.edgeID != currentDendrite.edgeID)
                            {
                                double distance = Utility::euclidean_distance(axonPoint, dendritePoint, 3, 1);
                                if((distance) < this->max_dist)
                                {
                                    
                                if(hasBouton == false)
                                        hasBouton = axonPoint[IS_BOUTON];
                                    
                                    
        // 			     
                                    bool thisIsABouton = axonPoint[IS_BOUTON];
                                    numPointsAway = 0;
                                    isMatched = true;
                                    
                                    point_list->push_back(axonPoint);
                                
                                    if(/*thisIsABouton &&*/ firstLoop)
                                        bouton_list->push_back(axonPoint);
                                    
                                    point_list->push_back(dendritePoint);
                                    firstLoop = false;
        //                             
                                    isfirstone = false;

                                }
                            }
                            
                        }// if end dendrite point in section
                    }
                    
                    if(point_list->size() > 0)                                
                    {
                        //if just exited a proximity zone, OR we're at the end of the current axon
                        if((isMatched == false && numPointsAway > this->numberOfPointsAway) || currentAxon.hasNext() == false)            
                        {                         
                                Edge * currentEdge = currentAxon.getEdge();
                                Edge * currentDedriteEdge = currentDendrite.getEdge();
                                int * IDs = new int[2];
                                IDs[0] = AxonID;
                                IDs[1] = DendriteID;
                                
                                Proximity * prox = new Proximity(point_list, this->original_image, IDs);
                                
                                ////std::cout<<"pushing "<<AxonID<<" "<<DendriteID<<std::endl;
                                
                                ProxMasterArray[amira_graph->getAxonCellIndex(AxonID)][amira_graph->getDendriteCellIndex(DendriteID)]->push_back(prox);

                                /*if(prox->Axon_ID == AXON_ID_1 && prox->Dendrite_ID == DENDRITE_ID_2)
                                proximity_list_1_2.push_back(prox);
                                
                                if(prox->Axon_ID == AXON_ID_2 && prox->Dendrite_ID == DENDRITE_ID_1)
                                proximity_list_2_1.push_back(prox);
                                if(prox->Axon_ID == AXON_ID_1 && prox->Dendrite_ID == DENDRITE_ID_3 )
                                proximity_list_1_3.push_back(prox);
                                if(prox->Axon_ID == AXON_ID_3 && prox->Dendrite_ID == DENDRITE_ID_1)
                                proximity_list_3_1.push_back(prox);
                                if(prox->Axon_ID == AXON_ID_2 && prox->Dendrite_ID == DENDRITE_ID_3)
                                proximity_list_2_3.push_back(prox);
                                if(prox->Axon_ID == AXON_ID_3 && prox->Dendrite_ID == DENDRITE_ID_2)
                                proximity_list_3_2.push_back(prox);
                                
                                proximity_list.push_back(prox);*/
                                
                                //////std::cout<< "after push"<<std::endl;
                                
    //                             this->proximity_list.push_back(prox);
    //                             this->proximity_list_nondecon.push_back(prox_nondecon);
                                prox->setBouton(hasBouton);
                                prox_counter ++;
    //                             
                                double proximity_confidence_value;
                                if(notestrun)
                                proximity_confidence_value = calculateConfidenceValue(currentEdge, bouton_list, distance_list);
                                else 
                                proximity_confidence_value = 0;
    //                             prox->calculateConfidenceValue(confidence_list, distance_list);
                                prox->setConfidenceValue(proximity_confidence_value);
    //                             
                                delete point_list;
    // 			    delete distance_list;
                                delete bouton_list;
                                distance_list.clear();
    //                             confidence_list.clear();
                                point_list = new std::vector<double *>;
    // 			    distance_list = new std::vector<double >;
                                bouton_list = new std::vector<double *>;
                                hasBouton = false;
                                isfirstone = true;
                                
                                
                        }
                        
                    }
    #ifdef DEBUG                
                    //////std::cout<< "numPointsAway: "<<numPointsAway<<std::endl;
                    numPointsAway++;
    #endif		
                } //if end axon point in section
            }
            
        }
        
        }
        
    }
    
    std::string masterhxfilename(this->outputpath);
    
    ////std::cout<<masterhxfilename<<std::endl;
    
    masterhxfilename += "/Masterhxload.hx";
    
    ////std::cout<<masterhxfilename<<std::endl;
    std::ofstream MasterHxData(masterhxfilename.c_str());
    
    MasterHxData << "# Amira Script"<< std::endl;
    MasterHxData << " "<< std::endl;
    
    
    
    for(int i = 0; i < num_of_cells; i++)
    {
        for(int j = 0; j < num_of_cells; j++)
        {
            if((i!=j) && (!ProxMasterArray[i][j]->empty()))
            {
                
               ////std::cout<<"axon "<<i+1 <<"dendrite"<<j+1 <<std::endl;
               mergreProximites(ProxMasterArray[i][j]);
               
               // for each merged prox list write out the images and related amira files 
               
               std::string format(this->outputpath);
              format += "/CELL#";
              format += Utility::inttochar(i+1);
              //////std::cout<<"inttochar "<<i+1<<" "<<Utility::inttochar(i+1)<<std::endl;
              format += "-CELL#";
              format += Utility::inttochar(j+1);
              //////std::cout<<"inttochar "<<j+1<<" "<<Utility::inttochar(j+1)<<std::endl;
              //////std::cout<<format<<std::endl;
              //format += "_CELL#1-->CELL#2";
              for(int k=0; k < ProxMasterArray[i][j]->size(); k++)
              {
                  //std::string format(this->outputpath);
                  //format += "_CELL#1-->CELL#2";
                  if(!ProxMasterArray[i][j]->empty())
                  {
                      format += "_roi_";
                      format += Utility::inttochar(k);
                      //(ProxMasterArray[i][j]->at(k))->writeProximityImage(format.c_str(), this->transformation);
                      MasterHxData<<"[ load "<<format<<"/loadROI.hx ]"<<std::endl;
                     
                      
                    
                  }
              }
              
              
        
               
            }
        }
    }
    
    
    MasterHxData.close();
    

};





// void AxonDendriteProximityFinder::markBoutons()
// {
//         //////std::cout<< "In markBoutons" << std::endl;
//         
//         std::vector< Edge * > * edges = amira_graph->edgesPointer();
//         int numOfEdges = edges->size();
//         bool inBouton = false;
//         
//         
// //        int pixel_value = 0;
//         float avg_std_dev = getAverageStandardDeviation();
//         
//         ////////std::cout<< "avg_std_dev: " << avg_std_dev << std::endl;
//         
//         
//         for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
//         {                       
//                 Edge * currentEdge = edges->at(pos);
//                 unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
// //           SimpleAxon currentAxon = axon_list.at(axon_count);
// //           Edge * currentEdge = currentAxon.getEdge();
//           std::list< double * >::iterator it;
//           
//           for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
//           {
//                 double * coords = *it;
//                 
//                 if(coords[IS_BOUTON]){
// //             
//                 float brightness = coords[LOCAL_BRIGHTNESS];
//                 ////////std::cout<<"brightness : "<<brightness<<std::endl;
//                 float std_dev = coords[LOCAL_SIGMA];
//                 ////////std::cout<<"local sigma : "<<std_dev<<std::endl;
//                 float avg_bright =  averageLocalValues(currentEdge, it, LOCAL_BRIGHTNESS, 3);
//                 ////////std::cout<<"average brightness : "<<avg_bright<<std::endl;
//                 float avg_rad = averageLocalValues(currentEdge, it, SURFACE, 3);;
//                 ////////std::cout<<"average surface : "<<avg_rad<<std::endl;
//                  float radius = coords[SURFACE];
//                  ////////std::cout<<"radius : "<<radius<<std::endl;
// //                 
// //              //////std::cout << "Point at (x, y): ("<< coords[X_COORD] << ", " << coords[Y_COORD] << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
//                 
//                 
//              if( (brightness > (avg_bright + 10)) && (radius > avg_rad)/* && (radius >= 3.5)*/ ){ 
//                coords[IS_BOUTON] = 1;
// //                double * tmp = new double[2];
// //                tmp[X_COORD] = coords[X_COORD];/*/XYSAMPLING;*/
// //                tmp[Y_COORD] = coords[Y_COORD];/*/XYSAMPLING;*/
// //                tmp[Z_COORD] = coords[Z_COORD];/*/ZSAMPLING;*/
// //              this->all_bouton_list->push_back(tmp);
//              }
//              else
//                coords[IS_BOUTON] = 0;
//                } 
// /*                
//                 if( coords[IS_BOUTON] && radius > 0.5){ 
//                   coords[IS_BOUTON] = 1;
//                 inBouton = true;
//                 }
//                 else
//                   coords[IS_BOUTON] = 0;*/
// 
//           }
//         }
//           
//           
//           
//         
//   
// };

// float AxonDendriteProximityFinder::boutonConfidenceValue(Edge * currentEdge, double * coords ){
//     
//     //computes the likelihood of these coordinates being an actual bouton
//     float confidence_value, average_surface, average_local_bright, bright_ratio, surface_ratio;
//     std::list<double*>::iterator axon_it;
//     
//     
//     for(axon_it = currentEdge->edgePointCoordinates.begin(); axon_it != currentEdge->edgePointCoordinates.end(); axon_it ++){
//     if(*axon_it == coords){
//     //////std::cout << "Point at (x, y): ("<< coords[X_COORD] << ", " << coords[Y_COORD] << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
//     average_surface = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//     //////std::cout<<"average surface "<<average_surface<<std::endl;
//     average_local_bright = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
//     //////std::cout<<"average local brightness "<<average_local_bright<<std::endl;
//     }
//     }
//     if(average_local_bright != 0){
//     bright_ratio = coords[LOCAL_BRIGHTNESS]/average_local_bright;
//     //////std::cout<<"bright ratio : "<<bright_ratio<<std::endl;
//     }
//     else
//         bright_ratio = 1;
//     if(average_surface != 0){    
//     surface_ratio = coords[SURFACE]/average_surface;
//     //////std::cout<<"surface ratio : "<<surface_ratio<<std::endl;
//       
//     }
//     else
//         surface_ratio = 1;
//     
//     confidence_value =  (bright_ratio + surface_ratio)/2;
//     
//     return confidence_value;
// }

double AxonDendriteProximityFinder::calculateConfidenceValue(Edge * currentEdge/*, double * coords*/, std::vector<double*> * confidenceList, std::vector<double>  distances){ 
  //     computes the likelihood of these coordinates really being a bouton
  
  
    float confidence_value, average_surface, average_local_bright, bright_ratio, surface_ratio, average_confidence;
    std::list<double*>::iterator axon_it;
    std::list<float> confidence_values;
    int numOfBoutons = confidenceList->size();
    ////////std::cout<<"number of boutons : "<<numOfBoutons<<std::endl;
  if(numOfBoutons > 0){  
    for(int i = 0; i< numOfBoutons; i++){
    double * coords = confidenceList->at(i);  
    for(axon_it = currentEdge->edgePointCoordinates.begin(); axon_it != currentEdge->edgePointCoordinates.end(); axon_it ++){
    if(*axon_it == coords){
   // //////std::cout << "Point at (x, y): ("<< coords[X_COORD] << ", " << coords[Y_COORD] << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
    average_surface = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
    ////////std::cout<<"average surface "<<average_surface<<std::endl;
    average_local_bright = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
    ////////std::cout<<"average local brightness "<<average_local_bright<<std::endl;
    }
    }
    if(average_local_bright != 0){
    bright_ratio = coords[LOCAL_BRIGHTNESS]/average_local_bright;
    ////////std::cout<<"bright ratio : "<<bright_ratio<<std::endl;
    }
    else
        bright_ratio = 1;
    if(average_surface != 0){    
    surface_ratio = coords[SURFACE]/average_surface;
    ////////std::cout<<"surface ratio : "<<surface_ratio<<std::endl;
      
    }
    else
        surface_ratio = 1;
    
    confidence_value =  (bright_ratio + surface_ratio)/2;
    confidence_values.push_back(confidence_value);
    average_confidence += confidence_value;
  
    }
    
    average_confidence = average_confidence/confidence_values.size();
    ////////std::cout<<"average confidence : "<<average_confidence<<std::endl;
    confidence_values.clear();
    
    
  }
  if(numOfBoutons == 0)
    average_confidence = 0;
  
    
//     double confidence_val;
//     
//     double proximity_avg_bright;
//     double proximity_avg_std_dev;
//     double proximity_avg_radius;
    double proximity_avg_distance;
    double min_distance = 100;
    int counter = 0;
    int counter2 = 0;
//     double bouton_confidence_val;
//     int count = 0;
// //     int boutons_detected = confidenceList.size();
//     double confidence;
//     
   
    
    for(int i = 0; i < distances.size(); i++){
        
            proximity_avg_distance += distances.at(i);
            counter2++;
            if(distances.at(i) < min_distance)
            min_distance = distances.at(i);
            
//             if(min_distance < 0)
//                 min_distance = 0;
            
            
            
            
        }
        
        ////////std::cout<<"minimum distance: "<<min_distance<<std::endl;
        
        
        
        
        
        proximity_avg_distance = proximity_avg_distance/counter2;
//        //////std::cout<<"average_distance:  "<<proximity_avg_distance<<std::endl;
//        //////std::cout<<"min_distance:  "<<min_distance<<std::endl;
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
    
//    confidence_val = proximity_avg_bright + proximity_avg_std_dev + proximity_avg_radius /*+ 2* min_distance*/;
//      //////std::cout<<"confidence list size: "<<confidenceList->size()<<std::endl;
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
    
//     //////std::cout<<"bouton_confidence_value: "<<bouton_confidence_val<<std::endl;
//     //////std::cout<<"min_distance:  "<<min_distance<<std::endl;
    //this->confidence_value = (bouton_confidence_val + 1/(0.1 + min_distance))/2;
    //this->confidence_value = bouton_confidence_val;
   // //////std::cout<<"confidence value: "<<this->confidence_value<<std::endl;
    
    ////////std::cout<<"average confidence : "<<average_confidence<<std::endl;
//     return (0.6*average_confidence + 1/(1 + min_distance));
    return min_distance;
    
    

    
    
}


// void AxonDendriteProximityFinder::markLocalMaximums()
// {
//     
//     //////std::cout<<"in MarkLocalMaximums"<<std::endl;
//     double buffer = MAX_DIST;
// //    SimpleAxon currentAxon = axon_list.at(axon_count);
// //    Edge * currentEdge = currentAxon.getEdge();
//     std::vector< Edge * > * edges = amira_graph->edgesPointer();
//     unsigned int numOfEdges = edges->size();
// //    double * currentPoint = currentAxon.getEdgePointCoordinates();
//     std::list< double * >::iterator axon_it;
// //    unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
//     float bright_offset = 5;
//     float radius_offset = 0;
//     
//     for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
//         {                       
//                 Edge * currentEdge = edges->at(pos);
//                 unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
//     
//     
//     
//     for(axon_it = currentEdge->edgePointCoordinates.begin() ; axon_it != currentEdge->edgePointCoordinates.end(); axon_it++){
//         
//                         double * current_coords = *axon_it;
//                         
//                         if(isInSectionBoundingBox(current_coords)){
//                         
//                         float diameter = current_coords[SURFACE];
//                         float sum = 0;
//                         float average = 0;
//                         long counter = 1;
//                         bool isRadiusMax = false;
// 			
// 			
// #ifdef DEBUG                        
//                         int x_pos, y_pos, z_pos;
//                                         
//                         x_pos = rint(current_coords[X_COORD]/XYSAMPLING);
//                         y_pos = rint(current_coords[Y_COORD]/XYSAMPLING);
//                         z_pos = rint(current_coords[Z_COORD]/ZSAMPLING);
//                                         
// 
// 			//////std::cout<<"x: "<<x_pos<<"  y: "<<y_pos<<"  z: "<<z_pos<<std::endl;
// #endif                        
//                         std::list< double * >::iterator next_it = axon_it;              
//                         std::list< double * >::iterator prev_it = axon_it;
//                         ++next_it;
//                         --prev_it;
//                         
//                         std::list< double * >::iterator next_next_it = next_it;         
//                         std::list< double * >::iterator prev_prev_it = prev_it;
//                         ++next_next_it;
//                         --prev_prev_it; 
//          
// 
//     
//                         if(axon_it == currentEdge->edgePointCoordinates.begin())     //if at beginning of edge
//                         {
//                                 double * next_coords = *next_it;
//                                 
// //                                //////std::cout<<"at the beginning!"<<std::endl;
//                                 
//                                 
//                         }
//                         else if(next_it == currentEdge->edgePointCoordinates.end())       //if at end of edge
//                         {
//                                 double * prev_coords = *prev_it;
//                                 
// //                                //////std::cout<<"at the end!"<<std::endl;
//                           
//                                 
//                         }
//                         else if(prev_it == currentEdge->edgePointCoordinates.begin() && edgeSize >= 4)  //second bouton from beginning
//                         {
//                                 
//                                 double * next_coords = *next_it;
//                                 double * prev_coords = *prev_it;
//                                 double * next_next_coords = *next_next_it;
//                                 
//                                 if(isInSectionBoundingBox(next_coords) && isInSectionBoundingBox(prev_coords) && isInSectionBoundingBox(next_next_coords)){
// #ifdef DEBUG                                
//                                //////std::cout<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                //////std::cout<<current_coords[SURFACE]<<std::endl;
// #endif                                
//                           
//                                 if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
//                                   current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
//                                 {
//                                     
//                                     
//                                     average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//                                     current_coords[AVG_SURF] = average;
//                                     
//                                     if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
//                                         current_coords[SURFACE] >= next_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     else if(next_coords[SURFACE] >= next_next_coords[SURFACE] && 
//                                         next_coords[SURFACE] >= current_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     
// //                                  if(isRadiusMax)
//                                     if(current_coords[SURFACE] > average+radius_offset)
//                                     {
//                                         
//                                         average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
//                                         current_coords[AVG_BRIGHT] = average;
//                                         
//                                        
//                                         if(current_coords[LOCAL_BRIGHTNESS] > (average + bright_offset)){                                            
//                                           current_coords[IS_BOUTON] = 1;
// #ifdef DEBUG					  
//                                         //////std::cout<<"just marked a bouton!"<<std::endl;
// #endif
//                                         }
//                                         else{
//                                           current_coords[IS_BOUTON] = 0;
// #ifdef DEBUG
//                                           //////std::cout<<"just unmarked a bouton!"<<std::endl;
// #endif
//                                         }
//                                     }
//                                 }
//                                 else
//                                   current_coords[IS_BOUTON] = 0;
// #ifdef DEBUG				
//                                 //////std::cout<<"no Bouton!"<<std::endl;
// #endif				
//                         }  //end if  
//                         }
// 
//                         else if(next_next_it == currentEdge->edgePointCoordinates.end() && edgeSize >= 4)  //second bouton from end
//                         {
//                                 double * next_coords = *next_it;
//                                 double * prev_coords = *prev_it;
//                                 double * prev_prev_coords = *prev_prev_it;
//                                  if(isInSectionBoundingBox(next_coords) && isInSectionBoundingBox(prev_coords) && isInSectionBoundingBox(prev_prev_coords)){
//                                 
//                                 
//                           
//                                 if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
//                                   current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
//                                 {
//                                     
//                                     average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//                                     current_coords[AVG_SURF] = average;
// //                                    //////std::cout<<"Brightness:    "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                     
// //                                    //////std::cout<<"average:  "<<average<<std::endl;
// //                                   //////std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                     
//                                     
//                                     if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
//                                         current_coords[SURFACE] >= next_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     else if(prev_coords[SURFACE] >= prev_prev_coords[SURFACE] && 
//                                         prev_coords[SURFACE] >= current_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     
// //                                  if(isRadiusMax)
//                                     if(current_coords[SURFACE] > average+radius_offset)
//                                     {
// 
// //                                        //////std::cout<<"Surface:   "<<current_coords[SURFACE]<<std::endl;
//                                         average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
//                                         current_coords[AVG_BRIGHT] = average;
//                                         
// //                                        //////std::cout<<"average:  "<<average<<std::endl;
// //                                        //////std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                       
//                                         if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset)){
//                                           current_coords[IS_BOUTON] = 1;
// //                                            //////std::cout<<"just marked a bouton!"<<std::endl;
//                                         }
//                                         else{
//                                           current_coords[IS_BOUTON] = 0;
// //                                          //////std::cout<<"just unmarked a bouton!"<<std::endl;
//                                             
//                                         }  
//                                     }
//                                 }
//                                 else{
//                                   current_coords[IS_BOUTON] = 0;
// //                                  //////std::cout<<"just unmarked a bouton!"<<std::endl;
//                                    
//                                 }
//                         }  //end if 
//                         }
//                         else if(edgeSize >= 5)                                                          //in middle of edge
//                         {
//                                 double * next_coords = *next_it;
//                                 double * prev_coords = *prev_it;
//                                 double * next_next_coords = *next_next_it;
//                                 double * prev_prev_coords = *prev_prev_it;
//                                 
//                                  if(isInSectionBoundingBox(next_coords) && isInSectionBoundingBox(prev_coords) && isInSectionBoundingBox(next_next_coords) && isInSectionBoundingBox(prev_prev_coords)){
//                                 
//                           
//                                 if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
//                                   current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS])
//                                 {
// 
//                                     average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//                                     current_coords[AVG_SURF] = average;
// //                                    //////std::cout<<"Brightness:     "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                     
// //                                    //////std::cout<<"average:  "<<average<<std::endl;
// //                                    //////std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                     
//                                     if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
//                                         current_coords[SURFACE] >= next_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     else if(next_coords[SURFACE] >= next_next_coords[SURFACE] && 
//                                         next_coords[SURFACE] >= current_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     else if(prev_coords[SURFACE] >= prev_prev_coords[SURFACE] && 
//                                         prev_coords[SURFACE] >= current_coords[SURFACE] )
//                                       isRadiusMax = true;
//                                     
// //                                  if(isRadiusMax)
//                                     if(current_coords[SURFACE] > average+radius_offset)
//                                     {
// 
//                                           average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
//                                           current_coords[AVG_BRIGHT] = average;
// //                                          //////std::cout<<"Surface:   "<<current_coords[SURFACE]<<std::endl;
// //                                          //////std::cout<<"average:                          "<<average<<std::endl;
// //                                          //////std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                           
//                                           if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset)){
//                                             current_coords[IS_BOUTON] = 1;
//  //                                           //////std::cout<<"just marked a bouton!"<<std::endl;
//                                           }
//                                           else{
//                                             current_coords[IS_BOUTON] = 0;
// //                                            //////std::cout<<"just unmarked a bouton!"<<std::endl;
//                                               
//                                         }
//                                     }
//                                 }
//                                 else
//                                   current_coords[IS_BOUTON] = 0;{
// //                                  //////std::cout<<"just unmarked a bouton!"<<std::endl;
//                                   }
//                         }   //end if
//                         }
//                         }// end if axon is in section
//                 }
// 
//         
//         }
//     
//             
//     //std::vector< Edge * > * edges = amira_graph->edgesPointer();
// 
// };


float AxonDendriteProximityFinder::averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection)
{
        
#ifdef DEBUG    
        //////std::cout<<"in averageLocalValues"<<std::endl;
#endif        
        
        float sum = 0;
        int counter = 0, i=0;
                
        unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
        std::list< double * >::iterator edge_it;

        for(edge_it = centerPoint, i=0; edge_it != currentEdge->edgePointCoordinates.end() && i<=numToAverageInEachDirection; edge_it++, i++) 
        {
                double * coords = *edge_it;
                
                if(isInSectionBoundingBox(coords)){
                sum += coords[property];
                counter ++;
                }

        }
        
        if(centerPoint != currentEdge->edgePointCoordinates.begin())
        {
          for(edge_it = (--centerPoint), i=0; edge_it != currentEdge->edgePointCoordinates.begin() && i<numToAverageInEachDirection; edge_it--, i++) 
          {
                  double * coords = *edge_it;
                  if(isInSectionBoundingBox(coords)){
                  sum += coords[property];
                  counter++;
                  }

          }
        }

  return sum/(counter);
};






float AxonDendriteProximityFinder::calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        //////std::cout<< "In Calculate Local Standard Deviation" << std::endl;
#endif
        Image2DType::RegionType local_region;
        Image2DType::IndexType local_index;
        Image2DType::SizeType local_size;

        local_index[0] = x0-radius;
        local_index[1] = y0-radius;

        local_size[0] = 2*radius+1;
        local_size[1] = 2*radius+1;

        if(local_index[0]<0)
        {       
                local_size[0] = local_size[0] + local_index[0];
                local_index[0]=0;
        }
        if(local_index[1]<0)
        {       
                local_size[1] = local_size[1] + local_index[1];
                local_index[1]=0;
        }
        if(local_index[0] + local_size[0] >= image_plane->GetLargestPossibleRegion().GetSize(0))
        {       
            if(local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0))
                local_size[0] = 0;
            else
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0)-1-local_index[0];
                local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0);
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {
            if(local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1))
                local_size[1] = 0;
            else
                local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1)-1-local_index[1];
                local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1);
        }

// #ifdef DEBUG
//         //////std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
// #endif

        local_region.SetIndex(local_index);
        local_region.SetSize(local_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIterator2DType histo_it(image_plane, local_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }
        
        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();
        
        CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();

        covariance_algorithm->SetInputSample( histo_sample );
        covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
        covariance_algorithm->Update();

        CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
        float variance = *(covariance.operator[](0));
        float standard_deviation = std::sqrt(variance);

        return standard_deviation;
};


float AxonDendriteProximityFinder::calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        //////std::cout<< "In Calculate Local Brightness" << std::endl;
#endif
        Image2DType::RegionType local_region;
        Image2DType::IndexType local_index;
        Image2DType::SizeType local_size;

        local_index[0] = x0-radius;
        local_index[1] = y0-radius;

        local_size[0] = 2*radius+1;
        local_size[1] = 2*radius+1;
        
//      //////std::cout<<local_size[0]<<std::endl;
//      //////std::cout<<local_size[1]<<std::endl;        
//      //////std::cout<<local_index[0]<<std::endl;
        
        
        

        if(local_index[0]<0)
        {       
                local_size[0] = local_size[0] + local_index[0];
                local_index[0]=0;
        }
        if(local_index[1]<0)
        {       
                local_size[1] = local_size[1] + local_index[1];
                local_index[1]=0;
        }
        if(local_index[0] + local_size[0] >= image_plane->GetLargestPossibleRegion().GetSize(0))
        {
            if(local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0))
                local_size[0] = 0;
            else
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0) - local_index[0] - 1;
 //               local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0);
                
//             else    
//             local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0) - local_index[0] - 1;
//             //local_size[0] = local_index[0] - image_plane->GetLargestPossibleRegion().GetSize(0)-1;
//                      
                            
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {   
            if(local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1))
                local_size[1] = 0;
            else
            local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1) - local_index[1] - 1;
//            local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1);
//             float diff = local_index[1] + 1;
//             if(image_plane->GetLargestPossibleRegion().GetSize(0) <= diff)
//                 local_size[1] = 0;
//            
//             else
//             local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1) - local_index[1] - 1;
//             //local_size[1] = local_index[1] - image_plane->GetLargestPossibleRegion().GetSize(1)-1;
        }
        
        
        
        

#ifdef DEBUG
        //////std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
#endif

        local_region.SetIndex(local_index);
        local_region.SetSize(local_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIterator2DType histo_it(image_plane, local_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }

        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();

        float mean = mean_algorithm->GetOutput()->GetElement(0);

        return mean;
};



float AxonDendriteProximityFinder::calculateThreshold(ImageType::Pointer image, float x0, float y0, float z0)
{
#ifdef DEBUG    
        //////std::cout<< "In Calculate Threshold" << std::endl;
#endif
        ImageType::RegionType threshold_region;
        ImageType::IndexType threshold_index;
        ImageType::SizeType threshold_size;

        threshold_index[0] = x0-6;
        threshold_index[1] = y0-6;
        threshold_index[2] = z0;

        threshold_size[0] = 12;
        threshold_size[1] = 12;
        threshold_size[2] = 1;

        if(threshold_index[0]<0)
        {       
                threshold_size[0] = threshold_size[0] + threshold_index[0];
                threshold_index[0]=0;
        }
        if(threshold_index[1]<0)
        {       
                threshold_size[1] = threshold_size[1] + threshold_index[1];
                threshold_index[1]=0;
        }
        if(threshold_index[0] + threshold_size[0] >= image->GetLargestPossibleRegion().GetSize(0))
        {
                threshold_size[0] = image->GetLargestPossibleRegion().GetSize(0)-1 -threshold_index[0];
                
        }
        if(threshold_index[1] + threshold_size[1] >= image->GetLargestPossibleRegion().GetSize(1))      
        {
                threshold_size[1] = image->GetLargestPossibleRegion().GetSize(1)-1 -threshold_index[1];
        }
        if(threshold_index[0] < image->GetLargestPossibleRegion().GetSize(0) && threshold_index[1] < image->GetLargestPossibleRegion().GetSize(1) && threshold_index[2] < image->GetLargestPossibleRegion().GetSize(2))
        {

#ifdef DEBUG
        //////std::cout<< "x= " << x0<< " y= " << y0 << " index= " << threshold_index << " size= " << threshold_size << std::endl << std::flush;
#endif

        threshold_region.SetIndex(threshold_index);
        threshold_region.SetSize(threshold_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIteratorType histo_it(image, threshold_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }

        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();

        CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();

        covariance_algorithm->SetInputSample( histo_sample );
        covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
        covariance_algorithm->Update();

        float mean = mean_algorithm->GetOutput()->GetElement(0);
        CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
        float variance = *(covariance.operator[](0));
        float standard_deviation = std::sqrt(variance);

#ifdef DEBUG    
        //////std::cout<< "Calculated Threshold" << mean << std::endl << std::flush;
#endif
        return (mean );
        }
        else
        return 0;    
};

float AxonDendriteProximityFinder::getAverageStandardDeviation()
{
#ifdef DEBUG
        //////std::cout<< "In getAverageStandardDeviation" << std::endl;
#endif        
  
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();        
        
        float sum = 0;
        int counter = 0;
        

        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();

                //for every point along edge
                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                {
                        double * coords = *edge_it;
                        
                        if(isInSectionBoundingBox(coords)){
                        
                        if(coords[IS_BOUTON])
                        {
                            sum += coords[LOCAL_SIGMA];
                            counter++;
                        }
                        }
                }
        }
  
  return sum/counter;
};

int AxonDendriteProximityFinder::setRadiusAndLocalBrightness(ImageType::Pointer image)
{

        //////std::cout<< "In SetRadiusAndLocalBrightness " << std::endl;
        
//        SimpleAxon currentAxon = axon_list.at(axon_count);
//        Edge * currentEdge = currentAxon.getEdge();

        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        

        unsigned int numOfEdges = edges->size();
        unsigned int end_x = image->GetLargestPossibleRegion().GetSize(0);
//        //////std::cout<<end_x<<std::endl;
        unsigned int end_y = image->GetLargestPossibleRegion().GetSize(1);
//        //////std::cout<<end_y<<std::endl;
        unsigned int end_z = image->GetLargestPossibleRegion().GetSize(2);
//        //////std::cout<<end_z<<std::endl;
        
//        image->Print(//////std::cout);
        
        int x_pos, y_pos, z_pos;

        Image2DType::Pointer image_plane_5;
        Image2DType::Pointer image_plane_3;
        
        int count = 1;

        for(int z = 0; z < end_z; z++)
        {
                image_plane_5 = getImagePlane(z, 5, image);
                image_plane_3 = getImagePlane(z,3,image);

//                image_plane_3->Print(//////std::cout);
//                image_plane_5->Print(//////std::cout);
        
                
                if(numOfEdges == 0)
                {
                        //////std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                        return -1;
                }
                else
                {

                        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                                
                        {      
                                Edge * currentEdge = edges->at(pos);
                                std::list< double * >::iterator edge_it;        
                                std::list< double * >::iterator next_it;
                                
//                                int counter = 0;


                                //for every point along edge
                                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                                {
                                        double * coords = *edge_it;
                                        if(isInSectionBoundingBox(coords)){                                        
//                                      counter ++;                                  
                                        x_pos = rint(coords[X_COORD]/XYSAMPLING);
                                        y_pos = rint(coords[Y_COORD]/XYSAMPLING);
                                        z_pos = rint(coords[Z_COORD]/ZSAMPLING);
                                        
//                                      //////std::cout<<"x: "<<x_pos<<"  y: "<<y_pos<<"  z: "<<z_pos<<std::endl;
                                        
                                        if(z_pos == z)
                                        {
                                                /*****************************************************************/
                                                /*set radius info                                                */
                                                /*****************************************************************/
                                          
                                                float x0 = x_pos;
                                                float y0 = y_pos;
                                                float threshold = 0;

//                                              threshold = CalculateThreshold(image_plane_3, x0, y0);
//                                              //////std::cout<<"try to set threshold"<<std::endl;
                                                threshold = coords[THRESHOLD];
//                                              //////std::cout<<"threshold :"<<coords[THRESHOLD]<<std::endl;
                                                
                                                unsigned int nr_of_rays = 20;
                                                float ray_length = 0.5;
                                                int index = 0;

                                                std::vector<VECTOR *> vectors;
                                                VECTOR * adjustment_vect;
                                                float distance = 0;

                                                x0 = x0 + 0.5;
                                                y0 = y0 + 0.5;
        
                                                for(unsigned n=1; n<=nr_of_rays; n++)
                                                {
                                                        VECTOR * tmp;

                                                        tmp = sendRay(x0, y0, ray_length, nr_of_rays, threshold, n, image_plane_3);
                                                        

                                                        vectors.push_back(tmp);
//                                                      //////std::cout<< "Distance: " << tmp_distance << std::endl;
                                                }

                                                index = addUpOpposingRaysAndFindShortest(vectors, &adjustment_vect, &distance);

                                                coords[SURFACE] = distance;
                                                
                                                /*****************************************************************/
                                                /*adjust midline based on radius info                            */
                                                /*****************************************************************/
                                                next_it = edge_it;
                                                next_it++;
                                                
//                                                 if(edge_it != currentEdge->edgePointCoordinates.begin() && next_it != currentEdge->edgePointCoordinates.end())
//                                                 {
// //                                                  //////std::cout<<"old       "<<coords[X_COORD]<<"         "<<coords[Y_COORD]<<std::endl;  
//                                                   float new_x = coords[X_COORD] + adjustment_vect->coords[X_COORD] * XYSAMPLING;
//                                                   float new_y = coords[Y_COORD] + adjustment_vect->coords[Y_COORD] * XYSAMPLING;
//                                                   
//                                                   if((new_x/XYSAMPLING) < end_x-1 && (new_y/XYSAMPLING) < end_y-1)
//                                                   {
//                                                           coords[X_COORD] = new_x;
//                                                           coords[Y_COORD] = new_y;
//                                                           //////std::cout<<"new       "<<coords[X_COORD]<<"         "<<coords[Y_COORD]<<std::endl;  
//                                                   }
//                                                 }
                                                
                                                /*****************************************************************/
                                                /*set brightness info                                            */
                                                /*****************************************************************/
                                                x0 = rint(coords[X_COORD]/XYSAMPLING);
                                                y0 = rint(coords[Y_COORD]/XYSAMPLING);
                                                
//                                                //////std::cout<<x0<<"    "<<y0<<std::endl;
                                        
                                                coords[LOCAL_BRIGHTNESS] = calculateLocalBrightness(image_plane_5, x0, y0, 1);
                                                coords[LOCAL_SIGMA] = calculateLocalStandardDeviation(image_plane_5, x0, y0, 2);
//                                                    coords[LOCAL_SIGMA] = 1;
                                        }
                                        
                                        } // end if    
                                }
                        }
                }
        }       

        smoothRadii();
        return 0;
};

VECTOR * AxonDendriteProximityFinder::sendRay(float x0, float y0, float ray_length, unsigned int nr_of_rays, float threshold, unsigned int n, Image2DType::Pointer image_plane)
{
#ifdef DEBUG
        //////std::cout<< "In SendRay" << std::endl;
#endif
        VECTOR * tmp_vect = new VECTOR;

        float phi = n*(2*PI/nr_of_rays);
        float x_r = 0, y_r = 0, x_f = x0, y_f = y0;
        
//      //////std::cout << "Phi: " << n*360/nr_of_rays << std::endl;

        float grey_value = bilinearInterpolation(x_f, y_f, image_plane);
        float old_grey_value = 0;
        unsigned int counter = 0;
        float gradient, max_gradient = 0;
        bool stillbrighter;
        int grey_counter = 0;
        
        bool shouldPrint = false;
        
        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());
        Image2DType::IndexType corner_index;


        while(grey_value >= threshold /*grey_counter <= 2*/)
        {
                  
                x_f = x_f + (ray_length * cos(phi));
                y_f = y_f + (ray_length * sin(phi));
                if(x_f <= 1 || y_f <= 1 || x_f >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_f >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                        break;
                else
                {
                        corner_index[0] = x_f;
                        corner_index[1] = y_f;
                        it.SetIndex(corner_index);
                        
                        old_grey_value = grey_value;
                        grey_value = it.Get();
                        
                        gradient = abs(old_grey_value-grey_value);
                        
                        if(gradient > max_gradient){
                            gradient = max_gradient;
                            stillbrighter = true;
                            grey_counter = 0;
                                                    }
                        if(gradient <= max_gradient){
                            stillbrighter = false;
                            grey_counter ++;
                            
                        }                            
                  
                  
                        //grey_value = BilinearInterpolation(x_f, y_f, image_plane);
                        counter++;
                }
                  
        }

        x_f = x_f - x0;
        y_f = y_f - y0;
        
        tmp_vect->coords[X_COORD] = x_f;
        tmp_vect->coords[Y_COORD] = y_f;
        tmp_vect->magnitude = sqrt((x_f*x_f)+(y_f*y_f));
//      tmp_vect->magnitude = ray_length * counter;
        
        return tmp_vect;
};


float AxonDendriteProximityFinder::bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane)
{
#ifdef DEBUG
        //////std::cout<< "In BilinearInterpolation" << std::endl;
#endif
        /*unsigned*/ int x1 = (/*unsigned*/ int) (x-1), y1 = (/*unsigned*/ int)(y-1);         //lower corners
        /*unsigned */int x2 = (/*unsigned*/ int)(x+1), y2 = (/*unsigned*/ int)(y+1);  //upper corners

        /*unsigned*/ char Q11 = 0, Q21 = 0, Q22 = 0, Q12 = 0;// grey values of corners  
        
//          //////std::cout<<x<<"   "<<y<<std::endl;
//          //////std::cout<<x1<<"   "<<y1<<std::endl;
//          //////std::cout<<x2<<"   "<<y2<<std::endl;
        
        
//        image_plane->Print(//////std::cout);

//        //////std::cout<<"GetLargestPossibleRegion"<<std::endl;
        
        int y_size, x_size;
        
        x_size = image_plane->GetLargestPossibleRegion().GetSize(0);
        y_size = image_plane->GetLargestPossibleRegion().GetSize(1);
        
        if(x < x_size && x >= 0 && y < y_size && y >=0){
        
        
        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());
        
//        //////std::cout<<"Got the largest possible region"<<std::endl;

        Image2DType::IndexType corner_index;
        corner_index[0] = x1;
//        //////std::cout<<corner_index[0]<<std::endl;
        corner_index[1] = y1;
//        //////std::cout<<corner_index[1]<<std::endl;
        it.SetIndex(corner_index);
//        //////std::cout<<"set index"<<std::endl;
        Q11 = it.Get();
//        //////std::cout<<"char gets index"<<std::endl;
//       //////std::cout<<" corner 1 1"<<std::endl;
        ;

        corner_index[0] = x2;
        corner_index[1] = y1;
        it.SetIndex(corner_index);
        Q21 = it.Get();
        
//        //////std::cout<<" corner 2 1"<<std::endl;

        corner_index[0] = x2;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q22 = it.Get();
        
//        //////std::cout<<" corner 2 2"<<std::endl;

        corner_index[0] = x1;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q12 = it.Get();
        
//        //////std::cout<<" corner 1 2"<<std::endl;

        float P = 0;
        float R1 = 0, R2 = 0;

        R1=((x2-x)/(x2-x1)*Q11)+((x-x1)/(x2-x1)*Q21);
        R2=((x2-x)/(x2-x1)*Q12)+((x-x1)/(x2-x1)*Q22);
        P =(((y2-y)/(y2-y1))*R1)+(((y-y1)/(y2-y1))*R2);

        return P;
        
        }
    else return 0;    
};

Image2DType::Pointer AxonDendriteProximityFinder::getImagePlane(int z, int depth, ImageType::Pointer input_image)
{
#ifdef DEBUG    
        //////std::cout<< "In GetImagePlane" << std::endl;
#endif
  
  
        Image2DType::Pointer image_plane = Image2DType::New();

        Image2DType::IndexType target_index;
        Image2DType::SizeType target_size;

        target_index[0] = 0;
        target_index[1] = 0;

        target_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
        target_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
        
        #ifdef DEBUG    
        //////std::cout<< target_size << std::endl;
        #endif

        Image2DType::RegionType target_region;

        target_region.SetSize(target_size);
        target_region.SetIndex(target_index) ;

        image_plane->SetRegions(target_region);

        image_plane->Allocate();
        image_plane->FillBuffer(0);

        Iterator2DType target_it(image_plane, image_plane->GetLargestPossibleRegion()); 

        int end = input_image->GetLargestPossibleRegion().GetSize(2);
        int displacement = (depth-1)/2;
        
        for(int margin = z-displacement; margin <= z+displacement && margin < end; margin++)
        {
                if(margin < 0)
                        margin = 0;

                ImageType::RegionType input_region;
                ImageType::IndexType input_index;
                ImageType::SizeType input_size;
        
                input_index[0] = 0;
                input_index[1] = 0;
                input_index[2] = margin;
                        
                input_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
                input_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
                input_size[2] = 1;
                

        
                input_region.SetIndex(input_index);
                input_region.SetSize(input_size);
                
                ConstIteratorType source_it(input_image, input_region);
        
                for(source_it.GoToBegin(), target_it.GoToBegin(); !source_it.IsAtEnd(); ++source_it, ++target_it)
                {
                        unsigned char projection_value = 0;
                        unsigned char value = 0;
                        
                        projection_value = target_it.Get();
                        value = source_it.Get();
                        
                        if(value > projection_value)
                                target_it.Set(value);
                }
        }

        return image_plane;
};


int AxonDendriteProximityFinder::addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance )
{
#ifdef DEBUG
        //////std::cout<< "AddUpOpposingRaysAndFindShortest" << std::endl;
#endif
        float min_distance_1 = 100000000, min_distance_2 = 100000000, min_distance_3 = 100000000;
        float min_distance = 0;
        std::vector<float> diameters;
        int index=0, index2=0, index3=0, end = vectors.size(), half = end/2;
        bool dist1set = false, dist2set= false, dist3set= false;

        for(int i = 0; i < half; i++)
                diameters.push_back(vectors[i]->magnitude + vectors[half+i]->magnitude);

        end = diameters.size();
        
        for(int j = 0; j < end; j++)
        {
                if(diameters[j]<min_distance_1 && dist1set == false)
                {
                        min_distance_1 = diameters[j];
                        index = j;
                        dist1set = true;
                }
                if(diameters[j]<min_distance && j != index && dist1set==true&& dist2set==false){
                    min_distance_2 = diameters[j];
                    index2 = j;
                    dist2set = true;
                    
                }
                if(diameters[j]<min_distance&& j != index2 && dist1set==true&& dist2set==true){
                    min_distance_3 = diameters[j];
                    index3 = j;
                    dist3set = true;
                    
                }
                
        }
        
        VECTOR * tmp_vect = new VECTOR;
        tmp_vect->coords[X_COORD] = ( vectors[index]->coords[X_COORD] + vectors[half+index]->coords[X_COORD] ) / 2;
        tmp_vect->coords[Y_COORD] = ( vectors[index]->coords[Y_COORD] + vectors[half+index]->coords[Y_COORD] ) / 2;
        *adjustment_vect = tmp_vect;
        min_distance = (min_distance_1 + min_distance_2 + min_distance_3)/3;
        *distance = min_distance;
        
        
#ifdef DEBUG
        ////////std::cout<< "min diam: " << min_distance << std::endl;
#endif
        
        return index;

};
float AxonDendriteProximityFinder::setAverageThreshold(int numOfPointsToAverage)
{

        //////std::cout<< "In SetAverageThreshold " << std::endl;
        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
        unsigned int numOfEdges = edges->size();                
        

        if(numOfEdges == 0)
        {
                //////std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                return -1;
        }
        else
        {
// #pragma omp parallel for schedule(dynamic,1)
                for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                {                       
                        int pointCount = 1;
                        float threshold_sum = 0, current_threshold=0;
                        Edge * currentEdge = edges->at(pos);
                        std::list< double * >::iterator primary_edge_it, secondary_edge_it;
                        std::queue<float> threshold_q;
                        bool isFirstLoop = true;
                        int minimum_threshold = 11;
                        
                        for(secondary_edge_it = currentEdge->edgePointCoordinates.begin(); secondary_edge_it != currentEdge->edgePointCoordinates.end();secondary_edge_it++)
                        {
                          double * secondaryPoint = *secondary_edge_it;
                          
                          if(isInSectionBoundingBox(secondaryPoint)){
                          
                          secondaryPoint[THRESHOLD] = minimum_threshold;
                          
                          ////////std::cout<<"point threshold : "<<secondaryPoint[THRESHOLD]<<std::endl;
                          
                          }
                          
                        }
                        
                        

                        //for every point along edge
                        for(primary_edge_it = currentEdge->edgePointCoordinates.begin(), secondary_edge_it = primary_edge_it;
                            primary_edge_it != currentEdge->edgePointCoordinates.end(); 
                            primary_edge_it++) 
                        {
                                double * primaryPoint = *primary_edge_it;
                                
                                if(isInSectionBoundingBox(primaryPoint)){
                                
                                int x_pos = rint( primaryPoint[X_COORD] / XYSAMPLING );
                                int y_pos = rint( primaryPoint[Y_COORD] / XYSAMPLING );
                                int z_pos = rint( primaryPoint[Z_COORD] / ZSAMPLING );
                                
//                              //////std::cout<< "Position "<< primaryPoint[X_COORD] << "  " << primaryPoint[Y_COORD]<< "  "<< primaryPoint[Z_COORD] << std::endl;
//                              //////std::cout<< "Position "<< x_pos << "  " << y_pos<< "  "<< z_pos << std::endl;
                                

                                if(pointCount < numOfPointsToAverage)
                                {
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  threshold_sum += current_threshold;
                                }
                                else if(isFirstLoop)
                                {
                                  isFirstLoop = false;
                                  double * secondaryPoint = *secondary_edge_it;
                                  
                                  if(isInSectionBoundingBox(secondaryPoint)){
                                  
                                  
                                  
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  
                                  threshold_sum += current_threshold;
                                  
                                  ////////std::cout<<"threshold sum first loop: "<<threshold_sum<<std::endl;
                                  
                                  float average_threshold = threshold_sum/numOfPointsToAverage;
                                  
                                  for(int i=0; i<=((numOfPointsToAverage/2)+1); i++, secondary_edge_it++)
                                  {
                                      
                                    secondaryPoint = *secondary_edge_it;
                                    if(isInSectionBoundingBox(secondaryPoint)){
                                    
                                    if(average_threshold > minimum_threshold)
                                      secondaryPoint[THRESHOLD] = average_threshold;
                                    else
                                      secondaryPoint[THRESHOLD] = minimum_threshold;
                                      }
                                  }
                                  }//end if isinsection secondarypoint
                                }
                                else
                                {
                                  double * secondaryPoint = *secondary_edge_it;
                                  
                                  if(isInSectionBoundingBox(secondaryPoint)){
                                  
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  threshold_sum -= threshold_q.front();
                                  threshold_q.pop();
                                  threshold_sum += current_threshold;
                                  
                                  ////////std::cout<<"threshold sum other loop: "<<threshold_sum<<std::endl;
                                  
                                  float average_threshold = threshold_sum/numOfPointsToAverage;
                                  
                                  if(average_threshold > minimum_threshold)
                                    secondaryPoint[THRESHOLD] = average_threshold;
                                  else
                                    secondaryPoint[THRESHOLD] = minimum_threshold;
                                  
                                } //end if isinsection secondary point
                                }
                                
                                pointCount++;
                                
                        } //end if isinsection
                        }
                        
                        float average_threshold = threshold_sum/numOfPointsToAverage;
                        secondary_edge_it++;
                        
                        for(;secondary_edge_it != primary_edge_it;secondary_edge_it++)
                        {
                          double * secondaryPoint = *secondary_edge_it;
                          
                          if(isInSectionBoundingBox(secondaryPoint)){
                          
                          if(average_threshold > minimum_threshold)
                            secondaryPoint[THRESHOLD] = average_threshold;
                          else
                            secondaryPoint[THRESHOLD] = minimum_threshold;
                          
                          ////////std::cout<<"point threshold : "<<secondaryPoint[THRESHOLD]<<std::endl;
                          
                          }
                          
                        }
                        
                }
        }       

        return 0;
};

// float AxonDendriteProximityFinder::setAverageThreshold()
// {
// 
//         //////std::cout<< "In Set Average threshold " << std::endl;
//         
//         std::vector< Edge * > * edges = amira_graph->edgesPointer();
//         
//         unsigned int numOfEdges = edges->size(); 
//         
//         int pointCount = 1;
//                         float threshold_sum = 0, current_threshold=0;
//                         
//                         std::list< double * >::iterator primary_edge_it, secondary_edge_it;
//                         std::queue<float> threshold_q;
//                         bool isFirstLoop = true;
//                         int minimum_threshold = 11;
//                         float average_threshold;
//         
// 
//         if(numOfEdges == 0)
//         {
//                 //////std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
//                 return -1;
//         }
//         else
//         {
//             for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
//                 {                       
//                         Edge * currentEdge = edges->at(pos);
//                         
//                         for(primary_edge_it = currentEdge->edgePointCoordinates.begin(); primary_edge_it != currentEdge->edgePointCoordinates.end(); primary_edge_it++) 
//                         {
//                                 double * primaryPoint = *primary_edge_it;
//                                 if(isInSectionBoundingBox(primaryPoint)){
//                                 
//                                 int x_pos = rint( primaryPoint[X_COORD] / XYSAMPLING );
//                                 int y_pos = rint( primaryPoint[Y_COORD] / XYSAMPLING );
//                                 int z_pos = rint( primaryPoint[Z_COORD] / ZSAMPLING );
//                                 
//                                 current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);    
//                                 //////std::cout<<"current_threshold : "<<current_threshold<<std::endl;
//                                 //threshold_q.push(current_threshold);
//                                 threshold_sum += current_threshold*0.00001;
//                                 //////std::cout<<"threshold sum : "<<threshold_sum<<std::endl;
//                                 
//                                 average_threshold = threshold_sum/(pointCount*0.00001);
//                                 //////std::cout<<average_threshold<<std::endl;
//                                 pointCount++;
//                                 }
//                         }
//                         
//                         
//                 }
//                 
//                 //////std::cout<<"average threshold : "<<average_threshold<<std::endl;
//                 
//                 
//                 for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
//                 {
//                 Edge * currentEdge = edges->at(pos);    
//                     
//                 for(secondary_edge_it = currentEdge->edgePointCoordinates.begin(); secondary_edge_it != currentEdge->edgePointCoordinates.end(); secondary_edge_it++) 
//                         {
//                             double * secondaryPoint;    
//                             secondaryPoint = *secondary_edge_it;
//                             if(isInSectionBoundingBox(secondaryPoint)){
//                                 if(average_threshold > minimum_threshold)
//                                 secondaryPoint[THRESHOLD] = average_threshold;
//                                 else
//                                 secondaryPoint[THRESHOLD] = minimum_threshold;
//                                 
//                                 ////////std::cout<<"coords threshold: "<<secondaryPoint[THRESHOLD]<<std::endl;
//                             }
//                         }
//                 }
// // // #pragma omp parallel for schedule(dynamic,1)
// //                 for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
// //                 {                       
// //                         int pointCount = 1;
// //                         float threshold_sum = 0, current_threshold=0;
// //                         Edge * currentEdge = edges->at(pos);
// //                         std::list< double * >::iterator primary_edge_it, secondary_edge_it;
// //                         std::queue<float> threshold_q;
// //                         bool isFirstLoop = true;
// //                         int minimum_threshold = 11;
// //                         
// // 
// //                         //for every point along edge
// //                         for(primary_edge_it = currentEdge->edgePointCoordinates.begin(), secondary_edge_it = primary_edge_it;
// //                             primary_edge_it != currentEdge->edgePointCoordinates.end(); 
// //                             primary_edge_it++) 
// //                         {
// //                                 double * primaryPoint = *primary_edge_it;
// //                                 
// //                                 if(isInSectionBoundingBox(primaryPoint)){
// //                                 
// //                                 int x_pos = rint( primaryPoint[X_COORD] / XYSAMPLING );
// //                                 int y_pos = rint( primaryPoint[Y_COORD] / XYSAMPLING );
// //                                 int z_pos = rint( primaryPoint[Z_COORD] / ZSAMPLING );
// //                                 
// // //                              //////std::cout<< "Position "<< primaryPoint[X_COORD] << "  " << primaryPoint[Y_COORD]<< "  "<< primaryPoint[Z_COORD] << std::endl;
// // //                              //////std::cout<< "Position "<< x_pos << "  " << y_pos<< "  "<< z_pos << std::endl;
// //                                 
// // 
// //                                 if(pointCount < numOfPointsToAverage)
// //                                 {
// //                                   current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);                                  
// //                                   threshold_q.push(current_threshold);
// //                                   threshold_sum += current_threshold;
// //                                   //////std::cout<<"threshold sum not firstloop: "<<threshold_sum<<std::endl;
// //                                 }
// //                                 else if(isFirstLoop)
// //                                 {
// //                                   isFirstLoop = false;
// //                                   double * secondaryPoint;
// //                                   
// //                                   if(isInSectionBoundingBox(secondaryPoint)){
// //                                   
// //                                   
// //                                   
// //                                   current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
// //                                   threshold_q.push(current_threshold);
// //                                   threshold_sum += current_threshold;
// //                                   //////std::cout<<"threshold sum fist loop: "<<threshold_sum<<std::endl;
// //                                   
// //                                   float average_threshold = threshold_sum/numOfPointsToAverage;
// //                                   //////std::cout<<"average threshold first loop: "<<average_threshold<<std::endl;
// //                                   
// //                                   for(int i=0; i<=((numOfPointsToAverage/2)+1); i++, secondary_edge_it++)
// //                                   {
// //                                       
// //                                     secondaryPoint = *secondary_edge_it;
// //                                     if(isInSectionBoundingBox(secondaryPoint)){
// //                                     
// //                                     if(average_threshold > minimum_threshold)
// //                                       secondaryPoint[THRESHOLD] = average_threshold;
// //                                     else
// //                                       secondaryPoint[THRESHOLD] = minimum_threshold;
// //                                       }
// //                                   }
// //                                   }//end if isinsection secondarypoint
// //                                 }
// //                                 else
// //                                 {
// //                                   double * secondaryPoint = *secondary_edge_it;
// //                                   
// //                                   if(isInSectionBoundingBox(secondaryPoint)){
// //                                   
// //                                   current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
// //                                   threshold_q.push(current_threshold);
// //                                   threshold_sum -= threshold_q.front();
// //                                   threshold_q.pop();
// //                                   threshold_sum += current_threshold;
// //                                   //////std::cout<<"threshold sum else: "<<threshold_sum<<std::endl;
// //                                   
// //                                   float average_threshold = threshold_sum/numOfPointsToAverage;
// //                                   
// //                                   if(average_threshold > minimum_threshold)
// //                                     secondaryPoint[THRESHOLD] = average_threshold;
// //                                   else
// //                                     secondaryPoint[THRESHOLD] = minimum_threshold;
// //                                   
// //                                 } //end if isinsection secondary point
// //                                 }
// //                                 
// //                                 pointCount++;
// //                                 //////std::cout<<"pointcount: "<<pointCount<<std::endl;
// //                                 
// //                         } //end if isinsection
// //                         }
// //                         
// //                         float average_threshold = threshold_sum/numOfPointsToAverage;
// //                         secondary_edge_it++;
// //                         
// //                         for(;secondary_edge_it != primary_edge_it;secondary_edge_it++)
// //                         {
// //                           double * secondaryPoint = *secondary_edge_it;
// //                           
// //                           if(isInSectionBoundingBox(secondaryPoint)){
// //                           
// //                           if(average_threshold > minimum_threshold)
// //                             secondaryPoint[THRESHOLD] = average_threshold;
// //                           else
// //                             secondaryPoint[THRESHOLD] = minimum_threshold;
// //                           
// //                           }
// //                           
// //                         }
// //                         
// //                 }
//         }       
// 
//         return 0;
// };




// ImageType::Pointer ProximityFinder::getOriginalImage(char* file_name, int start_index, int end_index)
// {
// 
//         char input_file[1024];
//         strcpy(input_file, file_name);
// 
// 
//         NameGeneratorType::Pointer name_gen = NameGeneratorType::New();         //DECLARE AND INITIALIZE NAME_GENERATOR
//         name_gen->SetSeriesFormat( input_file );
//         name_gen->SetStartIndex( start_index );
//         name_gen->SetEndIndex( end_index );
//         name_gen->SetIncrementIndex( 1 );
//         
//                         
//         SeriesReaderType::Pointer input_reader = SeriesReaderType::New();       //DECLARE AND INITIALIZE INPUT_READER
//         input_reader->SetImageIO( itk::PNGImageIO::New() );
//         input_reader->SetFileNames( name_gen->GetFileNames() );
// 
//         try
//         {
//                 input_reader->Update();
//         }
//         catch( itk::ExceptionObject & err )
//         {
//                 std::cerr << "ImageReaderExceptionObject caught !" << std::endl;
//                 std::cerr << err << std::endl;
//         }
// 
// 
//         ImageType::Pointer original_image;
// 
//         original_image = input_reader->GetOutput();
//         original_image->Update();
//         
//         return original_image;
// };


void AxonDendriteProximityFinder::smoothRadii()
{
#ifdef DEBUG
        //////std::cout<< "SmoothRadii!! " << std::endl;
#endif        

                
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();
        double length=0;
        
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
                length += currentEdge->physicalLength;
                

                //for every point along edge
                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                {

                        double * current_coords = *edge_it;
                        if(isInSectionBoundingBox(current_coords)){
                        float diameter = current_coords[SURFACE];
                        long counter = 1;
                        
                        std::list< double * >::iterator next_it = edge_it;              
                        std::list< double * >::iterator prev_it = edge_it;
                        ++next_it;
                        --prev_it;
                        
                        if(edge_it == currentEdge->edgePointCoordinates.begin())
                        {
                                double * next_coords = *next_it;
                                
                                
                                if(isWithinEuclideanRange(next_coords, current_coords, 5) && isInSectionBoundingBox(next_coords))
                                {
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                }
                                
                        }
                        else if(next_it == currentEdge->edgePointCoordinates.end())
                        {
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5) && isInSectionBoundingBox(prev_coords))
                                {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                }
                        }
                        else
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5) && isInSectionBoundingBox(prev_coords))
                                {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                }
                                
                                if(isWithinEuclideanRange(next_coords, current_coords, 5) && isInSectionBoundingBox(next_coords))
                                {
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                }       
                        }
                                
                        diameter = 0.6666667*(diameter/counter);
                        current_coords[SURFACE] = diameter;
//                      //////std::cout<< "Final Diam: " << diameter << std::endl;
                        }//end if
                }

        }
//        total_axon_length = length; //set the physical length
};

bool AxonDendriteProximityFinder::isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance)
{
#ifdef DEBUG
  //////std::cout<<"in isWithinEuclideanRange!"<<std::endl; 
#endif  
  float temp_x = 0, temp_y = 0, temp_z = 0, distance = 0;

  temp_x = (float)XYSAMPLING *(( tmp[X_COORD] - neighbor[X_COORD] ) * ( tmp[X_COORD] - neighbor[X_COORD] ));
  temp_y = (float)XYSAMPLING *(( tmp[Y_COORD] - neighbor[Y_COORD] ) * ( tmp[Y_COORD] - neighbor[Y_COORD] ));
  temp_z = (float)ZSAMPLING  *(( tmp[Z_COORD] - neighbor[Z_COORD] ) * ( tmp[Z_COORD] - neighbor[Z_COORD] ));

  distance = temp_x + temp_y + temp_z;

  if(distance < ldistance * ldistance)
    return true;
  else
    return false;
};




void AxonDendriteProximityFinder::mergreProximites(std::vector<Proximity *> *input_prox_list)
{
    
    ////std::cout<<"in merge"<<std::endl;
    std::vector<Proximity *> cleaned_prox_list;
    bool merge_proximity = false;
    bool proximities_are_close = true;
    int counter = 0;
    bool hasBouton;
    double * prox_box = new double[6];
    double * other_prox_box = new double[6];
    
    std::vector< Proximity * >::iterator proxitr;
    
    for(proxitr = input_prox_list->begin(); proxitr != input_prox_list->end(); ++proxitr)
    {
        ////std::cout<< ((*proxitr)->getCenterCoords())[0]<<" "<< ((*proxitr)->getCenterCoords())[1]<<" "<<((*proxitr)->getCenterCoords())[2]<<std::endl;
    }
    
    // the max distance between two adjecent points, beyond which they will be seperate prox zones
    // this the same as the diagonal of the prox box to avoid the overlap between prox zones and hence 
    // duplication of bouton and overlap stats
    const double inter_proxzone_distance = MAX_PROX_BOX_SIZE * sqrt(3);
    
    ////std::cout<<" "<<inter_proxzone_distance<<" "<<std::endl;
    //while(proximities_are_close)
    std::vector<Proximity *>::iterator first_it;
    std::vector<Proximity *>::iterator second_it;
    
        for(first_it = input_prox_list->begin(); first_it != input_prox_list->end(); first_it++)
        {    
            ////////std::cout<<"i ="<<i<<std::endl;
           
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = *first_it;//input_prox_list->at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(second_it = first_it; second_it != input_prox_list->end(); )
            {   
                //////std::cout<<"inside for loop "<<std::endl;
                ////////std::cout<<"j ="<<j<<std::endl;
                if(second_it != first_it)   
                { 
                    Proximity * other_prox = *second_it;//input_prox_list->at(j);
                    
                    double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                    
                    double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                    
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                    
                    if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                    {    
                        //////std::cout<< "dist = "<<distance_between_points<<std::endl;
                        point_list->push_back(center_coords);
                        point_list->push_back(other_center_coords);
                        merge_proximity = true;
                        hasBouton = prox->getBouton();
                        other_prox->getBoundingBoxCoords(other_prox_box);
                        ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                        prox->getBoundingBoxCoords(prox_box);
                        ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                        if(hasBouton == false)
                        hasBouton = other_prox->getBouton();  
                        double confidence_value_2 = other_prox->getConfidenceValue();
                        if(confidence_value_2 > confidence_value)
                            confidence_value = confidence_value_2;
                        
                        //////std::cout<<"before erase "<<std::endl;
                        //////std::cout<<"size = "<<input_prox_list->size()<<std::endl;
                        //////std::cout<<"removed = "<<other_center_coords[0]<<" "<<other_center_coords[1]<<" "<<other_center_coords[2]<<std::endl;
                        
                        second_it = input_prox_list->erase(second_it);
                        
                        //////std::cout<<"after erase "<<std::endl;
                        //////std::cout<<"size = "<<input_prox_list->size()<<std::endl;
                        
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                        double * min_coords = new double[3];
                        double * max_coords = new double[3];
                        double * other_min_coords = new double[3];
                        double * other_max_coords = new double[3];
                        
                        min_coords[0] = prox_box[3];
                        min_coords[1] = prox_box[4];
                        min_coords[2] = prox_box[5];
                        
                        max_coords[0] = prox_box[0];
                        max_coords[1] = prox_box[1];
                        max_coords[2] = prox_box[2];
                        
                        other_max_coords[0] = other_prox_box[0];
                        other_max_coords[1] = other_prox_box[1];
                        other_max_coords[2] = other_prox_box[2];
                        
                        other_min_coords[0] = other_prox_box[3];
                        other_min_coords[1] = other_prox_box[4];
                        other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                            point_list->push_back(other_min_coords);
                            point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                            point_list->push_back(min_coords);
                            point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                        
                        counter ++;
                        
                        {
                            Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
    /*                      Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                            mergedProximity->setConfidenceValue(confidence_value);
                            mergedProximity->setBouton(hasBouton);
                            cleaned_prox_list.push_back(mergedProximity);
            //                    isalreadymerged = true;
                            *first_it = mergedProximity;
                            
                            center_coords = mergedProximity->getCenterCoords();
                        }
                                                    
                    }
                    else
                    {
                       second_it++;   
                    }
                        
                                        
                }
                else
                {
                    second_it++;  
                }
                    
            }
            
            
            
            
//             if(point_list->size() > 1 && merge_proximity == true){
//                 Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
// /*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
//                 mergedProximity->setConfidenceValue(confidence_value);
//                 mergedProximity->setBouton(hasBouton);
//                 cleaned_prox_list.push_back(mergedProximity);
// //                    isalreadymerged = true;
//                 input_prox_list.at(i) = mergedProximity;
// //                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
//                 
//                 
//                 
//                 
//                 }
                
//                 if(merge_proximity == false){
//                 cleaned_prox_list.push_back(prox); 
//                 }
            
        } 
        
        ////std::cout<<"after merge"<<std::endl;
        for(proxitr = input_prox_list->begin(); proxitr != input_prox_list->end(); ++proxitr)
        {
                ////std::cout<< ((*proxitr)->getCenterCoords())[0]<<" "<< ((*proxitr)->getCenterCoords())[1]<<" "<<((*proxitr)->getCenterCoords())[2]<<std::endl;
        }
    
};


void AxonDendriteProximityFinder::removeDoubledProximites()
{
    // this function merges proximites that are too close and should therefore be one single proximity
    
#ifdef DEBUG    
    //////std::cout<<"in removeDoubledProximities"<<std::endl;
#endif
    
    std::vector<Proximity *> cleaned_prox_list;
    std::vector<Proximity *> cleaned_prox_list_1_2;
    std::vector<Proximity *> cleaned_prox_list_2_1;
    std::vector<Proximity *> cleaned_prox_list_1_3;
    std::vector<Proximity *> cleaned_prox_list_3_1;
    std::vector<Proximity *> cleaned_prox_list_2_3;
    std::vector<Proximity *> cleaned_prox_list_3_2;
    bool merge_proximity = false;
    bool proximities_are_close = true;
    int counter = 0;
    bool hasBouton;
    double * prox_box = new double[6];
    double * other_prox_box = new double[6];
    
    // the max distance between two adjecent points, beyond which they will be seperate prox zones
    // this the same as the diagonal of the prox box to avoid the overlap between prox zones and hence 
    // duplication of bouton and overlap stats
    const double inter_proxzone_distance = this->max_dist;//MAX_PROX_BOX_SIZE * sqrt(3);
    
    //////std::cout<<" "<<inter_proxzone_distance<<" "<<std::endl;
    while(proximities_are_close)
    {
        int counter = 0;
        for(int i = 0; i < this->proximity_list.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                        {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
			    hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
			    if(hasBouton == false)
			    hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list.erase(proximity_list.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
		    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list.push_back(prox); 
                    }
               
                }
                
                
                
                
                
       for(int i = 0; i < this->proximity_list_1_2.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_1_2.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_1_2.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_1_2.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_1_2.erase(proximity_list_1_2.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_1_2.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_1_2.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_1_2.push_back(prox); 
                    }
               
                }
                
                
                for(int i = 0; i < this->proximity_list_2_1.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_2_1.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_2_1.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_2_1.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_2_1.erase(proximity_list_2_1.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_2_1.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_2_1.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_2_1.push_back(prox); 
                    }
               
                }
                
                
                
                
                
                
                
                
                  for(int i = 0; i < this->proximity_list_1_3.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_1_3.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_1_3.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_1_3.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_1_3.erase(proximity_list_1_3.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_1_3.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_1_3.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_1_3.push_back(prox); 
                    }
               
                }  
                
                
                
                
                
                              for(int i = 0; i < this->proximity_list_3_1.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_3_1.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_3_1.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_3_1.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_3_1.erase(proximity_list_3_1.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_1_3.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_3_1.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_3_1.push_back(prox); 
                    }
               
                }
                
                
                
                
                
                        for(int i = 0; i < this->proximity_list_2_3.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_2_3.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_2_3.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_2_3.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_2_3.erase(proximity_list_2_3.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_2_3.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_2_3.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_2_3.push_back(prox); 
                    }
               
                }      
                
                
                
                
                
                              for(int i = 0; i < this->proximity_list_3_2.size(); i++)
        {    
           
            
                
            std::vector<double *> * point_list = new std::vector<double *>;

            Proximity * prox = this->proximity_list_3_2.at(i);
            double confidence_value = prox->getConfidenceValue();
            int * IDs = new int[1];
            IDs[0] = prox->Axon_ID;
            IDs[1] = prox->Dendrite_ID;
            
            
         
            double * center_coords = prox->getCenterCoords();

//          //////std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
          
            for(int j = i; j < this->proximity_list_3_2.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        Proximity * other_prox = this->proximity_list_3_2.at(j);
                        
                        

                        double * other_center_coords = other_prox->getCenterCoords();
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords<<std::endl;
//                       //////std::cout<<"Previous center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
                        
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                                              
                        
//                        //////std::cout<<"distance    :  "<<distance_between_points<<std::endl;
//                        ImageType::Pointer image = prox.getProximityImage();
//                        image->Print(//////std::cout);                       
                        
                        if(distance_between_points < inter_proxzone_distance && prox->Axon_ID == other_prox->Axon_ID && prox->Dendrite_ID == other_prox->Dendrite_ID)
                            {                                
                            point_list->push_back(center_coords);
                            point_list->push_back(other_center_coords);
                            merge_proximity = true;
                            hasBouton = prox->getBouton();
                            other_prox->getBoundingBoxCoords(other_prox_box);
                            ////////std::cout<<other_prox_box[0]<<"  "<<other_prox_box[1]<<"  "<<other_prox_box[2]<<"  "<<other_prox_box[3]<<"  "<<other_prox_box[4]<<"  "<<other_prox_box[5]<<"  "<<std::endl;
                            prox->getBoundingBoxCoords(prox_box);
                            ////////std::cout<<prox_box[0]<<"  "<<prox_box[1]<<"  "<<prox_box[2]<<"  "<<prox_box[3]<<"  "<<prox_box[4]<<"  "<<prox_box[5]<<"  "<<std::endl;
                            if(hasBouton == false)
                            hasBouton = other_prox->getBouton();  
                            double confidence_value_2 = other_prox->getConfidenceValue();
                            if(confidence_value_2 > confidence_value)
                                confidence_value = confidence_value_2;
                            proximity_list_1_3.erase(proximity_list_3_2.begin()+j);
//                             proximity_list_nondecon.erase(proximity_list_nondecon.begin()+j);
                            double * min_coords = new double[3];
                            double * max_coords = new double[3];
                            double * other_min_coords = new double[3];
                            double * other_max_coords = new double[3];
                            
                            min_coords[0] = prox_box[3];
                            min_coords[1] = prox_box[4];
                            min_coords[2] = prox_box[5];
                            
                            max_coords[0] = prox_box[0];
                            max_coords[1] = prox_box[1];
                            max_coords[2] = prox_box[2];
                            
                            other_max_coords[0] = other_prox_box[0];
                            other_max_coords[1] = other_prox_box[1];
                            other_max_coords[2] = other_prox_box[2];
                            
                            other_min_coords[0] = other_prox_box[3];
                            other_min_coords[1] = other_prox_box[4];
                            other_min_coords[2] = other_prox_box[5];
//                             for(int k = 0; k <6; k++){
// //                                 double * other_coords = new double[3];
// //                                 other_coords[0] = other_prox_box[k]
                                point_list->push_back(other_min_coords);
                                point_list->push_back(other_max_coords);
// //                                 
// //                                 //////std::cout<<other_prox_box[k]<<std::endl;
// //                                 double * coords = prox_box[k];
                                point_list->push_back(min_coords);
                                point_list->push_back(max_coords);
// //                                 //////std::cout<<prox_box[k]<<std::endl;
//                             }
                            
                            counter ++;
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list->size() > 1 && merge_proximity == true){
                    Proximity * mergedProximity = new Proximity(point_list, this->original_image, IDs);
/*                    Proximity * mergedProximity_nondecon = new Proximity(point_list, this->original_image_nondecon, IDs);*/                    
                    mergedProximity->setConfidenceValue(confidence_value);
                    mergedProximity->setBouton(hasBouton);
                    cleaned_prox_list_3_2.push_back(mergedProximity);
//                    isalreadymerged = true;
                    proximity_list_3_2.at(i) = mergedProximity;
//                     proximity_list_nondecon.at(i) = mergedProximity_nondecon;
                    
                    
                    
                    
                    }
                    
                    if(merge_proximity == false){
                    cleaned_prox_list_3_2.push_back(prox); 
                    }
               
                }
                
                
                
                
                
                
                
       if(counter == 0)
           proximities_are_close = false;
    }
                 
//     for(int k = 0; k< proximity_list.size(); k++){
//         double * prox_coords = proximity_list.at(k)->getCenterCoords();
//         if prox_coords[z]
//         
//     }
//           proximity_list.clear();
// 
//             for(int i = 0; i < cleaned_prox_list.size(); i++){
//             this->proximity_list.push_back(cleaned_prox_list.at(i));              
//           }
          
};

bool AxonDendriteProximityFinder::readAmiraTransformations()

{
    //reads Amira transformations from spatial graph file in order to transform the merged graph into it's original position in to the current section's image planes
  
  
  
//  std::ifstream inputStream1(offset);
//  double manual_z_scale;
//     double z_offset;
//     
//  if(!inputStream1.fail()){
// //  
//     
//     const char * numbers = "0123456789";
//     const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";  
//     const char * whitespace = "\n ";
//     const char * signs = "+-";
//     
//     std::string offset_parameters;
//     std::string point = ".";
//     
// while(!std::getline(inputStream1, offset_parameters).eof() /*&& line < 100*/)
//                 {    
//     
// if(true){
//     std::string::size_type loc1, loc2, loc5;
//     loc5= offset_parameters.find("z scale", 0);
//         if(loc5 ==0){
//         
//         loc1 = offset_parameters.find_first_of(numbers, 0);
//         //////std::cout<<"loc1 "<<loc1<<std::endl;
//         
//         loc2 = offset_parameters.find_first_of(whitespace, loc1);
//         //////std::cout<<"loc2 "<<loc2<<std::endl;
//         if(loc2 != std::string::npos)
//         if(loc2 < loc1)
//         loc1 = loc2;
// //         if(offset_parameters.find_first_of(numbers, loc2)){
// //            loc3 = offset_parameters.find("_Trans.am", loc2);
// //            //////std::cout<<"loc3 "<<loc3<<std::endl;
//            
//            char * tmp1 = new char[loc2 - loc1];
//            
//            offset_parameters.copy(tmp1, loc2- loc1 , loc1);
// //             for(int i = 0; i< sizeof(tmp1) -1; i++)
//                //////std::cout<<tmp1<<std::endl;
// 
// //            char * tmp2 = new char[loc3 - loc2];         
// //            
// //             offset_parameters.copy(tmp2, loc3 - loc2 -1, loc2 +1);
// //             for(int i = 0; i< sizeof(tmp2) -1; i++)
// //                //////std::cout<<tmp2[i]<<std::endl;
//             
//             
// //             std::string tmp1;
// //            std::string tmp2;
// 
// //            scale_name.copy(tmp2, loc3 - loc2, loc2);
// //            scale_name.copy(tmp1, loc2 - loc1, loc1);
//           std::string tmp3 = "";
//           tmp3 += tmp1;
//           //////std::cout<<tmp3<<std::endl;
// //           tmp3 += ".";
// //           tmp3 += tmp2;
//           //////std::cout<<tmp3<<std::endl;
//           const char * z_scale = new char[tmp3.size() +1];
//           z_scale = tmp3.c_str();
//           
// 
//             
// //            const char z_scale = new char [tmp3.size()];
// //            tmp3.copy(z_scale, tmp3.size(), 0);
//            
//            
//            
//            
// //            char * tmp2 = new char[loc3 - loc2];
// //             scale_name.copy(tmp2, loc3 - loc2, loc2);
// //             char * tmp1 = new char[loc2 - loc1];
// //             scale_name.copy(tmp1, loc2 - loc1, loc1);
// //             
//             double ftmp1 = atof(z_scale);
//             manual_z_scale = ftmp1;
//             //////std::cout<<"z_scale from string: "<<manual_z_scale<<std::endl;
//            } 
//            
//            
//          if(offset_parameters.find("z offset" , 0 ) /*!= std::string::npos*/){
//              std::string::size_type loc3, loc4, loc;
//              loc = offset_parameters.find("z offset", 0);
//          //////std::cout<<"loc "<<loc<<std::endl;
//             
//             loc3 = offset_parameters.find_first_of(numbers, loc);
//         loc4 = offset_parameters.find_first_of(signs, loc3);
//         if(loc4 != std::string::npos)
//                 if(loc4 < loc3)
//                         loc3 = loc4;
//         //////std::cout<<"loc3 "<<loc3<<std::endl;
//         loc4= offset_parameters.find_first_of(whitespace, loc3);
//         //////std::cout<<"loc4 "<<loc4<<std::endl;
//         if(loc4 != std::string::npos)
//         if(loc4 < loc3)
//         loc3 = loc4;
//         
//         char * tmp2 = new char[loc4- loc3];
//            
//         offset_parameters.copy(tmp2, loc4- loc3 , loc3);
//         
// //        //////std::cout<<tmp1<<std::endl;
//          /* const char * z_scale = new char[sizeof(tmp1) +1]*/;
//           z_offset = atof(tmp2);
//           //////std::cout<<"z_offset from string: "<<z_offset<<std::endl;
//         
//         
//             
//          } 
//         
//         
// //         loc = offset_parameters.find("z offset", 0);
// //         //////std::cout<<"loc "<<loc<<std::endl;
// //         if(true){
// //         std::string::size_type loc3, loc4;
//         
//           
//           
//           
//            
// 
//         //}
//         
// //        point.copy(tmp1, 1, 1);
//         
//         }
//         
//     }
// } 
// 
// inputStream1.close();
//     std::scanf(manual_z_scale, str);
//     std::printf(manual_z_scale);

// const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
//                 const char * numbers = "0123456789";
//                 const char * signs = "+-";
//                 const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
//                 const char * whitespace = "\t ";
//                 
//                 std::string currentLine;
//                 unsigned int line = 0;
//                 
//                  while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
//                 {
//                     std::string::size_type loc1, loc2, loc3;
//                     loc1 = currentLine.find("Zscaled, 0);
//                     loc2 = currentLine.find_first_of(signs, 0);
//                     char * tmp1 = new char[loc2 - loc1];
//                     currentLine.copy(tmp1, loc2 - loc1, loc1);
//                     double ftmp1 = atof(tmp1);
//                 }
//                 
//         }
//     //////std::cout<<transformFilename<<std::endl;
//     //////std::cout<<filename<<std::endl;
    

#ifdef DEBUG    
//////std::cout<<"in readAmiraTransformations"<<std::endl;
#endif
        //////std::cout<<"inputfilename: " << inputfilename <<std::endl;
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
//                                              //////std::cout << "found correct section transform parameters!" << std::endl;
                                                unsigned int count = 0;
						std::string::size_type loc1, loc2, loc3;
						loc1 = currentLine.find_first_of(numbers, 0);
						loc2 = currentLine.find_first_of(signs, 0);
						if(loc2 != std::string::npos)
							if(loc2 < loc1)
								loc1 = loc2;
						loc2 = currentLine.find_first_of(whitespace, loc1 + 1);	//ignores last value: is always 1 anyways
						while(loc2 != std::string::npos && count < 16)
						{
                                                        char * tmp1 = new char[20];
                                                
                                                        for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                        {
                                                            tmp1[i] = 'f';
                                                        }
							//char * tmp1 = new char[loc2 - loc1];
							currentLine.copy(tmp1, loc2 - loc1, loc1);
							double ftmp1 = atof(tmp1);
							sectionRotation[count%4][count/4]= ftmp1;	// amira files are columns after each other
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
                                                
//                                                 std::cout << "transformation matrix as read:" << std::endl;
//                                                 for(int ii = 0; ii < 4; ++ii)
//                                                 {
//                                                         std::cout << "[";
//                                                         for(int jj = 0; jj < 4; ++jj)
//                                                         {
//                                                                 if(jj < 3)
//                                                                         std::cout << sectionRotation[ii][jj] << ",\t";
//                                                                 else
//                                                                         std::cout << sectionRotation[ii][jj];
//                                                         }
//                                                         std::cout << "]" << std::endl;
//                                                         
//                                                 }
                                              
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
//                                                 //////std::cout<<"shift : "<<sectionTranslation[2][3]<<std::endl;
                                                
                                                //////std::cout << "translation matrix:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        //////std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                //if(jj < 3)
                                                                        //////std::cout << sectionTranslation[ii][jj] << ",\t";
                                                               // else;
                                                                        //////std::cout << sectionTranslation[ii][jj];
                                                        }
                                                        //////std::cout << "]" << std::endl;
                                                }
                                                
                                                //////std::cout << "rotation matrix with scaling:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        //////std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                //if(jj < 3)
                                                                        //////std::cout << sectionRotation[ii][jj] << ",\t";
                                                                //else;
                                                                        //////std::cout << sectionRotation[ii][jj];
                                                        }
                                                        //////std::cout << "]" << std::endl;
                                                }
                                                double square_scaling[2];
                                                double scaling[4];
                                                scaling[3] = 1;
                                                
                                                
                                                
                                                for(int jj = 0; jj< 3; ++jj){
                                                
                                                square_scaling[jj] = sectionRotation[0][jj]*sectionRotation[0][jj] + sectionRotation[1][jj]*sectionRotation[1][jj] + sectionRotation[2][jj]*sectionRotation[2][jj];
                                                scaling[jj] = sqrt(square_scaling[jj]);
                                                //////std::cout<<"scaling : "<<scaling[jj]<<std::endl;
                                                
                                                }
                                                //scaling[2] = manual_z_scale[0];
                                                ////////std::cout<<"scaling : "<<scaling[2]<<std::endl;
                                                
                                                for(int ii = 0; ii < 3; ++ii)
                                                {                                                        
                                                        for(int jj = 0; jj < 3; ++jj)                                                        
                                                            sectionRotation[ii][jj] = sectionRotation[ii][jj]/scaling[jj];                                                       
                                                } 
                                                
                                                
                                                //////std::cout << "rotation matrix without scaling:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        //////std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                //if(jj < 3)
                                                                        //////std::cout << sectionRotation[ii][jj] << ",\t";
                                                                //else;
                                                                        //////std::cout << sectionRotation[ii][jj];
                                                        }
                                                        //////std::cout << "]" << std::endl;
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
                                                                    
                                                                    ////////std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                                        mProduct[ii][jj] += sectionTranslation[ii][kk]*sectionRotation[kk][jj]*scaling[jj];
                                                                }
                                                                
                                                this->transformation = mProduct;
                                                
//                                                 std::cout << "transformation matrix:" << std::endl;
//                                                 for(int ii = 0; ii < 4; ++ii)
//                                                 {
//                                                         std::cout << "[";
//                                                         for(int jj = 0; jj < 4; ++jj)
//                                                         {
//                                                                 if(jj < 3)
//                                                                         std::cout << transformation[ii][jj] << ",\t";
//                                                                 else
//                                                                         std::cout << transformation[ii][jj];
//                                                         }
//                                                         std::cout << "]" << std::endl;
//                                                         
//                                                 }
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
                                                
                                                
                                                this->inverse_transformation = mInverse;
                                               
//                                                mInverse[0][3] = -1*(sectionTranslation[0][3]);
//                                                mInverse[1][3] = -1*(sectionTranslation[1][3]);
//                                                mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                mInverse[3][3] = 1;
                                               
                                               //this->inverse_transformation = mInverse;
                                               
                                               //////std::cout << "inverse transformation matrix:" << std::endl;
                                                for(int ii = 0; ii < 4; ++ii)
                                                {
                                                        //////std::cout << "[";
                                                        for(int jj = 0; jj < 4; ++jj)
                                                        {
                                                                //if(jj < 3)
                                                                        //////std::cout << inverse_transformation[ii][jj] << ",\t";
                                                                //else;
                                                                        //////std::cout << inverse_transformation[ii][jj];
                                                        }
                                                        //////std::cout << "]" << std::endl;
                                                        
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

bool AxonDendriteProximityFinder::manual_z_scaling(const char * offset, AmiraSpatialGraph * graph, bool inverse){
    
    // additional z-scaling and z-offset from a text file if a z-scale factor and/or offset is not saved in the transformation files
    std::ifstream inputStream1(offset);
double * manual_z_scale = new double[2];
double z_offset;


if(!inputStream1.fail()){
    
    const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t ";
                double ** transformation;
                
                transformation = new double *[4];
                for(int ii = 0; ii < 4; ++ii)
                {
                        transformation[ii] = new double[4];
                        for(int jj = 0; jj < 4; ++jj)
                        {
                                if(ii != jj)
                                        transformation[ii][jj] = 0;
                                else
                                        transformation[ii][jj] = 1;
                        }
                }
                
                
                
                
                std::string currentLine1;
                
                
                while(!std::getline(inputStream1, currentLine1).eof() /*&& line < 100*/)
                {
                    std::string::size_type loc1, loc2;
                    loc1 = currentLine1.find_first_of(whitespace, 0);
                    char * tmp1 = new char[loc1];
                    currentLine1.copy(tmp1, loc1, 0);
                    double ftmp1 = atof(tmp1);
                    manual_z_scale[0]= ftmp1;
                    //////std::cout<<"from file: "<<manual_z_scale[0]<<std::endl;
                    
                    loc2 = currentLine1.size();
                    char * tmp2 = new char[loc2];
                    currentLine1.copy(tmp2, loc2 - loc1, loc1);
                    double ftmp2 = atof(tmp2);
                    manual_z_scale[1]= ftmp2;
                    //////std::cout<<"from file: "<<manual_z_scale[1]<<std::endl;
        
                    
//                     if(currentLine1.size())
//                     if( != std::string::npos)
//                     unsigned int count = 0;
//                     std::string::size_type loc1, loc2, loc3;
//                     loc1 = currentLine1.find_first_of(numbers, 0);
//                     loc2 = currentLine1.find_first_of(signs, 0);
//                     if(loc2 != std::string::npos)
//                     if(loc2 < loc1)
//                     loc1 = loc2;
//                     loc2 = currentLine1.find_first_of(whitespace, loc1 );
//                     
//                     
//                     while(loc2 != std::string::npos && count < 2)
//                         {
//                                 char * tmp1 = new char[loc2 - loc1];
//                                 currentLine1.copy(tmp1, loc2 - loc1, loc1);
//                                 double ftmp1 = atof(tmp1);
//                                 manual_z_scale[count]= ftmp1;        // amira files are columns after each other
//                                 //////std::cout<<"from file: "<<manual_z_scale[count]<<std::endl;
//                                 loc3 = loc2;
//                                 loc1 = currentLine1.find_first_of(numbers, loc3);
//                                 loc2 = currentLine1.find_first_of(signs, loc3);
//                                 if(loc2 != std::string::npos)
//                                         if(loc2 < loc1)
//                                                 loc1 = loc2;
//                                 loc2 = currentLine1.find_first_of(whitespace, loc1 );
//                                 
//                                 delete [] tmp1;
//                                 count ++;
//                         }

                    
                }
                
                
                if(inverse){
                transformation[2][3] = manual_z_scale[1];                
                transformation[2][2] = manual_z_scale[0];
                }
                else{
                transformation[2][3] = -manual_z_scale[1];                
                transformation[2][2] = 1/manual_z_scale[0];
                    
                }
                graph->applyInverseTransformation(transformation);
                
                
    
    
    
}

inputStream1.close();


return 0;  
    
    
    
}

// double ** AxonDendriteProximityFinder::calculateTransAfterRotInverse()
// {
//         double ** mInverse = new double *[4];
//         for(int ii = 0; ii < 4; ++ii)
//         {
//                 mInverse[ii] = new double[4];
//                 for(int jj = 0; jj < 4; ++jj)
//                         mInverse[ii][jj] = 0;
//         }
//         for(int ii = 0; ii < 2; ++ii)
//                 for(int jj = 0; jj < 2; ++jj)
//                         mInverse[ii][jj] = sectionRotation[jj][ii];
//         mInverse[0][3] = -1*(sectionRotation[0][0]*sectionTranslation[0][3] + sectionRotation[1][0]*sectionTranslation[1][3]);
//         mInverse[1][3] = -1*(sectionRotation[0][1]*sectionTranslation[0][3] + sectionRotation[1][1]*sectionTranslation[1][3]);
//         mInverse[2][3] = -1*sectionTranslation[2][3];
//         mInverse[3][3] = 1;
//         
//         return mInverse;
// };




std::vector<double* > * AxonDendriteProximityFinder::transformCoordinates(std::vector<double *> * list)
{
    //retransformation of proximity coordinates from original position to their final position in the merged graph
    //////std::cout<<"in transformCoordinates"<<std::endl;    
    
    std::vector<double *> * new_list = list;
    
    
    
    double ** transform_matrix = this->transformation;

        
        for(int i =0;  i < list->size(); i++) 
                {

                                double tmpPoint[4];
                                double transPoint[4];
                                

                                        tmpPoint[0] = list->at(i)[X_COORD];
                                        tmpPoint[1] = list->at(i)[Y_COORD];
                                        tmpPoint[2] = list->at(i)[Z_COORD];
                                        tmpPoint[3] = 1;
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                transPoint[jj] = 0;
                                                for(int kk = 0; kk < 4; ++kk){
                                                        transPoint[jj] += transform_matrix[jj][kk]*tmpPoint[kk];
//                                                //////std::cout<<"Transform      "<<worldTransform[jj][kk]<<std::endl;
                                                }
                                        }
                                        
                                        
                                        new_list->at(i)[X_COORD] = transPoint[0];
                                        new_list->at(i)[Y_COORD] = transPoint[1];
                                        new_list->at(i)[Z_COORD] = transPoint[2];
//                                }
//                        }
//                }     
        }
        
        return new_list;
};

std::vector<double* > * AxonDendriteProximityFinder::transformCoordinatesZScale(std::vector<double *> * list, const char * offset)
{
    // additional z-scaling from a text file if a z-scale factor is not saved in the transformation files for retransformation
    //////std::cout<<"in transformCoordinates"<<std::endl;    
    
    std::vector<double *> * new_list = list;
        std::ifstream inputStream1(offset);
        
        double ** transformation;
        
        transformation = new double *[4];
                for(int ii = 0; ii < 4; ++ii)
                {
                        transformation[ii] = new double[4];
                        for(int jj = 0; jj < 4; ++jj)
                        {
                                if(ii != jj)
                                        transformation[ii][jj] = 0;
                                else
                                        transformation[ii][jj] = 1;
                        }
                }
                
        
        
double * manual_z_scale = new double[2];
double z_offset;


if(!inputStream1.fail()){
    
    const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t ";
                
                
                
                
//                 
                
                std::string currentLine1;
                
                
                while(!std::getline(inputStream1, currentLine1).eof() /*&& line < 100*/)
                {
                    std::string::size_type loc1, loc2;
                    loc1 = currentLine1.find_first_of(whitespace, 0);
                    char * tmp1 = new char[loc1];
                    currentLine1.copy(tmp1, loc1, 0);
                    double ftmp1 = atof(tmp1);
                    manual_z_scale[0]= ftmp1;
                    //////std::cout<<"from file: "<<manual_z_scale[0]<<std::endl;
                    
                    loc2 = currentLine1.size();
                    char * tmp2 = new char[loc2];
                    currentLine1.copy(tmp2, loc2 - loc1, loc1);
                    double ftmp2 = atof(tmp2);
                    manual_z_scale[1]= ftmp2;
                    //////std::cout<<"from file: "<<manual_z_scale[1]<<std::endl;
                    
                    transformation[2][3] = manual_z_scale[1];                
                    transformation[2][2] = manual_z_scale[0];
        
                    
//                     if(currentLine1.size())
//                     if( != std::string::npos)
//                     unsigned int count = 0;
//                     std::string::size_type loc1, loc2, loc3;
//                     loc1 = currentLine1.find_first_of(numbers, 0);
//                     loc2 = currentLine1.find_first_of(signs, 0);
//                     if(loc2 != std::string::npos)
//                     if(loc2 < loc1)
//                     loc1 = loc2;
//                     loc2 = currentLine1.find_first_of(whitespace, loc1 );
//                     
//                     
//                     while(loc2 != std::string::npos && count < 2)
//                         {
//                                 char * tmp1 = new char[loc2 - loc1];
//                                 currentLine1.copy(tmp1, loc2 - loc1, loc1);
//                                 double ftmp1 = atof(tmp1);
//                                 manual_z_scale[count]= ftmp1;        // amira files are columns after each other
//                                 //////std::cout<<"from file: "<<manual_z_scale[count]<<std::endl;
//                                 loc3 = loc2;
//                                 loc1 = currentLine1.find_first_of(numbers, loc3);
//                                 loc2 = currentLine1.find_first_of(signs, loc3);
//                                 if(loc2 != std::string::npos)
//                                         if(loc2 < loc1)
//                                                 loc1 = loc2;
//                                 loc2 = currentLine1.find_first_of(whitespace, loc1 );
//                                 
//                                 delete [] tmp1;
//                                 count ++;
//                         }

                    
                }
                
                
                
                
                
                
                
                
    
    
    
}

inputStream1.close();
    
    
    
    
    double ** transform_matrix = transformation;

        
        for(int i =0;  i < list->size(); i++) 
                {

                                double tmpPoint[4];
                                double transPoint[4];
                                

                                        tmpPoint[0] = list->at(i)[X_COORD];
                                        tmpPoint[1] = list->at(i)[Y_COORD];
                                        tmpPoint[2] = list->at(i)[Z_COORD];
                                        tmpPoint[3] = 1;
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                transPoint[jj] = 0;
                                                for(int kk = 0; kk < 4; ++kk){
                                                        transPoint[jj] += transform_matrix[jj][kk]*tmpPoint[kk];
//                                                //////std::cout<<"Transform      "<<worldTransform[jj][kk]<<std::endl;
                                                }
                                        }
                                        
                                        
                                        new_list->at(i)[X_COORD] = transPoint[0];
                                        new_list->at(i)[Y_COORD] = transPoint[1];
                                        new_list->at(i)[Z_COORD] = transPoint[2];
//                                }
//                        }
//                }     
        }
        
        return new_list;
};


// Edge * AxonDendriteProximityFinder::transformEdge(Edge * edge){
//   
//   double ** transform_matrix = this->inverse_transformation;
//   std::list< double * >::iterator edge_it;        
//   unsigned int edgeSize = edge->edgePointCoordinates.size();
//   for(edge_it = edge->edgePointCoordinates.begin();edge_it != edge->edgePointCoordinates.end(); edge_it++) 
//                 {
//                                         double * point = *edge_it;
//   
//   double tmpPoint[4];
//   double transPoint[4];
//                 
// 
//                         tmpPoint[0] = point[X_COORD];
//                         tmpPoint[1] = point[Y_COORD];
//                         tmpPoint[2] = point[Z_COORD];
//                         tmpPoint[3] = 1;
//                         for(int jj = 0; jj < 4; ++jj)
//                         {
//                                 transPoint[jj] = 0;
//                                 for(int kk = 0; kk < 4; ++kk){
//                                         transPoint[jj] += transform_matrix[jj][kk]*tmpPoint[kk];
// //                                                //////std::cout<<"Transform      "<<worldTransform[jj][kk]<<std::endl;
//                                 }
//                         }
//                         
//                         
//                         point[X_COORD] = transPoint[0];
//                         point[Y_COORD] = transPoint[1];
//                         point[Z_COORD] = transPoint[2];
//                                 }
//   return edge;
// }


double AxonDendriteProximityFinder::calculateGap(double * axonPt, double * dendritePt, ImageType::Pointer image, bool istest, int prox){
    
    ////////std::cout<<"In calculateGap"<<std::endl;

    float x_axon, y_axon, z_axon, x_dendrite, y_dendrite, z_dendrite;
    double axon_touch_x, axon_touch_y, dendrite_touch_x, dendrite_touch_y, min_x_dist, min_y_dist;
    Image2DType::Pointer image_plane_z;
    float z_focus;
    float threshold_axon = 0, threshold_dendrite = 0;
    double * axon_surface_point = new double[2];
    double * dendrite_surface_point = new double[2];
    
    std::vector<double> real_distance_x;
    std::vector<double> real_distance_y;
    std::vector<double*> * axon_surface_points = new std::vector<double *>;
    std::vector<double*> * dendrite_surface_points = new std::vector<double *>;
//    Image2DType::Pointer image_plane_3;
    
        
    
    
//    image_plane_3 = getImagePlane(z,3,image);
    
    
    x_axon = rint(axonPt[X_COORD]/XYSAMPLING);
    y_axon = rint(axonPt[Y_COORD]/XYSAMPLING);
    z_axon = rint(axonPt[Z_COORD]/ZSAMPLING);
    
//    //////std::cout<<x_axon<<"  "<<y_axon<<"  "<<z_axon<<std::endl;
    
    x_dendrite = rint(dendritePt[X_COORD]/XYSAMPLING);
    y_dendrite = rint(dendritePt[Y_COORD]/XYSAMPLING);
    z_dendrite = rint(dendritePt[Z_COORD]/ZSAMPLING);
    
//    //////std::cout<<x_dendrite<<"  "<<y_dendrite<<"  "<<z_dendrite<<std::endl;
    
    float z_dist = abs(z_axon - z_dendrite);
    if (z_axon > z_dendrite)
        z_focus = z_axon - z_dist/2;
    if (z_dendrite > z_axon)
        z_focus = z_dendrite - z_dist/2;
    
//     //////std::cout<<"z distance : "<<z_dist<<std::endl;
//     //////std::cout<<"z centre   : "<<z_focus<<std::endl;
    
    image_plane_z = getImagePlane(z_focus, z_dist, image);
    
  
    
    
    threshold_axon = axonPt[THRESHOLD];
    double diameter = axonPt[SURFACE];
//     //////std::cout<<"diameter : "<<diameter<<std::endl;
//     threshold_dendrite = dendritePt[THRESHOLD];
    
   // //////std::cout<<"axon threshold: "<<threshold_axon<<std::endl;
    ////////std::cout<<"dendrite threshold: "<<threshold_dendrite<<std::endl;
//                                                //////std::cout<<"threshold"<<coords[THRESHOLD]<<std::endl;
                                                
        unsigned int nr_of_rays = 30;
        float ray_length = 0.5;
        int index = 0;
        

        std::vector<VECTOR *> axon_vectors;
//        VECTOR * adjustment_vect;
        float distance = 0;

        x_axon = x_axon + 0.5;
        y_axon = y_axon + 0.5;

            for(unsigned n=1; n<=nr_of_rays; n++)
            {
            VECTOR * tmp_a;
//                                                        //////std::cout<<"try to sendRay"<<std::endl;
            tmp_a = sendRay(x_axon, y_axon, ray_length, nr_of_rays, threshold_axon, n, image_plane_z);
            ////////std::cout<<"tmp_a: "<<tmp_a->coords[X_COORD]<<"  "<<tmp_a->coords[Y_COORD]<<std::endl;
//                                                        //////std::cout<<"sent Ray"<<std::endl;
            axon_vectors.push_back(tmp_a);
            axon_surface_point[X_COORD] = tmp_a->coords[X_COORD] + double(x_axon);
            ////////std::cout<<"axon surface x "<<surface_point[X_COORD]<<std::endl;            
            axon_surface_point[Y_COORD] = tmp_a->coords[Y_COORD] + double(y_axon);
            axon_surface_point[Z_COORD] = double(z_axon);
            ////////std::cout<<"axon surface y "<<surface_point[Y_COORD]<<std::endl;
            double vector_length = tmp_a->magnitude;
            ////////std::cout<<"tmp_a length: "<<vector_length<<std::endl;
            
            axon_surface_points->push_back(axon_surface_point);
            delete tmp_a;
            
//                        //////std::cout<< "surface point push back success " << std::endl;
            }
            
//            for(int i = 0; i<axon_vectors.size(); i++){
//                //////std::cout<<"axon_vector_x"<<axon_vectors.at(i)->coords[X_COORD]<<"      "<<axon_vectors.at(i)->coords[Y_COORD]<<std::endl;}
            
        std::vector<VECTOR *> dendrite_vectors;
//         VECTOR * adjustment_vect;
//          float distance = 0;

        x_dendrite = x_dendrite + 0.5;
        y_dendrite = y_dendrite + 0.5;

            for(unsigned n=1; n<=nr_of_rays; n++)
            {
            VECTOR * tmp_d;
//                                                        //////std::cout<<"try to sendRay"<<std::endl;
            tmp_d = sendRay(x_dendrite, y_dendrite, ray_length, nr_of_rays, threshold_dendrite, n, image_plane_z);
//                                                        //////std::cout<<"sent Ray"<<std::endl;
            ////////std::cout<<"tmp_d: "<<tmp_d->coords[X_COORD]<<"  "<<tmp_d->coords[Y_COORD]<<std::endl;
            dendrite_vectors.push_back(tmp_d);
            
            dendrite_surface_point[X_COORD] = tmp_d->coords[X_COORD] + x_dendrite;
            
            dendrite_surface_point[Y_COORD] = tmp_d->coords[Y_COORD] + y_dendrite;
            
            dendrite_surface_point[Z_COORD] = z_dendrite;
            
            dendrite_surface_points->push_back(dendrite_surface_point);
            
            ////////std::cout<<"axon surface point : "<<axon_surface_point[0]<<"  "<<axon_surface_point[1]<<std::endl;
            ////////std::cout<<"dendrite surface point : "<<dendrite_surface_point[0]<<"  "<<dendrite_surface_point[1]<<std::endl;
            delete tmp_d;
            
            axon_surface_point = new double[2];
            dendrite_surface_point = new double[2];
            
            ////////std::cout<<"dendrite surface point : "<<surface_point[0]<<"  "<<surface_point[1]<<std::endl;
            
            
//             //////std::cout<< "Distance: " << tmp_distance << std::endl;
            }   
            
            
            
//             for( int i = 0; i < axon_vectors.size(); i++){
//             axon_touch_x = (axon_vectors.at(i)->coords[X_COORD]+ x_axon)*XYSAMPLING;
//             axon_touch_y = (axon_vectors.at(i)->coords[Y_COORD] + y_axon)*XYSAMPLING;
//             
//             
//             
// //             //////std::cout<<"axon_coords :"<<axon_vectors.at(i)->coords[X_COORD]<<std::endl;
// //              //////std::cout<<"axon_touch x :"<<axon_touch_x<<std::endl;
// //              //////std::cout<<"axon_touch y :"<<axon_touch_y<<std::endl;
// // //             
//                 for(int j = 0; j < dendrite_vectors.size(); j++){
//                     
// //                     //////std::cout<<"dendrite_coords :"<<dendrite_vectors.at(j)->coords[X_COORD]<<std::endl;
// //                     //////std::cout<<"x_dendrite :"<<x_dendrite<<std::endl;
//                     
//                     dendrite_touch_x = (dendrite_vectors.at(j)->coords[X_COORD] + x_dendrite)*XYSAMPLING;
//                     dendrite_touch_y = (dendrite_vectors.at(j)->coords[Y_COORD] + y_dendrite)*XYSAMPLING;
//                     
// /*                    //////std::cout<<"dendrite_touch x :"<<dendrite_touch_x<<std::endl;
//                     //////std::cout<<"dendrite_touch y :"<<dendrite_touch_y<<std::endl;
//      */               
//                     
//                    double dist_x = axon_touch_x - dendrite_touch_x;
//                    if( j == 0 && i == 0)
//                        min_x_dist = abs(dist_x);
//                     
//                    if(dist_x < min_x_dist)
//                    min_x_dist = abs(dist_x);
//                    
// //                   //////std::cout<<"min_x_dist: "<<min_x_dist<<std::endl;
//                    
//                    double dist_y = axon_touch_y - dendrite_touch_y;
//                    
//                    if( j == 0 && i == 0)
//                    min_y_dist = abs(dist_y);
//                     
//                    if(dist_y < min_y_dist)
//                    min_y_dist = abs(dist_y);
//                    
// //                   //////std::cout<<"min_y_dist: "<<min_y_dist<<std::endl;
//                    
//                    
// //                 dendrite_touch_x = dendritePt[X_COORD] + dendrite_vectors.at(i)->coords[X_COORD];
// //                 dendrite_touch_y = dendritePt[Y_COORD] + dendrite_vectors.at(i)->coords[Y_COORD];
// //                 
// //                 real_distance_x.push_back(abs( axon_touch_x - dendrite_touch_x));
// //                 real_distance_y.push_back(abs( axon_touch_y - dendrite_touch_y));
// //                 
//                     
//                 }
//             }

            double minimum_euclid_dist = 10000;
            int placeholder_axon, placeholder_dendrite = 0;
            
            for(int i = 0; i < axon_surface_points->size(); i++){
                
                
                
                
                for(int j = 0; j < dendrite_surface_points->size(); j++){
                    double * axon_point = new double[2];
                    double * dendrite_point = new double[2];
                    
                    axon_point = axon_surface_points->at(i);
                    dendrite_point = dendrite_surface_points->at(j);
                ////////std::cout<<"axon surface point = "<<axon_surface_points->at(i)[0]<<" "<<axon_surface_points->at(i)[1]<<"   dendrite surface point"<<dendrite_surface_points->at(j)[0]<<" "<<dendrite_surface_points->at(j)[1]<<std::endl;
                double dist = Utility::euclidean_distance(axon_point, dendrite_point, 3, 1);
                
//                 if(dist < 3)
//                //////std::cout<<"distance = "<<dist<<std::endl;
                //if(distance != 0)
                
                
                if(dist < minimum_euclid_dist /*&& dist != 0*/)
                    minimum_euclid_dist = dist;
                    placeholder_axon = i;
                    placeholder_dendrite = j;
//                     
                
                }
                
                
            }
            
            
                
            minimum_euclid_dist = minimum_euclid_dist + 0;
            
            if(istest){
            double * axonSRF = new double[2];
            double * dendriteSRF = new double[2];
            
            axonSRF = axon_surface_points->at(placeholder_axon);
            axonSRF[X_COORD] = axonSRF[X_COORD]*XYSAMPLING;
            axonSRF[Y_COORD] = axonSRF[Y_COORD]*XYSAMPLING;
            axonSRF[Z_COORD] = axonSRF[Z_COORD]*ZSAMPLING;
            dendriteSRF = dendrite_surface_points->at(placeholder_dendrite);
            dendriteSRF[X_COORD] = dendriteSRF[X_COORD]*XYSAMPLING;
            dendriteSRF[Y_COORD] = dendriteSRF[Y_COORD]*XYSAMPLING;
            dendriteSRF[Z_COORD] = dendriteSRF[Z_COORD]*ZSAMPLING;
            
            std::string landmarkname_axon = this->inputfilename;
            std::string landmarkname_dendrite = this->inputfilename;
            
            
            landmarkname_axon += "axon_surface_point";
            landmarkname_dendrite += "dendrite_surface_point";
            
//             writeSingleLandmarkFile(axonSRF, landmarkname_axon.c_str(), prox);
//             writeSingleLandmarkFile(dendriteSRF, landmarkname_dendrite.c_str(), prox);
            
            }
            
            ////////std::cout<<"min_distance = "<<minimum_euclid_dist<<std::endl;
            
            return minimum_euclid_dist;
                   
//            real_distance_x.push_back(minimum_euclid_dist);
    
    
    
}




// double ** AxonDendriteProximityFinder::readAndGetAmiraTransformation(const char * inputfilename)
// {
//     
//        
// //     //////std::cout<<transformFilename<<std::endl;
// //     //////std::cout<<filename<<std::endl;
//     
// //////std::cout<<"in readAmiraTransformations"<<std::endl;
//         double ** transformation = new double *[4];
//                 for(int ii = 0; ii < 4; ++ii)
//                 {
//                         transformation[ii] = new double[4];
//                         for(int jj = 0; jj < 4; ++jj)
//                         {
//                                 if(ii != jj)
//                                         transformation[ii][jj] = 0;
//                                 else
//                                         transformation[ii][jj] = 1;
//                         }
//                 }
// 
// 
//         std::ifstream inputStream(inputfilename);
//      
//         if(!inputStream.fail())
//         {
//                 const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
//                 const char * numbers = "0123456789";
//                 const char * signs = "+-";
//                 const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
//                 const char * whitespace = "\t ";
//                 
//                 std::string currentLine;
//                 unsigned int line = 0;
//                 
//                
//                 
//                 bool parameters = 1;
//                 bool transform = 0;
//                 bool correctSection = 1;
//                 bool correctPrevSection = 0;
//                 int sectionID = 0;
//                 unsigned int brackets = 0, transformBrackets = 0;
//                 unsigned int currentIndex = 0;
//                 
//                 while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
//                 {
//                         
//                         if(currentLine.size())
//                                 if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
//                                         {
// //                                              //////std::cout << "found correct section transform parameters!" << std::endl;
//                                                 unsigned int count = 0;
//                                                 std::string::size_type loc1, loc2, loc3;
//                                                 loc1 = currentLine.find_first_of(numbers, 0);
//                                                 loc2 = currentLine.find_first_of(signs, 0);
//                                                 if(loc2 != std::string::npos)
//                                                         if(loc2 < loc1)
//                                                                 loc1 = loc2;
//                                                 loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
//                                                 while(loc2 != std::string::npos && count < 16)
//                                                 {
//                                                         char * tmp1 = new char[loc2 - loc1];
//                                                         currentLine.copy(tmp1, loc2 - loc1, loc1);
//                                                         double ftmp1 = atof(tmp1);
//                                                         transformation[count%4][count/4]= ftmp1;        // amira files are columns after each other
//                                                         loc3 = loc2;
//                                                         loc1 = currentLine.find_first_of(numbers, loc3);
//                                                         loc2 = currentLine.find_first_of(signs, loc3);
//                                                         if(loc2 != std::string::npos)
//                                                                 if(loc2 < loc1)
//                                                                         loc1 = loc2;
//                                                         loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
//                                                         ++count;
//                                                         delete [] tmp1;
//                                                 }
//                                              //////std::cout << "transformation matrix:" << std::endl;
//                                              for(int ii = 0; ii < 4; ++ii)
//                                              {
//                                                      //////std::cout << "[";
//                                                      for(int jj = 0; jj < 4; ++jj)
//                                                      {
//                                                              if(jj < 3)
//                                                                      //////std::cout << transformation[ii][jj] << ",\t";
//                                                              else
//                                                                      //////std::cout << transformation[ii][jj];
//                                                      }
//                                                      //////std::cout << "]" << std::endl;
//                                              }
//                                                 //remove numeric artifacts from z-axis:
// //                                                 for(int ii = 0; ii < 2; ++ii)
// //                                                 {
// //                                                         transformation[2][ii] = 0;
// //                                                         transformation[ii][2] = 0;
// //                                                 }
// //                                                 transformation[2][2] = 1;
//                                         }
//                        
//                 }
//         }
//         
//         inputStream.close();
//         return transformation;
// }

// void AxonDendriteProximityFinder::getIDsFromMergedFile(AmiraSpatialGraph * merged_graph){
//     
//     std::vector< Edge * > * edges = amira_graph->edgesPointer();
//     std::vector< Edge * > * labeled_edges = merged_graph->edgesPointer();
//     
//     
//     
//     unsigned int numOfEdges = edges->size();
//     unsigned int numOfMergeEdges = labeled_edges->size();
// //    double * currentPoint = currentAxon.getEdgePointCoordinates();
//     std::list< double * >::iterator edge_it;
//     std::list<double * >::iterator merge_edge_it;
// //    unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
// 
//     
//     for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
//     {                       
//                 Edge * currentEdge = edges->at(pos);
//                 unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
//     
//     
//     
//         for(edge_it = currentEdge->edgePointCoordinates.begin() ; edge_it != currentEdge->edgePointCoordinates.end(); edge_it++)
//         {
//             
//                             double * current_coords = *edge_it;
//                             
//                          for(long merge_pos = numOfMergeEdges -1; merge_pos >= 0; merge_pos --)
//                          {
//                              bool isSameEdge = false;
//                              
//                              Edge * currentLabeledEdge = labeled_edges->at(merge_pos);
//                              unsigned int labeledEdgeSize = currentLabeledEdge->edgePointCoordinates.size();
//                              
//                              for(merge_edge_it = currentLabeledEdge->edgePointCoordinates.begin() ; merge_edge_it != currentLabeledEdge->edgePointCoordinates.end(); merge_edge_it++)
//                              {
//                                  double * current_labeled_coords = *merge_edge_it;
//                                  
//                                  if(current_coords[X_COORD] == current_labeled_coords[X_COORD] && current_coords[Y_COORD] == current_labeled_coords[Y_COORD] && current_coords[Z_COORD] == current_labeled_coords[Z_COORD])
//                                  {
//                                      isSameEdge = true;
//                                  }    
//                                 else
//                                 {                                    
//                                 isSameEdge = false;
//                                 break;
//                                 }
//                                  
//                                 
//                                  
//                               }
//                            
//                            currentEdge->label = currentLabeledEdge->label;
//                              
//                              
//                         }  
//         
//         }
//     
//     }
// }

// void AxonDendriteProximityFinder::cutOffGraph(AmiraSpatialGraph * whole_graph){
//     
// //    BoundingBox bounding_box;
// //    getBoundingBox();
//     
//     double * center_coords;
// //
//     AmiraSpatialGraph * merged_graph = whole_graph;
// //     AmiraSpatialGraph * transformed_graph = this->amira_graph;
//     AmiraSpatialGraph * cut_out_merge = new AmiraSpatialGraph();
// //    merged_graph->applyInverseTransformation();
//     
// //     transformation = new double *[4];
// //                 for(int ii = 0; ii < 4; ++ii)
// //                 {
// //                         transformation[ii] = new double[4];
// //                         for(int jj = 0; jj < 4; ++jj)
// //                         {
// //                                 if(ii != jj)
// //                                         transformation[ii][jj] = 0;
// //                                 else
// //                                         transformation[ii][jj] = 1;
// //                         }
// //                 }
// //                 
//     double boundingbox[6];
// //    transformed_graph->applyInverseTransformation();
//     amira_graph->getBoundingBox(boundingbox);
// //     center_coords = boundingbox.getCenterCoords();
// //     //////std::cout<<"X center: "<<center_coords[X_COORD]<<"  Y center:  "<<center_coords[Y_COORD]<<"  Z center: "<<center_coords[Z_COORD]<<std::endl;
//     //////std::cout<<"X center: "<<boundingbox[0]<<"  Y center:  "<<boundingbox[2]<<"  Z center: "<<boundingbox[4]<<std::endl;
//     //////std::cout<<"X center: "<<boundingbox[1]<<"  Y center:  "<<boundingbox[3]<<"  Z center: "<<boundingbox[5]<<std::endl;
//     
//     int x_dim, y_dim, z_dim;
// //     x_dim = graph_bounding_box.X_dim();
// //     y_dim = graph_bounding_box.Y_dim();
// //     z_dim = graph_bounding_box.Z_dim();
//     
// //     double x_min = center_coords[X_COORD] - x_dim/2;
// //     double y_min = center_coords[Y_COORD] - y_dim/2;
// //     double z_min = center_coords[Z_COORD] - z_dim/2;
// //     
// //     double x_max = center_coords[X_COORD] + x_dim/2;
// //     double y_max = center_coords[Y_COORD] + y_dim/2;
// //     double z_max = center_coords[Z_COORD] + z_dim/2;
//     
//     double x_min = boundingbox[0];
//     double y_min = boundingbox[2];
//     double z_min = boundingbox[4];
//     
//     double x_max = boundingbox[1];
//     double y_max = boundingbox[3];
//     double z_max = boundingbox[5];
//     
//     
// 
//     
//     
//     
//                 
//     
//      std::vector< Edge * > * edges = amira_graph->edgesPointer();
//      std::vector< Edge * > * labeled_edges = merged_graph->edgesPointer();
//      std::vector< Vertex *> * labeled_vertices = merged_graph->verticesPointer();
// //     
// //     
// //     
//      unsigned int numOfEdges = edges->size();
//      unsigned int numOfMergeEdges = labeled_edges->size();
//      
//      //////std::cout<<"merged edges size : "<<numOfMergeEdges<<std::endl;
// // //    unsigned int numOfMergedVertices = labeled_vertices->size();
// // //    double * currentPoint = currentAxon.getEdgePointCoordinates();
// //     std::list< double * >::iterator edge_it;
//      std::list< double * >::iterator medge_it;
// //     
// //    unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
// 
//     
// //     for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
// //     {           
// //                 std::vector<double *> * tmp_list = new std::vector<double *>;
// //                 Edge * currentEdge = edges->at(pos);
// //                 unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
// //     
// //                 double * axonPoint = currentAxon.getNext();
// //     
// //         for(edge_it = currentEdge->edgePointCoordinates.begin() ; edge_it != currentEdge->edgePointCoordinates.end(); edge_it++)
// //         {
// //             
// //                             double * current_coords = *edge_it;
// //                             
// //                            
// //                             //////std::cout<<"X: "<<current_coords[X_COORD]<<"  Y:  "<<current_coords[Y_COORD]<<"  Z: "<<current_coords[Z_COORD]<<std::endl;
// //                             
// //                              tmp_list->push_back(axonPoint);
// //                             
// //                             
// // //                         
// //         
// //         }
// //         for(int i = 0; i < tmp_list->size(); i++){
// //         //////std::cout<<tmp_list->at(i)[X_COORD]<<"  "<<tmp_list->at(i)[Y_COORD]<<"  "<<tmp_list->at(i)[Z_COORD]<<"  "<<std::endl;
// //         graph_bounding_box.expandBox_amiraCoords(tmp_list->at(i));
// //         //////std::cout<<"X_dim : "<<graph_bounding_box.X_dim()<<" Y_dim : "<<graph_bounding_box.Y_dim()<<" Z_dim : "<<graph_bounding_box.Z_dim()<<std::endl;
// //         }
// //     
// //     }
// 
// //     for(long pos = numOfMergedVertices -1; pos >= 0; pos--){
// //         
// //         Vertex * current_vertex = labeled_vertices->at(pos);
// //         
// //         if(current_vertex[X_COORD] < x_max && current_vertex[Y_COORD] < y_max && current_vertex[Z_COORD] < z_max && current_vertex[X_COORD] > x_min && current_vertex[Y_COORD] > y_min && current_vertex[Z_COORD] > z_min)
// //         {
// //             cut_out_merge->addVertex(current_vertex);
// //             
// //         }
// //         
// //         
// //     }
//     
//     
//     for(long pos = numOfMergeEdges -1; pos >= 0; pos--)      //for each edge in list
//     {                       
//                 Edge * currentEdge = labeled_edges->at(pos);
// //                Edge * transformedEdge = currentEdge;
// //                unsigned int edgeSize = transformedEdge->edgePointCoordinates.size();
//     
//     
//     
//         for(medge_it = currentEdge->edgePointCoordinates.begin() ; medge_it != currentEdge->edgePointCoordinates.end(); medge_it++)
//         {
//             
//                             double * current_coords = *medge_it;
// //                             //////std::cout<<"X: "<<current_coords[X_COORD]<<"  Y:  "<<current_coords[Y_COORD]<<"  Z: "<<current_coords[Z_COORD]<<" max x :"<<x_max<<" max y :"<<y_max<<" max z :"<<z_max<<std::endl;
//                             
//                             if(current_coords[X_COORD] < x_max && current_coords[Y_COORD] < y_max && current_coords[Z_COORD] < z_max && current_coords[X_COORD] > x_min && current_coords[Y_COORD] > y_min && current_coords[Z_COORD] > z_min)
//                             {
//                                 cut_out_merge->addEdge(currentEdge);
//                                 
// //                                //////std::cout<<"added Edge"<<std::endl;
//                             }
//                             
//     
//                             
//                             
//                             
//         }                    
//     }
//     std::vector< Edge * > * cutEdges = cut_out_merge->edgesPointer();
//     unsigned int numOfMergeEdgesAfterCut = cutEdges->size();
//      
//      //////std::cout<<"cut merged edges size : "<<numOfMergeEdgesAfterCut<<std::endl;
//      cut_out_merge->applyInverseTransformation();
//      this->amira_graph = cut_out_merge;
//      double merge_boundingbox[6];
//      cut_out_merge->getBoundingBox(merge_boundingbox);
//      
//      const char * output_name = "a";
//      
//      
//        Reader checkGraph(inputfilename, output_name);
//        checkGraph.setSpatialGraph(cut_out_merge);
//        checkGraph.writeSpatialGraphFile();
// 
//      
// //     center_coords = boundingbox.getCenterCoords();
// //     //////std::cout<<"X center: "<<center_coords[X_COORD]<<"  Y center:  "<<center_coords[Y_COORD]<<"  Z center: "<<center_coords[Z_COORD]<<std::endl;
//     //////std::cout<<"X merge: "<<merge_boundingbox[0]<<"  Y merge:  "<<merge_boundingbox[2]<<"  Z merge: "<<merge_boundingbox[4]<<std::endl;
//     //////std::cout<<"X merge: "<<merge_boundingbox[1]<<"  Y merge:  "<<merge_boundingbox[3]<<"  Z merge: "<<merge_boundingbox[5]<<std::endl;
// //     
// //     
// //     ;    for(long pos = numOfMergedVertices -1; pos >= 0; pos--){
// //         
// //         Vertex * current_vertex = labeled_vertices->at(pos);
// //         
// //         if(current_vertex[X_COORD] < x_max && current_vertex[Y_COORD] < y_max && current_vertex[Z_COORD] < z_max && current_vertex[X_COORD] > x_min && current_vertex[Y_COORD] > y_min && current_vertex[Z_COORD] > z_min)
// //         {
// //             cut_out_merge->addVertex(current_vertex);
// //             
// //         }
// //         
// //         
// //     }
// //     
// //     
// //     
// //     
// //     
// }

// void AxonDendriteProximityFinder::applyInverseTransformation(AmiraSpatialGraph * graph_to_be_transformed)
// {
//     AmiraSpatialGraph * transformed_graph = new AmiraSpatialGraph();
//     
//         if(/*!isIdentity && !inverseTransformationApplied*/true)
//         {
// //              printTransformation();
//                 graph_to_be_transformed->inverseTransformationApplied = 1;
// //              std::flush(//////std::cout << "Transforming " << vertices.size() << " vertices..." << std::endl);
//                 std::vector<Vertex *> * vertices = graph_to_be_transformed->verticesPointer();
//                 std::vector< Vertex * >::iterator vertexIt;
//                 
// //              int vertexCount = 1;
//                 for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
//                 {
// //                      //////std::cout << "Transforming Vertex " << vertexCount << " of " << vertices.size() << std::endl;
// //                      ++vertexCount;
//                         double oldCoords[4], newCoords[4];
//                         for(int ii = 0; ii < 3; ++ii)
//                         {
//                                 oldCoords[ii] = (*vertexIt)->coordinates[ii];
//                                 newCoords[ii] = 0;
//                         }
//                         oldCoords[3] = 1;
//                         newCoords[3] = 1;
//                         for(int ii = 0; ii < 3; ++ii)
//                                 for(int jj = 0; jj < 4; ++jj)
//                                         newCoords[ii] += inverse_transformation[ii][jj]*oldCoords[jj];
//                         
//                         for(int ii = 0; ii < 3; ++ii)
//                                 (*vertexIt)->coordinates[ii] = newCoords[ii];
//                 }
//                 
// //              std::flush(//////std::cout << "Transforming " << edges.size() << " edges..." << std::endl);
//                 std::vector< Edge * > * edges = graph_to_be_transformed->edgesPointer();
//                 std::vector< Edge * >::iterator edgeIt;
// //              int edgeCount = 1;
//                 for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
//                 {
// //                      //////std::cout << "Transforming Edge " << edgeCount << " of " << edges.size() << std::endl;
// //                      //////std::cout << (*edgeIt)->edgePointCoordinates.size() << " points" << std::endl;
// //                      ++edgeCount;
//                         std::list< double * >::iterator edgeListIt;
//                         for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
//                         {
//                                 double oldCoords[4], newCoords[4];
//                                 for(int ii = 0; ii < 3; ++ii)
//                                 {
//                                         oldCoords[ii] = (*edgeListIt)[ii];
//                                         newCoords[ii] = 0;
//                                 }
//                                 oldCoords[3] = 1;
//                                 newCoords[3] = 1;
//                                 for(int ii = 0; ii < 3; ++ii)
//                                         for(int jj = 0; jj < 4; ++jj)
//                                                 newCoords[ii] += inverse_transformation[ii][jj]*oldCoords[jj];
//                                 
//                                 for(int ii = 0; ii < 3; ++ii)
//                                         (*edgeListIt)[ii] = newCoords[ii];
//                         }
//                 }
// //              std::flush(//////std::cout << "done!" << std::endl);
//         }
//         
//         return transformed_graph;
// };

bool AxonDendriteProximityFinder::isInSectionBoundingBox(double * coords){
//ensures that the point of the transformed whole cell graph is in the current section's bounding box    
    
#ifdef DEBUG    
    ////////std::cout<<"in isInSectionBoundingBox!"<<std::endl;
#endif    
    
    bool isInSection = false;
//    double boundingbox[6];
//    transformed_graph->applyInverseTransformation();
    
//     center_coords = boundingbox.getCenterCoords();
//     //////std::cout<<"X center: "<<center_coords[X_COORD]<<"  Y center:  "<<center_coords[Y_COORD]<<"  Z center: "<<center_coords[Z_COORD]<<std::endl;
/*    //////std::cout<<"X center: "<<boundingbox[0]<<"  Y center:  "<<boundingbox[2]<<"  Z center: "<<boundingbox[4]<<std::endl;
    //////std::cout<<"X center: "<<boundingbox[1]<<"  Y center:  "<<boundingbox[3]<<"  Z center: "<<boundingbox[5]<<std::endl;
  */  
//    int x_dim, y_dim, z_dim;
//     x_dim = graph_bounding_box.X_dim();
//     y_dim = graph_bounding_box.Y_dim();
//     z_dim = graph_bounding_box.Z_dim();
    
//     double x_min = center_coords[X_COORD] - x_dim/2;
//     double y_min = center_coords[Y_COORD] - y_dim/2;
//     double z_min = center_coords[Z_COORD] - z_dim/2;
//     
//     double x_max = center_coords[X_COORD] + x_dim/2;
//     double y_max = center_coords[Y_COORD] + y_dim/2;
//     double z_max = center_coords[Z_COORD] + z_dim/2;
    
    double x_min = sectionboundingbox[0];
    double y_min = sectionboundingbox[2];
    double z_min = sectionboundingbox[4];
    
    double x_max = sectionboundingbox[1];
    double y_max = sectionboundingbox[3];
    double z_max = sectionboundingbox[5];
    
    //////std::cout<<"X: "<<coords[X_COORD]<<"  Y:  "<<coords[Y_COORD]<<"  Z: "<<coords[Z_COORD]<<" max x :"<<x_max<<" max y :"<<y_max<<" max z :"<<z_max<<std::endl;
    
    
    if(coords[X_COORD] < x_max && coords[Y_COORD] < y_max && coords[Z_COORD] < z_max && coords[X_COORD] > x_min && coords[Y_COORD] > y_min  && coords[Z_COORD] > z_min)
                            {
                                isInSection = true;
                                
//                                //////std::cout<<"added Edge"<<std::endl;
                            }
    
    return isInSection;
}


bool AxonDendriteProximityFinder::isTheSameCell(int AxonID, int DendriteID)
{
    
    if(amira_graph->getAxonCellIndex(AxonID) == amira_graph->getDendriteCellIndex(DendriteID))
    {
        return true;
    }
    else
    {
        return false;
    }
    
};

bool AxonDendriteProximityFinder::isAxon(int ID)
{
    std::vector<int >* id_array = amira_graph->getAxonIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
        }
        
    }
    
    return false;
    
    
};

bool AxonDendriteProximityFinder::isDendrite(int ID)
{
    std::vector<int >* id_array = amira_graph->getDendriteIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
        }
        
    }
    
    return false;
    
    
};



int AxonDendriteProximityFinder::getCellIndex(std::vector<int> *IDArray, int ID)
{
    int id_index = 0xFFFF;
    
    for(int i = 0; i < IDArray->size(); i++ )
    {
        if(IDArray->at(i) == ID)
        {
            return i;
        }
    }
    return id_index;
    
};

/*
void AxonDendriteProximityFinder::writeProximityImages(const char * outputFilename, std::vector<Proximity *>* proxlistlist[][])
{
    //std::vector<Proximity *>* proxlistlist = *proxlistlistptr[][];

    for(int i = 0; i < amira_graph->getNumberofCells(); i++)
    {
        for(int j = 0; j < amira_graph->getNumberofCells(); j++)
        {
            if(i!=j)
            {
              
              std::string format(outputFilename);
              format += "_CELL#";
              //format += itoa(i+1);
              format += "-->CELL#";
              //format += itoa(j+1);
              
              //format += "_CELL#1-->CELL#2";
              for(int k=0; k < proxlistlist[i][j]->size(); k++)
              {
                  std::string format(outputFilename);
                  //format += "_CELL#1-->CELL#2";
                  (proxlistlist[i][j]->at(k))->writeProximityImage(format.c_str(), k);
              }
              //proxlistlist[i][j]->writeProximityImage(format.c_str(), i);
            }
        }
    }
    
    for(int i=0; i<proximity_list_1_2.size(); i++)
    {
        std::string format(outputFilename);
        format += "_CELL#1-->CELL#2";
        proximity_list_1_2.at(i)->writeProximityImage(format.c_str(), i);
    }
    for(int i=0; i<proximity_list_2_1.size(); i++)
    {
        std::string format1(outputFilename);
        format1 += "_CELL#2-->CELL#1";
        proximity_list_2_1.at(i)->writeProximityImage(format1.c_str(), i);
    }
    for(int i=0; i<proximity_list_1_3.size(); i++)
    {
        std::string format2(outputFilename);
        format2 += "_CELL#1-->CELL#3";
        proximity_list_1_3.at(i)->writeProximityImage(format2.c_str(), i);
    }
    for(int i=0; i<proximity_list_3_1.size(); i++)
    {
        std::string format3(outputFilename);
        format3 += "_CELL#3-->CELL#1";
        proximity_list_3_1.at(i)->writeProximityImage(format3.c_str(), i);
    }
    for(int i=0; i<proximity_list_2_3.size(); i++)
    {
        std::string format4(outputFilename);
        format4 += "_CELL#2-->CELL#3";
        proximity_list_2_3.at(i)->writeProximityImage(format4.c_str(), i);
    }
    for(int i=0; i<proximity_list_3_2.size(); i++)
    {
        std::string format5(outputFilename);
        format5 += "_CELL#3-->CELL#2";
        proximity_list_3_2.at(i)->writeProximityImage(format5.c_str(), i);
    }
    
};*/

void AxonDendriteProximityFinder::set_edge_lists(AmiraSpatialGraph * merged_input_graph, bool normalmode)
{
        int count = 0; 
        std::vector< Edge * > * edges = merged_input_graph->edgesPointer();
        int numOfEdges = edges->size();
        bool isInSection = false;
        
        //if(normalmode)
        {

            
            for(int i=0; i<numOfEdges; i++)  //for each edge
            {
                Edge * currentEdge = edges->at(i);
                std::list< double * >::iterator edge_it;
                
                for(edge_it = currentEdge->edgePointCoordinates.begin() ; edge_it != currentEdge->edgePointCoordinates.end(); edge_it++)
                {
                    double * coords = *edge_it;
                    //if(isInSectionBoundingBox(coords))
                    {
                            
                        //if(currentEdge->label == DENDRITE_ID_1 || currentEdge->label == DENDRITE_ID_2 || currentEdge->label == DENDRITE_ID_3)
                        if(isDendrite(currentEdge->label))
                        {
                            dendrite_list.push_back(SimpleDendrite(currentEdge, count));
                        }
                        /*else if(currentEdge->label == UNKNOWN_ID && include_unknown)
                        {
                            axon_list.push_back(SimpleAxon(currentEdge, count));
                            dendrite_list.push_back(SimpleDendrite(currentEdge, count));
                        }*/
                        //else if(currentEdge->label == AXON_ID_1 || currentEdge->label == AXON_ID_2 || currentEdge->label == AXON_ID_3)
                        if(isAxon(currentEdge->label))
                        {
                            axon_list.push_back(SimpleAxon(currentEdge, count));
                        }
                        count++;
                       
                    }
                    
                }
                    
            }
            
        }
     
     
       
        
};

#if 0
/********************************************************************************************************************
Method: writeProximityLandmarks()

Gets the centers of mass for each proximity and then prints them into a landmark file
*********************************************************************************************************************/
void AxonDendriteProximityFinder::writeProximityLandmarks(const char * outputFilename)
{
    
    //////std::cout<<"In writeProximityLandmarks"<<std::endl;
    std::vector<double *> * center_list_1_2 = new std::vector<double *>;
    std::vector<double *> * center_list_1_3 = new std::vector<double *>;
    std::vector<double *> * center_list_2_3 = new std::vector<double *>;
    std::vector<double *> * center_list_2_1 = new std::vector<double *>;
    std::vector<double *> * center_list_3_1 = new std::vector<double *>;
    std::vector<double *> * center_list_3_2 = new std::vector<double *>;
    
//    std::string format(outputFilename);
//    format += "/";
//     format += outputFilename;
    
    //////std::cout<<outputFilename<<std::endl;
    
    std::string boutonname(outputFilename);
    boutonname += "_boutons";
    
    
    for(int i=0; i<this->proximity_list.size(); i++)
    {
        Proximity * prox = proximity_list.at(i);
        if(prox->Axon_ID == AXON_ID_1 && prox->Dendrite_ID == DENDRITE_ID_2 )
        center_list_1_2->push_back(prox->getCenterCoords());

        if( prox->Axon_ID == AXON_ID_2 && prox->Dendrite_ID == DENDRITE_ID_1)
        center_list_2_1->push_back(prox->getCenterCoords());
        
        if(prox->Axon_ID == AXON_ID_1 && prox->Dendrite_ID == DENDRITE_ID_3 )
        center_list_1_3->push_back(prox->getCenterCoords());
        
        if( prox->Axon_ID == AXON_ID_3 && prox->Dendrite_ID == DENDRITE_ID_1)
        center_list_3_1->push_back(prox->getCenterCoords());
        
        if(prox->Axon_ID == AXON_ID_2 && prox->Dendrite_ID == DENDRITE_ID_3)
        center_list_2_3->push_back(prox->getCenterCoords());
        
        if(prox->Axon_ID == AXON_ID_3 && prox->Dendrite_ID == DENDRITE_ID_2)
        center_list_3_2->push_back(prox->getCenterCoords());

    }
    
    if(center_list_1_2->size() > 0){
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        std::string outputFilename_1_2(outputFilename);
        outputFilename_1_2 += "_CELL#1-->CELL#2";
        //////std::cout<<outputFilename_1_2<<std::endl;
        writeInfoFile(center_list_1_2, outputFilename_1_2.c_str(), this->proximity_list_1_2);
        writeHXFile(center_list_1_2, outputFilename_1_2.c_str(), this->proximity_list_1_2);

            
        center_list_1_2 = transformCoordinates(center_list_1_2);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_1_2, outputFilename_1_2.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_1_2, outputFilename_1_2.c_str());
        

    
        for(int i=0; i<center_list_1_2->size(); i++)
        {
            double * coords = center_list_1_2->at(i);
            writeSingleLandmarkFile(coords, outputFilename_1_2.c_str(), i);
        }
    
    }
    
    if(center_list_2_1->size() > 0){
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        std::string outputFilename_2_1(outputFilename);
        outputFilename_2_1 += "_CELL#2-->CELL#1";
        //////std::cout<<outputFilename_2_1<<std::endl;
        writeInfoFile(center_list_2_1, outputFilename_2_1.c_str(), this->proximity_list_2_1);
        writeHXFile(center_list_2_1, outputFilename_2_1.c_str(), this->proximity_list_2_1);

            
        center_list_2_1 = transformCoordinates(center_list_2_1);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_2_1, outputFilename_2_1.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_2_1, outputFilename_2_1.c_str());
        

    
        for(int i=0; i<center_list_2_1->size(); i++)
        {
            double * coords = center_list_2_1->at(i);
            writeSingleLandmarkFile(coords, outputFilename_2_1.c_str(), i);
        }
    
    }
    

    
        if(center_list_1_3->size() > 0){
        std::string outputFilename_1_3(outputFilename);
        outputFilename_1_3 += "_CELL#1-->CELL#3";    
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        writeInfoFile(center_list_1_3, outputFilename_1_3.c_str(), this->proximity_list_1_3);
        writeHXFile(center_list_1_3, outputFilename_1_3.c_str(), this->proximity_list_1_3);

        center_list_1_3 = transformCoordinates(center_list_1_3);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_1_3, outputFilename_1_3.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_1_3, outputFilename_1_3.c_str());
        

    
        for(int i=0; i<center_list_1_3->size(); i++)
        {
            double * coords = center_list_1_3->at(i);
            writeSingleLandmarkFile(coords, outputFilename_1_3.c_str(), i);
        }
    
    }    
    
    if(center_list_3_1->size() > 0){
        std::string outputFilename_3_1(outputFilename);
        outputFilename_3_1 += "_CELL#3-->CELL#1";    
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        writeInfoFile(center_list_3_1, outputFilename_3_1.c_str(), this->proximity_list_3_1);
        writeHXFile(center_list_3_1, outputFilename_3_1.c_str(), this->proximity_list_3_1);

        center_list_3_1 = transformCoordinates(center_list_3_1);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_3_1, outputFilename_3_1.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_3_1, outputFilename_3_1.c_str());
        

    
        for(int i=0; i<center_list_3_1->size(); i++)
        {
            double * coords = center_list_3_1->at(i);
            writeSingleLandmarkFile(coords, outputFilename_3_1.c_str(), i);
        }
    
    }   

//     if(center_list_3_1->size() > 0){
//         std::string outputFilename_3_1(outputFilename);
//         outputFilename_3_1 += "_CELL#3-->CELL#1";
//         //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
//         writeInfoFile(center_list_3_1, outputFilename_3_1.c_str(), this->proximity_list_3_1);
//         writeHXFile(center_list_3_1, outputFilename_3_1.c_str(), this->proximity_list_3_1);
// 
//         center_list_3_1 = transformCoordinates(center_list_3_1);
//         //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
//         
//         writeLandmarkFile(center_list_3_1, outputFilename_3_1.c_str());
//         //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
//         writeListFile(center_list_3_1, outputFilename_3_1.c_str());
//         
// 
//     
//         for(int i=0; i<center_list_3_1->size(); i++)
//         {
//             double * coords = center_list_3_1->at(i);
//             writeSingleLandmarkFile(coords, outputFilename_3_1.c_str(), i);
//         }
//     
//     }
    
    if(center_list_2_3->size() > 0){
        std::string outputFilename_2_3(outputFilename);
        outputFilename_2_3 += "_CELL#2-->CELL#3";
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        writeInfoFile(center_list_2_3, outputFilename_2_3.c_str(), this->proximity_list_2_3);
        writeHXFile(center_list_2_3, outputFilename_2_3.c_str(), this->proximity_list_2_3);

            
        center_list_2_3 = transformCoordinates(center_list_2_3);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_2_3, outputFilename_2_3.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_2_3, outputFilename_2_3.c_str());
        

    
        for(int i=0; i<center_list_2_3->size(); i++)
        {
            double * coords = center_list_2_3->at(i);
            writeSingleLandmarkFile(coords, outputFilename_2_3.c_str(), i);
        }
    
    }    
    
    if(center_list_3_2->size() > 0){
        std::string outputFilename_3_2(outputFilename);
        outputFilename_3_2 += "_CELL#3-->CELL#2";
        //center_list = transformCoordinatesZScale(center_list, this->offsetfile);
        writeInfoFile(center_list_3_2, outputFilename_3_2.c_str(), this->proximity_list_3_2);
        writeHXFile(center_list_3_2, outputFilename_3_2.c_str(), this->proximity_list_3_2);
//         std::string outputname_nondecon = outputFilename_3_2;
//         outputname_nondecon += "_nondecon";
//         //////std::cout<<"outputname_nondecon : "<<outputname_nondecon<<std::endl;
//         char * outputname = new char[outputname_nondecon.size()];
//         outputname_nondecon.copy(outputname, outputname_nondecon.size(), 0);
    
//         writeInfoFile(center_list_3_2, outputname_nondecon.c_str());
//         writeHXFile(center_list_3_2, outputname_nondecon.c_str());
            
        center_list_3_2 = transformCoordinates(center_list_3_2);
        //this->all_bouton_list = transformCoordinates(this->all_bouton_list);
        
        writeLandmarkFile(center_list_3_2, outputFilename_3_2.c_str());
        //writeLandmarkFile(this->all_bouton_list, boutonname.c_str());
        writeListFile(center_list_3_2, outputFilename_3_2.c_str());
        

    
        for(int i=0; i<center_list_3_2->size(); i++)
        {
            double * coords = center_list_3_2->at(i);
            writeSingleLandmarkFile(coords, outputFilename_3_2.c_str(), i);
        }
        
    
    }       
};


/*******************************************************************************/
/*                                                                             */
/*writeLandmarkFile()                                                          */  
/*      writes a list into a valid landmarkAscii file                          */
/*                                                                             */
/*******************************************************************************/
void AxonDendriteProximityFinder::writeLandmarkFile(std::vector<double*> * list, const char* outputFilename)
{
  //////std::cout<< "in writeLandmarkFile" << std::endl;
  
// transformToWorldCoordinates(list);
  
  std::string format = outputFilename;
  format += "_locations.landmarkAscii";
  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << list->size()     << std::endl;
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
  
  for(int i=0; i < list->size(); i++)
  {
    LandMarkData << list->at(i)[X_COORD] << " " << list->at(i)[Y_COORD] << " " << list->at(i)[Z_COORD] << std::endl;
  } 
  
  LandMarkData.close();
};

/********************************************************************************/
/*writeListFile()                                                               */
/*      writes a list into a text file                                          */
/********************************************************************************/


void AxonDendriteProximityFinder::writeListFile(std::vector<double*> * list, const char* outputFilename)
{
//      writes a list-file with coordinates, bouton presence and confidence value
    std::string format(outputFilename);
    format += "_locations_list.txt";
    std::ofstream ListData( format.c_str() );
    ImageType::SizeType region_size;   
    
    ListData << "Proximity location's coordinates"<< std::endl;
    ListData << " "<< std::endl;
    ListData << " "<< std::endl;
    ListData << "     Bouton detected     Confidence Value       Amira coordinates               center coordinates"<< std::endl;
    
    for(int i=0; i < list->size(); i++)
    {
        Proximity * prox = this->proximity_list.at(i);
        bool has_bouton = prox->getBouton();

        double confidence_value = prox->getConfidenceValue();
 
        
        ImageType::RegionType output_region = prox->getOutputRegion();
        region_size = output_region.GetSize();
   
        
        ListData <<std::setfill(' ')<< std::setw(2)<< i <<std::setfill(' ')<< std::setw(7)<<std::boolalpha<<has_bouton<<std::setfill(' ')<< std::setw(23)<<confidence_value<<std::setfill(' ')<< std::setw(23)<< list->at(i)[X_COORD] - region_size[0]/2*XYSAMPLING << " " << std::setfill(' ')<< std::setw(8) << list->at(i)[Y_COORD] - region_size[1]/2*XYSAMPLING <<" "<< std::setfill(' ')<< std::setw(8) << list->at(i)[Z_COORD] - region_size[2]/2*ZSAMPLING <<" "<<std::setfill(' ')<< std::setw(13)<< list->at(i)[X_COORD] <<" "<< std::setfill(' ')<< std::setw(8)<< list->at(i)[Y_COORD]<<" " << std::setfill(' ')<< std::setw(8) << list->at(i)[Z_COORD]<<" "<<std::endl;
    }
    
    ListData.close();
};


void AxonDendriteProximityFinder::writeHXFile(std::vector<double*> * list, const char* outputFilename, std::vector<Proximity *> proximity_list)
{
  //writes an HXFile to ensure right positioning of the proximity images in Amira  
    
    
    for(int i=0; i < list->size(); i++)
    {

         
        
        Proximity * prox = proximity_list.at(i);
        
        ImageType::Pointer prox_image = prox->getProximityImage();        
        ImageType::RegionType output_region = prox->getOutputRegion();
        ImageType::RegionType iterator_region = prox->getIteratorRegion();
        ImageType::SizeType region_size;
        
        int Axon_ID = prox->Axon_ID;
        int Dendrite_ID = prox->Dendrite_ID;
        
           
    if(Axon_ID == AXON_ID_1)
        Axon_ID = 1;
    if(Axon_ID == AXON_ID_2)
        Axon_ID = 2;
    if(Axon_ID == AXON_ID_3)
        Axon_ID = 3;
    
    if(Dendrite_ID == DENDRITE_ID_1)
        Dendrite_ID = 1;
    if(Dendrite_ID == DENDRITE_ID_2)
        Dendrite_ID = 2;
    if(Dendrite_ID == DENDRITE_ID_3)
        Dendrite_ID = 3;
        

     
        region_size = output_region.GetSize();
        ImageType::IndexType start = iterator_region.GetIndex();
        

    
    
        std::string format2 = outputFilename;
        format2 += "_roi";
        format2 += Utility::makeString(i);
        format2 += "/loadROI.hx";
        //////std::cout<<"HX Filename: "<<format2<<std::endl;
        std::ofstream ListData( format2.c_str() );
        
        std::string proximity_label = "_Axon_";
        proximity_label += Utility::makeString(Axon_ID);
        proximity_label += "_Dendrite_";
        proximity_label += Utility::makeString(Dendrite_ID);
        proximity_label += "_";
        

        
        ListData << "# Amira Script"<< std::endl;
        ListData << " "<< std::endl;
        ListData << "[ load ${SCRIPTDIR}/proximity"<<proximity_label<<i<<".info ] setLabel {proximity"<<proximity_label<<i<<".info}"<< std::endl;
        ListData << "       "<<std::endl;
        
        ListData << "proximity"<<proximity_label<<i<<".info setTranslation " <<start[X_COORD]*XYSAMPLING<<" "<<start[Y_COORD]*XYSAMPLING<<" "<<start[Z_COORD]*ZSAMPLING<<" "<<std::endl;
        //ListData << "proximity"<<proximity_label<<i<<".info setTranslation " <<list->at(i)[X_COORD] - region_size[0]/2*XYSAMPLING<<" "<<list->at(i)[Y_COORD] - region_size[1]/2*XYSAMPLING<<" "<<list->at(i)[Z_COORD] - region_size[2]/2*ZSAMPLING<<std::endl;
        ListData << "proximity"<<proximity_label<<i<<".info applyTransform "<<  std::endl;
        ListData << "proximity"<<proximity_label<<i<<".info setTransform " <<transformation[0][0]<< " "<<transformation[1][0]<< " "<<transformation[2][0]<< " "<<transformation[3][0]<< " "<<transformation[0][1]<< " "<<transformation[1][1]<< " "<<transformation[2][1]<< " "<<transformation[3][1]<< " "<<transformation[0][2]<< " "<<transformation[1][2]<< " "<<transformation[2][2]<< " "<<transformation[3][2]<< " "<<transformation[0][3]<< " "<<transformation[1][3]<< " "<<transformation[2][3]<< " "<<transformation[3][3]<< " "<< std::endl;
        
        
            
        ListData.close(); 
   
        
        }
    
    
}

void AxonDendriteProximityFinder::writeMasterHxFile(const char* outputFilename, const char * amiraMergeName, int cellIDs){
  //writes an HXFile to ensure right positioning of the proximity images in Amira  
    std::string format1 = outputFilename;
            
    format1 += "_loadAll.hx";
    std::string AmiraMergedTotalString(amiraMergeName);
    unsigned int folderPosition = AmiraMergedTotalString.find_last_of("/");
    std::string AmiraMergedFileString = AmiraMergedTotalString.substr(folderPosition+1);
    int AxonID, DendriteID;
    
    std::ofstream ListData(format1.c_str());
    
    ListData << "# Amira Script"<< std::endl;
    ListData << " "<< std::endl;
    

        if(cellIDs == 12){        
            
            for(int i=0; i < proximity_list_1_2.size(); i++)
        {         
            Proximity * prox = proximity_list_1_2.at(i);
         AxonID = prox->Axon_ID;
         DendriteID = prox->Dendrite_ID;
          
         std::string format2 = "[ load ${SCRIPTDIR}/proximity";
         
         if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_2 || AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#2";

        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#3";

        if(AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_2)
        format2 += "_CELL#2-->CELL#3";

        
        format2 += "_roi";
        format2 += Utility::makeString(i);
        format2 += "/loadROI.hx";    

        ListData <<format2<<" ]"<< std::endl;         
   
        
        } 
        

            
        }
        
                if(cellIDs == 13){        
            
            for(int i=0; i < proximity_list_1_3.size(); i++)
        {         
            Proximity * prox = proximity_list_1_3.at(i);
         AxonID = prox->Axon_ID;
         DendriteID = prox->Dendrite_ID;
          
         std::string format2 = "[ load ${SCRIPTDIR}/proximity";
         
        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_2 || AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#2";

        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#3";

        if(AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_2)
        format2 += "_CELL#2-->CELL#3";

         
        
        format2 += "_roi";
        format2 += Utility::makeString(i);
        format2 += "/loadROI.hx";    

        ListData <<format2<<" ]"<< std::endl;         
   
        
        } 
        

            
        }
        
                if(cellIDs == 23){        
            
            for(int i=0; i < proximity_list_2_3.size(); i++)
        {         
            Proximity * prox = proximity_list_2_3.at(i);
         AxonID = prox->Axon_ID;
         DendriteID = prox->Dendrite_ID;
          
         std::string format2 = "[ load ${SCRIPTDIR}/proximity";
         
        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_2 || AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#2";

        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_1)
        format2 += "_CELL#1-->CELL#3";

        if(AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_2)
        format2 += "_CELL#2-->CELL#3";

         
        
        format2 += "_roi";
        format2 += Utility::makeString(i);
        format2 += "/loadROI.hx";    

        ListData <<format2<<" ]"<< std::endl;         
   
        
        }
        

            
        }
               
        
        
      
       
       
     std::string format3 = "[ load ${SCRIPTDIR}/proximity";   
        
     if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_2 || AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_1)
        format3 += "_CELL#1-->CELL#2";
//         if(AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_1)
//         format3 += "_CELL#2-->CELL#1";
        if(AxonID == AXON_ID_1 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_1)
        format3 += "_CELL#1-->CELL#3";
//         if(AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_1)
//         format3 += "_CELL#3-->CELL#1";
        if(AxonID == AXON_ID_2 && DendriteID == DENDRITE_ID_3 || AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_2)
        format3 += "_CELL#2-->CELL#3";
/*        if(AxonID == AXON_ID_3 && DendriteID == DENDRITE_ID_2)
        format3 += "_CELL#3-->CELL#2";  */   
        
     ListData <<format3<<"_locations.landmarkAscii ]"<< std::endl;
     ListData <<"[ load ${SCRIPTDIR}/../../" << AmiraMergedFileString.c_str() << " ]"<< std::endl;
     if(cellIDs == 12){
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_1.am ]"<< std::endl;
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_2.am ]"<< std::endl;
     }
          if(cellIDs == 13){
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_1.am ]"<< std::endl;
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_3.am ]"<< std::endl;
     }
          if(cellIDs == 23){
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_2.am ]"<< std::endl;
     ListData <<"[ load ${SCRIPTDIR}/proximitycell_3.am ]"<< std::endl;
     }
     //ListData <<"[ load "<< <<" ]"<< std::endl;
    
     ListData << "       "<<std::endl;
        
        
        
        
     ListData.close();
    
    
}

void AxonDendriteProximityFinder::writeSingleLandmarkFile(double * coord, const char* outputFilename, int number)
{
//  //////std::cout<< "in writeLandmarkFile" << std::endl;
  
  std::string format = outputFilename;
  format += "_roi";
  format += Utility::makeString(number);
  format += ".landmarkAscii";
  

  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << 1     << std::endl;
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
  
  
    LandMarkData << coord[X_COORD] << " " << coord[Y_COORD] << " " << coord[Z_COORD] << std::endl;
  
  
  LandMarkData.close();
};


void AxonDendriteProximityFinder::writeInfoFile( const char* outputFilename, Proximity * prox, int i)
{
    //writes an InfoFile to ensure right positioning of the proximity images in Amira
    std::string format = outputFilename;
    format += Utility::makeString(i);
    format += ".info";
    
    ////std::cout<<"info Filename: "<<format<<std::endl;
    
    std::ofstream ListData( format.c_str() );
    
    ImageType::SizeType region_size;
    ImageType::Pointer input_image = prox->getProximityImage();
    ImageType::RegionType region = input_image->GetLargestPossibleRegion();
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
        
        std::string imagename = "proximity";
           
           
        
        ListData <<" "<<imagename<<std::setfill('0')<< std::setw(3)<<i<<".tif"<<" "<<first_slice + 0.5*i <<std::endl;
    }
    
    ListData.close();  
    
//     for(int i=0; i < list->size(); i++)
//     {
// 
//          
//         
//     Proximity * prox = proximity_list.at(i);
//     int Axon_ID = prox->Axon_ID;
//     if(Axon_ID == AXON_ID_1)
//         Axon_ID = 1;
//     if(Axon_ID == AXON_ID_2)
//         Axon_ID = 2;
//     if(Axon_ID == AXON_ID_3)
//         Axon_ID = 3;
//     int Dendrite_ID = prox->Dendrite_ID;
//     if(Dendrite_ID == DENDRITE_ID_1)
//         Dendrite_ID = 1;
//     if(Dendrite_ID == DENDRITE_ID_2)
//         Dendrite_ID = 2;
//     if(Dendrite_ID == DENDRITE_ID_3)
//         Dendrite_ID = 3;
//     
//     std::string format = outputFilename;
//     format += "_roi";
//     format += Utility::makeString(i);
//     format += "/proximity";
//     format += "_Axon_";
//     format +=  Utility::makeString(Axon_ID);
//     format += "_Dendrite_";
//     format += Utility::makeString(Dendrite_ID);
//     format += "_";
//     format += Utility::makeString(i);
//     format += ".info";
//     
//     //////std::cout<<"info Filename: "<<format<<std::endl;
//     
//     std::ofstream ListData( format.c_str() );
//     
//     ImageType::SizeType region_size;
//     ImageType::Pointer input_image = prox->getProximityImage();
//     ImageType::RegionType region = input_image->GetLargestPossibleRegion();
//     double pixel_size = XYSAMPLING;
//     
//     region_size = region.GetSize();
//     ImageType::IndexType start = region.GetIndex();
//     
// 
//     
//     const unsigned int first_slice = 0;
//     const unsigned int last_slice  = region_size[2] -1;
//     
//     ListData << "# Amira Stacked Slices"<< std::endl;
//     ListData << " "<< std::endl;
//     ListData << "pixelsize "<<pixel_size<<"  "<<pixel_size<< std::endl;
//     ListData << "       "<<std::endl;
//     
//     for(int i = 0 ; i <= last_slice; i++)
//     {
//         
//         std::string imagename = "proximity";
//            
//            
//         
//         ListData <<" "<<imagename<<std::setfill('0')<< std::setw(3)<<i<<".tif"<<" "<<first_slice + 0.5*i <<std::endl;
//     }
//     
//     ListData.close();  
//     
//     }
    
};

#endif
