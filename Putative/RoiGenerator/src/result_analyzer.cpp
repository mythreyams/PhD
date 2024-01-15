/****************************************************************************/
/*                                                                          */
/* File:      result_analyzer.cpp 					    */
/*                                                                          */
/* Purpose:                                                                 */
/*	                                            			    */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail: marcel.oberlaender@mpimf-heidelberg.mpg.de                        */
/*                                                                          */
/* History:   24.06.2013                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/
#include "result_analyzer.h"



ResultAnalyzer::ResultAnalyzer(AmiraSpatialGraph* processed_graph)
{
  
  amira_graph = processed_graph;

  
};

/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::analyzeProcessedGraph()
{

  std::vector<double*> * detectedBoutons = getSimpleBoutonList();
  writeLandmarkFile(detectedBoutons, "detected_boutons");
  
  std::vector<double*> * actualBoutons = new std::vector<double*>;
 // readLandmarkfile(landmark_filename, actualBoutons);
  
  compareBoutons(detectedBoutons, actualBoutons);
  
  //writeTestImage();

};

/*******************************************************************************************/
/*getSimpleBoutonList()  comverts the amira graph points to a simple coordinate list       */
/*******************************************************************************************/
std::vector<double*> * ResultAnalyzer::getSimpleBoutonList()
{
  std::cout<< "In getSimpleBoutonList" << std::endl;
  
  std::vector<double*> * landmarkList = new std::vector<double *>;
  
  std::vector< Edge * > * edges = amira_graph->edgesPointer();
	unsigned int numOfEdges = edges->size();	
	
	
// #pragma omp parallel for schedule(dynamic,1)
	for(long pos = numOfEdges -1; pos >= 0; pos--)	//for each edge in list
	{			
		Edge * currentEdge = edges->at(pos);
		std::list< double * >::iterator edge_it;	
		unsigned int edgeSize = currentEdge->edgePointCoordinates.size();

		//for every point along edge
		for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
		{
			double * coords = *edge_it;
			
			if(coords[IS_BOUTON])
			{
			  double * tmp = new double[ARRAY_LENGTH];
			  for(int i=0;i<ARRAY_LENGTH;i++)
			    tmp[i] = coords[i];
			  
			  landmarkList->push_back(tmp);
			}
	
		}
	}
  
  return landmarkList;
};

/*******************************************************************************************/
/*writeLandmarkFile()  writes a bouton list into a valid landmarkAscii file                */
/*******************************************************************************************/
void ResultAnalyzer::writeLandmarkFile(std::vector<double*> * list, char* outputFilename)
{
  std::cout<< "in writeLandmarkFile" << std::endl;
  
  std::string format = outputFilename;
  format += ".landmarkAscii";
  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0" 		<< std::endl;
  LandMarkData << "" 					<< std::endl;
  LandMarkData << "" 					<< std::endl;
  LandMarkData << "define Markers " << list->size() 	<< std::endl;
  LandMarkData << "" 					<< std::endl;
  LandMarkData << "Parameters {" 			<< std::endl;
  LandMarkData << "    NumSets 1," 			<< std::endl;
  LandMarkData << "    ContentType \"LandmarkSet\"" 	<< std::endl;
  LandMarkData << "}" 					<< std::endl;
  LandMarkData << "" 					<< std::endl;
  LandMarkData << "Markers { float[3] Coordinates } @1" << std::endl;
  LandMarkData << "" 					<< std::endl;
  LandMarkData << "# Data section follows" 		<< std::endl;
  LandMarkData << "@1" 					<< std::endl;
  
  for(int i=0; i < list->size(); i++)
  {
    LandMarkData << list->at(i)[X_COORD] << " " << list->at(i)[Y_COORD] << " " << list->at(i)[Z_COORD] << std::endl;
  } 
  
};

/*******************************************************************************************/
/*readLandmarkfile()  populates a list of bouton coordinates                               */
/*******************************************************************************************/
void ResultAnalyzer::readLandmarkfile(char* filename, std::vector<double*> * landmarkList)
{
  
  std::ifstream inputStream(filename);
  std::string currentLine = "";
  int count = 0;
  
  if(!inputStream.fail())
  {
    while(!std::getline(inputStream, currentLine).eof() && currentLine.find("@1") != 0 ){}
    while(!std::getline(inputStream, currentLine).eof() && currentLine.size() > 0 )
    { 
      double * tmp = new double[ARRAY_LENGTH];
      
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[X_COORD], &tmp[Y_COORD], &tmp[Z_COORD]);
      
      landmarkList->push_back(tmp);
      
    }
  }

};

/*******************************************************************************************/
/*compareBoutons()  compares the results of program with "ground truth" and gives stats    */
/*******************************************************************************************/
void ResultAnalyzer::compareBoutons(std::vector<double*> * detectedBoutons, std::vector<double*> * actualBoutons)
{
  int detectedIndex = -1, actualIndex = -1;
  double distance = 0, total=0, false_pos=0, false_neg=0, minDist = 100000;
  
  std::vector<int> * correctIndicesDetected = new std::vector<int>;
  std::vector<int> * correctIndicesActual = new std::vector<int>;
  
  for(int i=0; i< detectedBoutons->size();i++)
  {
    minDist = 100000;
    detectedIndex = -1;
    actualIndex = -1;
    
    for(int j = 0; j < actualBoutons->size(); j++)
    {
      distance = BoutonFinder::getEuclideandistance(detectedBoutons->at(i), actualBoutons->at(j));
      
      if(distance < minDist)
      {
	minDist = distance;
	detectedIndex = i;
	actualIndex = j;
      }      
    }
    
    if(minDist <= MAX_DIST)
    {		     
      if(!contains(correctIndicesDetected, detectedIndex) && !contains(correctIndicesActual, actualIndex))
      {
	correctIndicesDetected->push_back(detectedIndex);
	correctIndicesActual->push_back(actualIndex);
      }
      
    }
     
  }
  
  writeBoutonStatistics(detectedBoutons, actualBoutons, correctIndicesDetected, correctIndicesActual);
  
};

/*******************************************************************************************/
/*writeBoutonStatistics()  write the important statistics regarding the boutons to file    */
/*******************************************************************************************/
void ResultAnalyzer::writeBoutonStatistics(std::vector<double*> * detectedBoutons, std::vector<double*> * actualBoutons, std::vector<int> * correctIndicesDetected, std::vector<int> * correctIndicesActual)
{
  float false_pos=0, false_neg=0, total=0, correct=0;
  float correct_sum=0, false_pos_sum=0, false_neg_sum=0;
  float correct_mean=0, false_pos_mean=0, false_neg_mean=0;
  float correct_sigma=0, false_pos_sigma=0, false_neg_sigma=0;
  float error=0, brightness=0, std_dev=0;
  
  //if(correctIndicesDetected->size() != correctIndicesactual->size())
    //std::cout << "ERROR: INDEX LISTS ARE NOT THE SAME SIZE";
    
  
  correct = correctIndicesDetected->size();
  false_pos = detectedBoutons->size() - correct;
  false_neg = actualBoutons->size() - correct;
  total = correct + false_neg + false_pos;
  
  
  std::cout <<"Total boutons found: " << total << std::endl;
  std::cout <<"Correctly identified: " << correct << " = " << (correct/total)*100 << "%" << std::endl;
  std::cout <<"False positives (overcount): " << false_pos << " = " << (false_pos/total)*100 << "%" << std::endl;
  std::cout <<"False negatives (undercount): " << false_neg << " = " << (false_neg/total)*100 << "%" << std::endl<< std::endl;
  
  
  std::cout <<"#true\t#det\t#fpos\t#fneg\tlen" << std::endl;
  //std::cout <<actualBoutons->size()<<"\t"<<detectedBoutons->size()<<"\t"<<false_pos<<"\t"<<false_neg<<"\t"<<total_axon_length << std::endl;
  
  
  std::ofstream BoutonStats( "bouton_statistics.txt" );
  std::ofstream OverallStats("/home/neuromorph/Data_NeuroMorph/BoutonProject/bouton_ROIs/overall_bouton_stats.txt", std::ofstream::out | std::ofstream::app);
  
  BoutonStats <<"Total boutons found: " << total << std::endl;
  BoutonStats <<"Correctly identified: " << correct << " = " << (correct/total)*100 << "%" << std::endl;
  BoutonStats <<"False positives (overcount): " << false_pos << " = " << (false_pos/total)*100 << "%" << std::endl;
  BoutonStats <<"False negatives (undercount): " << false_neg << " = " << (false_neg/total)*100 << "%" << std::endl<<std::endl;
  
  BoutonStats <<"#true\t#det\t#fpos\t#fneg\tlen" << std::endl;
  BoutonStats <<actualBoutons->size()<<"\t"<<detectedBoutons->size()<<"\t"<<false_pos<<"\t"<<false_neg<<"\t"<<total_axon_length << std::endl<<std::endl;
  
  OverallStats <<actualBoutons->size()<<"\t"<<detectedBoutons->size()<<"\t"<<false_pos<<"\t"<<false_neg<<"\t"<<total_axon_length << std::endl;
  
  //*****************Calculate means*******************************
  for(int i=0; i<detectedBoutons->size();i++)
  { 
    if(contains(correctIndicesDetected, i))
      correct_sum += detectedBoutons->at(i)[LOCAL_BRIGHTNESS];
    else
      false_pos_sum += detectedBoutons->at(i)[LOCAL_BRIGHTNESS];    
  }
  false_pos_mean = false_pos_sum / false_pos;
  correct_mean = correct_sum / correct;
  
  for(int i=0; i<actualBoutons->size();i++)
  { 
    if(contains(correctIndicesActual, i))
    {
      
    }
    else
      false_neg_sum += actualBoutons->at(i)[LOCAL_BRIGHTNESS]; 
  }
  false_neg_mean = false_neg_sum / false_neg;
  
  //****************Calculate standard devs************************
//   correct_sum=0;
//   for(int i=0; i<correctDetectedBoutons->size();i++)
//   { 
//     brightness = correctDetectedBoutons->at(i)[LOCAL_BRIGHTNESS];
//     std_dev = correctDetectedBoutons->at(i)[LOCAL_SIGMA];
//     error = brightness - correct_mean;
//     correct_sum += error*error;    
//     
//     BoutonStats << i << " - local_brightness: " << brightness << "\tlocal_sigma: " << std_dev << std::endl;
//   }
//   correct_sigma = sqrt( correct_sum / correct );
//   BoutonStats << "Correctly detected boutons -> mean: " << correct_mean << " standard deviation: " << correct_sigma << std::endl<<std::endl;
  
  
  false_pos_sum=0;
  for(int i=0; i<detectedBoutons->size();i++)
  { 
    brightness = detectedBoutons->at(i)[LOCAL_BRIGHTNESS];
    std_dev = detectedBoutons->at(i)[LOCAL_SIGMA];
    error = brightness - false_pos_mean;
    false_pos_sum += error*error;    
    
    BoutonStats << i << " - local_brightness: " << brightness << "\tlocal_sigma: " << std_dev << std::endl;
  }
  false_pos_sigma = sqrt( false_pos_sum / false_pos );
  BoutonStats << "False positive boutons -> mean: " << false_pos_mean << " standard deviation: " << false_pos_sigma << std::endl<<std::endl;
  
  
  false_neg_sum=0;
  for(int i=0; i<actualBoutons->size();i++)
  { 
    brightness = actualBoutons->at(i)[LOCAL_BRIGHTNESS];
    std_dev = actualBoutons->at(i)[LOCAL_SIGMA];
    error = brightness - false_neg_mean;
    false_neg_sum += error*error;    
    
    BoutonStats << i << " - local_brightness: " << brightness << "\tlocal_sigma: " << std_dev << std::endl;
  }
  false_neg_sigma = sqrt( false_neg_sum / false_neg );
  BoutonStats << "False negative boutons -> mean: " << false_neg_mean << " standard deviation: " << false_neg_sigma << std::endl<<std::endl;


};








