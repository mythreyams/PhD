/****************************************************************************/
/*                                                                          */
/* File:      result_analyzer.h 					    */
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
#include "bouton_finder.h"


#ifndef RESULTANALYZER
#define RESULTANALYZER


//#define DEBUG

//****************************

class ResultAnalyzer
{
 public:
   ResultAnalyzer();
   ResultAnalyzer(AmiraSpatialGraph* processed_graph);

  
 private:
  AmiraSpatialGraph* amira_graph;
  
  
  std::vector<double*> * getSimpleBoutonList();
  void writeBoutonStatistics(std::vector<double*> * detectedBoutons, std::vector<double*> * actualBoutons, std::vector<int> * correctIndicesDetected, std::vector<int> * correctIndicesactual);
  void writeLandmarkFile(std::vector<double*> * list, char* outputFilename);
  void readLandmarkfile(char* filename, std::vector<double*> * landmarkList);
  void compareBoutons(std::vector<double*> * detectedBoutons, std::vector<double*> * actualBoutons);
  void analyzeProcessedGraph();
  
  
  bool contains(std::vector<int>* list, int index)
  {
    for(int i=0;i<list->size();i++)
    {
	if(list->at(i) == index)
	  return true;
    }
    
    return false;
  };

};

#endif

