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
#include "bouton_cluster.h"


//extern bouton_params BoutonParamsArray[];

ResultAnalyzer::ResultAnalyzer( int section_num, int roi_num, int axon_num, int dend_num, const char* bouton_landmark_file,const char* overlap_landmark_file, const char* outputErrorLandmarkPathname, const char* outputErrorStatsPathname, std::list<double *>*autoBoutons, std::list<double *>*autoOverlap)
{
  ////////////std:://cout<<"in result main"<<std::endl;
  //set the manual landmarks;
  std::list<double *>manualBoutons;
  std::list<double *>manualOverlap;
  
  readLandmarkfile(bouton_landmark_file,  &manualBoutons );

  readLandmarkfile(overlap_landmark_file, &manualOverlap );   
  
  //bouton_cluster * boutonclusterobj = new bouton_cluster();
  // compare this with the automatic landmarks
  compareBoutonLandmarks(autoBoutons, &manualBoutons);
  
  compareOverlapLandmarks(autoOverlap, &manualOverlap);
  
  writeBoutonErrorLandmarks(outputErrorLandmarkPathname);
  
  writeOverlapErrorLandmarks(outputErrorLandmarkPathname);
  
  writeStats(outputErrorStatsPathname, section_num, roi_num, axon_num, dend_num);
  
  
};

/*******************************************************************************************/
/*readLandmarkfile()  populates a list of bouton coordinates                               */
/*******************************************************************************************/
void ResultAnalyzer::readLandmarkfile(const char* filepath, std::list<double *>*manual_landmarks)
{
    
  ////////////std:://cout<<"reading"<<std::endl;
    ////////////std:://cout<<filepath<<std::endl;
//   
//   std::string format = filepath;
//   format += "/Landmarks.landmarkAscii";
//   ////////////std:://cout<<format<<std::endl;
  std::ifstream inputStream(filepath);
  std::string currentLine = "";
  int count = 0;
  
  //////////////std:://cout<<format<<std::endl;
  if(!inputStream.fail())
  {
    ////////////std:://cout<<"valid path"<<std::endl;
    while(!std::getline(inputStream, currentLine).eof() && currentLine.find("@1") != 0 ){}
    while(!std::getline(inputStream, currentLine).eof() && currentLine.size() > 0 )
    { 
      double * tmp = new double[ARRAY_LENGTH];
      
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[X_COORD], &tmp[Y_COORD], &tmp[Z_COORD]);
      
      manual_landmarks->push_back(tmp);
      
    }
  }

};

/*******************************************************************************************/
/*readLandmarkfile()  populates a list of bouton coordinates                               */
/*******************************************************************************************/
void ResultAnalyzer::compareBoutonLandmarks(std::list<double*>*autoLandmarks, std::list<double*>*manualLandmarks)
{
    ////////////std:://cout<<"comparing"<<std::endl;
    std::list<double*>::iterator manual_it;
    std::list<double*>::iterator auto_it;
    

    // copy the manual and auto landmarks first to error lists.. then remove the matching entries
    for(manual_it = manual_landmarks.begin(); manual_it != manual_landmarks.end(); manual_it++)
    {
        ////////////std:://cout<<"in initializing loop manual"<<std::endl;
        ////////////std:://cout<< *manual_it<<std::endl;
       
        fn_landmarks_1.push_back(*manual_it);
        fn_landmarks_2.push_back(*manual_it);
        fn_landmarks_3.push_back(*manual_it);
        fn_landmarks_4.push_back(*manual_it);
        fn_landmarks_5.push_back(*manual_it);
        
    }
    
    for(auto_it = autoLandmarks->begin(); auto_it != autoLandmarks->end(); auto_it++)
    {
        ////////////std:://cout<<"in initializing loop auto"<<std::endl;
         ////////////std:://cout<< *auto_it<<std::endl;
        fp_landmarks_1.push_back(*auto_it);
        fp_landmarks_2.push_back(*auto_it);
        fp_landmarks_3.push_back(*auto_it);
        fp_landmarks_4.push_back(*auto_it);
        fp_landmarks_5.push_back(*auto_it);
    }
    
    
    // iterate through the list of manul landmarks and see if there is a auto bouton nearby 
    // for a distance measure of 1... 5 microns.. if so remove the matching ones from the error lists
   
   //for(int errorresol = 1; errorresol <= 5; errorresol++ )
   {
      removeMatchingLandmarks(&fn_landmarks_1, &fp_landmarks_1, 1);
      removeMatchingLandmarks(&fn_landmarks_2, &fp_landmarks_2, 2);
      removeMatchingLandmarks(&fn_landmarks_3, &fp_landmarks_3, 3);
      removeMatchingLandmarks(&fn_landmarks_4, &fp_landmarks_4, 4);
      removeMatchingLandmarks(&fn_landmarks_5, &fp_landmarks_5, 5);
      
     
   }
   
    
};

/*******************************************************************************************/
/*readLandmarkfile()  populates a list of bouton coordinates                               */
/*******************************************************************************************/
void ResultAnalyzer::compareOverlapLandmarks(std::list<double*>*autoLandmarks, std::list<double*>*manualLandmarks)
{
    ////////////std:://cout<<"comparing"<<std::endl;
    std::list<double*>::iterator manual_it;
    std::list<double*>::iterator auto_it;
    

    // copy the manual and auto landmarks first to error lists.. then remove the matching entries
    for(manual_it = manual_landmarks.begin(); manual_it != manual_landmarks.end(); manual_it++)
    {
        ////////////std:://cout<<"in initializing loop manual"<<std::endl;
        ////////////std:://cout<< *manual_it<<std::endl;
       
        overlap_fn_landmarks_1.push_back(*manual_it);
        overlap_fn_landmarks_2.push_back(*manual_it);
        overlap_fn_landmarks_3.push_back(*manual_it);
        overlap_fn_landmarks_4.push_back(*manual_it);
        overlap_fn_landmarks_5.push_back(*manual_it);
        
    }
    
    for(auto_it = autoLandmarks->begin(); auto_it != autoLandmarks->end(); auto_it++)
    {
        ////////////std:://cout<<"in initializing loop auto"<<std::endl;
         ////////////std:://cout<< *auto_it<<std::endl;
        overlap_fp_landmarks_1.push_back(*auto_it);
        overlap_fp_landmarks_2.push_back(*auto_it);
        overlap_fp_landmarks_3.push_back(*auto_it);
        overlap_fp_landmarks_4.push_back(*auto_it);
        overlap_fp_landmarks_5.push_back(*auto_it);
    }
    
    
    // iterate through the list of manul landmarks and see if there is a auto bouton nearby 
    // for a distance measure of 1... 5 microns.. if so remove the matching ones from the error lists
   
   //for(int errorresol = 1; errorresol <= 5; errorresol++ )
   {
      removeMatchingLandmarks(&overlap_fn_landmarks_1, &overlap_fp_landmarks_1, 1);
      removeMatchingLandmarks(&overlap_fn_landmarks_2, &overlap_fp_landmarks_2, 2);
      removeMatchingLandmarks(&overlap_fn_landmarks_3, &overlap_fp_landmarks_3, 3);
      removeMatchingLandmarks(&overlap_fn_landmarks_4, &overlap_fp_landmarks_4, 4);
      removeMatchingLandmarks(&overlap_fn_landmarks_5, &overlap_fp_landmarks_5, 5);
      
     
   }
   
    
};

/*******************************************************************************************/
/*readLandmarkfile()  populates a list of bouton coordinates                               */
/*******************************************************************************************/
void ResultAnalyzer::removeMatchingLandmarks(std::list<double*>*manual_landmarks, std::list<double*>*autoLandmarks, int distMeasure)
{
    std::list<double*>::iterator manual_it;
    std::list<double*>::iterator auto_it;
    
    for(manual_it = manual_landmarks->begin(); manual_it != manual_landmarks->end(); /*manual_it++*/)
    {
        double min_dist = 0xffff;
        
	std::list<double*>::iterator min_auto_it;
        for(auto_it = autoLandmarks->begin(); auto_it != autoLandmarks->end(); auto_it++)
        {
            
            // look for a match under each error resolution
            double dist = Utility::euclidean_distance(*manual_it, *auto_it,3,1);
            
            if(dist < min_dist)
            {
                min_dist = dist;
		min_auto_it = auto_it;
		
            }
            
        }
        
        ////////////std:://cout<<"initialized"<<std::endl;
        
        if(min_dist < distMeasure)
        {
            manual_it = manual_landmarks->erase(manual_it);
            autoLandmarks->erase(min_auto_it);
            
        }
        else
	{
	    manual_it++;
	}
        
       
    }
  
  
};




/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::writeBoutonErrorLandmarks(const char* outputPathname)
{
    ////////////std:://cout<<"in writeErrorStatsandLandmarks"<<std::endl;
    
    for(int i= 1; i <= 5; i++)
    {
        if(i == 1)
        {
            writeErrorLandmarks(outputPathname, "bouton_fn_landmarks_1", fn_landmarks_1);
            writeErrorLandmarks(outputPathname, "bouton_fp_landmarks_1", fp_landmarks_1); 
        }
        if(i == 2)
        {
            writeErrorLandmarks(outputPathname, "bouton_fn_landmarks_2", fn_landmarks_2);
            writeErrorLandmarks(outputPathname, "bouton_fp_landmarks_2", fp_landmarks_2); 
        }
        if(i == 3)
        {
            writeErrorLandmarks(outputPathname, "bouton_fn_landmarks_3", fn_landmarks_3);
            writeErrorLandmarks(outputPathname, "bouton_fp_landmarks_3", fp_landmarks_3);
        }
        if(i == 4)
        {
            writeErrorLandmarks(outputPathname, "bouton_fn_landmarks_4", fn_landmarks_4);
            writeErrorLandmarks(outputPathname, "bouton_fp_landmarks_4", fp_landmarks_4);
        }
        if(i == 5)
        {
            writeErrorLandmarks(outputPathname, "bouton_fn_landmarks_5", fn_landmarks_5);
            writeErrorLandmarks(outputPathname, "bouton_fp_landmarks_5", fp_landmarks_5);
        }
        
    }
    
    
};

/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::writeOverlapErrorLandmarks(const char* outputPathname)
{
    ////////////std:://cout<<"in writeErrorStatsandLandmarks"<<std::endl;
    
    for(int i= 1; i <= 5; i++)
    {
        if(i == 1)
        {
            writeErrorLandmarks(outputPathname, "overlap_fn_landmarks_1", overlap_fn_landmarks_1);
            writeErrorLandmarks(outputPathname, "overlap_fp_landmarks_1", overlap_fp_landmarks_1); 
        }
        if(i == 2)
        {
            writeErrorLandmarks(outputPathname, "overlap_fn_landmarks_2", overlap_fn_landmarks_2);
            writeErrorLandmarks(outputPathname, "overlap_fp_landmarks_2", overlap_fp_landmarks_2); 
        }
        if(i == 3)
        {
            writeErrorLandmarks(outputPathname, "overlap_fn_landmarks_3", overlap_fn_landmarks_3);
            writeErrorLandmarks(outputPathname, "overlap_fp_landmarks_3", overlap_fp_landmarks_3);
        }
        if(i == 4)
        {
            writeErrorLandmarks(outputPathname, "overlap_fn_landmarks_4", overlap_fn_landmarks_4);
            writeErrorLandmarks(outputPathname, "overlap_fp_landmarks_4", overlap_fp_landmarks_4);
        }
        if(i == 5)
        {
            writeErrorLandmarks(outputPathname, "overlap_fn_landmarks_5", overlap_fn_landmarks_5);
            writeErrorLandmarks(outputPathname, "overlap_fp_landmarks_5", overlap_fp_landmarks_5);
        }
        
    }
    
};

/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::writeErrorLandmarks(const char* outputPathname, const char* outputfilename, std::list<double*> error_list)
{
    ////////////std:://cout<<"in writeErrorLandmarks"<<std::endl;
  ////////////////std:://cout<< "in writeLandmarkFile" << std::endl;
  ////////////////std:://cout<< dia_list->size() << std::endl;
  ////////////////std:://cout<< outputPathname << std::endl;
  
// transformToWorldCoordinates(list);
   if(chdir(outputPathname) != 0)
            perror("Couldn't open image drirectory!");
  std::string format = outputfilename;
  format += "_locations.landmarkAscii";
  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << error_list.size()<< std::endl;
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
  
  std::list< double * >::const_iterator diaIt = error_list.begin();
        ////////////////std:://cout<< "in before loop" << std::endl;
        //////////////std:://cout<< imageTranslation[0]<<" " <<imageTranslation[1]<<" " <<imageTranslation[2]<<" " << std::endl;
        while(diaIt != error_list.end())
        {
            //////////////////std:://cout<< "in loop" << std::endl;
                double * nextPt = *diaIt;
                
                {
                    LandMarkData << nextPt[X_COORD] /*+ imageTranslation[0] */ << " " << nextPt[Y_COORD]/*+imageTranslation[1]*/ << " " << nextPt[Z_COORD]/*+imageTranslation[2]*/ << std::endl;
                }
                //////////////std:://cout << nextPt[X_COORD]  << " " << nextPt[Y_COORD] << " " << nextPt[Z_COORD] << std::endl;
                diaIt++;
                
        }
        
  
  LandMarkData.close();

    
};

/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::writeStats(const char* outputPathname, int section_num, int roi_num, int axon_num, int dend_num)
{
    
    ////////////std:://cout<<"in writeStats"<<std::endl;
    
    bool file_exsists = 0;
    std::string s1("Errorstats.xls");
    
    if(chdir(outputPathname) != 0)
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
  
    std::fstream xlsParamWriter;
   
    xlsParamWriter.open ("Errorstats.xls", std::fstream::in | std::fstream::out | std::fstream::app);
    
    // write header if it does not exists
    if(!file_exsists)
    {
        xlsParamWriter<<"Section Number"<<'\t'<<"ROI Number"<<'\t'
        <<"From Cell #"<<'\t'<<"To Cell #"<<'\t'
        <<"Bouton FN Resolution 1"<<'\t'<<"Bouton FP Resolution 1"<<'\t'<<"Bouton FN Resolution 2"<<'\t'<<"Bouton FP Resolution 2"<<'\t'<<"Bouton FN Resolution 3"<<'\t'
        <<"Bouton FP Resolution 3"<<'\t'<<"Bouton FN Resolution 4"<<'\t'<<"Bouton FP Resolution 4"<<'\t'<<"Bouton FN Resolution 5"<<'\t'<<"Bouton FP Resolution 5"<<'\t'
        <<"Putative FN Resolution 1"<<'\t'<<"Putative FP Resolution 1"<<'\t'<<"Putative FN Resolution 2"<<'\t'<<"Putative FP Resolution 2"<<'\t'<<"Putative FN Resolution 3"<<'\t'
        <<"Putative FP Resolution 3"<<'\t'<<"Putative FN Resolution 4"<<'\t'<<"Putative FP Resolution 4"<<'\t'<<"Putative FN Resolution 5"<<'\t'<<"Putative FP Resolution 5"<<'\n';
    }
       
    xlsParamWriter<<section_num<<'\t'<<roi_num<<'\t'
    <<axon_num<<'\t'<<dend_num<<'\t'
    <<fn_landmarks_1.size()<<'\t'<<fp_landmarks_1.size()<<'\t'<<fn_landmarks_2.size()<<'\t'<<fp_landmarks_2.size()<<'\t'
    <<fn_landmarks_3.size()<<'\t'<<fp_landmarks_3.size()<<'\t'<<fn_landmarks_4.size()<<'\t'<<fp_landmarks_4.size()<<'\t'
    <<fn_landmarks_5.size()<<'\t'<<fp_landmarks_5.size()<<'\t'
    <<overlap_fn_landmarks_1.size()<<'\t'<<overlap_fp_landmarks_1.size()<<'\t'<<overlap_fn_landmarks_2.size()<<'\t'<<overlap_fp_landmarks_2.size()<<'\t'
    <<overlap_fn_landmarks_3.size()<<'\t'<<overlap_fp_landmarks_3.size()<<'\t'<<overlap_fn_landmarks_4.size()<<'\t'<<overlap_fp_landmarks_4.size()<<'\t'
    <<overlap_fn_landmarks_5.size()<<'\t'<<overlap_fp_landmarks_5.size()<<'\n';
    
    xlsParamWriter.close();


};









