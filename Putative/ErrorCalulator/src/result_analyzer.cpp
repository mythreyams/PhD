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



//extern bouton_params BoutonParamsArray[];

ResultAnalyzer::ResultAnalyzer( std::list<double *>*manual_landmarks, std::list<double *>*auto_landmarks, const char * outputErrorLandmarkPathname, const char * outputErrorStatsPathname, const char * outputfilename)
{
  compareBoutonLandmarks(auto_landmarks, manual_landmarks);
  
  writeAllErrorLandmarks(outputErrorLandmarkPathname, outputfilename);
  
  writeStats(outputErrorStatsPathname);
  
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
    for(manual_it = manualLandmarks->begin(); manual_it != manualLandmarks->end(); manual_it++)
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
//         cout<<"in initializing loop auto"<<std::endl;
//         cout<< (*auto_it)[0]<<" "<<(*auto_it)[1]<<" "<<(*auto_it)[2]<<endl;
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
            double dist = euclidean_distance(*manual_it, *auto_it,3,1);
            
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
void ResultAnalyzer::writeAllErrorLandmarks(const char* outputPathname, const char* outputfilename)
{
    ////////////std:://cout<<"in writeErrorStatsandLandmarks"<<std::endl;
    
    for(int i= 1; i <= 5; i++)
    {
      std::string str = outputfilename;
      str += "_";
      str += inttochar(i);
      str += "_";
      
      writeErrorLandmarks(outputPathname, (str+"fn").c_str(), fn_landmarks_1);
      writeErrorLandmarks(outputPathname, (str+"fp").c_str(), fp_landmarks_1);
        /*if(i == 1)
        {
            writeErrorLandmarks(outputPathname, "fn_landmarks_1", fn_landmarks_1);
            writeErrorLandmarks(outputPathname, "fp_landmarks_1", fp_landmarks_1); 
        }
        if(i == 2)
        {
            writeErrorLandmarks(outputPathname, "fn_landmarks_2", fn_landmarks_2);
            writeErrorLandmarks(outputPathname, "fp_landmarks_2", fp_landmarks_2); 
        }
        if(i == 3)
        {
            writeErrorLandmarks(outputPathname, "fn_landmarks_3", fn_landmarks_3);
            writeErrorLandmarks(outputPathname, "fp_landmarks_3", fp_landmarks_3);
        }
        if(i == 4)
        {
            writeErrorLandmarks(outputPathname, "fn_landmarks_4", fn_landmarks_4);
            writeErrorLandmarks(outputPathname, "fp_landmarks_4", fp_landmarks_4);
        }
        if(i == 5)
        {
            writeErrorLandmarks(outputPathname, "fn_landmarks_5", fn_landmarks_5);
            writeErrorLandmarks(outputPathname, "fp_landmarks_5", fp_landmarks_5);
        }*/
        
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
                    LandMarkData << nextPt[0] /*+ imageTranslation[0] */ << " " << nextPt[1]/*+imageTranslation[1]*/ << " " << nextPt[2]/*+imageTranslation[2]*/ << std::endl;
                }
                //////////////std:://cout << nextPt[X_COORD]  << " " << nextPt[Y_COORD] << " " << nextPt[Z_COORD] << std::endl;
                diaIt++;
                
        }
        
  
  LandMarkData.close();

    
};

/****************************************************************************/
/*analyzeProcessedGraph()  primary function for ResultAnalyzer class        */
/****************************************************************************/
void ResultAnalyzer::writeStats(const char* outputPathname )
{
    
    ////////////std:://cout<<"in writeStats"<<std::endl;
    
    bool file_exsists = 0;
    std::string s1("Errorstats.xls");
    
    if(chdir(outputPathname) != 0)
            perror("Couldn't open image drirectory!");
    
   
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
//     if(!file_exsists)
//     {
//         xlsParamWriter<<"FN Resolution 1"<<'\t'<<"FP Resolution 1"<<'\t'<<"FN Resolution 2"<<'\t'<<"FP Resolution 2"<<'\t'<<"FN Resolution 3"<<'\t'
//         <<"FP Resolution 3"<<'\t'<<"FN Resolution 4"<<'\t'<<"FP Resolution 4"<<'\t'<<"FN Resolution 5"<<'\t'<<"FP Resolution 5"<<'\n';
//     }
       
    xlsParamWriter<<fn_landmarks_1.size()<<'\t'<<fp_landmarks_1.size()<<'\t'<<fn_landmarks_2.size()<<'\t'<<fp_landmarks_2.size()<<'\t'
    <<fn_landmarks_3.size()<<'\t'<<fp_landmarks_3.size()<<'\t'<<fn_landmarks_4.size()<<'\t'<<fp_landmarks_4.size()<<'\t'
    <<fn_landmarks_5.size()<<'\t'<<fp_landmarks_5.size()<<'\n';
    
    xlsParamWriter.close();


};


double ResultAnalyzer::euclidean_distance(double * a, double * b, int dim, float z_scalar)
{
    double dist = 0;
    double diff = 0;
    
    for(int i=0; i< dim; i++)
    {
	diff = a[i]-b[i];
	if(i == 2)
	  dist += z_scalar*diff*diff;
	else 
	  dist += diff*diff;
    }
    
    cout<< "dist: "<< sqrt(dist)<< endl;
    
    return sqrt(dist);
};

char * ResultAnalyzer::inttochar(int i)
{
  /* Room for INT_DIGITS digits, - and '\0' */
  static char buf[3 + 2];
  char *p = buf + 3 + 1;       /* points to terminating '\0' */
  if (i >= 0) {
    do {
      *--p = '0' + (i % 10);
      i /= 10;
    } while (i != 0);
    return p;
  }
  else {                        /* i < 0 */
    do {
      *--p = '0' - (i % 10);
      i /= 10;
    } while (i != 0);
    *--p = '-';
  }
  return p;
};






