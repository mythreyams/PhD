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



#ifndef RESULTANALYZER
#define RESULTANALYZER
#include "string"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <math.h>


using namespace std;
//#define DEBUG

//****************************

class ResultAnalyzer
{
 public:
   
   ResultAnalyzer( list<double *>*manual_landmarks, list<double *>*auto_landmarks, const char * outputErrorLandmarkPathname, const char * outputErrorStatsPathname, const char* outputfilename);
  
 private:
  
  //std::list<double*> manual_landmarks;
  std::list<double*> fp_landmarks_1;
  std::list<double*> fp_landmarks_2;
  std::list<double*> fp_landmarks_3;
  std::list<double*> fp_landmarks_4;
  std::list<double*> fp_landmarks_5;
  std::list<double*> fn_landmarks_1;
  std::list<double*> fn_landmarks_2;
  std::list<double*> fn_landmarks_3;
  std::list<double*> fn_landmarks_4;
  std::list<double*> fn_landmarks_5;
  
  
  
  void compareBoutonLandmarks(std::list<double*>*autoLandmarks,std::list<double*>*manualLandmarks);
  
  void writeAllErrorLandmarks(const char* outputPathname, const char* outputfilename);
  
  void writeErrorLandmarks(const char* outputPathname, const char* outputfilename, std::list<double*> error_list);
  void writeStats(const char* outputPathname);
  
  void removeMatchingLandmarks(std::list<double*>*manual_landmarks, std::list<double*>*autoLandmarks, int distMeasure);
  double euclidean_distance(double * a, double * b, int dim, float z_scalar);
  char * inttochar(int i);
  
};

#endif

