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
#include "bouton_params.h"


#ifndef RESULTANALYZER
#define RESULTANALYZER


//#define DEBUG

//****************************

class ResultAnalyzer
{
 public:
   
   ResultAnalyzer::ResultAnalyzer( int sectin_num, int roi_num, int axon_num, int dend_num, const char* bouton_landmark_file,const char* overlap_landmark_file, const char* outputErrorLandmarkPathname, const char* outputErrorStatsPathname, std::list<double *>*autoBoutons, std::list<double *>*autoOverlap);
  
 private:
  
  std::list<double*> manual_landmarks;
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
  
  std::list<double*> overlap_fp_landmarks_1;
  std::list<double*> overlap_fp_landmarks_2;
  std::list<double*> overlap_fp_landmarks_3;
  std::list<double*> overlap_fp_landmarks_4;
  std::list<double*> overlap_fp_landmarks_5;
  std::list<double*> overlap_fn_landmarks_1;
  std::list<double*> overlap_fn_landmarks_2;
  std::list<double*> overlap_fn_landmarks_3;
  std::list<double*> overlap_fn_landmarks_4;
  std::list<double*> overlap_fn_landmarks_5;
  
  void compareBoutonLandmarks(std::list<double*>*autoLandmarks,std::list<double*>*manualLandmarks);
  void compareOverlapLandmarks(std::list<double*>*autoLandmarks,std::list<double*>*manualLandmarks);
  void writeBoutonErrorLandmarks(const char* outputPathname);
  void writeOverlapErrorLandmarks(const char* outputPathname);
  void writeErrorLandmarks(const char* outputPathname, const char* outputfilename, std::list<double*> error_list);
  void writeStats(const char* outputPathname, int,int,int,int);
  void readLandmarkfile(const char* filename, std::list<double*>*);
  void removeMatchingLandmarks(std::list<double*>*manual_landmarks, std::list<double*>*autoLandmarks, int distMeasure);
  
};

#endif

