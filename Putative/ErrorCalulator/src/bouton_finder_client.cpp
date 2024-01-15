
//STL includes
#include "string"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <complex>
#include <utility>
#include <ctime>

//GSL includes
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>

#include "result_analyzer.h"

using namespace std;

double euclidean_distance(double * a, double * b, int dim)
{
    double dist = 0;
    double diff = 0;
    
    for(int i=0; i< dim; i++)
    {
	diff = a[i]-b[i];
	dist += diff*diff;
    }
    return sqrt(dist);
};

void readLandmarkfile(const char* filename, list<double *>  *landmarks)
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
    //////\b//std::cout\b<<"valid path"<<std::endl;
    while(!std::getline(inputStream, currentLine).eof() && currentLine.find("@1") != 0 ){}
    while(!std::getline(inputStream, currentLine).eof() && currentLine.size() > 0 )
    { 
      double * tmp = new double[3];
      
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[0], &tmp[1], &tmp[2]);
      
      landmarks->push_back(tmp);
      
      
    }
    //////\b//std::cout\b<<"done reading prox landmarks"<<std::endl;
  }

};

void writeLandmarkListFile( const char* outputFilename, list<double *> *listinImageCoords )
{
 
  std::ofstream LandMarkData;
  
  LandMarkData.open (outputFilename, std::ofstream::out);
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << listinImageCoords->size()    << std::endl;
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
  
  std::list<double*>::iterator list_it1;
  for(list_it1 = listinImageCoords->begin() ; list_it1 != listinImageCoords->end(); list_it1++)
  {
    
    
    LandMarkData <<  (*list_it1)[0] << " " <<  (*list_it1)[1] << " " <<  (*list_it1)[2] << std::endl;
    //cout <<  (*list_it1)[0] << " " <<  (*list_it1)[1] << " " <<  (*list_it1)[2] << std::endl;  
  }

  LandMarkData.close();
};




int main( int argc , char * argv[])
{
  
  char* manual_landmark_file = argv[1];
  char* auto_landmark_file = argv[2];
  char* error_landmark_path = argv[3];
  char* error_stat_path = argv[4];
  char* outputFilename = argv[5];
  
  list<double *> manual_landmarks;
  list<double *> auto_landmarks;
  
  cout<< "got in "<< endl;
  readLandmarkfile(manual_landmark_file, &manual_landmarks);
  readLandmarkfile(auto_landmark_file, &auto_landmarks);
  
  std::list<double *>::iterator it;
  cout << "manual "<< endl;
  for(it = manual_landmarks.begin() ; it != manual_landmarks.end(); it++)
  {
    cout << (*it)[0] << " "<< (*it)[1] << " "<< (*it)[2] << endl;
  }
  
  cout << "auto "<< endl;
  for(it = auto_landmarks.begin() ; it != auto_landmarks.end(); it++)
  {
    cout << (*it)[0] << " "<< (*it)[1] << " "<< (*it)[2] << endl;
  }
    
 ResultAnalyzer * result_analyzer = new ResultAnalyzer(&manual_landmarks, &auto_landmarks, error_landmark_path, error_stat_path, outputFilename);
        
 delete result_analyzer;
  
  
  
  return 0;

//#endif
};
