/****************************************************************************/
/*                                                                          */
/* File:      bouton_finder.h                                               */
/*                                                                          */
/* Purpose:                                                                 */
/*                                                                          */
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
#include "typedefs.h"
#include "utility.h"

//#include "typedefs_two.h"
#include "amiraReader.h"
#include "bouton_params.h"


#ifndef BOUTONFINDER
#define BOUTONFINDER


#define DEBUG

//****************************
//#define ARRAY_LENGTH 16
#define CORRESPONDANCE 3
#define MAX_DIST 30

#define Z_SCALER 1

//extern std::list< bouton_params * > BoutonParamsList;

class BoutonFinder 
{
 public:
  BoutonFinder(const char * input, bool type);
  BoutonFinder(char* amira_file, char* image_folder_name, int startindex, int endindex, char* loadhxfile, char* outputname);
  //double total_axon_length;

  void findBoutons(void);
  void findSpatialGraphPoint(void);
  void findOrientationOfAxon(const double * landmarkPoint, int numNieghbours, double*);
  void findBrightness(Image2DType::Pointer image_plane, double * ContourList,  std::list<double*>cont_min_rad, double *  BrightnessBB, double *  BrightnessContour, double *  BrightnessMinContour, std::list<double *> *brightness_points_list);
  //void  findRadius(Image2DType::Pointer, double * landmarkPoint,  double *orientationlist,double *, double *);  
  //std::list <double *>  findCountour(Image2DType::Pointer, double * landmarkPoint, std::list< double * > orientationlist );  
  Image2DType::Pointer getImagePlane(int z, int depth, ImageType::Pointer input_image);
  void transformGraphandLandmark(void);
  void getNonBoutonPoints(ImageType::Pointer image);
  static double getEuclideandistance(double* tmp, double* neighbor)
  {
    double temp_x = 0, temp_y = 0, temp_z = 0, distance = 0;

    temp_x = tmp[X_COORD] - neighbor[X_COORD];
    temp_y = tmp[Y_COORD] - neighbor[Y_COORD];
    temp_z = tmp[Z_COORD] - neighbor[Z_COORD];
    
    temp_x *= temp_x;
    temp_y *= temp_y;
    temp_z *= temp_z;

    distance = sqrt(temp_x + temp_y + Z_SCALER*temp_z);

    return distance;
  };
  
  ImageType::Pointer original_image;

  
 //private:
  Reader * amiraReader;
  AmiraSpatialGraph* amira_graph;
  
  
  double  sectionTranslation[4][4];// = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double  sectionRotation[4][4];// =  {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double  transformation[4][4];//{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double  inverse_transformation[4][4];//= {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  double imageTranslation[3];// = {0,0,0};
  PixelType OtsuThreshold; 
  const char * outputFilename;
  char * loadhxfile;
  
   
  
  
  ImageType::Pointer getOriginalImage(char* file_name, int startindex, int endindex);
  bool getTransformationandInverseTransformation(char* inputfilename);
  bool getTranslation(char* loadhxfile);
  
  void writeTestImage();
  //static void writeImagePlanes(ImageType::Pointer input_image, char output_file[1024]);
  float calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius);
  float calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius);
  //float setAverageThreshold(int numOfPointsToAverage);
  float getAverageStandardDeviation();
  void markBoutons();
  //void markLocalMaximums();
  int addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance );
  void writeBoutonStatistics(std::vector<double*> * detectedBoutons, std::vector<double*> * actualBoutons, std::vector<int> * correctIndicesDetected, std::vector<int> * correctIndicesactual);
  float averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection);
  void writeBoutonLandmarkfile(void);
  void UndoTransformation(void);
  char* landmark_filename;
  double total_axon_length;
  
  
  //int setRadiusAndLocalBrightness(ImageType::Pointer image);
  //Image2DType::Pointer getImagePlane(int z, int depth, ImageType::Pointer input_image);
  float calculateThreshold(ImageType::Pointer image_plane, float x0, float y0, float z0);
  void sendRay(double x0, double y0,  Image2DType::Pointer image_plane, double *, double *, std::list<double*>* ,double);
  float bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane);
  void smoothRadii();
  bool isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance);
  PixelType calculateOtsuThreshold(Image2DType::Pointer image, float x, float y, float z);
  
  
  
  

};

#endif

