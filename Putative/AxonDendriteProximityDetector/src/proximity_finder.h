/****************************************************************************/
/*                                                                          */
/* File:      proximity_finder.h                                            */
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
/* History:   26.03.2014                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


#ifndef PROXIMITYFINDER
#define PROXIMITYFINDER

#include "edge_wrapper.h"
#include "proximity.h"
#include "amiraReader.h"
#include "brick_image_reader.h"
#include "typedefs.h"
#include "typedefs_two.h"


//#define DEBUG

#define MAX_DIST        3

#define AXON_ID         5
#define DENDRITE_ID_1   6
#define UNKNOWN_ID      7
#define DENDRITE_ID_2   8

typedef EdgeWrapper SimpleAxon;
typedef EdgeWrapper SimpleDendrite;

//****************************

class ProximityFinder
{
 public:
     ProximityFinder(AmiraSpatialGraph * input_graph, ImageType::Pointer input_image, const char * filename, int removeDoubledProximites, int numOfPointsAway, double max_dist )
     {
         set_edge_lists(input_graph);
         amira_graph = input_graph;
         original_image = input_image;
         original_image->Update();
         this->max_dist = max_dist;
         this->inputfilename = filename;
         this->eliminateDoubledProximites = int(removeDoubledProximites);
         this->numberOfPointsAway = numOfPointsAway;
         
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
         
         //std::cout<<"maxdist: "<<max_dist<<std::endl;
     };
     ProximityFinder(AmiraSpatialGraph * input_graph, ImageType::Pointer input_image, const char * filename, int removeDoubledProximites, int numOfPointsAway)
     {
         set_edge_lists(input_graph);
         original_image = input_image;
         original_image->Update();
         amira_graph = input_graph;
         this->inputfilename = filename;
         this->eliminateDoubledProximites = int(removeDoubledProximites);
         this->numberOfPointsAway = numOfPointsAway;
        
         this->max_dist = MAX_DIST;
         
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
     };
       
     void findProximities();
     void writeProximityLandmarks(char * outputFilename);
     
     bool readAmiraTransformations();

     std::vector<double *> * transformCoordinates(std::vector<double *> * list);
     
     void writeProximityImages(char * outputFilename)
     {
            for(int i=0; i<proximity_list.size(); i++)
            {
                    proximity_list.at(i)->writeProximityImage(outputFilename, i);
            }
     };
     
     void markLocalMaximums();
     float averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection);
     bool markBoutons(double * axon_point);     
     float calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius);
     float calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius);
     float calculateThreshold(ImageType::Pointer image, float x0, float y0, float z0);
     float setAverageThreshold(int numOfPointsToAverage);
     int setRadiusAndLocalBrightness(ImageType::Pointer image);
      int addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance );
      Image2DType::Pointer getImagePlane(int z, int depth, ImageType::Pointer input_image);
      VECTOR * sendRay(float x0, float y0, float ray_length, unsigned int nr_of_rays, float threshold, unsigned int n, Image2DType::Pointer image_plane);
      float bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane);
      void smoothRadii();
      bool isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance);
      float getAverageStandardDeviation();
//      ImageType::Pointer getOriginalImage(char* file_name, int start_index, int end_index);
//      
     
;
     


  
 private:
     std::vector<SimpleAxon> axon_list;
     std::vector<SimpleDendrite> dendrite_list;
     std::vector<Proximity *> proximity_list;
     ImageType::Pointer original_image;
     float max_dist;
     AmiraSpatialGraph * amira_graph;
     float ** sectionTranslation;
     float ** sectionRotation;
     double ** transformation;
     const char * inputfilename;
     bool eliminateDoubledProximites;
     int numberOfPointsAway;
     


     
    
     void writeLandmarkFile(std::vector<double*> * list, char* outputFilename);
     void writeListFile(std::vector<double*> * list, char* outputFilename);

     void set_edge_lists(AmiraSpatialGraph * input_graph)
     {
        int count = 0; 
        std::vector< Edge * > * edges = input_graph->edgesPointer();
        int numOfEdges = edges->size();
                
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
                Edge * currentEdge = edges->at(i);
                
                if(currentEdge->label == DENDRITE_ID_1 || currentEdge->label == DENDRITE_ID_2)
                {
                    dendrite_list.push_back(SimpleDendrite(currentEdge, count));
                }
                else if(currentEdge->label == UNKNOWN_ID)
                {
                    axon_list.push_back(SimpleAxon(currentEdge, count));
                    dendrite_list.push_back(SimpleDendrite(currentEdge, count));
                }
                else if(currentEdge->label == AXON_ID)
                {
                    axon_list.push_back(SimpleAxon(currentEdge, count));
                }
                count++;
        }
     };
     
};

  


#endif

