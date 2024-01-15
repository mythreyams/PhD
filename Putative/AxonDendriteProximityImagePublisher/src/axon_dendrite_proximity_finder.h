/****************************************************************************/
/*                                                                          */
/* File:      proximity_finder.h                                            */
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




#ifndef AXONDENDRITEPROXIMITYFINDER
#define AXONDENDRITEPROXIMITYFINDER

#include "edge_wrapper.h"
#include "proximity.h"
#include "amiraReader.h"
#include "brick_image_reader.h"
#include "typedefs.h"
#include "typedefs_two.h"


//#define DEBUG

#define MAX_DIST        4

// prox regoins with centroids closer than this are merged (in microns)
#define MERGE_PROX_DISTANCE   15



typedef EdgeWrapper SimpleAxon;
typedef EdgeWrapper SimpleDendrite;

//****************************

class AxonDendriteProximityFinder
{
 public:
     
     AxonDendriteProximityFinder(const char* landmarkfilename, const char* xlsinputfilename, const char* xlsoutputfilename, int section_number, const char*deconvolvedImageFileName, int imagestartindex, int imageendindex,  const char * sectiongraphname, const char * mergedgraphname, int include_unknown,const char *outputFilename/*, ImageType::Pointer nondecon_image*/);
     
     AxonDendriteProximityFinder(Reader * ameraReaderPointer, AmiraSpatialGraph * input_graph, AmiraSpatialGraph * merged_graph, ImageType::Pointer input_image, const char * inputFilename, int include_unknown, double max_dist, const char * manual_z_scale,const char* outputFilename);
     
     AxonDendriteProximityFinder(Reader * ameraReaderPointer, AmiraSpatialGraph * input_graph, ImageType::Pointer input_image, const char * inputFilename, int include_unknown, int numOfPointsAway, double max_dist,const char* outputFilename);
     
     AxonDendriteProximityFinder(Reader * ameraReaderPointer, AmiraSpatialGraph * input_graph, ImageType::Pointer input_image, const char * inputFilename,const char* outputFilename);
     
       
     void findProximities(bool regard_boutons);
     void findProximities1(bool regard_boutons);
     void sortProxDist( std::vector<Proximity *> *prox_list_unsorted, std::vector<Proximity *> *prox_list_sorted );
     void findValleys(std::vector<Proximity*>*inputList, std::vector<Proximity*>*outputProxVector,double threshold, double dist_to_merge);
     void writeProxXls(const char* filename, std::vector<Proximity*> *sorted_prox_vector);
     void writeProxLandmarks( const char* outputFilename, std::vector<Proximity*> *sorted_prox_vector);
     void readLandmarkfile(const char* filename, std::vector<double *>*prox_landmarks);
     void applyInverseTransformationToPoint(double * point, double **inverse_transformation);
     //void addSectionNumberToXls(const char * xlsfilename, int row_number, int column_num, int section_number);
     //void writeProximityLandmarks(const char * outputFilename);
     
     bool readAmiraTransformations();
     void removeDoubledProximites();
     void mergreProximites(std::vector<Proximity *>*);
    
     bool manual_z_scaling(const char * offset, AmiraSpatialGraph * graph, bool inverse);
     
//     double ** readAndGetAmiraTransformation();
     

     std::vector<double *> * transformCoordinates(std::vector<double *> * list);
     
     std::vector<double *> ProximityLandmarks;
     
     std::vector<double *> inSectionLandmarks;
     
      std::vector<int> sectionLineNumbers;
      
      
     
     //void writeProximityImages(const char * outputFilename, std::vector<Proximity *>* proxlistlist[][]);
     
     
     void markLocalMaximums();
     float averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection);
     void markBoutons();     
     float calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius);
     float calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius);
     float calculateThreshold(ImageType::Pointer image, float x0, float y0, float z0);
     float setAverageThreshold();
     float setAverageThreshold(int numOfPointsToAverage);
     int setRadiusAndLocalBrightness(ImageType::Pointer image);
      int addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance );
      Image2DType::Pointer getImagePlane(int z, int depth, ImageType::Pointer input_image);
      VECTOR * sendRay(float x0, float y0, float ray_length, unsigned int nr_of_rays, float threshold, unsigned int n, Image2DType::Pointer image_plane);
      float bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane);
      void smoothRadii();
      bool isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance);
      float getAverageStandardDeviation();
      void getIDsFromMergedFile(AmiraSpatialGraph * merged_graph);
      std::vector<double* > * transformCoordinatesZScale(std::vector<double *> * list, const char * offset);
      double calculateConfidenceValue(Edge * currentEdge, /*double * coords , */std::vector<double*> * confidenceList, std::vector<double>  distances);
      void writeLandmarks( const char* outputFilename, std::vector<double*> *landmark_vector);
      void getMacthingGraphPoints(AmiraSpatialGraph * sectiongraph_to_be_matched, AmiraSpatialGraph * main_graph_to_be_matched);
      //void writeMasterHxFile(const char* outputFilename, const char * mergedAmiraFileName, int cellIDs);
//      float boutonConfidenceValue(Edge * currentEdge, double * coords);
//      void cutOffGraph(AmiraSpatialGraph * merged_graph);
      
      
//      ImageType::Pointer getOriginalImage(char* file_name, int start_index, int end_index);
//      
     
;
     


  
 private:
     Reader * graphReaderPointer;
     std::vector<SimpleAxon> axon_list;
     std::vector<SimpleDendrite> dendrite_list;
     
     //std::vector<std::vector<Proximity *>*> proximity_master_list;
     
     std::vector<Proximity *> proximity_list_1_2;
     std::vector<Proximity *> proximity_list_2_1;
     std::vector<Proximity *> proximity_list_1_3;
     std::vector<Proximity *> proximity_list_3_1;
     std::vector<Proximity *> proximity_list_2_3;
     std::vector<Proximity *> proximity_list_3_2;
     std::vector<Proximity *> proximity_list;
//      std::vector<Proximity *> proximity_list_nondecon;
     std::vector<double *> * all_bouton_list;
     ImageType::Pointer original_image;
#ifdef NONDECON
     ImageType::Pointer original_image_non_decon;
#endif
//      ImageType::Pointer original_image_nondecon;
     float max_dist;
     AmiraSpatialGraph * amira_graph;
     AmiraSpatialGraph * merged_graph;
     AmiraSpatialGraph * untransformed_section_graph;
     AmiraSpatialGraph * transformed_section_graph;
     double ** sectionTranslation;
     double ** sectionRotation;
     double ** transformation;
     double ** inverse_transformation;
     const char * inputfilename;
     bool eliminateDoubledProximites;
     bool include_unknown;
     int numberOfPointsAway;
     BoundingBox graph_bounding_box;
     BoundingBox bounding_box_merged;
     double * sectionboundingbox;
     double * bound_box;
     const char * offsetfile;
     const char * outputpath;

    
     double calculateGap(double * axonPt, double * dendritePt, ImageType::Pointer image, bool istest,int prox);
     bool isTheSameCell(int AxonID, int DendriteID);
     bool isAxon(int ID);
     bool isDendrite(int ID);
     //bool getCellIndex(int AxonID, int DendriteID);
     int getCellIndex(std::vector<int> *IDArray, int ID);
     //void AxonDendriteProximityFinder::InitProxMasterList(void);
     
     
     bool isInSectionBoundingBox(double * coords);
     
//     void getBoundingBox();
     
//     Edge * transformEdge(Edge * edge);
     

     void set_edge_lists(AmiraSpatialGraph * merged_input_graph, bool normalmode);
     
     
};

  


#endif
