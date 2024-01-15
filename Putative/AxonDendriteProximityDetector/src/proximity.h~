/****************************************************************************/
/*                                                                          */
/* File:      proximity.h                                            */
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


#ifndef PROXIMITY
#define PROXIMITY


#include "bounding_box.h"
#include "utility.h"
#include "basics.h"
#include "amiraReader.h"
#include "typedefs.h"
//#include "axon_dendrite_proximity_finder.h"




//****************************
#define XY_BUFFER          20
#define Z_BUFFER           15
#define MAX_PROX_BOX_SIZE  5
#define MAX_PROX_BOX_SIZE_IMAGE_COORD_XY   MAX_PROX_BOX_SIZE/XYSAMPLING
#define MAX_PROX_BOX_SIZE_IMAGE_COORD_Z    MAX_PROX_BOX_SIZE/ZSAMPLING

/****************************************************************************/
/*                                                                          */
/* Class: Proximity                                                         */
/*                                                                          */
/*      Defines a proximity zone where a dendrite and axon are near to      */
/*      each other. Contains the BoundingBox marking that zone as well as   */
/*      methods for writing the portion of the image represented by the     */
/*      Proximity object into an imagestack contained within its own file   */
/*                                                                          */
/*                                                                          */
/****************************************************************************/
class Proximity
{
 public:
    //Proximity(){}
     
    Proximity(std::vector<double *> * pointList, ImageType::Pointer og_image, int * IDs)
    {
        original_image = og_image;
        Axon_ID = IDs[0];
        Dendrite_ID = IDs[1];
        
        for(int i=0; i<pointList->size();i++)
        {
            bounding_box.expandBox_amiraCoords(pointList->at(i));
            
            ////std::cout<<pointList->at(i)[X_COORD]<<"  "<<pointList->at(i)[Y_COORD]<<"  "<<pointList->at(i)[Z_COORD]<<std::endl;
            
        }
        
        this->contains_bouton = false;
	
       axonDistanceFromSoma = 0.0;
       dendriteDistanceFromSoma = 0.0;
    }
    
    Proximity(double * axon, double * dendrite, double dist, int  axonID, int dendID, ImageType::Pointer og_image, int linearid)
    {
        original_image = og_image;
        
       axon_point[0] = axon[0];
       axon_point[1] = axon[1];
       axon_point[2] = axon[2];
       
       
       dendrite_point[0] = dendrite[0];
       dendrite_point[1] = dendrite[1];
       dendrite_point[2] = dendrite[2];
       
       
       
       centriod[0] = (axon_point[0]+dendrite_point[0]) / 2;
       centriod[1] = (axon_point[1]+dendrite_point[1]) / 2;
       centriod[2] = (axon_point[2]+dendrite_point[2]) / 2;
        
       
       distance = dist;
       
       Axon_ID = axonID;
       Dendrite_ID = dendID;
       
       linear_id = linearid;
       
       axonDistanceFromSoma = 0.0;
       dendriteDistanceFromSoma = 0.0;
        
    }
    
    Proximity(double * axon, double * dendrite, double dist, int  axonID, int dendID, ImageType::Pointer og_image, int linearid, double _axonDistanceFromSoma, double _dendriteDistanceFromSoma)
    {
        original_image = og_image;
        
       axon_point[0] = axon[0];
       axon_point[1] = axon[1];
       axon_point[2] = axon[2];
       
       
       dendrite_point[0] = dendrite[0];
       dendrite_point[1] = dendrite[1];
       dendrite_point[2] = dendrite[2];
       
       
       
       centriod[0] = (axon_point[0]+dendrite_point[0]) / 2;
       centriod[1] = (axon_point[1]+dendrite_point[1]) / 2;
       centriod[2] = (axon_point[2]+dendrite_point[2]) / 2;
        
       
       distance = dist;
       
       Axon_ID = axonID;
       Dendrite_ID = dendID;
       
       linear_id = linearid;
       
       axonDistanceFromSoma = _axonDistanceFromSoma;
       dendriteDistanceFromSoma = _dendriteDistanceFromSoma;
        
    }
    
    Proximity(Proximity * input_prox, ImageType::Pointer og_image, int linearid)
    {
        original_image =  og_image;
        
       axon_point[0] = input_prox->getAxonPoint()[0];
       axon_point[1] = input_prox->getAxonPoint()[1];
       axon_point[2] = input_prox->getAxonPoint()[2];
       
       
       dendrite_point[0] = input_prox->getDendritePoint()[0];
       dendrite_point[1] = input_prox->getDendritePoint()[1];
       dendrite_point[2] = input_prox->getDendritePoint()[2];
       
       
       
       centriod[0] = (axon_point[0]+dendrite_point[0]) / 2;
       centriod[1] = (axon_point[1]+dendrite_point[1]) / 2;
       centriod[2] = (axon_point[2]+dendrite_point[2]) / 2;
        
       
       distance = input_prox->getDistance();
       
       Axon_ID = input_prox->Axon_ID;
       Dendrite_ID = input_prox->Dendrite_ID;
       
       linear_id = linearid;
       
       axonDistanceFromSoma = input_prox->axonDistanceFromSoma;
       dendriteDistanceFromSoma = input_prox->dendriteDistanceFromSoma;
        
    }
    
    Proximity(Proximity * input_prox, ImageType::Pointer og_image)
    {
        original_image =  og_image;
        
       axon_point[0] = input_prox->getAxonPoint()[0];
       axon_point[1] = input_prox->getAxonPoint()[1];
       axon_point[2] = input_prox->getAxonPoint()[2];
       
       
       dendrite_point[0] = input_prox->getDendritePoint()[0];
       dendrite_point[1] = input_prox->getDendritePoint()[1];
       dendrite_point[2] = input_prox->getDendritePoint()[2];
       
       
       
       centriod[0] = (axon_point[0]+dendrite_point[0]) / 2;
       centriod[1] = (axon_point[1]+dendrite_point[1]) / 2;
       centriod[2] = (axon_point[2]+dendrite_point[2]) / 2;
        
       
       distance = input_prox->getDistance();
       
       Axon_ID = input_prox->Axon_ID;
       Dendrite_ID = input_prox->Dendrite_ID;
       
       linear_id = input_prox->getLinearID();
       
       axonDistanceFromSoma = input_prox->axonDistanceFromSoma;
       dendriteDistanceFromSoma = input_prox->dendriteDistanceFromSoma;
        
    }
    
    double getDistance()
    {
        return distance;
    }
    
    double * getAxonPoint()
    {
        return axon_point;
    }
    
    double * getDendritePoint()
    {
        return dendrite_point;
    }
    
    double * getRealCentroid()
    {
        return centriod;
    }
    
    double * getRealCentroidInImageCoords()
    {
        double * centre_point = new double[3];
        centre_point[0] = centriod[0]/XYSAMPLING;
        centre_point[1] = centriod[1]/XYSAMPLING;
        centre_point[2] = centriod[2]/ZSAMPLING;
        return centre_point;
    }
     
    double * getCenterCoords()
    {
        return bounding_box.getCenter_amiraCoords();
    }
   
    int getLinearID()
    {
        return linear_id;
    }
    
    int getAxonID()
    {
        return Axon_ID;
    }
    
    int getDendriteID()
    {
        return Dendrite_ID;
    }
    
    ImageType::RegionType getOutputRegion(); 
    
    ImageType::Pointer getProximityImage();
    void writeProximityImage( const char * outputFileName, double** transformation);
    void removeDoubledProximities();
    
    void setBouton(bool in_bouton){
        this->contains_bouton = in_bouton;
        
    }
    
    bool getBouton(){
        return this->contains_bouton;
    }
    
    double getConfidenceValue(){
        double confidence = this->confidence_value;
        return confidence;        
    };
    
    void calculateConfidenceValue(std::vector<double* > * confidenceList, std::vector<double>   distances) ;
    void setConfidenceValue(double value){
        this->confidence_value = value;
    };
    
    
    void getBoundingBoxCoords(double box_coords[6]){
//         double * box_coords = new double[6];
        this->bounding_box.getBoundingBox(box_coords);
        ////std::cout<<box_coords[0]<<"  "<<box_coords[1]<<"  "<<box_coords[2]<<"  "<<box_coords[3]<<"  "<<box_coords[4]<<"  "<<box_coords[5]<<"  "<<std::endl;
//         return box_coords;
        
    }
    ImageType::RegionType getIteratorRegion();
    
    int Axon_ID, Dendrite_ID;
    double axonDistanceFromSoma;
    double dendriteDistanceFromSoma;
  
 private:
    double axon_point[3];
    double dendrite_point[3];
    double distance;
    double centriod[3];
    
    BoundingBox bounding_box;
    ImageType::Pointer original_image;
    AmiraSpatialGraph* amira_graph;
    bool contains_bouton;
    double confidence_value;
    int linear_id;
    
    ImageType::Pointer getProximityImage(int x_dim, int y_dim);
    
    void writeImagePlanes(ImageType::Pointer input_image, const char * output_file);
    
    void writeLandmarkFile( const char* outputFilename);
     //void writeListFile(std::vector<double*> * list, const char* outputFilename);
        
     //void writeSingleLandmarkFile(double * coord, const char* outputFilename, int number);
    void writeInfoFile(const char* outputFilename);
     
    void writeHXFile( const char* outputFilename, double** transformation);   
    
        
};




#endif

