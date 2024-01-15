/*
    Copyright 2015 NeuroMorph <email>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/
#include "typedefs.h"
#include "utility.h"
#include "typedefs_two.h"
#include "basics.h"

#ifndef BOUTON_PARAMS_H
#define BOUTON_PARAMS_H

#define RAY_LEN_PER_DIRECTION   0.75

#define NUM_RAYS                360

#define Axon_1_ID  15
#define Axon_2_ID  16
#define Axon_3_ID  17

// Distance from the boundary to ignore in microns
// overlap between prox regions are 5 microns
#define X_OFFSET   1
#define Y_OFFSET   1
#define Z_OFFSET   1

#define NUMBER_OF_PIXELS_FOR_OVERLAP_RAY   2
class bouton_params
{
    
public:
    
    // Orientation at different points in xy, xz and yz planes
    // +/- 1 pt, 2 pts, 3 pts, 5 pts, 7 pts, 10 pts
    //bool isBouton;
    double landmark[3];
    //double nearestPoint[4];//4th coor gives the edge num
    double transformed_landmark[4];
    //double distanceFromLandmark;
    double orientation[15];
    double Axon_1_XZ_Angle;
    double Axon_2_XZ_Angle;
    double Axon_3_XZ_Angle;
    //double orientation_global[3];
    double radius_1_plane_XY[2];// one for oriented other for min radius in 1 z plane image
    double radius_1_plane_XY_perpendicular[2];
    double radius_1_plane_XY_otsu[2];
    double radius_1_plane_XY_perpendicular_otsu[2];
    double radius_3_plane_XY[2];
    double radius_5_plane_XY[2];
    double contour_1_plane_XY[NUM_RAYS * 4];
    std::list<double*> contour_1_plane_XY_min_rad;
    double perpendicular_otsu_front[3];
    double perpendicular_otsu_back[3];
    double contour_3_plane_XY[NUM_RAYS * 4];
    std::list<double*> contour_3_plane_XY_min_rad;
    double contour_5_plane_XY[NUM_RAYS * 4];
    std::list<double*> contour_5_plane_XY_min_rad;
    double Brightness_1_plane_Total;
    double Brightness_1_plane_XY_BB[1];
    double Brightness_3_plane_XY_BB[1];
    double Brightness_5_plane_XY_BB[1];
    double Brightness_1_plane_XY_Contour[1];
    double Brightness_3_plane_XY_Contour[1];
    double Brightness_5_plane_XY_Contour[1];
    double Brightness_1_plane_XY_Contour_Min[1];
    double Brightness_3_plane_XY_Contour_Min[1];
    double Brightness_5_plane_XY_Contour_Min[1];
    std::list<double *> brightness_points_1;
    std::list<double *> brightness_points_3;
    std::list<double *> brightness_points_5;
    unsigned int which_cluster;
    unsigned int prev_which_cluster;
    bool  overlap[2];
    

public:
    
    bouton_params()
    {
    }
    
    bouton_params(double * landmark_coord, bool isbouton)
    {
        
        
        //transformed_landmark = new double [3];
        transformed_landmark[0] = 0;
        transformed_landmark[1] = 0;
        transformed_landmark[2] = 0;
        transformed_landmark[3] = 0;
        
        
        
        //orientation = new double[15];
        //////////////////std:://cout<<" "<<orientation<<std::endl;
        /*
        radius_1_plane = new std::list< double * >;
        ////////////////std:://cout<<" "<<radius_1_plane<<std::endl;
        radius_3_plane = new std::list< double * >;
        radius_5_plane = new std::list< double * >;
        contour_1_plane = new std::list< double * >;
        contour_3_plane = new std::list< double * >;
        contour_5_plane = new std::list< double * >;
        Brightness_1_plane = new std::list< double * >;*/
        
        //////////////////std:://cout<<landmark[0]<<'\t'<<landmark[1]<<'\t'<<landmark[2]<<'\n';
         
    }
    
    // Read the landmark locations
    //void getLandmarks(char *);
    
    void readLandmarkCoordinates(std::ifstream* inputStream, bool isBouton);
        
    // Find the nearest point of axon to the landmoark location along with the distance
    void findNearestPoint();
    
    // Get the orientation of each point taking different neighbors
    void findOrientationOfAxon();
    
    // Find the radius at each point in 3 planes
    void findRadius();
    
     
    
    
    
};

#endif // BOUTON_PARAMS_H
