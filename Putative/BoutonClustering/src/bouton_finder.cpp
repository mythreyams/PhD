/****************************************************************************/
/*                                                                          */
/* File:      bouton_finder.cpp                                                     */
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
#include "bouton_finder.h"

#include "math.h"

//#define PI 3.14159265



//extern std::list< bouton_params * > BoutonParamsList;
extern bouton_params BoutonParamsArray[20];
extern int NumOfLandmarks;
extern int AxonID;
/****************************************************************************/
/* This constructor reads the given amira file and copies the image stack   */
/****************************************************************************/
BoutonFinder::BoutonFinder(char* amira_file, char* image_folder_name, int start_index,int end_index, char* hxfile, char* outputpath )
{
    
  std::cout << "Begin reading" << std::endl;
  
  outputFilename = outputpath;
  // Read cell's amira graph
  amiraReader = new Reader(amira_file, outputpath);//"amira_file_withRad_centerline");  
  amiraReader->readSpatialGraphFile(false);
  
  
  amira_graph = amiraReader->getSpatialGraph();
  
  std::cout << "Begin reading image" << std::endl;
  
  // read the image stack of the proximity region
  original_image = getOriginalImage(image_folder_name, start_index, end_index);
  
  std::cout << "end reading image" << std::endl;
  
  loadhxfile = new char[50];
  
  loadhxfile = hxfile;
  
 
//#if 0
  // read the loadHx file so that translation and transformation can be read

  if(!getTransformationandInverseTransformation(loadhxfile))
  {
      std::cout << "Error getting inverse transform" << std::endl;
  }
  
  
  
  // transform the spatial graph to image coordinates

          // first inverse transform so that the whole cell graph is inline with
          // the section's graph
          amira_graph->applyInverseTransformation(inverse_transformation);
          
          //std::cout << "ok here" << std::endl;

          // now apply the inverse translation to match the image section
          getTranslation(loadhxfile);

          //std::cout << "ok here is it" << std::endl;
        // now apply the inverse scaling
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        std::vector< Vertex * > * vertices = amira_graph->verticesPointer();
        int numOfEdges = edges->size();
        int numofVertices = vertices->size();
        
        //int pixel_value = 0;
        for(int i=0; i<numofVertices; i++)  //for each edge
        {
            
            Vertex * currentVertex = vertices->at(i);
            
            double * coords = (currentVertex)->coordinates;
            coords[X_COORD] = coords[X_COORD] - imageTranslation[X_COORD];
            coords[Y_COORD] = coords[Y_COORD] - imageTranslation[Y_COORD];
            coords[Z_COORD] = coords[Z_COORD] - imageTranslation[Z_COORD];
            
        }
        
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
          
          Edge * currentEdge = edges->at(i);
          std::list< double * >::iterator it;
          
          for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
          {
            double * coords = *it;
                
            //if(coords[IS_BOUTON])
            {
                ImageType::IndexType pixelIndex;
                
                // Apply image translation before scaling
                coords[X_COORD] = coords[X_COORD] - imageTranslation[X_COORD];
                coords[Y_COORD] = coords[Y_COORD] - imageTranslation[Y_COORD];
                coords[Z_COORD] = coords[Z_COORD] - imageTranslation[Z_COORD];
              /*
                int x_pos = rint( coords[X_COORD] / XYSAMPLING );
                int y_pos = rint( coords[Y_COORD] / XYSAMPLING );
                int z_pos = rint( coords[Z_COORD] / ZSAMPLING);
                
                coords[X_COORD] = x_pos;
                coords[Y_COORD] = y_pos;
                coords[Z_COORD] = z_pos;*/
                
                //pixel_value = 255;
                    
                //output_image->SetPixel(pixelIndex, pixel_value);
                
//              std::cout << "Bouton at (x, y): ("<< x_pos << ", " << y_pos << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
            }
          }
          
        }
        
        // Get the axon points
        getAxonPoints(original_image);
        
//#endif
        
        // write it so I can see it
        
        //amiraReader->writeSpatialGraphFile2();
      
      // Undo the translation and transformation before writing the final output
#if 0
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
          
          Edge * currentEdge = edges->at(i);
          std::list< double * >::iterator it;
          
          for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
          {
            double * coords = *it;
                
            //if(coords[IS_BOUTON])
            {
                ImageType::IndexType pixelIndex;
                
                // Apply image translation before scaling
                coords[X_COORD] = coords[X_COORD] + imageTranslation[X_COORD];
                coords[Y_COORD] = coords[Y_COORD] + imageTranslation[Y_COORD];
                coords[Z_COORD] = coords[Z_COORD] + imageTranslation[Z_COORD];
              /*
                int x_pos = rint( coords[X_COORD] / XYSAMPLING );
                int y_pos = rint( coords[Y_COORD] / XYSAMPLING );
                int z_pos = rint( coords[Z_COORD] / ZSAMPLING);
                
                coords[X_COORD] = x_pos;
                coords[Y_COORD] = y_pos;
                coords[Z_COORD] = z_pos;*/
                
                pixel_value = 255;
                    
                //output_image->SetPixel(pixelIndex, pixel_value);
                
//              std::cout << "Bouton at (x, y): ("<< x_pos << ", " << y_pos << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
            }
          }
          
        }       
        
        amira_graph->applyTransformation();
        
        amiraReader->writeSpatialGraphFile2();
        
#endif


};

/****************************************************************************/
/*findSpatialGraphPoint()  primary function for BoutonFinder class                    */
/****************************************************************************/
#if 0
void BoutonFinder::findSpatialGraphPoint(/*std::list< bouton_params * > *boutonParamsList*/)
{
    
    
    //std::cout<<landmark[0]<<'\t'<<landmark[1]<<'\t'<<landmark[2]<<std::endl;
    
    //double * min_coord = NULL;
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    
    unsigned int numOfEdges = edges->size();
    
    if(numOfEdges == 0)
    {
            std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
            //return -1;
    }
    else
    {
        
        //std::list<bouton_params*>::iterator paramIterator;
                                
        for(int i = 0; i < NumOfLandmarks;  i++)
        {
            double min_distance = 0xffffffffffffffff;
                
        
// #pragma omp parallel for schedule(dynamic,1)
            for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
            {                       
                    Edge * currentEdge = edges->at(pos);
                    std::list< double * >::iterator edge_it;        
                    std::list< double * >::iterator next_it;

                    //for every point along edge
                    for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                    {                 
                        double * coords = *edge_it;
                        
                        
                        double dist = Utility::euclidean_distance(BoutonParamsArray[i].landmark,coords,3,1);
                        
                        if(dist < min_distance)
                        {
                            
                            min_distance = dist;
                            BoutonParamsArray[i].nearestPoint[0] = coords[0];
                            BoutonParamsArray[i].nearestPoint[1] = coords[1];
                            BoutonParamsArray[i].nearestPoint[2] = coords[2];
                            //min_coord = coords;
                            BoutonParamsArray[i].nearestPoint[3] = pos;
                            
                            AxonID = currentEdge->label;
                            
                            
                            //std::cout<<dist<<'\t'<<min_coord[0]<<'\t'<<min_coord[1]<<'\t'<<min_coord[2]<<std::endl;
                        }
                            
                    }
                    
                    // save the axon id we are working with
                    
                    
                    
                    
            }
            
            std::cout<<'\t'<<BoutonParamsArray[i].nearestPoint[0]<<'\t'<<BoutonParamsArray[i].nearestPoint[1]<<'\t'<<BoutonParamsArray[i].nearestPoint[2]<<std::endl;
        }
        
        // Get all the points on the axon within 10 um
        
        
        // Arrange the points as per their distance from the point
        
        // 
        
        
    }
   
    //return min_coord;
                            
}
#endif
/****************************************************************************/
/*findOrientationOfAxon()  primary function for BoutonFinder class                    */
/****************************************************************************/
void BoutonFinder::findOrientationOfAxon(const double * landmarkPoint, int numNieghbours,  double* angleslist)
{
    // we want to take 5 - 10 neiboring points in either dierction
    // the problem is that the current edge may not ahve all the points needed
    // and next edge and prev edge may not be the neiboring ones in
    
    // lets search for the nearest points and make sure that they belong to the same axon
    // Direction : x,y,z greater than then +ve direction.. otherwise -ve direction
    //std::cout<<landmarkPoint[0]<<'\t'<<landmarkPoint[1]<<'\t'<<landmarkPoint[2]<<std::endl;
        
    int ind = 0;
    std::list< double * > points;
    //std::list< double * > angleslist;
    //double * distCoord = new double[4];
    //std::list< double * > distances;  
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    unsigned int numOfEdges = edges->size();
    if(numOfEdges == 0)
    {
            std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
            //return -1;
    }
    else
    {
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                std::list< double * >::iterator next_it;

                if(currentEdge->label == edges->at(landmarkPoint[3])->label)
                {
                    // 
                
                    //for every point along edge
                    for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                    {    
                        
                        double * coords = *edge_it;
                        double * distCoord = new double[4];
                        double * landmarkCoord = new double[3];
                        landmarkCoord[0] = landmarkPoint[0];
                        landmarkCoord[1] = landmarkPoint[1];
                        landmarkCoord[2] = landmarkPoint[2];
                        
                        double dist = 0;
                        dist = Utility::euclidean_distance(landmarkCoord,coords,3,1);
                        
                        if(dist < 5)
                        {
                            distCoord[0] = dist;
                            distCoord[1] = coords[0];
                            distCoord[2] = coords[1];
                            distCoord[3] = coords[2];
                            points.push_back(distCoord);
                            //points.push_back(&dist);
                            
                            //min_coord[3] = pos;
                            
                            //std::cout<<dist<<'\t'<<coords[0]<<'\t'<<coords[1]<<'\t'<<coords[2]<<std::endl;
                        }
                            
                    }
                }
                
          }
          // Does not work
          //points.sort();
          std::list< double * >::iterator point_it; 
          
          std::list< double * > SortedPoints;
          
          while(!points.empty())
          {
              double mindist = 0xffffffffffffffff;
              double * mindistCoord = new double(4); 
             // Optimize this sorting algorithm
             for(point_it = points.begin(); point_it !=  points.end(); point_it++)
             {
                
                 
                 if((*point_it)[0] < mindist)
                 {
                      mindist = (*point_it)[0];
                      mindistCoord = *point_it;
                     
                 }
                 
                 // sort the points according to their distance from the landmark point
                 //std::cout<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<(*point_it)[3]<<std::endl;
                
                
             }
             SortedPoints.push_back(mindistCoord);
             points.remove(mindistCoord);
          }
          
          // print
          for(point_it = SortedPoints.begin(); point_it !=  SortedPoints.end(); point_it++)
          {
         
            // sort the points according to their distance from the landmark point
            //std::cout<<'\t'<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<(*point_it)[3]<<std::endl;
        
          }
          
          std::list< double * > PositivePoints;
          std::list< double * > NegativePoints;
          // select the nearest points along +ve and -ve directions
          for(point_it = (++SortedPoints.begin()); point_it !=  SortedPoints.end(); point_it++)
          {
            double * SortedCoord = new double(3);
            SortedCoord[0] = (*point_it)[1];
            SortedCoord[1] = (*point_it)[2];
            SortedCoord[2] = (*point_it)[3];
            if((SortedCoord[0] >= (*SortedPoints.begin())[1]) && (SortedCoord[1] >= (*SortedPoints.begin())[2]) && (SortedCoord[2] >= (*SortedPoints.begin())[3])) 
            {
                PositivePoints.push_back(SortedCoord);
                // sort the points according to their distance from the landmark point
                //std::cout<<'\t'<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<(*point_it)[3]<<std::endl;
            }
            else
            {
                NegativePoints.push_back(SortedCoord);
            }
        
          }
          
          
          /*for(point_it = PositivePoints.begin(); point_it !=  PositivePoints.end(); point_it++)
          {
         
            // sort the points according to their distance from the landmark point
            std::cout<<'\t'<<'\t'<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<std::endl;
        
          }
          
          for(point_it = NegativePoints.begin(); point_it !=  NegativePoints.end(); point_it++)
          {
         
            // sort the points according to their distance from the landmark point
            std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<std::endl;
        
          }*/
          std::list< double * >::iterator positive_it;
          std::list< double * >::iterator negative_it;
          //std::list< double * >::iterator angles_it;
          
          int i = 0;
          // Find the slope between points
          for(positive_it = PositivePoints.begin(),negative_it = NegativePoints.begin(); i< numNieghbours; positive_it++,negative_it++,i++)
          {
              double * angles = new double[3];
              
              
              //(*positive_it)[0] = (*positive_it)[0] - (*SortedPoints.begin())[1];
              //(*positive_it)[1] = (*positive_it)[1] - (*SortedPoints.begin())[2];
              //(*positive_it)[2] = (*positive_it)[2] - (*SortedPoints.begin())[3];
              
              //(*negative_it)[0] = (*negative_it)[0] - (*SortedPoints.begin())[1];
              //(*negative_it)[1] = (*negative_it)[1] - (*SortedPoints.begin())[2];
              //(*negative_it)[2] = (*negative_it)[2] - (*SortedPoints.begin())[3];
              // xy angle atan(y2-y1/x2-x1)
              //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<(*positive_it)[0]<<'\t'<<(*positive_it)[1]<<'\t'<<(*positive_it)[2]<<'\t'<<std::endl;
              //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<(*negative_it)[0]<<'\t'<<(*negative_it)[1]<<'\t'<<(*negative_it)[2]<<'\t'<<std::endl;
              
              if((*positive_it)[0] != (*negative_it)[0])
              {
                  //angles[0] = asin(((*positive_it)[1] - (*negative_it)[1]) / ((*positive_it)[0] - (*negative_it)[0]) ) * 180 / PI; 
                  //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<angles[0]<<std::endl;
                  //angles[0] = acos(((*positive_it)[1] - (*negative_it)[1]) / ((*positive_it)[0] - (*negative_it)[0]) ) * 180 / PI; 
                  //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<angles[0]<<std::endl;
                  angles[0] = atan2(((*positive_it)[1] - (*negative_it)[1]) , ((*positive_it)[0] - (*negative_it)[0]) ) * 180 / PI; 
                  //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<angles[0]<<std::endl;
              }
              else
              {
                  angles[0] = 90;
              }
              
              // yz angle
              if((*positive_it)[1] != (*negative_it)[1])
              {
                  angles[1] = atan2(((*positive_it)[2] - (*negative_it)[2]) , ((*positive_it)[1] - (*negative_it)[1]) ) * 180 / PI; 
              }
              else
              {
                  angles[1] = 90;
              }
              // xz angle
              if((*positive_it)[0] != (*negative_it)[0])
              {
                  angles[2] = atan2(((*positive_it)[2] - (*negative_it)[2]) , ((*positive_it)[0] - (*negative_it)[0]) ) * 180 / PI; 
              }
              else
              {
                  angles[2] = 90;
              }
              //std::cout<<'\t'<<'\t'<<'\t'<<'\t'<<angles[0]<<'\t'<<angles[1]<<'\t'<<angles[2]<<'\t'<<std::endl;
              
              
              //angleslist.push_back(angles);
              angleslist[ind]= 180 - abs(angles[0]);
              ind++;
              angleslist[ind]= 180 - abs(angles[1]);
              ind++;
              angleslist[ind]= 180 - abs(angles[2]);
              ind++;
              
           }
          
          /*
          // print angles
          for(point_it = angleslist.begin(); point_it !=  angleslist.end(); point_it++)
          {
         
            // sort the points according to their distance from the landmark point
            std::cout<<'\t'<<'\t'<<'\t'<<(*point_it)[0]<<'\t'<<(*point_it)[1]<<'\t'<<(*point_it)[2]<<'\t'<<std::endl;
        
          }*/
    
    }
    
    //return &orientation;
}

/****************************************************************************/
/*findRadius()  primary function for BoutonFinder class                    */
/****************************************************************************/
void  BoutonFinder::findBrightness(Image2DType::Pointer image_plane, double * ContourList,  std::list<double*>cont_min_rad, double *  BrightnessBB, double *  BrightnessContour, double *  BrightnessMinContour, std::list<double *> *brightness_points_list)
{
    double min_x = 0xffffffffffffffff; 
    double min_y = 0xffffffffffffffff;
    double max_x = 0;
    double max_y = 0;
    
    // find min and max x,y coor  
    
    
    
    for(int i = 0; i < NUM_RAYS*4; i = i+2)
    {
        
        if (ContourList[i] <= min_x)
            min_x = ContourList[i];
        
        if (ContourList[i+1] <= min_y)
            min_y = ContourList[i+1];
        
        if (ContourList[i] >= max_x)
            max_x = ContourList[i];
        
        if (ContourList[i+1] >= max_y)
            max_y = ContourList[i+1];
        
    }
    
    //std::cout<<min_x<<'\t'<<min_y<<'\t'<<max_x<<'\t'<<max_y<<'\n';
    
     Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());

     Image2DType::IndexType center_index;
       
       
     // Get the brightness of the bounding box
    
    for(int ix = min_x; ix<= max_x; ix++)
    {
        for(int iy = min_y; iy <= max_y; iy++)
        {
            
            if ((ix >= min_x)&&(ix <= max_x) && (iy >= min_y)&&(iy <= max_y) )
            {
//                 double * tmp = new double[3];
//                 int grey = 0;
                center_index[0] = ix;
                center_index[1] = iy;
                it.SetIndex(center_index);
//                 grey = it.Get();
                //double grey_value = it.Get();//bilinearInterpolation(x_f, y_f, image_plane);//
//                 tmp[0] = ix;
//                 tmp[1] = iy;
//                 tmp[2] = grey;
//                 brightness_points_list->push_back(tmp);
                
                BrightnessBB[0] += it.Get();
            }
                
                
        }
    }
    
//     bool found = false;
    
    
    
    // for each point within the Box need to check whether it is within the contour
    // fr each column
    //          check whether  contXmin<x<contXmax contYmin<y<contYmax
    
    
    // Brightness of the contour
    for(int ix = min_x; ix<= max_x; ix++)
    {
        int mincol = 0xffff;
        int maxcol = 0;
        
        std::list<double*>::iterator cont_min_rad_it;
        for(cont_min_rad_it = cont_min_rad.begin(); cont_min_rad_it != cont_min_rad.end(); cont_min_rad_it++)
        {
            //std::cout<<"     "<<ContourList[i]<<std::endl;
            if ( ix == (*cont_min_rad_it)[0]) /*&& (iy == ContourList[i+1] && iy <= ContourList[i+1]*/
            {
                //std::cout<<"before"<<ix<<" "<<mincol<<" "<<maxcol<<std::endl;
                if((*cont_min_rad_it)[1] >=  maxcol)
                {
                    maxcol = (*cont_min_rad_it)[1];
                    
                }
                    
                if((*cont_min_rad_it)[1] <= mincol)
                {
                    mincol = (*cont_min_rad_it)[1];
                        
                }
                
                //std::cout<<"after"<<ix<<" "<<mincol<<" "<<maxcol<<std::endl;
                
            }
            
        }
        
        //std::cout<<"--------------"<<std::endl;
        
        //std::cout<<ix<<std::endl;
        //std::cout<<min_y<<" "<<max_y<<std::endl;
        
        //std::cout<<mincol<<" "<<maxcol<<std::endl;
        
        // go through all the y's for this x and find the cont max and min ys for this column
        for(int iy = min_y; iy <= max_y; iy++)
        {
            // see if the pixel falls within the bouton contour
            if((iy >= mincol) && (iy<= maxcol))
            {
                double * tmp = new double[3];
                int grey = 0;
                center_index[0] = ix;
                center_index[1] = iy;
                it.SetIndex(center_index);
                grey = it.Get();
                //double grey_value = it.Get();//bilinearInterpolation(x_f, y_f, image_plane);//
                tmp[0] = ix;
                tmp[1] = iy;
                tmp[2] = grey;
                brightness_points_list->push_back(tmp);
                
                BrightnessMinContour[0] += it.Get();
            }
                
        }
    }
    
    
    // Brightness of the min cont
    for(int ix = min_x; ix<= max_x; ix++)
    {
    
        int mincol = 0xffff;
        int maxcol = 0;
        
        
        for(int i = 0; i < NUM_RAYS*4; i = i+2)
        {
            //std::cout<<"     "<<ContourList[i]<<std::endl;
            if ( ix == ContourList[i]) /*&& (iy == ContourList[i+1] && iy <= ContourList[i+1]*/
            {
               // std::cout<<"before"<<ix<<" "<<mincol<<" "<<maxcol<<std::endl;
                if(ContourList[i+1] >=  maxcol)
                {
                    maxcol = ContourList[i+1];
                    
                }
                    
                if(ContourList[i+1] <= mincol)
                {
                    mincol = ContourList[i+1];
                        
                }
                
               // std::cout<<"after"<<ix<<" "<<mincol<<" "<<maxcol<<std::endl;
                
            }
            
        }
        
        //std::cout<<"--------------"<<std::endl;
        
        //std::cout<<ix<<std::endl;
        //std::cout<<min_y<<" "<<max_y<<std::endl;
        
        //std::cout<<mincol<<" "<<maxcol<<std::endl;
        
        // go through all the y's for this x and find the cont max and min ys for this column
        for(int iy = min_y; iy <= max_y; iy++)
        {
            // see if the pixel falls within the bouton contour
            if((iy >= mincol) && (iy<= maxcol))
            {
                double * tmp = new double[3];
                int grey = 0;
                center_index[0] = ix;
                center_index[1] = iy;
                it.SetIndex(center_index);
                grey = it.Get();
                //double grey_value = it.Get();//bilinearInterpolation(x_f, y_f, image_plane);//
                tmp[0] = ix;
                tmp[1] = iy;
                tmp[2] = grey;
                //brightness_points_list->push_back(tmp);
                
                BrightnessContour[0] += it.Get();
            }
                
        }
    }
    
    
    /*std::list<double*>::iterator itr;
            
    for(itr= brightness_points_list.begin(); itr != brightness_points_list.end(); itr++)
    {
        //xlsParamWriter<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
        std::cout<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
        
    }*/
    
}


#if 0
/****************************************************************************/
/*findRadius()  primary function for BoutonFinder class                    */
/****************************************************************************/
void  BoutonFinder::findRadius(Image2DType::Pointer image_plane, double * landmarkPoint,  double *  orientationList, double *rad_list, double  *contour_list)
{
     //std::list <double *> * rad_contour_list =  new std::list <double *>;
     //std::list< double * > rad_list;
     std::list< double * > rad_list_local;
     //std::list< double * > contour_list;
     double min_rad = 0xffffffffffffffff;
    
    // match image with graph by tranforming graph
    
    // read the loadHx file so that translation and transformation can be read
    
        
        //int z_pos = rint(transformed_landmark[Z_COORD]/ZSAMPLING);
        //Image2DType::Pointer image_plane_3 = getImagePlane(z_pos,3,original_image);
    
    // from the projection image get the histogram along vertical, horizontal and oriented lines
    // get the pixels at half the histogram max or the min if half max is not found
    // which is the case when a dendrite is near by
        
        // Take the xy 1 pt angle from the list
        //std::list< double * >::iterator orientationList_it;
        //orientationList_it = orientationList->begin();
        std::list<double *> return_list; // has rad as first entry and 2 contour pts as the next two entries
        std::list< double * >::iterator return_list_it;
        
        //std::cout<< " " << &rad_list << '\t' <<&contour_list << '\t'<<std::endl;
        
        double phi = (orientationList[0]+90) * PI/180;
        double * radius = new double[1];
        return_list = sendRay(landmarkPoint[0],landmarkPoint[1],RAY_LEN_PER_DIRECTION,phi,image_plane);
        return_list_it = return_list.begin();
         *radius = (*return_list_it)[0];
            return_list.pop_front();
            contour_list->push_back( *(return_list.begin()));
            return_list.pop_front();
            contour_list->push_back( *(return_list.begin()));
        //std::cout<< "radius: " << *     radius << "angle: " << phi*180/PI<<std::endl;
        
        rad_list->push_back(radius);
        
        for(int i = 0; i < NUM_RAYS; i++)
        {
            double * radius = new double;
            double phi = i*(PI/NUM_RAYS);
            return_list = sendRay(landmarkPoint[0],landmarkPoint[1],RAY_LEN_PER_DIRECTION,phi,image_plane,rad_list,contour_list);
            return_list_it = return_list.begin();
            *radius = (*return_list_it)[0];
            return_list.pop_front();
            contour_list->push_back( *(return_list.begin()));
            return_list.pop_front();
            contour_list->push_back( *(return_list.begin()));
            //std::cout<< "radius: " << *radius << '\t'  << '\t'<<std::endl;
            rad_list_local.push_back(radius);
            
            
            if(*radius < min_rad )
            {
                min_rad = *radius;
            }
            
        }
        
        rad_list->push_back(&min_rad);
        
        std::list<double *>::iterator rad_it;
        
        for(rad_it = rad_list->begin(); rad_it != rad_list->end(); rad_it++)
        {
            //std::cout<< "radius: " << (*rad_it)[0] << '\t'  <<*rad_it<< '\t'<<std::endl;
        }
    
        for(rad_it = contour_list->begin(); rad_it != contour_list->end(); rad_it++)
        {
            //std::cout<< "radius: " << (*rad_it)[0] << '\t' <<(*rad_it)[1] << '\t'<<std::endl;
        }
    // for xy
    
    // take 3 or 4 projection images
    
    // for yz
    
    // for xz
        //rad_contour_list->push_back(rad_list);
        //rad_contour_list->push_back(contour_list);
        
        
        
        //return rad_list;
    
}


/****************************************************************************/
/*findRadius()  primary function for BoutonFinder class                    */
/****************************************************************************/
std::list <double *> BoutonFinder::findCountour(Image2DType::Pointer image_plane, double * landmarkPoint, std::list< double * > orientationList)
{
     //std::list <double *> * rad_contour_list =  new std::list <double *>;
     std::list< double * > rad_list;
     std::list< double * > rad_list_local;
     std::list< double * > contour_list;
     double min_rad = 0xffffffffffffffff;
    // TODO : extend to multiple rays
    
    // match image with graph by tranforming graph
    
    // read the loadHx file so that translation and transformation can be read
    
        
        //int z_pos = rint(transformed_landmark[Z_COORD]/ZSAMPLING);
        //Image2DType::Pointer image_plane_3 = getImagePlane(z_pos,3,original_image);
    
    // from the projection image get the histogram along vertical, horizontal and oriented lines
    // get the pixels at half the histogram max or the min if half max is not found
    // which is the case when a dendrite is near by
        
        // Take the xy 1 pt angle from the list
        std::list< double * >::iterator orientationList_it;
        orientationList_it = orientationList.begin();
        std::list<double *> return_list; // has rad as first entry and 2 contour pts as the next two entries
        std::list< double * >::iterator return_list_it;
        
        float phi = ((*orientationList_it)[0]+90) * PI/180;
        double * radius = new double[1];
        return_list = sendRay(landmarkPoint[0],landmarkPoint[1],RAY_LEN_PER_DIRECTION,phi,image_plane);
        return_list_it = return_list.begin();
         *radius = (*return_list_it)[0];
            return_list.pop_front();
            contour_list.push_back( *(return_list.begin()));
            return_list.pop_front();
            contour_list.push_back( *(return_list.begin()));
        std::cout<< "radius: " << *     radius << "angle: " << phi*180/PI<<std::endl;
        
        
        
        for(int i = 0; i < NUM_RAYS; i++)
        {
            double * radius = new double;
            double phi = i*(PI/NUM_RAYS);
            return_list = sendRay(landmarkPoint[0],landmarkPoint[1],RAY_LEN_PER_DIRECTION,phi,image_plane);
            return_list_it = return_list.begin();
            *radius = (*return_list_it)[0];
            return_list.pop_front();
            contour_list.push_back( *(return_list.begin()));
            return_list.pop_front();
            contour_list.push_back( *(return_list.begin()));
            std::cout<< "radius: " << *radius << '\t'  << '\t'<<std::endl;
            rad_list_local.push_back(radius);
            
            
            if(*radius < min_rad )
            {
                min_rad = *radius;
            }
            
        }
        
        rad_list.push_back(&min_rad);
        
        std::list<double *>::iterator rad_it;
        
        for(rad_it = rad_list.begin(); rad_it != rad_list.end(); rad_it++)
        {
            std::cout<< "radius: " << (*rad_it)[0] << '\t'  << '\t'<<std::endl;
        }
    
        for(rad_it = contour_list.begin(); rad_it != contour_list.end(); rad_it++)
        {
            std::cout<< "radius: " << (*rad_it)[0] << '\t' <<(*rad_it)[1] << '\t'<<std::endl;
        }
    // for xy
    
    // take 3 or 4 projection images
    
    // for yz
    
    // for xz
        //rad_contour_list->push_back(rad_list);
        //rad_contour_list->push_back(contour_list);
        
        
        
        return contour_list;
    
}

#endif

/****************************************************************************/
/*findBoutons()  primary function for BoutonFinder class                    */
/****************************************************************************/
void BoutonFinder::findBoutons()
{
  
    
      
  //if(setRadiusAndLocalBrightness(original_image) != 0)  
   // exit(-1);
   
  
  //writeImagePlanes(original_image, "raster_amira_slice");
  
  //markLocalMaximums();
  
  //markBoutons();
  
  UndoTransformation();
  
  //writeBoutonLandmarkfile();
  
  amiraReader->writeSpatialGraphFile();
  


};




/****************************************************************************/
/*getOriginalImage()  reads in the cropped deconvolved image stack from the given path */
/****************************************************************************/
ImageType::Pointer BoutonFinder::getOriginalImage(char* path_name, int start_index, int end_index)
{
        ImageType::Pointer output_image = ImageType::New();
        char input_file[1024];
        strcpy(input_file, "roi%03d.png");
        
        
        
        std::string format = "roi";
          //format += "_seg";
          format += "%03d.";
          format += "png";
        
          
          
        // Get into the given folder
        if(chdir(path_name) != 0)
                perror("Couldn't open image drirectory!");

        NameGeneratorType::Pointer name_gen = NameGeneratorType::New();         //DECLARE AND INITIALIZE NAME_GENERATOR
        name_gen->SetSeriesFormat( format.c_str() );
        name_gen->SetStartIndex( start_index );
        name_gen->SetEndIndex( end_index );
        name_gen->SetIncrementIndex( 1 );
        
                        
        
        SeriesReaderType::Pointer input_reader = SeriesReaderType::New();       //DECLARE AND INITIALIZE INPUT_READER
        input_reader->SetImageIO( itk::PNGImageIO::New() );
        input_reader->SetFileNames( name_gen->GetFileNames() );

        try
        {
                input_reader->Update();
        }
        catch( itk::ExceptionObject & err )
        {
                std::cerr << "ImageReaderExceptionObject caught !" << std::endl;
                std::cerr << err << std::endl;
        }
        
        ImageType::Pointer input_image;

        input_image = input_reader->GetOutput();
        input_image->Update();
        
#if 0
          
          
         ImageType::IndexType input_index;
         ImageType::SizeType input_size;
                
        input_index[0] = 0;
        input_index[1] = 0;
        input_index[2] = 0;

        input_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
        input_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
        input_size[2] = end_index - start_index +1;
        
        
        std::cout << "input_size: " << input_size[0] << " " << input_size[1] << " " << input_size[2] << std::endl;
        
        ImageType::RegionType input_region;// = input_image->GetLargestPossibleRegion();             //DECLARE REGION
        input_region.SetSize( input_size );
        input_region.SetIndex( input_index );

        output_image->SetRegions( input_region );        //ALLOCATE REGION
        output_image->Allocate();
        
        PasteFilterType::Pointer paster = PasteFilterType::New();

        paster->SetSourceImage( input_image );
        paster->SetSourceRegion( input_region );
        paster->SetDestinationImage( output_image );
        paster->SetDestinationIndex( input_index );
        paster->Update();

        output_image = paster->GetOutput();

        output_image->Update();

#endif
        
        return input_image;
};

bool BoutonFinder::getTransformationandInverseTransformation(char* inputfilename)
{
    std::cout<<"In getTranformfromfile"<<std::endl;
    // intialize transformation matrices
    
    //this->sectionTranslation = new double *[4];
    //this->sectionRotation = new double *[4];

    
            for(int ii = 0; ii < 4; ++ii)
            {
                    //this->sectionTranslation[ii] = new double[4];
                    //this->sectionRotation[ii] = new double[4];
                    for(int jj = 0; jj < 4; ++jj)
                    {
                            this->sectionTranslation[ii][jj] = 0;
                            this->sectionRotation[ii][jj] = 0;
                    }
                    this->sectionTranslation[ii][ii] = 1;
                    this->sectionRotation[ii][ii] = 1;
            }
        
            //this->transformation = new double *[4];
                    for(int ii = 0; ii < 4; ++ii)
                    {
                            //this->transformation[ii] = new double[4];
                            for(int jj = 0; jj < 4; ++jj)
                            {
                                    if(ii != jj)
                                            this->transformation[ii][jj] = 0;
                                    else
                                            this->transformation[ii][jj] = 1;
                            }
                    }
                
            //this->inverse_transformation = new double *[4];
                    for(int ii = 0; ii < 4; ++ii)
                    {
                            //this->inverse_transformation[ii] = new double[4];
                            for(int jj = 0; jj < 4; ++jj)
                            {
                                    if(ii != jj)
                                            this->inverse_transformation[ii][jj] = 0;
                                    else
                                            this->inverse_transformation[ii][jj] = 1;
                            }
                    }
                    
    //std::cout<<"inputfilename: " << inputfilename <<std::endl;
        std::ifstream inputStream(inputfilename);
     
        if(!inputStream.fail())
        {
                //const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                //const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t\n ";
                
                std::string currentLine;
                //unsigned int line = 0;
                
               
                
                //bool parameters = 1;
                //bool transform = 0;
 //               bool correctSection = 1;
 //               bool correctPrevSection = 0;
 //               int sectionID = 0;
                //unsigned int brackets = 0, transformBrackets = 0;
                //unsigned int currentIndex = 0;
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
                        
                        if(currentLine.size())
                            
                            
                                if(currentLine.find("setTransform ", 0) != std::string::npos)
                                {
                                    
                                              //std::cout << "found correct section transform parameters!" << std::endl;
                                        unsigned int count = 0;
                                        std::string::size_type loc1, loc2, loc3;
                                        loc1 = currentLine.find("setTransform ", 0);
                                        loc1 += 13;// Lenth of setTransform
                                        loc2 = currentLine.find_first_of(signs, 0);
                                        if(loc2 != std::string::npos)
                                                if(loc2 < loc1)
                                                        loc1 = loc2;
                                        loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                        while(loc2 != std::string::npos && count < 16)
                                        {
//                                             std::cout << loc1 << '\t'<<loc2<<std::endl;
                                            
                                            //std::cout<<"May crash here if the compiler is different.. seee comments in code"<<std::endl;
                                            
//##################WARNING#############################################################
                                            // eventhough we are creating the buffer at run time
                                            // the variable still has 4-5 bytes allocated and when we do atof the 
                                            // junk value (if a number) will get converted and hence the read values will be different
                                            // Therefore for this version of the compiler, I am initializing 5 bytes eventhough the variable is supposed to have 
                                            // as many number of bytes as required at runtime.. may cause issues if another compiler does not work like this
                                            
                                                char * tmp = new char[20];
                                                
                                                for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                {
                                                    tmp[i] = 'f';
                                                }
                                                
//                                                 std::cout << tmp<<std::endl;
                                                
                                                currentLine.copy(tmp, loc2 - loc1, loc1);
                                                
//                                                 std::cout << tmp<<std::endl;
                                                
                                                
                                                double ftmp1 = atof(tmp);
                                                
                                                //std::cout << "before!" << std::endl;
                                                //sectionRotation[0][0]= 1;        // amira files are columns after each other
                                                this->sectionRotation[count%4][count/4]= ftmp1;        // amira files are columns after each other
//                                                 std::cout << ftmp1 << '\t'<<tmp<<std::endl;
                                                loc3 = loc2;
                                                loc1 = currentLine.find_first_of(numbers, loc3);
                                                loc2 = currentLine.find_first_of(signs, loc3);
                                                if(loc2 != std::string::npos)
                                                        if(loc2 < loc1)
                                                                loc1 = loc2;
                                                loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                                                ++count;
                                                delete [] tmp;
                                        }
                                        
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        
                                                               
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
                                        //sectionRotation = transformation;
                                        //remove numeric artifacts from z-axis:
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         sectionRotation[2][ii] = 0;
//                                                         sectionRotation[ii][2] = 0;
//                                                 }
//                                                 sectionRotation[2][2] = 1;

                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                this->sectionTranslation[ii][3] = this->sectionRotation[ii][3];
                                                this->sectionRotation[ii][3] = 0;
                                        }
                                        this->sectionRotation[3][3] = 1;
                                        
//                                                 sectionTranslation[2][3] += manual_z_scale[1]/**manual_z_scale[0]*/;
//                                                 std::cout<<"shift : "<<sectionTranslation[2][3]<<std::endl;
                                        
                                        std::cout << "translation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionTranslation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionTranslation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
                                        std::cout << "rotation matrix with scaling:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionRotation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        double square_scaling[3];//myth increased to 3 from 2
                                        double scaling[4];
                                        scaling[3] = 1;
                                        
                                        
                                        
                                        for(int jj = 0; jj< 3; ++jj){
                                        
                                        square_scaling[jj] = this->sectionRotation[0][jj]*this->sectionRotation[0][jj] + this->sectionRotation[1][jj]*this->sectionRotation[1][jj] + this->sectionRotation[2][jj]*this->sectionRotation[2][jj];
                                        scaling[jj] = sqrt(square_scaling[jj]);
                                        std::cout<<"scaling : "<<scaling[jj]<<std::endl;
                                        
                                        }
                                        //scaling[2] = manual_z_scale[0];
                                        //std::cout<<"scaling : "<<scaling[2]<<std::endl;
                                        
                                        for(int ii = 0; ii < 3; ++ii)
                                        {                                                        
                                                for(int jj = 0; jj < 3; ++jj)                                                        
                                                    this->sectionRotation[ii][jj] = this->sectionRotation[ii][jj]/scaling[jj];                                                       
                                        } 
                                        
                                        
                                        std::cout << "rotation matrix without scaling:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->sectionRotation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->sectionRotation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                        }
                                        
//                                                 double ** mInverse = new double *[4];
//                                                 for(int ii = 0; ii < 4; ++ii)
//                                                 {
//                                                         mInverse[ii] = new double[4];
//                                                         for(int jj = 0; jj < 4; ++jj)
//                                                                 mInverse[ii][jj] = 0;
//                                                 }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
//                                                 mInverse[0][3] = -1*(sectionRotation[0][0]*sectionTranslation[0][3] + sectionRotation[1][0]*sectionTranslation[1][3]);
//                                                 mInverse[1][3] = -1*(sectionRotation[0][1]*sectionTranslation[0][3] + sectionRotation[1][1]*sectionTranslation[1][3]);
//                                                 mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                 mInverse[3][3] = 1;
//                                                 
//                                                 this->inverse_transformation = mInverse;
                                        
                                        
                                        
                                        
                                        double ** mProduct = new double *[4];
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                mProduct[ii] = new double[4];
                                                for(int jj = 0; jj < 4; ++jj)
                                                        mProduct[ii][jj] = 0;
                                        }
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                                for(int jj = 0; jj < 4; ++jj)
                                                        for(int kk = 0; kk < 4; ++kk){
                                                            
                                                            //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                                mProduct[ii][jj] += this->sectionTranslation[ii][kk]*this->sectionRotation[kk][jj]*scaling[jj];
                                                        }
                                               
                                        /*for(int ii = 0; ii < 4; ++ii)
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                    
                                                    //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                        this->transformation[ii][jj] = mProduct[ii][jj];
                                        }*/
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                //std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        this->transformation[ii][jj] = mProduct[ii][jj];
                                                }
                                                
                                                
                                        }
                                        //this->transformation = mProduct;
                                        std::cout << "transformation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->transformation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->transformation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                                
                                        }
//                                                 for(int ii = 0; ii < 2; ++ii)
//                                                 {
//                                                         transformation[2][ii] = 0;
//                                                         transformation[ii][2] = 0;
//                                                 }
//                                                 transformation[2][2] = 1;

                                        double ** mInverse = new double *[4];
                                        
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                mInverse[ii] = new double[4];
                                                for(int jj = 0; jj < 4; ++jj)
                                                        mInverse[ii][jj] = 0;
                                        }
                                        
                                        double a1 = this->sectionRotation[0][0]*this->sectionRotation[1][1]*this->sectionRotation[2][2]; 
                                        double a2 = this->sectionRotation[0][1]*this->sectionRotation[1][2]*this->sectionRotation[2][0];  
                                        double a3 = this->sectionRotation[0][2]*this->sectionRotation[1][0]*this->sectionRotation[2][1];
                                        double a4 = this->sectionRotation[2][0]*this->sectionRotation[1][1]*this->sectionRotation[0][2];
                                        double a5 = this->sectionRotation[2][1]*this->sectionRotation[1][2]*this->sectionRotation[0][0];
                                        double a6 = this->sectionRotation[2][2]*this->sectionRotation[1][0]*this->sectionRotation[0][1];
                                        
                                        double det = (a1 + a2 + a3 -a4 -a5 -a6)*scaling[0]*scaling[1]*scaling[2];
                                        
                                        mInverse[0][0] = (this->sectionRotation[1][1]*this->sectionRotation[2][2] - this->sectionRotation[1][2]*this->sectionRotation[2][1])*scaling[1]*scaling[2]/det;
                                        mInverse[0][1] = (this->sectionRotation[0][2]*this->sectionRotation[2][1] - this->sectionRotation[0][1]*this->sectionRotation[2][2])*scaling[1]*scaling[2]/det;
                                        mInverse[0][2] = (this->sectionRotation[0][1]*this->sectionRotation[1][2] - this->sectionRotation[0][2]*this->sectionRotation[1][1])*scaling[1]*scaling[2]/det;
                                        mInverse[1][0] = (this->sectionRotation[1][2]*this->sectionRotation[2][0] - this->sectionRotation[1][0]*this->sectionRotation[2][2])*scaling[0]*scaling[2]/det;
                                        mInverse[1][1] = (this->sectionRotation[0][0]*this->sectionRotation[2][2] - this->sectionRotation[0][2]*this->sectionRotation[2][0])*scaling[0]*scaling[2]/det;
                                        mInverse[1][2] = (this->sectionRotation[0][2]*this->sectionRotation[1][0] - this->sectionRotation[0][0]*this->sectionRotation[1][2])*scaling[0]*scaling[2]/det;
                                        mInverse[2][0] = (this->sectionRotation[1][0]*this->sectionRotation[2][1] - this->sectionRotation[1][1]*this->sectionRotation[2][0])*scaling[1]*scaling[0]/det;
                                        mInverse[2][1] = (this->sectionRotation[0][1]*this->sectionRotation[2][0] - this->sectionRotation[0][0]*this->sectionRotation[2][1])*scaling[1]*scaling[0]/det;
                                        mInverse[2][2] = (this->sectionRotation[0][0]*this->sectionRotation[1][1] - this->sectionRotation[0][1]*this->sectionRotation[1][0])*scaling[1]*scaling[0]/det;
                                        
                                        
                                        
//                                                mInverse[0][0] = sectionRotation[1][1]/scaling[1]*sectionRotation[2][2]/scaling[2] - sectionRotation[1][2]/scaling[2]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[0][1] = sectionRotation[0][2]/scaling[2]*sectionRotation[2][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[0][2] = sectionRotation[0][1]/scaling[1]*sectionRotation[1][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[1][1]/scaling[1];
//                                                mInverse[1][0] = sectionRotation[1][2]/scaling[2]*sectionRotation[2][0]/scaling[0] - sectionRotation[1][0]/scaling[0]*sectionRotation[2][2]/scaling[2];
//                                                mInverse[1][1] = sectionRotation[0][0]/scaling[0]*sectionRotation[2][2]/scaling[2] - sectionRotation[0][2]/scaling[2]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[1][2] = sectionRotation[0][2]/scaling[2]*sectionRotation[1][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[1][2]/scaling[2];
//                                                mInverse[2][0] = sectionRotation[1][0]/scaling[0]*sectionRotation[2][1]/scaling[1] - sectionRotation[1][1]/scaling[1]*sectionRotation[2][0]/scaling[0];
//                                                mInverse[2][1] = sectionRotation[0][1]/scaling[1]*sectionRotation[2][0]/scaling[0] - sectionRotation[0][0]/scaling[0]*sectionRotation[2][1]/scaling[1];
//                                                mInverse[2][2] = sectionRotation[0][0]/scaling[0]*sectionRotation[1][1]/scaling[1] - sectionRotation[0][1]/scaling[1]*sectionRotation[1][0]/scaling[0];
//                                                
                                        
//                                                for(int ii = 0; ii < 2; ++ii)
//                                                         for(int jj = 0; jj < 2; ++jj)
//                                                                 mInverse[ii][jj] = sectionRotation[jj][ii];
                                        //scaling[2] = 1;
                                        mInverse[0][3] = -1*(this->sectionRotation[0][0]/scaling[0]*this->sectionTranslation[0][3] + this->sectionRotation[1][0]/scaling[0]*this->sectionTranslation[1][3] + this->sectionRotation[2][0]/scaling[0]*this->sectionTranslation[2][3]);
                                        mInverse[1][3] = -1*(this->sectionRotation[0][1]/scaling[1]*this->sectionTranslation[0][3] + this->sectionRotation[1][1]/scaling[1]*this->sectionTranslation[1][3] + this->sectionRotation[2][1]/scaling[1]*this->sectionTranslation[2][3]);
                                        mInverse[2][3] = -1*(this->sectionRotation[0][2]/scaling[2]*this->sectionTranslation[0][3] + this->sectionRotation[1][2]/scaling[2]*this->sectionTranslation[1][3] + this->sectionRotation[2][2]/scaling[2]*this->sectionTranslation[2][3]);
                                        mInverse[3][3] = 1;
                                        /*
                                        for(int ii = 0; ii < 4; ++ii)
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                    
                                                    //std::cout<<sectionTranslation[ii][kk]<<" * "<<sectionRotation[kk][jj]<<" * "<<scaling[jj]<<std::endl;
                                                        this->inverse_transformation[ii][jj] = mInverse[ii][jj];
                                        }*/
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                //std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        this->inverse_transformation[ii][jj] = mInverse[ii][jj];
                                                }
                                                
                                                
                                        }
                                        //this->inverse_transformation = mInverse;
                                        
//                                                mInverse[0][3] = -1*(sectionTranslation[0][3]);
//                                                mInverse[1][3] = -1*(sectionTranslation[1][3]);
//                                                mInverse[2][3] = -1*sectionTranslation[2][3];
//                                                mInverse[3][3] = 1;
                                        
                                        //this->inverse_transformation = mInverse;
                                        
                                        std::cout << "inverse transformation matrix:" << std::endl;
                                        for(int ii = 0; ii < 4; ++ii)
                                        {
                                                std::cout << "[";
                                                for(int jj = 0; jj < 4; ++jj)
                                                {
                                                        if(jj < 3)
                                                                std::cout << this->inverse_transformation[ii][jj] << ",\t";
                                                        else
                                                                std::cout << this->inverse_transformation[ii][jj];
                                                }
                                                std::cout << "]" << std::endl;
                                                
                                        } 
//                                                for(int ii = 0; ii < 3; ++ii)
//                                                 {
//                                                         
//                                                         for(int jj = 0; jj < 3; ++jj){
//                                                             
//                                                             mIn 
//                                                             
//                                                         }
//                                                 }
                                }
                       
                }
        }
        std::cout << "out of transformReading" << std::endl;
        
        inputStream.close();
        return 1;
    
}

#if 0
void BoutonFinder::transformGraphandLandmark(/*std::list< bouton_params * > *boutonParamsList*/)
{
    std::cout<<"In tarnsformGraphLandmark"<<std::endl;
    //double * transformed_landmark = new double[3];

        if(!getTransformationandInverseTransformation(loadhxfile))
        {
            std::cout << "Error getting inverse transform" << std::endl;
        }
        
        
        
        // transform the spatial graph to image coordinates

        // first inverse transform so that the whole cell graph is inline with
        // the section's graph
        //std::list<bouton_params*>::iterator paramIterator;
                                
        /*for(int i = 0; i < NumOfLandmarks; i++)
        {
            std::cout<<" "<<BoutonParamsArray[i].landmark[0]<<" "<<BoutonParamsArray[i].landmark[1]<<" "<<BoutonParamsArray[i].landmark[2]<<std::endl;
            std::cout<<" "<<BoutonParamsArray[i].nearestPoint[0]<<" "<<BoutonParamsArray[i].nearestPoint[1]<<" "<<BoutonParamsArray[i].nearestPoint[2]<<std::endl;
            std::cout<<" "<<BoutonParamsArray[i].transformed_landmark[0]<<" "<<BoutonParamsArray[i].transformed_landmark[1]<<" "<<BoutonParamsArray[i].transformed_landmark[2]<<std::endl;
        }*/
        amira_graph->applyInverseTransformation(inverse_transformation/*, boutonParamsList*/);
        
        /*for(int i = 0; i < NumOfLandmarks; i++)
        {
            std::cout<<" "<<BoutonParamsArray[i].landmark[0]<<" "<<BoutonParamsArray[i].landmark[1]<<" "<<BoutonParamsArray[i].landmark[2]<<std::endl;
            std::cout<<" "<<BoutonParamsArray[i].nearestPoint[0]<<" "<<BoutonParamsArray[i].nearestPoint[1]<<" "<<BoutonParamsArray[i].nearestPoint[2]<<std::endl;
            std::cout<<" "<<BoutonParamsArray[i].transformed_landmark[0]<<" "<<BoutonParamsArray[i].transformed_landmark[1]<<" "<<BoutonParamsArray[i].transformed_landmark[2]<<std::endl;
        }*/
        
        //std::cout << "ok here" << std::endl;

        // now apply the inverse translation to match the image section
        getTranslation(loadhxfile);

        //std::cout << "ok here is it" << std::endl;
        // now apply the inverse scaling
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        std::vector< Vertex * > * vertices = amira_graph->verticesPointer();
        int numOfEdges = edges->size();
        int numofVertices = vertices->size();
        
        //int pixel_value = 0;
        for(int i=0; i<numofVertices; i++)  //for each edge
        {
            
            Vertex * currentVertex = vertices->at(i);
            
            double * coords = (currentVertex)->coordinates;
            coords[X_COORD] = coords[X_COORD] - imageTranslation[X_COORD];
            coords[Y_COORD] = coords[Y_COORD] - imageTranslation[Y_COORD];
            coords[Z_COORD] = coords[Z_COORD] - imageTranslation[Z_COORD];
            
        }
        
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
          
          Edge * currentEdge = edges->at(i);
          std::list< double * >::iterator it;
          
          for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
          {
            double * coords = *it;
                
            //if(coords[IS_BOUTON])
            {
                ImageType::IndexType pixelIndex;
                
                // transform the landmarks as well
                //std::list<bouton_params*>::iterator paramIterator;
                                
                for(int i = 0; i < NumOfLandmarks; i++)
                {
                        //std::cout<<'\t'<<(*paramIterator)->transformed_landmark[0]<<'\t'<<(*paramIterator)->transformed_landmark[1]<<'\t'<<(*paramIterator)->transformed_landmark[2]<<std::endl;
        
                    if((BoutonParamsArray[i].transformed_landmark[0] == coords[0]) && (BoutonParamsArray[i].transformed_landmark[1] == coords[1]) && 
                        (BoutonParamsArray[i].transformed_landmark[2] == coords[2] ))
                    {
                        //std::cout<<'\t'<<(*paramIterator)->transformed_landmark[0]<<'\t'<<(*paramIterator)->transformed_landmark[1]<<'\t'<<(*paramIterator)->transformed_landmark[2]<<std::endl;
        
                        BoutonParamsArray[i].transformed_landmark[0] = coords[X_COORD] - imageTranslation[X_COORD];
                        BoutonParamsArray[i].transformed_landmark[1] = coords[Y_COORD] - imageTranslation[Y_COORD];
                        BoutonParamsArray[i].transformed_landmark[2] = coords[Z_COORD] - imageTranslation[Z_COORD];
                        
                        //std::cout<<'\t'<<BoutonParamsArray[i].transformed_landmark[0]<<'\t'<<BoutonParamsArray[i].transformed_landmark[1]<<'\t'<<BoutonParamsArray[i].transformed_landmark[2]<<std::endl;
        
                    }
                
                    
                }
                
                // Apply image translation before scaling
                coords[X_COORD] = coords[X_COORD] - imageTranslation[X_COORD];
                coords[Y_COORD] = coords[Y_COORD] - imageTranslation[Y_COORD];
                coords[Z_COORD] = coords[Z_COORD] - imageTranslation[Z_COORD];
                
              /*
                int x_pos = rint( coords[X_COORD] / XYSAMPLING );
                int y_pos = rint( coords[Y_COORD] / XYSAMPLING );
                int z_pos = rint( coords[Z_COORD] / ZSAMPLING);
                
                coords[X_COORD] = x_pos;
                coords[Y_COORD] = y_pos;
                coords[Z_COORD] = z_pos;*/
                
                //pixel_value = 255;
                    
                //output_image->SetPixel(pixelIndex, pixel_value);
                
//              std::cout << "Bouton at (x, y): ("<< x_pos << ", " << y_pos << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
            }
          }
          
        }
        
        /*
        transformed_landmark[X_COORD] = transformed_landmark[X_COORD] - imageTranslation[X_COORD];
        transformed_landmark[Y_COORD] = transformed_landmark[Y_COORD] - imageTranslation[Y_COORD];
        transformed_landmark[Z_COORD] = transformed_landmark[Z_COORD] - imageTranslation[Z_COORD];*/
        
        
        // take projections x , y , z using different number of z planes
        
        // get the image plane that this point lies on
        std::cout<<"out of TransformGraphandLandmark"<<std::endl;
        
        //return transformed_landmark;
    
    
    
}
#endif
bool BoutonFinder::getTranslation(char* inputfilename)
{
    std::cout<<"inputfilename: " << inputfilename <<std::endl;
        std::ifstream inputStream(inputfilename);
     
        if(!inputStream.fail())
        {
                //const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                //const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t\n\r\f ";
                
                std::string currentLine;
                //unsigned int line = 0;
                
                /*this->imageTranslation = new double [3];
                for(int ii = 0; ii < 3; ++ii)
                {
                        this->imageTranslation[ii] = 0;
                        
                }*/
               
                
               // bool parameters = 1;
               // bool transform = 0;
 //               bool correctSection = 1;
 //               bool correctPrevSection = 0;
 //               int sectionID = 0;
                //unsigned int brackets = 0, transformBrackets = 0;
                //unsigned int currentIndex = 0;
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
                        
                        if(currentLine.size())
                                if(currentLine.find("setTranslation ", 0) != std::string::npos)
                                {
                                    std::cout << "found correct section translation parameters!" << std::endl;
                                        unsigned int count = 0;
                                        std::string::size_type loc1, loc2, loc3;
                                        loc1 = currentLine.find("setTranslation ", 0);
                                        loc1 += 15;// Lenth of setTransform
                                        loc2 = currentLine.find_first_of(signs, 0);
                                        if(loc2 != std::string::npos)
                                                if(loc2 < loc1)
                                                        loc1 = loc2;
                                        loc2 = currentLine.find_first_of(whitespace, loc1); //ignores last value: is always 1 anyways
                                        while(/*loc2 != std::string::npos &&*/ count < 3)
                                        {
                                            std::cout<<loc1<<'\t'<<loc2<<std::endl;
                                                char * tmp1 = new char[20];
                                                for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                {
                                                    // 
                                                    tmp1[i] = 'f';
                                                }
                                                currentLine.copy(tmp1, loc2 - loc1, loc1);
                                                double ftmp1 = atof(tmp1);
                                                //std::cout << "before!" << std::endl;
                                                //sectionRotation[0][0]= 1;        // amira files are columns after each other
                                                imageTranslation[count%3]= ftmp1;        // amira files are columns after each other
                                                std::cout << imageTranslation[count%3] << std::endl;
                                                std::cout<<ftmp1<<'\t'<<tmp1<<std::endl;
                                                loc3 = loc2;
                                                loc1 = currentLine.find_first_of(numbers, loc3);
                                                loc2 = currentLine.find_first_of(signs, loc3);
                                                if(loc2 != std::string::npos)
                                                        if(loc2 < loc1)
                                                                loc1 = loc2;
                                                        
                                                loc2 = currentLine.find_first_of(whitespace, loc1);
                                                ++count;
                                                delete [] tmp1;
                                        }
                                        
                                }
                }
        }
        inputStream.close();
        
        return true;
    
}

void BoutonFinder::getAxonPoints(ImageType::Pointer image)
{
    double axon1_minz = 0xffffffffffffffff;
    double axon1_maxz = 0;
    double axon2_minz = 0xffffffffffffffff;
    double axon2_maxz = 0;
    double axon3_minz = 0xffffffffffffffff;
    double axon3_maxz = 0;
    
    double axon1_min[3] = {0xffff,0xffff,0xffff};
    double axon1_max[3] = {0xffff,0xffff,0xffff};
    double axon2_min[3] = {0xffff,0xffff,0xffff};
    double axon2_max[3] = {0xffff,0xffff,0xffff};
    double axon3_min[3] = {0xffff,0xffff,0xffff};
    double axon3_max[3] = {0xffff,0xffff,0xffff};
    // 
    //int num_of_boutons = NumOfLandmarks;
    
    std::vector< Edge * > * edges = amira_graph->edgesPointer();

    unsigned int numOfEdges = edges->size();
    unsigned int end_x = image->GetLargestPossibleRegion().GetSize(0)*XYSAMPLING;
    unsigned int end_y = image->GetLargestPossibleRegion().GetSize(1)*XYSAMPLING;
    unsigned int end_z = image->GetLargestPossibleRegion().GetSize(2)*ZSAMPLING;

    std::cout<< end_x << " " << end_y << " " << end_z<<std::endl;
    
    double x_pos, y_pos, z_pos;
    
        if(numOfEdges == 0)
        {
            std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
            //return -1;
        }
        else
        {
// #pragma omp parallel for schedule(dynamic,1)
            for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
            {
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;
                std::list< double * >::iterator next_it;

//                               std::cout<< "1 "  << std::endl;
                //if( (currentEdge->label == Axon_1_ID) || (currentEdge->label == Axon_2_ID) || (currentEdge->label == Axon_3_ID ) ) 
                if( ( currentEdge->label == Axon_1_ID) || ( currentEdge->label == Axon_2_ID) || (currentEdge->label == Axon_3_ID ))
		{
		  
		  

                    //for every point along edge
                    for(edge_it = currentEdge->edgePointCoordinates.begin(); edge_it != currentEdge->edgePointCoordinates.end(); edge_it++)
                    {
    //                                   std::cout<< "2 " << OtsuThreshold << std::endl;

                        double * coords = *edge_it;

                        x_pos = (coords[X_COORD]);
                        y_pos = (coords[Y_COORD]);
                        z_pos = (coords[Z_COORD]);

                        

                        //if(z_pos == z)
                        if(((x_pos >= (0+X_OFFSET)) &&  (x_pos < (end_x-X_OFFSET))) && ((y_pos >= (0+Y_OFFSET)) &&  (y_pos < (end_y-Y_OFFSET))) && ((z_pos >= 0) &&  (z_pos < (end_z-Z_OFFSET))))
                        {
                            std::cout<< "x " << x_pos << " y " << y_pos << std::endl;
                                
				if((currentEdge->label == Axon_1_ID) && (z_pos < axon1_minz))
				{
				  axon1_minz = z_pos;
				  axon1_min[0] = x_pos;
				  axon1_min[1] = y_pos;
				  axon1_min[2] = z_pos;
				  
				}
				if((currentEdge->label == Axon_1_ID) && (z_pos > axon1_maxz))
				{
				  axon1_maxz = z_pos;
				  axon1_max[0] = x_pos;
				  axon1_max[1] = y_pos;
				  axon1_max[2] = z_pos;
				}
				
				if((currentEdge->label == Axon_2_ID) && (z_pos < axon2_minz))
				{
				  axon2_minz = z_pos;
				  axon2_min[0] = x_pos;
				  axon2_min[1] = y_pos;
				  axon2_min[2] = z_pos;
				  
				}
				if((currentEdge->label == Axon_2_ID) && (z_pos > axon2_maxz))
				{
				  axon2_maxz = z_pos;
				  axon2_max[0] = x_pos;
				  axon2_max[1] = y_pos;
				  axon2_max[2] = z_pos;
				}
				
				if((currentEdge->label == Axon_3_ID) && (z_pos < axon3_minz))
				{
				  axon3_minz = z_pos;
				  axon3_min[0] = x_pos;
				  axon3_min[1] = y_pos;
				  axon3_min[2] = z_pos;
				  
				}
				if((currentEdge->label == Axon_3_ID) && (z_pos > axon3_maxz))
				{
				  axon3_maxz = z_pos;
				  axon3_max[0] = x_pos;
				  axon3_max[1] = y_pos;
				  axon3_max[2] = z_pos;
				}
			  
			  std::cout<< "done with axon stuff "  << std::endl;
                            
                                BoutonParamsArray[NumOfLandmarks].transformed_landmark[0] = coords[0];
                                BoutonParamsArray[NumOfLandmarks].transformed_landmark[1] = coords[1];
                                BoutonParamsArray[NumOfLandmarks].transformed_landmark[2] = coords[2];
                                BoutonParamsArray[NumOfLandmarks].transformed_landmark[3] = pos;
                                //BoutonParamsArray[NumOfLandmarks].isBouton = 0;
                                
                                
				// Apply image translation before scaling
                                BoutonParamsArray[NumOfLandmarks].landmark[0] = BoutonParamsArray[NumOfLandmarks].transformed_landmark[0] + imageTranslation[X_COORD];
                                BoutonParamsArray[NumOfLandmarks].landmark[1] = BoutonParamsArray[NumOfLandmarks].transformed_landmark[1] + imageTranslation[Y_COORD];
                                BoutonParamsArray[NumOfLandmarks].landmark[2] = BoutonParamsArray[NumOfLandmarks].transformed_landmark[2] + imageTranslation[Z_COORD];
                    
                                
                                // transform to get the original point
                                double oldCoords[4], newCoords[4];
                                    for(int ii = 0; ii < 3; ++ii)
                                    {
                                            oldCoords[ii] = BoutonParamsArray[NumOfLandmarks].landmark[ii];
                                            newCoords[ii] = 0;
                                    }
                                    oldCoords[3] = 1;
                                    newCoords[3] = 1;
                                    for(int ii = 0; ii < 3; ++ii)
                                            for(int jj = 0; jj < 4; ++jj)
                                                    newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
                                    
                                    for(int ii = 0; ii < 3; ++ii)
                                    {
                                            BoutonParamsArray[NumOfLandmarks].landmark[ii] = newCoords[ii];
                                            //BoutonParamsArray[NumOfLandmarks].nearestPoint[ii] = newCoords[ii];
                                    }
                                    
                                //BoutonParamsArray[NumOfLandmarks].nearestPoint[3] = pos;
				
				
				
                                NumOfLandmarks++;
                                
                                
                           
                            
                            
                            
                        }
                        
                    }
                    
                }
                
            }
            
            //write the min and max z coors for all pts as it is global and will be the same fr the entire roi
            /*for (int i = 0 ; i < NumOfLandmarks; i++)
	    {
	      BoutonParamsArray[i].Axon1MinZ[0] = axon1_min[0] ;
	      BoutonParamsArray[i].Axon1MinZ[1] = axon1_min[1] ;
	      BoutonParamsArray[i].Axon1MinZ[2] = axon1_min[2];
	      
	      BoutonParamsArray[i].Axon1MaxZ[0] = axon1_max[0] ;
	      BoutonParamsArray[i].Axon1MaxZ[1] = axon1_max[1] ;
	      BoutonParamsArray[i].Axon1MaxZ[2] = axon1_max[2];
	      
	      BoutonParamsArray[i].Axon2MinZ[0] = axon2_min[0] ;
	      BoutonParamsArray[i].Axon2MinZ[1] = axon2_min[1] ;
	      BoutonParamsArray[i].Axon2MinZ[2] = axon2_min[2];
	      
	      BoutonParamsArray[i].Axon2MaxZ[0] = axon2_max[0] ;
	      BoutonParamsArray[i].Axon2MaxZ[1] = axon2_max[1] ;
	      BoutonParamsArray[i].Axon2MaxZ[2] = axon2_max[2];
	      
	      BoutonParamsArray[i].Axon3MinZ[0] = axon3_min[0] ;
	      BoutonParamsArray[i].Axon3MinZ[1] = axon3_min[1] ;
	      BoutonParamsArray[i].Axon3MinZ[2] = axon3_min[2];
	      
	      BoutonParamsArray[i].Axon3MaxZ[0] = axon3_max[0] ;
	      BoutonParamsArray[i].Axon3MaxZ[1] = axon3_max[1] ;
	      BoutonParamsArray[i].Axon3MaxZ[2] = axon3_max[2];
	      
	    }*/
	    
	    for (int i = 0 ; i < NumOfLandmarks; i++)
	    {
	      // xz angle
	      if(axon1_max[0] != axon1_min[0])
	      {
		  BoutonParamsArray[i].Axon_1_XZ_Angle = 180 - abs(atan2((axon1_max[2] - axon1_min[2]) , (axon1_max[0] - axon1_min[0]) ) * 180 / PI); 
	      }
	      else
	      {
		if((axon1_max[0] == 0xffff) && (axon1_min[0] == 0xffff))
		  BoutonParamsArray[i].Axon_1_XZ_Angle = 0xffff;
		else
		  BoutonParamsArray[i].Axon_1_XZ_Angle = 90;
	      }
	      
	      if(axon2_max[0] != axon2_min[0])
	      {
		  BoutonParamsArray[i].Axon_2_XZ_Angle = 180 - abs(atan2((axon2_max[2] - axon2_min[2]) , (axon2_max[0] - axon2_min[0]) ) * 180 / PI); 
	      }
	      else
	      {
		if((axon2_max[0] == 0xffff) && (axon2_min[0] == 0xffff))
		  BoutonParamsArray[i].Axon_2_XZ_Angle = 0xffff;
		else
		  BoutonParamsArray[i].Axon_2_XZ_Angle = 90;
	      }
	      
	      if(axon3_max[0] != axon3_min[0])
	      {
		  BoutonParamsArray[i].Axon_3_XZ_Angle = 180 - abs(atan2((axon3_max[2] - axon3_min[2]) , (axon3_max[0] - axon3_min[0]) ) * 180 / PI); 
	      }
	      else
	      {
		if((axon3_max[0] == 0xffff) && (axon3_min[0] == 0xffff))
		  BoutonParamsArray[i].Axon_3_XZ_Angle = 0xffff;
		else
		  BoutonParamsArray[i].Axon_3_XZ_Angle = 90;
	      }
	    }
           
            
        }       
//           
    
    
}

/****************************************************************************/
/*writeTestImage()  writes 2D planes with bouton pixels identified          */
/****************************************************************************/
void BoutonFinder::writeTestImage()
{
        ImageType::Pointer output_image = ImageType::New(); 
        
        ImageType::IndexType start;
        start[0] =   0;  // first index on X
        start[1] =   0;  // first index on Y
        start[2] =   0;  // first index on Z
        
        ImageType::SizeType  size;
        size[0]  = original_image->GetLargestPossibleRegion().GetSize(0);  // size along X
        size[1]  = original_image->GetLargestPossibleRegion().GetSize(1);  // size along Y
        size[2]  = original_image->GetLargestPossibleRegion().GetSize(2);  // size along Z
        
        ImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );
        
        output_image->SetRegions( region );
        output_image->Allocate(); 
        
        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        int numOfEdges = edges->size();
        
        int pixel_value = 0;
        
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
          
          Edge * currentEdge = edges->at(i);
          std::list< double * >::iterator it;
          
          for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
          {
            double * coords = *it;
                
            if(coords[IS_BOUTON])
            {
                ImageType::IndexType pixelIndex;
              
                int x_pos = rint( coords[X_COORD] / XYSAMPLING );
                int y_pos = rint( coords[Y_COORD] / XYSAMPLING );
                int z_pos = rint( coords[Z_COORD] / ZSAMPLING);
                
                pixelIndex[0] = x_pos;
                pixelIndex[1] = y_pos;
                pixelIndex[2] = z_pos;
                
                pixel_value = 255;
                    
                output_image->SetPixel(pixelIndex, pixel_value);
                
//              std::cout << "Bouton at (x, y): ("<< x_pos << ", " << y_pos << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
            }
          }
          
        }
        
        //DILATE POINTS FOR NICER LOOK
//      DilateFilterType::Pointer grayscale_dilate = DilateFilterType::New();
//      
//      StructuringElementType  structuring_element;
//      structuring_element.SetRadius( 2 );
//      structuring_element.CreateStructuringElement();
// 
//      grayscale_dilate->SetKernel( structuring_element );
//      grayscale_dilate->SetInput( output_image );
//      output_image = grayscale_dilate->GetOutput();
//      output_image->Update();


        
        //std::cout << "Begin writing" << std::endl;
        //writeImagePlanes(output_image, "raster_amira_slice");
  
};

#if 0
/****************************************************************************/
/*WriteImagePlanes() evokes writing of 2D png planes         */
/****************************************************************************/

void BoutonFinder::writeImagePlanes(ImageType::Pointer input_image, char output_file[1024])
{
        #ifdef DEBUG
          std::cout<< "Generate Writer !!" << std::endl;
        #endif

          Writer2DType::Pointer writer = Writer2DType::New();
          writer->SetInput( input_image );

          NameGeneratorType::Pointer writer_name_gen = NameGeneratorType::New();

          std::string format = output_file;
          //format += "_seg";
          format += "%03d.";
          format += "tif";

          writer_name_gen->SetSeriesFormat( format.c_str() );

          ImageType::RegionType region = input_image->GetLargestPossibleRegion();
          ImageType::IndexType start = region.GetIndex();
          ImageType::SizeType size = region.GetSize();

          const unsigned int first_slice = start[2];
          const unsigned int last_slice  = start[2] + size[2] -1;

          writer_name_gen->SetStartIndex( first_slice );
          writer_name_gen->SetEndIndex( last_slice );
          writer_name_gen->SetIncrementIndex( 1 );

        #ifdef DEBUG
          std::cout << "-------------------------------------------------------------------" << std::endl;
          std::cout<< "Start Writing !!" << std::endl;
          std::cout << "-------------------------------------------------------------------" << std::endl;
        #endif

          writer->SetFileNames( writer_name_gen->GetFileNames() );

          try
            {
              writer->Update();
            }
          catch (itk::ExceptionObject & err )
            {
              std::cerr << "SegWriterExceptionObject caught !" << std::endl;
              std::cerr << err << std::endl;
            }

        #ifdef DEBUG
          std::cout<< "Writing Done !!" << std::endl;
        #endif
}; 
#endif

/************************************************************************************/
/*markBoutons() refines the conditions for what constitutes a bouton                */
/************************************************************************************/
void BoutonFinder::markBoutons()
{
        std::cout<< "In markBoutons" << std::endl;
  
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        int numOfEdges = edges->size();
        //bool inBouton = false;
        
        //int pixel_value = 0;
        float avg_std_dev = getAverageStandardDeviation();
        
        std::cout<< "avg_std_dev: " << avg_std_dev << std::endl;
        
        for(int i=0; i<numOfEdges; i++)  //for each edge
        {
          
          Edge * currentEdge = edges->at(i);
          std::list< double * >::iterator it;
          
          for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
          {
                double * coords = *it;
            
                //float brightness = coords[LOCAL_BRIGHTNESS];
                float std_dev = coords[LOCAL_SIGMA];
                //float avg_bright = coords[AVERAGE_BRIGHTNESS];
                //float avg_rad = coords[AVERAGE_SURFACE];
                float radius = coords[SURFACE];
                
//              std::cout << "Point at (x, y): ("<< coords[X_COORD] << ", " << coords[Y_COORD] << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
                
                
//              if( (brightness > (avg_bright + 10)) && (radius > avg_rad) && (radius >= 3.5) ) 
//                coords[IS_BOUTON] = 1;
//              else
//                coords[IS_BOUTON] = 0;
                
                
                if( coords[IS_BOUTON] && radius > 0.5  && std_dev > avg_std_dev*0.66 ) 
                  coords[IS_BOUTON] = 1;
                else
                  coords[IS_BOUTON] = 0;

          }
          
        }
  
};

#if 0
/************************************************************************************/
/*MarkLocalMaximums()  marks every local brightness max as a potential bouton       */
/************************************************************************************/
void BoutonFinder::markLocalMaximums()
{

        std::cout<< "In MarkLocalMaximums!! " << std::endl;

                
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();
        
        float bright_offset = 5;
        float radius_offset = 0;
        
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
                

                //for every point along edge
                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                {

                        double * current_coords = *edge_it;
                        //float diameter = current_coords[SURFACE];
                        //float sum = 0;
                        float average = 0;
                        //long counter = 1;
                        bool isRadiusMax = false;
                        
                        std::list< double * >::iterator next_it = edge_it;              
                        std::list< double * >::iterator prev_it = edge_it;
                        ++next_it;
                        --prev_it;
                        
                        std::list< double * >::iterator next_next_it = next_it;         
                        std::list< double * >::iterator prev_prev_it = prev_it;
                        ++next_next_it;
                        --prev_prev_it;
                        
                        if(edge_it == currentEdge->edgePointCoordinates.begin())         //if at beginning of edge
                        {
                                //double * next_coords = *next_it;
                                
                                
                        }
                        else if(next_it == currentEdge->edgePointCoordinates.end())       //if at end of edge
                        {
                                //double * prev_coords = *prev_it;
                          
                                
                        }
                        else if(prev_it == currentEdge->edgePointCoordinates.begin() && edgeSize >= 4)  //second bouton from beginning
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                                double * next_next_coords = *next_next_it;
                                
                          
                                if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
                                  current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
                                {
                                    
                                    average = averageLocalValues(currentEdge, edge_it, SURFACE, 3);
                                    
                                    if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
                                        current_coords[SURFACE] >= next_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(next_coords[SURFACE] >= next_next_coords[SURFACE] && 
                                        next_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    
//                                  if(isRadiusMax)
                                    if(current_coords[SURFACE] > average+radius_offset)
                                    {
                                        
                                        average = averageLocalValues(currentEdge, edge_it, LOCAL_BRIGHTNESS, 1);
                                      
                                        if(current_coords[LOCAL_BRIGHTNESS] > (average + bright_offset)) 
                                          current_coords[IS_BOUTON] = 1;
                                        else
                                          current_coords[IS_BOUTON] = 0;
                                    }
                                }
                                else
                                  current_coords[IS_BOUTON] = 0;
                        }
                        else if(next_next_it == currentEdge->edgePointCoordinates.end() && edgeSize >= 4)  //second bouton from end
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                                double * prev_prev_coords = *prev_prev_it;
                                
                                
                          
                                if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
                                  current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
                                {
                                    
                                    average = averageLocalValues(currentEdge, edge_it, SURFACE, 3);
                                    
                                    
                                    if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
                                        current_coords[SURFACE] >= next_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(prev_coords[SURFACE] >= prev_prev_coords[SURFACE] && 
                                        prev_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    
//                                  if(isRadiusMax)
                                    if(current_coords[SURFACE] > average+radius_offset)
                                    {

                                        
                                        average = averageLocalValues(currentEdge, edge_it, LOCAL_BRIGHTNESS, 1);
                                      
                                        if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset))
                                          current_coords[IS_BOUTON] = 1;
                                        else
                                          current_coords[IS_BOUTON] = 0;
                                    }
                                }
                                else
                                  current_coords[IS_BOUTON] = 0;
                                
                        }
                        else if(edgeSize >= 5)                                                          //in middle of edge
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                                double * next_next_coords = *next_next_it;
                                double * prev_prev_coords = *prev_prev_it;
                                
                          
                                if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
                                  current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS])
                                {

                                    average = averageLocalValues(currentEdge, edge_it, SURFACE, 3);
                                    
                                    if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
                                        current_coords[SURFACE] >= next_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(next_coords[SURFACE] >= next_next_coords[SURFACE] && 
                                        next_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(prev_coords[SURFACE] >= prev_prev_coords[SURFACE] && 
                                        prev_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    
//                                  if(isRadiusMax)
                                    if(current_coords[SURFACE] > average+radius_offset)
                                    {

                                          average = averageLocalValues(currentEdge, edge_it, LOCAL_BRIGHTNESS, 1);
                                          
                                          if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset))
                                            current_coords[IS_BOUTON] = 1;
                                          else
                                            current_coords[IS_BOUTON] = 0;
                                    }
                                }
                                else
                                  current_coords[IS_BOUTON] = 0;
                                
                        }
        
                }

        }
};




/************************************************************************************/
/*SetRadius() measures the local brightness and radius at every point               */
/************************************************************************************/
int BoutonFinder::setRadiusAndLocalBrightness(ImageType::Pointer image)
{

        std::cout<< "In SetRadiusAndLocalBrightness " << std::endl;

        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
        unsigned int numOfEdges = edges->size();
        unsigned int end_x = image->GetLargestPossibleRegion().GetSize(0);
        unsigned int end_y = image->GetLargestPossibleRegion().GetSize(1);
        unsigned int end_z = image->GetLargestPossibleRegion().GetSize(2);
        
        std::cout<< "xend " << end_x << "yend " << end_y << "zend" << end_z << std::endl;
        
        int x_pos, y_pos, z_pos;

        Image2DType::Pointer image_plane_5;
        Image2DType::Pointer image_plane_3;
        
        //int count = 1;

        for(unsigned int z = 0; z < end_z; z++)
        {
                image_plane_5 = getImagePlane(z, 5, image);
                image_plane_3 = getImagePlane(z,3,image);
        
                
                if(numOfEdges == 0)
                {
                        std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                        return -1;
                }
                else
                {
// #pragma omp parallel for schedule(dynamic,1)
                        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                        {                       
                                Edge * currentEdge = edges->at(pos);
                                std::list< double * >::iterator edge_it;        
                                std::list< double * >::iterator next_it;

 //                               std::cout<< "1 "  << std::endl;
                                
                                //for every point along edge
                                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                                {
 //                                   std::cout<< "2 " << OtsuThreshold << std::endl;
                                    
                                        double * coords = *edge_it;
                                        
                                        x_pos = rint(coords[X_COORD]/XYSAMPLING);
                                        y_pos = rint(coords[Y_COORD]/XYSAMPLING);
                                        z_pos = rint(coords[Z_COORD]/ZSAMPLING);
                                        
                                        //std::cout<< "x " << x_pos << "y " << y_pos << std::endl;
                                        
                                        //if(z_pos == z)
                                        if(((x_pos >= 0) &&  (x_pos < end_x)) && ((y_pos >= 0) &&  (y_pos < end_y)) && (z_pos == z))
                                        {
 //                                               std::cout<< "3 " << OtsuThreshold << std::endl;
                                                /*****************************************************************/
                                                /*set radius info                                                */
                                                /*****************************************************************/
                                          
                                                float x0 = x_pos;
                                                float y0 = y_pos;
                                                PixelType threshold = 0;

                                                threshold = calculateOtsuThreshold(image_plane_3, x0,y0,z_pos);
                                                //threshold = coords[THRESHOLD];
                                                
                                                unsigned int nr_of_rays = 20;
                                                //float ray_length = 0.5;
                                                int index = 0;

                                                std::vector<VECTOR *> vectors;
                                                VECTOR * adjustment_vect;
                                                float distance = 0;

                                                //x0 = x0 + 0.5;
                                                //y0 = y0 + 0.5;
                                                std::cout<< "x0" << x0 << "y0" << y0 << "z0" << z_pos << std::endl;

                                                for(unsigned n=1; n<=nr_of_rays; n++)
                                                {
                                                        VECTOR * tmp;
                                                        std::cout<< "z pos " << z_pos << std::endl;
                                                        //tmp = sendRay(x0, y0, ray_length, nr_of_rays, threshold, n, image_plane_3);
                                                        vectors.push_back(tmp);
                                                       
                                                }
                                
                                                index = addUpOpposingRaysAndFindShortest(vectors, &adjustment_vect, &distance); 
#ifdef DEBUG
                                                std::cout<< "Distance: " << distance <<std::endl;
#endif
                                                coords[SURFACE] = distance/**XYSAMPLING*/;
                                                
                                                /*****************************************************************/
                                                /*adjust midline based on radius info                            */
                                                /*****************************************************************/
                                                next_it = edge_it;
                                                next_it++;
                                                
                                                if(edge_it != currentEdge->edgePointCoordinates.begin() && next_it != currentEdge->edgePointCoordinates.end())
                                                {
                                                  float new_x = coords[X_COORD] + adjustment_vect->coords[X_COORD] * XYSAMPLING;
                                                  float new_y = coords[Y_COORD] + adjustment_vect->coords[Y_COORD] * XYSAMPLING;
                                                  
                                                  if((new_x/XYSAMPLING) < end_x-1 && (new_y/XYSAMPLING) < end_y-1)
                                                  {
                                                          coords[X_COORD] = new_x;
                                                          coords[Y_COORD] = new_y;
                                                  }
                                                }
                                                
                                                /*****************************************************************/
                                                /*set brightness info                                            */
                                                /*****************************************************************/
                                                x0 = rint(coords[X_COORD]/XYSAMPLING);
                                                y0 = rint(coords[Y_COORD]/XYSAMPLING);
                                        
                                                coords[LOCAL_BRIGHTNESS] = calculateLocalBrightness(image_plane_5, x0, y0, 1);
                                                coords[LOCAL_SIGMA] = calculateLocalStandardDeviation(image_plane_5, x0, y0, 2);
                                        }
                                        
                                        
                                }
                        }
                }
        }       

        smoothRadii();
        return 0;
};

#endif

/************************************************************************************/
/*SetAverageThreshold() calculate the threshold to use for a certain # of points    */
/************************************************************************************/
PixelType BoutonFinder::calculateOtsuThreshold(Image2DType::Pointer image, float x0, float y0, float z0 )
{
        std::cout << "In calculateThreshold" << std::endl;
    
        Image2DType::RegionType threshold_region;
        Image2DType::IndexType threshold_index;
        Image2DType::SizeType threshold_size;

        threshold_index[0] = x0-6*7;
        threshold_index[1] = y0-6*7;
        threshold_index[2] = z0;

        threshold_size[0] = 12*7;
        threshold_size[1] = 12*7;
        threshold_size[2] = 1;

        if(threshold_index[0]<=0)
        {       
                threshold_size[0] = threshold_size[0] + threshold_index[0];
                threshold_index[0]=0;
        }
        if(threshold_index[1]<=0)
        {       
                threshold_size[1] = threshold_size[1] + threshold_index[1];
                threshold_index[1]=0;
        }
        if(threshold_index[0] + threshold_size[0] >= image->GetLargestPossibleRegion().GetSize(0))
        {
                threshold_size[0] = image->GetLargestPossibleRegion().GetSize(0)-1-threshold_index[0];
        }
        if(threshold_index[1] + threshold_size[1] >= image->GetLargestPossibleRegion().GetSize(1))      
        {
                threshold_size[1] = image->GetLargestPossibleRegion().GetSize(1)-1-threshold_index[1];
        }

#ifdef DEBUG
        std::cout<< "x0= " << x0 << " y0= " << y0 << " z0= " << z0 << " index= " << threshold_index << " size= " << threshold_size << std::endl << std::flush;
        std::cout<< "xsize= " << image->GetLargestPossibleRegion().GetSize(0)<< " ysize= " << image->GetLargestPossibleRegion().GetSize(1) << std::endl << std::flush;
#endif

        threshold_region.SetIndex(threshold_index);
        threshold_region.SetSize(threshold_size);
        
        Image2DType::RegionType iterator_region = image->GetLargestPossibleRegion();
        Image2DType::RegionType output_region = threshold_region;
            
        Image2DType::Pointer output_image = Image2DType::New(); 
        output_image->SetRegions( output_region );
        output_image->Allocate(); 
        
        Iterator2DType og_iterator(image, iterator_region);
        Iterator2DType output_iterator(output_image, output_region);
        
        for (og_iterator.GoToBegin(), output_iterator.GoToBegin();!output_iterator.IsAtEnd(); ++og_iterator, ++output_iterator ) 
        {
            unsigned char value = og_iterator.Get();      
            output_iterator.Set((unsigned char)value);
        }
    
        OtsuThresholdType::Pointer OtsuThreshold = OtsuThresholdType::New(); 
        
        OtsuThreshold->SetImage(output_image);
       
        OtsuThreshold->Compute();
        
        
        return (OtsuThreshold->GetThreshold());
        

}

#if 0
/************************************************************************************/
/*SetAverageThreshold() calculate the threshold to use for a certain # of points    */
/************************************************************************************/
float BoutonFinder::setAverageThreshold(int numOfPointsToAverage)
{

        std::cout<< "In Set Average threshold " << std::endl;
        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
        unsigned int numOfEdges = edges->size();
        
        unsigned int end_x = original_image->GetLargestPossibleRegion().GetSize(0);
        unsigned int end_y = original_image->GetLargestPossibleRegion().GetSize(1);
        unsigned int end_z = original_image->GetLargestPossibleRegion().GetSize(2);
        
#ifdef DEBUG
        std::cout<< "end_x "<< end_x << " end_y " << end_y << " end_z "<< end_z << std::endl;
#endif     
        if(numOfEdges == 0)
        {
                std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                return -1;
        }
        else
        {
// #pragma omp parallel for schedule(dynamic,1)
                for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                {                       
                        int pointCount = 1;
                        float threshold_sum = 0, current_threshold=0;
                        Edge * currentEdge = edges->at(pos);
                        std::list< double * >::iterator primary_edge_it, secondary_edge_it;
                        std::queue<float> threshold_q;
                        bool isFirstLoop = true;
                        int minimum_threshold = 11;
                        
                        

                        //for every point along edge
                        for(primary_edge_it = currentEdge->edgePointCoordinates.begin(), secondary_edge_it = primary_edge_it;
                            primary_edge_it != currentEdge->edgePointCoordinates.end(); 
                            primary_edge_it++) 
                        {
                            double * primaryPoint = *primary_edge_it;
                                
                            int x_pos = rint( primaryPoint[X_COORD] / XYSAMPLING );
                            int y_pos = rint( primaryPoint[Y_COORD] / XYSAMPLING );
                            int z_pos = rint( primaryPoint[Z_COORD] / ZSAMPLING );
                            // if the point is on the z planes
                            if(((x_pos >= 0) &&  (x_pos < end_x)) && ((y_pos >= 0) &&  (y_pos < end_y)) && ((z_pos >= 0) &&  (z_pos < end_z)))
                            {
                                
                                
//                              std::cout<< "Position "<< primaryPoint[X_COORD] << "  " << primaryPoint[Y_COORD]<< "  "<< primaryPoint[Z_COORD] << std::endl;
//                              std::cout<< "Position "<< x_pos << "  " << y_pos<< "  "<< z_pos << std::endl;
                                

                                if(pointCount < numOfPointsToAverage)
                                {
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  threshold_sum += current_threshold;
                                }
                                else if(isFirstLoop)
                                {
                                  isFirstLoop = false;
                                  double * secondaryPoint;
                                  
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  threshold_sum += current_threshold;
                                  
                                  float average_threshold = threshold_sum/numOfPointsToAverage;
                                  
                                  for(int i=0; i<=((numOfPointsToAverage/2)+1); i++, secondary_edge_it++)
                                  {
                                    secondaryPoint = *secondary_edge_it;
                                    
                                    if(average_threshold > minimum_threshold)
                                      secondaryPoint[THRESHOLD] = average_threshold;
                                    else
                                      secondaryPoint[THRESHOLD] = minimum_threshold;
                                  }
                                  
                                }
                                else
                                {
                                  double * secondaryPoint = *secondary_edge_it;
                                  
                                  current_threshold = calculateThreshold(original_image, x_pos, y_pos, z_pos);
                                  threshold_q.push(current_threshold);
                                  threshold_sum -= threshold_q.front();
                                  threshold_q.pop();
                                  threshold_sum += current_threshold;
                                  
                                  float average_threshold = threshold_sum/numOfPointsToAverage;
                                  
                                  if(average_threshold > minimum_threshold)
                                    secondaryPoint[THRESHOLD] = average_threshold;
                                  else
                                    secondaryPoint[THRESHOLD] = minimum_threshold;
                                }
                                
                                pointCount++;
                            }
                        }
                        
                        float average_threshold = threshold_sum/numOfPointsToAverage;
                        secondary_edge_it++;
                        
                        for(;secondary_edge_it != primary_edge_it;secondary_edge_it++)
                        {
                          double * secondaryPoint = *secondary_edge_it;
                          
                          if(average_threshold > minimum_threshold)
                            secondaryPoint[THRESHOLD] = average_threshold;
                          else
                            secondaryPoint[THRESHOLD] = minimum_threshold;
                          
                        }
                        
                }
        }       

        return 0;
};

#endif

/************************************************************************************/
/*image_plane GetImagePlane(z, image) creates a projection image from z local planes*/
/************************************************************************************/
Image2DType::Pointer BoutonFinder::getImagePlane(int z, int depth, ImageType::Pointer input_image)
{
#ifdef DEBUG    
        std::cout<< "In GetImagePlane" << std::endl;
#endif
  
  
        Image2DType::Pointer image_plane = Image2DType::New();

        Image2DType::IndexType target_index;
        Image2DType::SizeType target_size;

        target_index[0] = 0;
        target_index[1] = 0;

        
               
        
        target_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
        target_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
        
        
        
        #ifdef DEBUG    
        //std::cout<< target_size[0] << std::endl;
        #endif

        Image2DType::RegionType target_region;

        target_region.SetSize(target_size);
        target_region.SetIndex(target_index) ;

        image_plane->SetRegions(target_region);

        image_plane->Allocate();
        image_plane->FillBuffer(0);

        Iterator2DType target_it(image_plane, image_plane->GetLargestPossibleRegion()); 

        int end = input_image->GetLargestPossibleRegion().GetSize(2);
        int displacement = (depth-1)/2;
        
        if(depth == 1)
        {
            //
            ImageType::RegionType input_region;
            ImageType::IndexType input_index;
            ImageType::SizeType input_size;
    
            input_index[0] = 0;
            input_index[1] = 0;
            input_index[2] = z;
                    
            input_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
            input_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
            input_size[2] = 1;
    
            input_region.SetIndex(input_index);
            input_region.SetSize(input_size);
            
            ConstIteratorType source_it(input_image, input_region);
    
            for(source_it.GoToBegin(), target_it.GoToBegin(); !source_it.IsAtEnd(); ++source_it, ++target_it)
            {
                   // unsigned char projection_value = 0;
                    unsigned char value = 0;
                    
                    //projection_value = target_it.Get();
                    value = source_it.Get();
                    
                    //if(value > projection_value)
                            target_it.Set(value);
            }
            
        }
        
        else
        {
            for(int margin = z-displacement; margin <= z+displacement && margin < end; margin++)
            {
                    if(margin < 0)
                            margin = 0;

                    ImageType::RegionType input_region;
                    ImageType::IndexType input_index;
                    ImageType::SizeType input_size;
            
                    input_index[0] = 0;
                    input_index[1] = 0;
                    input_index[2] = margin;
                            
                    input_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
                    input_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
                    input_size[2] = 1;
            
                    input_region.SetIndex(input_index);
                    input_region.SetSize(input_size);
                    
                    ConstIteratorType source_it(input_image, input_region);
            
                    for(source_it.GoToBegin(), target_it.GoToBegin(); !source_it.IsAtEnd(); ++source_it, ++target_it)
                    {
                            unsigned char projection_value = 0;
                            unsigned char value = 0;
                            
                            projection_value = target_it.Get();
                            value = source_it.Get();
                            
                            if(value > projection_value)
                                    target_it.Set(value);
                    }
            }
        }
        
        

        return image_plane;
};
/************************************************************************************/
/*threshold CalculateThreshold(Image2DType::Pointer image_plane, float x0, float y0)*/
/************************************************************************************/
float BoutonFinder::calculateThreshold(ImageType::Pointer image, float x0, float y0, float z0)
{
    
#ifdef DEBUG    
        std::cout<< "In Calculate Threshold" << std::endl;
#endif
        ImageType::RegionType threshold_region;
        ImageType::IndexType threshold_index;
        ImageType::SizeType threshold_size;

        threshold_index[0] = x0-6;
        threshold_index[1] = y0-6;
        threshold_index[2] = z0;

        threshold_size[0] = 12;
        threshold_size[1] = 12;
        threshold_size[2] = 1;

        if(threshold_index[0]<=0)
        {       
                threshold_size[0] = threshold_size[0] + threshold_index[0];
                threshold_index[0]=0;
        }
        if(threshold_index[1]<=0)
        {       
                threshold_size[1] = threshold_size[1] + threshold_index[1];
                threshold_index[1]=0;
        }
        if(threshold_index[0] + threshold_size[0] >= image->GetLargestPossibleRegion().GetSize(0))
        {
                threshold_size[0] = image->GetLargestPossibleRegion().GetSize(0)-1-threshold_index[0];
        }
        if(threshold_index[1] + threshold_size[1] >= image->GetLargestPossibleRegion().GetSize(1))      
        {
                threshold_size[1] = image->GetLargestPossibleRegion().GetSize(1)-1-threshold_index[1];
        }

#ifdef DEBUG
        std::cout<< "x= " << x0<< " y= " << y0 << " z= " << z0 << " index= " << threshold_index << " size= " << threshold_size << std::endl << std::flush;
        std::cout<< "x= " << image->GetLargestPossibleRegion().GetSize(0)<< " y= " << image->GetLargestPossibleRegion().GetSize(1) << std::endl << std::flush;
#endif

        threshold_region.SetIndex(threshold_index);
        threshold_region.SetSize(threshold_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIteratorType histo_it(image, threshold_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }

        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();

        CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();

        covariance_algorithm->SetInputSample( histo_sample );
        covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
        covariance_algorithm->Update();

        float mean = mean_algorithm->GetOutput()->GetElement(0);
        CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
        //float variance = *(covariance.operator[](0));
       // float standard_deviation = std::sqrt(variance);

#ifdef DEBUG    
        std::cout<< "In Calculate Threshold" << mean << std::endl << std::flush;
#endif
        return (mean);
};

/************************************************************************************/
/*CalculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0)    */
/************************************************************************************/
float BoutonFinder::calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        std::cout<< "In Calculate Local Brightness" << std::endl;
#endif
        Image2DType::RegionType local_region;
        Image2DType::IndexType local_index;
        Image2DType::SizeType local_size;

        local_index[0] = x0-radius;
        local_index[1] = y0-radius;

        local_size[0] = 2*radius+1;
        local_size[1] = 2*radius+1;

        if(local_index[0]<0)
        {       
                local_size[0] = local_size[0] + local_index[0];
                local_index[0]=0;
        }
        if(local_index[1]<0)
        {       
                local_size[1] = local_size[1] + local_index[1];
                local_index[1]=0;
        }
        if(local_index[0] + local_size[0] >= image_plane->GetLargestPossibleRegion().GetSize(0))
        {
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0)-1-local_index[0];
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {
                local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1)-1-local_index[1];
        }

#ifdef DEBUG
        std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
#endif

        local_region.SetIndex(local_index);
        local_region.SetSize(local_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIterator2DType histo_it(image_plane, local_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }

        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();

        float mean = mean_algorithm->GetOutput()->GetElement(0);

        return mean;
};

/************************************************************************************/
/*CalculateLocalStandardDeviation                                                   */
/************************************************************************************/
float BoutonFinder::calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        std::cout<< "In Calculate Local Standard Deviation" << std::endl;
#endif
        Image2DType::RegionType local_region;
        Image2DType::IndexType local_index;
        Image2DType::SizeType local_size;

        local_index[0] = x0-radius;
        local_index[1] = y0-radius;

        local_size[0] = 2*radius+1;
        local_size[1] = 2*radius+1;

        if(local_index[0]<0)
        {       
                local_size[0] = local_size[0] + local_index[0];
                local_index[0]=0;
        }
        if(local_index[1]<0)
        {       
                local_size[1] = local_size[1] + local_index[1];
                local_index[1]=0;
        }
        if(local_index[0] + local_size[0] >= image_plane->GetLargestPossibleRegion().GetSize(0))
        {
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0)-1-local_index[0];
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {
                local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1)-1-local_index[1];
        }

#ifdef DEBUG
        std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
#endif

        local_region.SetIndex(local_index);
        local_region.SetSize(local_size);
        
        HistogramSampleType::Pointer histo_sample = HistogramSampleType::New();
        histo_sample->SetMeasurementVectorSize(1);

        ConstIterator2DType histo_it(image_plane, local_region);

        for(histo_it.GoToBegin(); !histo_it.IsAtEnd(); ++histo_it)
        {
                MeasurementHistogramType greyvalue;
                greyvalue = histo_it.Get();
                histo_sample->PushBack(greyvalue);
        }
        
        MeanAlgorithmType::Pointer mean_algorithm = MeanAlgorithmType::New();

        mean_algorithm->SetInputSample( histo_sample );
        mean_algorithm->Update();
        
        CovarianceAlgorithmType::Pointer covariance_algorithm = CovarianceAlgorithmType::New();

        covariance_algorithm->SetInputSample( histo_sample );
        covariance_algorithm->SetMean( mean_algorithm->GetOutput() );
        covariance_algorithm->Update();

        CovarianceAlgorithmType::OutputType covariance = *(covariance_algorithm->GetOutput());
        float variance = *(covariance.operator[](0));
        float standard_deviation = std::sqrt(variance);

        return standard_deviation;
};


/************************************************************************************/
/*SendRay(x0, y0, ray_legth, nr_of_rays, threshold, n, image_plane) sends a ray     */
/************************************************************************************/
void BoutonFinder::sendRay(double x0, double y0, Image2DType::Pointer image_plane, double *rad_array, double *contour_array, std::list<double*>* minCont, double angle)
{
//std::cout<<"entry "<<angle<<std::endl;
    double min_rad = 0xffffffffffffffff;
    int coutour_index = 0;
    x0 = rint(x0/XYSAMPLING);
    y0 = rint(y0/XYSAMPLING);


    // list of radius 
    std::list<double *> radius;
    std::list<double > skipped_angle;
    
    for(int j = 0; j < NUM_RAYS; j++)
    {
        double phi = j*(PI/NUM_RAYS) ;
        
        //std::cout<<"delta "<<PI/NUM_RAYS<<std::endl;
       // std::cout<<"converted "<<abs(angle)*PI/180<<std::endl;
       // std::cout<<"lower "<<abs(angle)*PI/180 - (PI/NUM_RAYS)<<"Upper "<<abs(angle)*PI/180 + (PI/NUM_RAYS)<<std::endl;
    
        //std::cout<<phi<<'\t'<<y0<<std::endl;
     
        if ( (phi >= ((abs(angle))*PI/180 - 2*(PI/NUM_RAYS) ) &&  (phi <=  ((abs(angle))*PI/180 + 2*(PI/NUM_RAYS)))))
        {
            double  localPhi = phi;
            // this is a direction +/- 10 degree for 18 rays along the axon
            // therefore break
            
            skipped_angle.push_back(localPhi);
            
            //std::cout<<"skipped"<<phi*180/PI<<std::endl;
            
            continue;
            
        }
            
            
        
        
        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());

        Image2DType::IndexType center_index;
        center_index[0] = x0;
        center_index[1] = y0;
        it.SetIndex(center_index);
        center_index = it.GetIndex();
        
        double grey_value = it.Get();//bilinearInterpolation(x_f, y_f, image_plane);//
        
        Image2DType::IndexType corner_index;


        std::list< double * > linear_profile;
        std::list< double * > linear_profile_front;
        std::list< double * > linear_profile_back;
        
        double x_f = (x0);
        double y_f = (y0);
        
            
        double * front_coord = new double[2];
        double * back_coord = new double[2];
        // forward
        for(int k = 0; k < rint(RAY_LEN_PER_DIRECTION/XYSAMPLING); k++)
        {
            double * profile_entry = new double[3];
            double slope = 0;
            double y_new = 0;
            double x_new = 0;
            
            x_f = x_f + 1;
            //y_f = y_f + 1;
            
            y_f = y_f - y0;
            x_f = x_f - x0;
            
            y_new = rint(y_f*cos(phi) - x_f*sin(phi));
            x_new = rint(y_f*sin(phi) + x_f*cos(phi));
            
            y_new = y_new + y0;
            x_new = x_new + x0;
            
            y_f = y_f + y0;
            x_f = x_f + x0;
            
           // std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
                   
            
//             if(phi == PI/2)
//             {
//                 //slope = 1;
//                 y_f = y_f +1;
//                 //std::cout<<x_f<<"    " <<y_f<<std::endl; 
//                 
//             }
//             else
//             {
//                 slope = tan(phi);
//                 x_f = x_f + 1;
//                 y_f = rint(slope* x_f - slope * x0 + y0);
//                 //y_f = rint(y_f);
//             }
            
            
                    //std::cout<<x_f<<"    " <<y_f<<std::endl; 
            if(x_new <= 1 || y_new <= 1 || x_new >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_new >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                        break;
            
            else
            {
                corner_index[0] = (x_new);
                corner_index[1] = (y_new);
                it.SetIndex(corner_index);
                
                profile_entry[0] = (x_new);
                profile_entry[1] = (y_new);
                
                
                grey_value = it.Get();
                profile_entry[2] = grey_value;
                
                //std::cout<<'\t'<<x_f<<'\t'<<y_f<<'\t'<<grey_value<<'\t'<<phi<<std::endl;
                
                linear_profile.push_back(profile_entry);
                linear_profile_front.push_back(profile_entry);
                
            }
            
        }
        
        x_f = (x0);
        y_f = (y0);
        //back
        for(int l = 0; l < rint(RAY_LEN_PER_DIRECTION/XYSAMPLING); l++)
        {
            double * profile_entry = new double[3];
            double slope = 0;
            
            double y_new = 0;
            double x_new = 0;
            
            x_f = x_f -1;
            //y_f = y_f -1;
            
            // conjugate rotation: move the point to centre instead of x0, y0 then rotate and then move it back again
            
            y_f = y_f - y0;
            x_f = x_f - x0;
            
            y_new = rint(y_f*cos(phi) - x_f*sin(phi));
            x_new = rint(y_f*sin(phi) + x_f*cos(phi));
            
            y_new = y_new + y0;
            x_new = x_new + x0;
            
            y_f = y_f + y0;
            x_f = x_f + x0;
            
           // std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
            
//             if(phi == PI/2)
//             {
//                 slope = 1;
//                 y_f = y_f -1;
//                 
//             }
//             else
//             {
//                 slope = tan(phi);
//                 x_f = x_f - 1;
//                 y_f = rint(slope* x_f - slope * x0 + y0);
//                 //y_f = rint(y_f);
//             }
                    
            if(x_new <= 1 || y_new <= 1 || x_new >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_new >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                        break;
            
            else
            {
                corner_index[0] = (x_new);
                corner_index[1] = (y_new);
                it.SetIndex(corner_index);
                
                profile_entry[0] = (x_new);
                profile_entry[1] = (y_new);
                
                
                grey_value = it.Get();
                profile_entry[2] = grey_value;
                
//                         std::cout<<'\t'<<x_f<<'\t'<<y_f<<'\t'<<grey_value<<'\t'<<phi<<std::endl;
                
                linear_profile.push_back(profile_entry);
                linear_profile_back.push_back(profile_entry);    
            }
            
        }
        
        
        
        
        std::list< double * >::iterator profile_it; 
        std::list< double * >::iterator next_it; 
        double max_gray = 0;
        
                // for printing only
        
                for(profile_it = linear_profile_front.begin(); profile_it != linear_profile_front.end(); profile_it++)
                {
                    //std::cout<<'\t'<<(*profile_it)[0]<<'\t'<<(*profile_it)[1]<<'\t'<<(*profile_it)[2]<<std::endl;
                   
                }
                
                for(profile_it = linear_profile_back.begin(); profile_it != linear_profile_back.end(); profile_it++)
                {
                    
                   //std::cout<<'\t'<<(*profile_it)[0]<<'\t'<<(*profile_it)[1]<<'\t'<<(*profile_it)[2]<<std::endl;
                   
                }
                
        
        
        
        // find the max gray value
        for(profile_it = linear_profile.begin(); profile_it != linear_profile.end(); profile_it++)
        {
            //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<max_gray<<'\t'<<std::endl;
            if((*profile_it)[2] > max_gray)
            {
                max_gray = (*profile_it)[2];
            
            }
        }
        
        
        //std::cout<<'\t'<<'\t'<<max_gray<<'\t'<<std::endl;
        // find the half max pixels or the min on either direction
        bool found_threshold = false;
        
        profile_it = linear_profile_front.begin();
        next_it = profile_it;
        next_it++;
        for(; next_it != linear_profile_front.end(); profile_it++,next_it++)
        {
            //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<(*next_it)[2]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
            
            if(((*profile_it)[2] >= max_gray/2) && ((*next_it)[2] <= max_gray/2))
            {
                found_threshold = true;
                front_coord[0] = (*next_it)[0];
                front_coord[1] = (*next_it)[1];
                
                if(coutour_index<(NUM_RAYS*4))
                {
                    contour_array[coutour_index] = front_coord[0];
                    
                    coutour_index++;
                    
                    //if(coutour_index<(NUM_RAYS*4))
                    contour_array[coutour_index] = front_coord[1];
                    coutour_index++;
                    
                    //std::cout<<std::endl;
                   // std::cout<<'\t'<<front_coord[0]<<'\t'<<front_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    //std::cout<<'\t'<<coutour_index<<std::endl;
                    
                }
                
                break;
                //rad_contour.push_back(front_coord);
                
                // Lets find the point between the 2 pixels for accuracy using interpolation
                // y = y2-y1 / x2-x1 * x-x1 + y1.. where y is the pixel number, x is gray value
                
//                         std::cout<<'\t'<<front_coord[0]<<'\t'<<front_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
//                         std::cout<<'\t'<<coutour_index<<std::endl;
            }
            
                
        }
        if(found_threshold == false)
        {
            // may be there is no half max.. may be there is another guy close to it
            // find a minimum of this profile_entry
            double min_gray = 0xffffffffffffffff;
            for(profile_it = linear_profile_front.begin(); profile_it != linear_profile_front.end(); profile_it++)
            {
                if((*profile_it)[2] < min_gray)
                {
                    min_gray = (*profile_it)[2];
                    front_coord[0] = (*profile_it)[0];
                    front_coord[1] = (*profile_it)[1];
                    
                    
                    
                    //rad_contour.push_back(front_coord);
//                             std::cout<<'\t'<<front_coord[0]<<'\t'<<front_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    
                }
                    
            }
            
            if(coutour_index<(NUM_RAYS*4))
            {
                contour_array[coutour_index] = front_coord[0];
                coutour_index++;
                contour_array[coutour_index] = front_coord[1];
                coutour_index++;
                
                //std::cout<<std::endl;
                    //std::cout<<'\t'<<front_coord[0]<<'\t'<<front_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                   // std::cout<<'\t'<<coutour_index<<std::endl;
                    
            }
            
            //else
                //std::cout<<"countour arr out of bound"<<std::endl;
            
            
        }
        
        
        
        found_threshold = false;
        
        profile_it = linear_profile_back.begin();
        next_it = profile_it;
        next_it++;
        for(; next_it != linear_profile_back.end(); profile_it++,next_it++)
        {
            //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<(*next_it)[2]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
            if(((*profile_it)[2] >= max_gray/2) && ((*next_it)[2] <= max_gray/2))
            {
                found_threshold = true;
                // find the point between the two
                
                back_coord[0] = (*next_it)[0];
                back_coord[1] = (*next_it)[1];
                
                if(coutour_index<(NUM_RAYS*4))
                {
                    contour_array[coutour_index] = back_coord[0];
                    coutour_index++;
                    contour_array[coutour_index] = back_coord[1];
                    coutour_index++;
                    
                    //std::cout<<std::endl;
                    //std::cout<<'\t'<<back_coord[0]<<'\t'<<back_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    //std::cout<<'\t'<<coutour_index<<std::endl;
                    
                }
                //else
                    //std::cout<<"countour arr out of bound"<<std::endl;
            
                
                break;
                //rad_contour.push_back(back_coord);
//                         std::cout<<'\t'<<(back_coord)[0]<<'\t'<<(back_coord)[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
//                         std::cout<<'\t'<<coutour_index<<std::endl;
                
            }
        }
        
        
        if(found_threshold == false)
        {
            // may be there is no half max.. may be there is another guy close to it
            // find a minimum of this profile_entry
            double min_gray = 0xffffffffffffffff;
            int cnt = 0;
            for(profile_it = linear_profile_back.begin(); profile_it != linear_profile_back.end(); profile_it++,cnt++)
            {
                
                
            
                if((*profile_it)[2] < min_gray)
                {
                        
                    min_gray = (*profile_it)[2];
                    back_coord[0] = (*profile_it)[0];
                    back_coord[1] = (*profile_it)[1];
                    
                    //rad_contour.push_back(back_coord);
                }
                    
            } 
            
            if(coutour_index<(NUM_RAYS*4))
            {
                contour_array[coutour_index] = back_coord[0];
                
                coutour_index++;
                
                contour_array[coutour_index] = back_coord[1];
                coutour_index++;
                
                 //std::cout<<std::endl;
                   // std::cout<<'\t'<<back_coord[0]<<'\t'<<back_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    //std::cout<<'\t'<<coutour_index<<std::endl;
            }
            //else
                //std::cout<<"countour arr out of bound"<<std::endl;
            
        }
        
        
        // find the radius
        double rad_ray = XYSAMPLING * Utility::euclidean_distance(back_coord,front_coord,2,1);
        //std::cout<<rad_ray<< std::endl;
//                 contour_array[coutour_index] = back_coord[0];
//                 coutour_index++;
//                 contour_array[coutour_index] = back_coord[1];
//                 coutour_index++;
               
        radius.push_back(&rad_ray);
        
        if(rad_ray < min_rad)
        {
            min_rad = rad_ray;
        }
           
    }
    
    rad_array[0] = min_rad;
    
    // Now take the avg of radii to fill in the missing countour along the axon
    std::list<double *>::iterator radit; 
    
    unsigned int cnt = 0;
    double sum = 0;
    for(radit = radius.begin(); radit != radius.end(); radit++)
    {
        cnt++;
        sum = sum + (*radit)[0];
        
    }
    
    double avg_rad = sum/cnt;
    //std::cout<<"avg= "<<avg_rad<< std::endl;
    
    rad_array[1] = avg_rad;
    
     // Now go avg_rad length along the direction of axon and mark those points as contour points
     std::list<double>::iterator skipped_angle_it;
     
     
     
     for(skipped_angle_it = skipped_angle.begin(); skipped_angle_it != skipped_angle.end(); skipped_angle_it++)
     {
         double x_f = x0;
         double y_f = y0;
         double x_new = 0;
         double y_new = 0;
         
        // start with x0,y0 and add the avg length to it.. then rotate it
         x_f = rint(x_f + (avg_rad/2)/XYSAMPLING);
         
         //std::cout<<'\t'<<(*skipped_angle_it)<<std::endl;
        // std::cout<<'\t'<<x_f<<std::endl;
         
         y_f = y_f - y0;
         x_f = x_f - x0;
         
         
         
         y_new = rint(y_f*cos((*skipped_angle_it)) - x_f*sin((*skipped_angle_it)));
         x_new = rint(y_f*sin((*skipped_angle_it)) + x_f*cos((*skipped_angle_it)));
         
         y_new = y_new + y0;
         x_new = x_new + x0;
         
         if(coutour_index<(NUM_RAYS*4))
         {
             contour_array[coutour_index] = x_new;
             
             coutour_index++;
             
             contour_array[coutour_index] = y_new;
             coutour_index++;
             
                 // std::cout<<std::endl;
                 // std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
                // std::cout<<'\t'<<coutour_index<<std::endl;
        }
         else
         {
             //std::cout<<"countour arr out of bound"<<std::endl;
         }
         
         x_f = x0;
         y_f = y0;
         
        x_f =  rint(x_f - (avg_rad/2)/XYSAMPLING);
         
         
         y_f = y_f - y0;
         x_f = x_f - x0;
        
         
         
         y_new = rint(y_f*cos((*skipped_angle_it)) - x_f*sin((*skipped_angle_it)));
         x_new = rint(y_f*sin((*skipped_angle_it)) + x_f*cos((*skipped_angle_it)));
         
         y_new = y_new + y0;
         x_new = x_new + x0;
         
         if(coutour_index<(NUM_RAYS*4))
        {
             contour_array[coutour_index] = x_new;
             
             coutour_index++;
             
             contour_array[coutour_index] = y_new;
             coutour_index++;
             
                 //std::cout<<std::endl;
                 //std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
                 // std::cout<<'\t'<<coutour_index<<std::endl;
         }
         else
        {
            // std::cout<<"countour arr out of bound"<<std::endl;
        }
         
        
         
     }
     
    

    // Min rad contour
    
    for(unsigned int m = 0; m < NUM_RAYS; m++ )
    {
        
        double theta = m*(PI/NUM_RAYS) ;
        
        double x_f = x0;
         double y_f = y0;
         double x_new = 0;
         double y_new = 0;
         
        // start with x0,y0 and add the avg length to it.. then rotate it
         x_f = rint(x_f + (min_rad/2)/XYSAMPLING);
         
        // std::cout<<'\t'<<(theta)<<std::endl;
        // std::cout<<'\t'<<x_f<<std::endl;
         
         y_f = y_f - y0;
         x_f = x_f - x0;
         
         
         
         y_new = rint(y_f*cos((theta)) - x_f*sin((theta)));
         x_new = rint(y_f*sin((theta)) + x_f*cos((theta)));
         
         y_new = y_new + y0;
         x_new = x_new + x0;
         
         double * coord = new double[2];
         
         coord[0] = x_new;
         coord[1] = y_new;
         
         //std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
         minCont->push_back(coord);
         
         x_f = x0;
         y_f = y0;
         
        x_f =  rint(x_f - (min_rad/2)/XYSAMPLING);
         
         
         y_f = y_f - y0;
         x_f = x_f - x0;
        
         
         
         y_new = rint(y_f*cos((theta)) - x_f*sin((theta)));
         x_new = rint(y_f*sin((theta)) + x_f*cos((theta)));
         
         y_new = y_new + y0;
         x_new = x_new + x0;
         
         double * coord1 = new double[2];
         coord1[0] = x_new;
         coord1[1] = y_new;
         
         //std::cout<<'\t'<<x_new<<'\t'<<y_new<<'\t'<<std::endl;
         minCont->push_back(coord1);
         
        
    }
    
    
    
};




#if 0
/************************************************************************************/
/*SendRay(x0, y0, ray_legth, nr_of_rays, threshold, n, image_plane) sends a ray     */
/************************************************************************************/
std::list<double *>  BoutonFinder::sendRay(double x0, double y0, float ray_length_per_direction, double phi, Image2DType::Pointer image_plane, double *rad_array, double *contour_array)
{

    
    {
    
            {
                
                std::list<double *> rad_contour;
                VECTOR * tmp_vect = new VECTOR;

                //float phi = n*(2*PI/nr_of_rays);
                //float x_r = 0, y_r = 0, x_f = x0, y_f = y0;
                
                
        //      std::cout << "Phi: " << n*360/nr_of_rays << std::endl;
                
                Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());

                Image2DType::IndexType center_index;
                center_index[0] = x0;
                center_index[1] = y0;
                it.SetIndex(center_index);
                center_index = it.GetIndex();
                
                double grey_value = it.Get();//bilinearInterpolation(x_f, y_f, image_plane);//
                
        #ifdef DEBUG
                //std::cout << "greyvalue " << grey_value << std::endl;
        #endif
                
                float old_grey_value = 0;
                unsigned int counter = 0;
                
                bool shouldPrint = false;
                
                //Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());
                Image2DType::IndexType corner_index;

        #ifdef DEBUG
                
                //std::cout<<"threshold "<< threshold<< std::endl;
        #endif
                
                // collect the gray values along this ray in both directions
                
                std::list< double * > linear_profile;
                std::list< double * > linear_profile_front;
                std::list< double * > linear_profile_back;
                x0 = rint(x0/XYSAMPLING);
                y0 = rint(y0/XYSAMPLING);
                double x_f = (x0);
                double y_f = (y0);
                
                //std::cout<<x_f<<"    " <<y_f<<std::endl; 
                
                double * front_coord = new double[2];
                double * back_coord = new double[2];
                // forward
                for(int i = 0; i < rint(ray_length_per_direction/XYSAMPLING); i++)
                {
                    double * profile_entry = new double[3];
                    double slope = 0;
                    if(phi == PI/2)
                    {
                        //slope = 1;
                        y_f = y_f +1;
                        //std::cout<<x_f<<"    " <<y_f<<std::endl; 
                        
                    }
                    else
                    {
                        slope = tan(phi);
                        x_f = x_f + 1;
                        y_f = rint(slope* x_f - slope * x0 + y0);
                        //y_f = rint(y_f);
                    }
                    
                    
                            //std::cout<<x_f<<"    " <<y_f<<std::endl; 
                    if(x_f <= 1 || y_f <= 1 || x_f >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_f >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                                break;
                    
                    else
                    {
                        corner_index[0] = (x_f);
                        corner_index[1] = (y_f);
                        it.SetIndex(corner_index);
                        
                        profile_entry[0] = (x_f);
                        profile_entry[1] = (y_f);
                        
                        
                        grey_value = it.Get();
                        profile_entry[2] = grey_value;
                        
                        //std::cout<<'\t'<<x_f<<'\t'<<y_f<<'\t'<<grey_value<<'\t'<<phi<<std::endl;
                        
                        linear_profile.push_back(profile_entry);
                        linear_profile_front.push_back(profile_entry);
                        
                    }
                    
                }
                
                x_f = (x0);
                y_f = (y0);
                //back
                for(int i = 0; i < rint(ray_length_per_direction/XYSAMPLING); i++)
                {
                    double * profile_entry = new double[3];
                    double slope = 0;
                    if(phi == PI/2)
                    {
                        slope = 1;
                        y_f = y_f -1;
                        
                    }
                    else
                    {
                        slope = tan(phi);
                        x_f = x_f - 1;
                        y_f = rint(slope* x_f - slope * x0 + y0);
                        //y_f = rint(y_f);
                    }
                            
                    if(x_f <= 1 || y_f <= 1 || x_f >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_f >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                                break;
                    
                    else
                    {
                        corner_index[0] = (x_f);
                        corner_index[1] = (y_f);
                        it.SetIndex(corner_index);
                        
                        profile_entry[0] = (x_f);
                        profile_entry[1] = (y_f);
                        
                        
                        grey_value = it.Get();
                        profile_entry[2] = grey_value;
                        
                        //std::cout<<'\t'<<x_f<<'\t'<<y_f<<'\t'<<grey_value<<'\t'<<phi<<std::endl;
                        
                        linear_profile.push_back(profile_entry);
                        linear_profile_back.push_back(profile_entry);    
                    }
                    
                }
                
                
                std::list< double * >::iterator profile_it; 
                std::list< double * >::iterator next_it; 
                double max_gray = 0;
                
                for(profile_it = linear_profile.begin(); profile_it != linear_profile.end(); profile_it++)
                {
                    //std::cout<<'\t'<<(*profile_it)[0]<<'\t'<<(*profile_it)[1]<<'\t'<<(*profile_it)[2]<<'\t'<<std::endl;
                }
                
                
                // find the max gray value
                for(profile_it = linear_profile.begin(); profile_it != linear_profile.end(); profile_it++)
                {
                    //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<max_gray<<'\t'<<std::endl;
                    if((*profile_it)[2] > max_gray)
                    {
                        max_gray = (*profile_it)[2];
                    
                    }
                }
                //std::cout<<'\t'<<'\t'<<max_gray<<'\t'<<std::endl;
                // find the half max pixels or the min on either direction
                bool found_threshold = false;
                
                profile_it = linear_profile_front.begin();
                next_it = profile_it;
                next_it++;
                for(; next_it != linear_profile_front.end(); profile_it++,next_it++)
                {
                    //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<(*next_it)[2]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    
                    if(((*profile_it)[2] >= max_gray/2) && ((*next_it)[2] <= max_gray/2))
                    {
                        found_threshold = true;
                        front_coord[0] = (*next_it)[0];
                        front_coord[1] = (*next_it)[1];
                        
                        rad_contour.push_back(front_coord);
                        
                        // Lets find the point between the 2 pixels for accuracy using interpolation
                        // y = y2-y1 / x2-x1 * x-x1 + y1.. where y is the pixel number, x is gray value
                        
                        
                    }
                    
                        
                }
                if(found_threshold == false)
                {
                    // may be there is no half max.. may be there is another guy close to it
                    // find a minimum of this profile_entry
                    double min_gray = 0xffffffffffffffff;
                    for(profile_it = linear_profile_front.begin(); profile_it != linear_profile_front.end(); profile_it++)
                    {
                        if((*profile_it)[2] < min_gray)
                        {
                            min_gray = (*profile_it)[2];
                            front_coord[0] = (*profile_it)[0];
                            front_coord[1] = (*profile_it)[1];
                            
                            rad_contour.push_back(front_coord);
                            
                        }
                            
                    }
                    
                    //std::cout<<'\t'<<front_coord[0]<<'\t'<<front_coord[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    
                }
                
                found_threshold = false;
                
                profile_it = linear_profile_back.begin();
                next_it = profile_it;
                next_it++;
                for(; next_it != linear_profile_back.end(); profile_it++,next_it++)
                {
                    //std::cout<<'\t'<<(*profile_it)[2]<<'\t'<<(*next_it)[2]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    if(((*profile_it)[2] >= max_gray/2) && ((*next_it)[2] <= max_gray/2))
                    {
                        found_threshold = true;
                        // find the point between the two
                        
                        back_coord[0] = (*next_it)[0];
                        back_coord[1] = (*next_it)[1];
                        rad_contour.push_back(back_coord);
                        //std::cout<<'\t'<<(back_coord)[0]<<'\t'<<(back_coord)[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    }
                }
                
                if(found_threshold == false)
                {
                    // may be there is no half max.. may be there is another guy close to it
                    // find a minimum of this profile_entry
                    double min_gray = 0xffffffffffffffff;
                    for(profile_it = linear_profile_back.begin(); profile_it != linear_profile_back.end(); profile_it++)
                    {
                        if((*profile_it)[2] < min_gray)
                        {
                            min_gray = (*profile_it)[2];
                            back_coord[0] = (*profile_it)[0];
                            back_coord[1] = (*profile_it)[1];
                            rad_contour.push_back(back_coord);
                        }
                            
                    }
                    //std::cout<<'\t'<<(back_coord)[0]<<'\t'<<(back_coord)[1]<<'\t'<<max_gray/2<<'\t'<<found_threshold<<std::endl;
                    
                }
                
                
                // find the radius
                double radius = XYSAMPLING * Utility::euclidean_distance(back_coord,front_coord,2,1);
                rad_contour.push_front(&radius);   
                
                
            }
            
    }

        
        //return rad_contour;
        //return tmp_vect;
};
#endif

/************************************************************************************/
/*grey_value BilinearInterpolation(x_f, y_f, image_plane)        see Wikipedia      */
/************************************************************************************/
float BoutonFinder::bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane)
{
#ifdef DEBUG
        std::cout<< "In BilinearInterpolation" << std::endl;
#endif
        unsigned int x1 = (unsigned int) (x-1), y1 = (unsigned int)(y-1);         //lower corners
        unsigned int x2 = (unsigned int)(x+1), y2 = (unsigned int)(y+1);  //upper corners

        unsigned char Q11 = 0, Q21 = 0, Q22 = 0, Q12 = 0;                 // grey values of corners

        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());

        Image2DType::IndexType corner_index;
        corner_index[0] = x1;
        corner_index[1] = y1;
        it.SetIndex(corner_index);
        Q11 = it.Get();

        corner_index[0] = x2;
        corner_index[1] = y1;
        it.SetIndex(corner_index);
        Q21 = it.Get();

        corner_index[0] = x2;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q22 = it.Get();

        corner_index[0] = x1;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q12 = it.Get();

        float P = 0;
        float R1 = 0, R2 = 0;

        R1=((x2-x)/(x2-x1)*Q11)+((x-x1)/(x2-x1)*Q21);
        R2=((x2-x)/(x2-x1)*Q12)+((x-x1)/(x2-x1)*Q22);
        P =(((y2-y)/(y2-y1))*R1)+(((y-y1)/(y2-y1))*R2);
        
#ifdef DEBUG
        std::cout<< "Interpolated grey value" << P << std::endl;
#endif

        return P;
};
/************************************************************************************/
/*AddUpOpposingRaysAndFindShortest(distances);                                      */
/************************************************************************************/
int BoutonFinder::addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance )
{
#ifdef DEBUG
        std::cout<< "AddUpOpposingRaysAndFindShortest" << std::endl;
#endif
        float min_distance = 100000000;
        std::vector<float> diameters;
        int index=0, end = vectors.size(), half = end/2;

        for(int i = 0; i < half; i++)
                diameters.push_back(vectors[i]->magnitude + vectors[half+i]->magnitude);

        end = diameters.size();
        
        for(int j = 0; j < end; j++)
        {
                if(diameters[j]<min_distance)
                {
                        min_distance = diameters[j];
                        index = j;
                }
        }
        
        VECTOR * tmp_vect = new VECTOR;
        tmp_vect->coords[X_COORD] = ( vectors[index]->coords[X_COORD] + vectors[half+index]->coords[X_COORD] ) / 2;
        tmp_vect->coords[Y_COORD] = ( vectors[index]->coords[Y_COORD] + vectors[half+index]->coords[Y_COORD] ) / 2;
        *adjustment_vect = tmp_vect;
        
        *distance = min_distance;
        
        
#ifdef DEBUG
        std::cout<< "min diam: " << min_distance << std::endl;
#endif
        
        return index;

};
/************************************************************************************/
/*SmoothRadii()  reduces drastic jumps in radius between adjacent points            */
/************************************************************************************/
void BoutonFinder::smoothRadii()
{

        std::cout<< "SmoothRadii!! " << std::endl;

                
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();
        double length=0;
        
        // create the edge list of the points within the proximity region
        
        
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
               // unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
                length += currentEdge->physicalLength;
                

                //for every point along edge
                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                {

                        double * current_coords = *edge_it;
                        float diameter = current_coords[SURFACE];
                        long counter = 1;
                        
                        std::list< double * >::iterator next_it = edge_it;              
                        std::list< double * >::iterator prev_it = edge_it;
                        ++next_it;
                        --prev_it;
                        
                        if(edge_it == currentEdge->edgePointCoordinates.begin())
                        {
                                double * next_coords = *next_it;
                                
                                if(isWithinEuclideanRange(next_coords, current_coords, 5))
                                {
                                    if(next_coords[SURFACE])
                                    {
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                    }
                                }
                                
                        }
                        else if(next_it == currentEdge->edgePointCoordinates.end())
                        {
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5))
                                {
                                    if(prev_coords[SURFACE])
                                    {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                    }
                                }
                        }
                        else
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5))
                                {
                                    if(prev_coords[SURFACE])
                                    {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                    }
                                }
                                
                                if(isWithinEuclideanRange(next_coords, current_coords, 5))
                                {
                                    if(next_coords[SURFACE])
                                    {
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                    }
                                }       
                        }
                                
                        diameter = XYSAMPLING *(diameter/counter);
                        current_coords[SURFACE] = diameter;
                        //std::cout<< "Final Diam: " << diameter << std::endl;
        
                }

        }
        total_axon_length = length; //set the physical length
};


/*******************************************************************************************/
/*IsWithinEuclideanRange(a,b,range) returns true if within the specified distance          */
/*******************************************************************************************/

bool BoutonFinder::isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance)
{
  float temp_x = 0, temp_y = 0, temp_z = 0, distance = 0;

  temp_x = (float)XYSAMPLING *(( tmp[X_COORD] - neighbor[X_COORD] ) * ( tmp[X_COORD] - neighbor[X_COORD] ));
  temp_y = (float)XYSAMPLING *(( tmp[Y_COORD] - neighbor[Y_COORD] ) * ( tmp[Y_COORD] - neighbor[Y_COORD] ));
  temp_z = (float)ZSAMPLING  *(( tmp[Z_COORD] - neighbor[Z_COORD] ) * ( tmp[Z_COORD] - neighbor[Z_COORD] ));

  distance = temp_x + temp_y + temp_z;

  if(distance < ldistance * ldistance)
    return true;
  else
    return false;
};


/*******************************************************************************************/
/*getAverageStandardDeviation()  returns the average std_dev of all bouton canidates       */
// /*******************************************************************************************/
float BoutonFinder::getAverageStandardDeviation()
{
        std::cout<< "In getAverageStandardDeviation" << std::endl;
  
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();        
        
        float sum = 0;
        int counter = 0;
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                //unsigned int edgeSize = currentEdge->edgePointCoordinates.size();

                //for every point along edge
                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                {
                        double * coords = *edge_it;
                        
                        if(coords[IS_BOUTON])
                        {
                            sum += coords[LOCAL_SIGMA];
                            counter++;
                        }
        
                }
        }
  
  return sum/counter;
};

/**********************************************************************************************/
/*averageLocalValues()  returns the local average of the desired stat in the area of the point*/
/**********************************************************************************************/
float BoutonFinder::averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection)
{
        
        float sum = 0;
        int counter = 0, i=0;
                
       // unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
        std::list< double * >::iterator edge_it;

        for(edge_it = centerPoint, i=0; edge_it != currentEdge->edgePointCoordinates.end() && i<=numToAverageInEachDirection; edge_it++, i++) 
        {
            double * coords = *edge_it;
            // Take average of only non zero values
            if(coords[property] == 0)
            {
                i--;
            }
            else
            {
                sum += coords[property];
                counter ++;
            }
        }
        
        if(centerPoint != currentEdge->edgePointCoordinates.begin())
        {
          for(edge_it = (--centerPoint), i=0; edge_it != currentEdge->edgePointCoordinates.begin() && i<numToAverageInEachDirection; edge_it--, i++) 
          {
                double * coords = *edge_it;
                
                if(coords[property] == 0)
                {
                    i--;
                }
                else
                {
                    sum += coords[property];
                    counter ++;
                }
          }
        }

  return sum/(counter);
};

/**********************************************************************************************/
/*averageLocalValues()  returns the local average of the desired stat in the area of the point*/
/**********************************************************************************************/

#if 0
void BoutonFinder::writeBoutonLandmarkfile(void)
{
#if 0
    std::vector<double *> * Boutonlist = new std::vector<double *>;
    
    // loop through the edge list and add bouton coordinates to the bouton list
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
    unsigned int numOfEdges = edges->size();
    
    for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
    {                       
        Edge * currentEdge = edges->at(pos);
        std::list< double * >::iterator edge_it;        
        
        
        //for every point along edge
        for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
        {
            double * coords = *edge_it;
            
            if(coords[IS_BOUTON])
            {
                Boutonlist->push_back(coords);
                
            }
        }
    }
#endif
    // write the landmark files
    
    std::string format = outputFilename;
    format += "_Contour.landmarkAscii";
    std::ofstream LandMarkData( format.c_str() );
    
    LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
    LandMarkData << ""                                    << std::endl;
    LandMarkData << ""                                    << std::endl;
    LandMarkData << "define Markers " << NUM_RAYS*2     << std::endl;
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
   
    for(unsigned int i=0; i < (NUM_RAYS) * 4; i=i+2)
    {
        //
        double * con_landmarks = new double[3];
        
        con_landmarks[0] = BoutonParamsArray[14].contour_1_plane_XY[i] * XYSAMPLING + imageTranslation[X_COORD];
        con_landmarks[1] = BoutonParamsArray[14].contour_1_plane_XY[i+1] * XYSAMPLING + imageTranslation[Y_COORD];
        

        
        // transform to get the original point
        double oldCoords[4], newCoords[4];
            for(int ii = 0; ii < 3; ++ii)
            {
                    oldCoords[ii] = con_landmarks[ii];
                    newCoords[ii] = 0;
            }
            oldCoords[3] = 1;
            newCoords[3] = 1;
            for(int ii = 0; ii < 3; ++ii)
                    for(int jj = 0; jj < 4; ++jj)
                            newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
            
            for(int ii = 0; ii < 3; ++ii)
            {
                    con_landmarks[ii] = newCoords[ii];
                    //BoutonParamsArray[NumOfLandmarks].nearestPoint[ii] = newCoords[ii];
            }
            
            //LandMarkData << con_landmarks[0] << " " << con_landmarks[1] << " " << BoutonParamsArray[14].transformed_landmark[2] << std::endl;
        
        

     }
    
    
    LandMarkData.close();


    std::string format1 = outputFilename;
    format1 += "_MinRadContour.landmarkAscii";
    std::ofstream LandMarkData1( format1.c_str() );
    
    LandMarkData1 << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
    LandMarkData1 << ""                                    << std::endl;
    LandMarkData1 << ""                                    << std::endl;
    LandMarkData1 << "define Markers " << NUM_RAYS*2     << std::endl;
    LandMarkData1 << ""                                    << std::endl;
    LandMarkData1 << "Parameters {"                        << std::endl;
    LandMarkData1 << "    NumSets 1,"                      << std::endl;
    LandMarkData1 << "    ContentType \"LandmarkSet\""     << std::endl;
    LandMarkData1  << "}"                                   << std::endl;
    LandMarkData1  << ""                                    << std::endl;
    LandMarkData1  << "Markers { float[3] Coordinates } @1" << std::endl;
    LandMarkData1  << ""                                    << std::endl;
    LandMarkData1  << "# Data section follows"              << std::endl;
    LandMarkData1  << "@1"                                  << std::endl;
   
    std::list<double*>::iterator min_cont_it;
    
    for(min_cont_it = BoutonParamsArray[14].contour_1_plane_XY_min_rad.begin(); min_cont_it != BoutonParamsArray[14].contour_1_plane_XY_min_rad.end(); min_cont_it++)
    {
        //
        double * con_landmarks = new double[3];
        
        con_landmarks[0] = (*min_cont_it)[0] * XYSAMPLING + imageTranslation[X_COORD];
        con_landmarks[1] = (*min_cont_it)[1] * XYSAMPLING + imageTranslation[Y_COORD];
        

        
        // transform to get the original point
        double oldCoords[4], newCoords[4];
            for(int ii = 0; ii < 3; ++ii)
            {
                    oldCoords[ii] = con_landmarks[ii];
                    newCoords[ii] = 0;
            }
            oldCoords[3] = 1;
            newCoords[3] = 1;
            for(int ii = 0; ii < 3; ++ii)
                    for(int jj = 0; jj < 4; ++jj)
                            newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
            
            for(int ii = 0; ii < 3; ++ii)
            {
                    con_landmarks[ii] = newCoords[ii];
                    //BoutonParamsArray[NumOfLandmarks].nearestPoint[ii] = newCoords[ii];
            }
            
            LandMarkData1 << con_landmarks[0] << " " << con_landmarks[1] << " " << BoutonParamsArray[14].landmark[2] << std::endl;
        
        

     }
    
    
    LandMarkData1.close();



    
    
}
#endif
void BoutonFinder::UndoTransformation(void)
{
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    std::vector< Vertex * > * vertices = amira_graph->verticesPointer();
    int numOfEdges = edges->size();
    int numofVertices = vertices->size();
    
    // Undo translation for the vertices
    for(int i=0; i<numofVertices; i++)  //for each edge
    {
        
        Vertex * currentVertex = vertices->at(i);
        
        double * coords = (currentVertex)->coordinates;
        coords[X_COORD] = coords[X_COORD] + imageTranslation[X_COORD];
        coords[Y_COORD] = coords[Y_COORD] + imageTranslation[Y_COORD];
        coords[Z_COORD] = coords[Z_COORD] + imageTranslation[Z_COORD];
        
    }
        
    // Undo translation for the edges
    for(int i=0; i<numOfEdges; i++)  //for each edge
    {
        
        Edge * currentEdge = edges->at(i);
        std::list< double * >::iterator it;
        
        for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
        {
            double * coords = *it;
                
            //if(coords[IS_BOUTON])
            {
                ImageType::IndexType pixelIndex;
                
                // Apply image translation before scaling
                coords[X_COORD] = coords[X_COORD] + imageTranslation[X_COORD];
                coords[Y_COORD] = coords[Y_COORD] + imageTranslation[Y_COORD];
                coords[Z_COORD] = coords[Z_COORD] + imageTranslation[Z_COORD];
                
            }
        }
        
    }       
    
    amira_graph->applyTransformation();
    
}







