/****************************************************************************/
/*                                                                          */
/* File:      proximity_finder.cpp                                          */
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
/* History:   30.03.2014                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



#include "proximity_finder.h"
#include "typedefs.h"




/********************************************************************************************************************
Method: findProximities()

Performs the bulk of the work for the ProximityFinder class.
Cycles through every axon-dendrite pairing and creates local regions called Proximities where the axon 
and dendrite are close to each other. This 'closeness' is determined by the global constant MAX_DIST
*********************************************************************************************************************/

void ProximityFinder::findProximities()
{   
    //std::cout<<"In findProximities()"<<std::endl;
    
    //std::cout<<"axonlist size: "<<axon_list.size()<<" dendritelist size: "<< dendrite_list.size()<<std::endl;
    
      
                
    
    

    if(setRadiusAndLocalBrightness(original_image)  != 0)  
       exit(-1);
    
    markLocalMaximums();
    
    
    
    
    // For each axon-dendrite pairing
//    for(int axonCount = 0; axonCount < axon_list.size(); axonCount++)
    for(int axonCount = 0; axonCount < axon_list.size(); axonCount++)    
    {
        SimpleAxon currentAxon = axon_list.at(axonCount);
        
        
        for(int dendriteCount = 0; dendriteCount < dendrite_list.size(); dendriteCount++)
        {
            SimpleDendrite currentDendrite = dendrite_list.at(dendriteCount);
            std::vector<double *> * point_list = new std::vector<double *>;
            
            // This variable ensures that edge segments whose distances hover around max_dist 
            // are not accidentally grouped into separate proximities when they should be one
            int numPointsAway = 0;
            
            currentAxon.setToBegin();
            while(currentAxon.hasNext())
            {
                
                double * axonPoint = currentAxon.getNext();
                
                bool isMatched = false;
                
                currentDendrite.setToBegin();
                while(currentDendrite.hasNext())
                {
                    double * dendritePoint = currentDendrite.getNext();
//                     //std::cout<<"axonID: "<<currentAxon.edgeID<< "  dendriteID: "<< currentDendrite.edgeID<<std::endl;
                    if(currentAxon.edgeID != currentDendrite.edgeID)
                    {
                        double distance = Utility::euclidean_distance(axonPoint, dendritePoint, 3, 1);
                        if(distance < this->max_dist)
                        {
//                          //std::cout<<"pointDist: "<<distance<<std::endl;
                           
                            
/*
                                //std::cout<<"axon_count:         "<<axonCount<<std::endl;
                                //std::cout<<"axonPoint:          "<<axonPoint[0]<<"      "<<axonPoint[1]<<"      "<<axonPoint[2]<<std::endl;
                                //std::cout<<"dendritePoint:          "<<dendritePoint[0]<<"      "<<dendritePoint[1]<<"      "<<dendritePoint[2]<<std::endl;
                                //std::cout<<"dendrite_count:     "<<dendriteCount<<std::endl;*/
//                                //std::cout<<"edge_counter:     "<<edge_counter<<std::endl;

//                            
//  /*                          SimpleAxon prox_axon = axon_list.at(axonCount);
//                            for(int i = 0; i <= 6; i++)
//                            {
//                                double * prox_point = prox_axon.getNext();
//    */                         
                             
//                             bool isBouton = markBoutons(axonPoint);
                             bool isBouton = true;
// //                            //std::cout<<"isBouton = "<<isBouton<<std::endl;
                             if(isBouton == true)
                             {
                            numPointsAway = 0;
                            isMatched = true;
                            point_list->push_back(axonPoint);
  //                          //std::cout<<"axonPoint:          "<<axonPoint[0]<<"      "<<axonPoint[1]<<"      "<<axonPoint[2]<<std::endl;
                            
                            point_list->push_back(dendritePoint);
//                            //std::cout<<"dendritePoint:          "<<dendritePoint[0]<<"      "<<dendritePoint[1]<<"      "<<dendritePoint[2]<<std::endl;
                            
                            }
                            else;
//                            }
                        }
                    }
                    
                }
                
                if(point_list->size() > 0)                                
                {
                    //if just exited a proximity zone, OR we're at the end of the current axon
                    if((isMatched == false && numPointsAway > this->numberOfPointsAway) || currentAxon.hasNext() == false)            
                    {
//                         //std::cout<<"size: "<<point_list->size()<<std::endl;
//                         //std::cout<<"axonID: "<<currentAxon.edgeID<< "  dendriteID: "<< currentDendrite.edgeID<<std::endl;
                            Proximity * prox = new Proximity(point_list, this->original_image);
                            this->proximity_list.push_back(prox);
                            ////std::cout<< "numPointsAway: "<<numPointsAway<<std::endl;
                            ////std::cout<<"Axon: "<<currentAxon.edgeID<<"Dendrite: "<<currentDendrite.edgeID<<std::endl;
                            
                            
                            delete point_list;
                            point_list = new std::vector<double *>;
                    }
                    
                }
//                //std::cout<< "numPointsAway: "<<numPointsAway<<std::endl;
                numPointsAway++;
            }
            
            
        }
        
    }
    
     //std::cout<<"Total proximity locations found before: "<< proximity_list.size()<<std::endl;
     
     
    if(this->eliminateDoubledProximites)
    removeDoubledProximities();
    
    //std::cout<<"Total proximity locations found: "<< proximity_list.size()<<std::endl;
    
};

bool ProximityFinder::markBoutons(double * axon_point)
{
//        //std::cout<< "In markBoutons" << std::endl;
        
        //std::vector< Edge * > * edges = amira_graph->edgesPointer();
        //int numOfEdges = edges->size();
        bool inBouton = false;
        
        
        int pixel_value = 0;
        float avg_std_dev = getAverageStandardDeviation();
        
        ////std::cout<< "avg_std_dev: " << avg_std_dev << std::endl;
        
        
        
//           SimpleAxon currentAxon = axon_list.at(axon_count);
//           Edge * currentEdge = currentAxon.getEdge();
//           std::list< double * >::iterator it;
//           
//           for(it = currentEdge->edgePointCoordinates.begin();it != currentEdge->edgePointCoordinates.end(); it++) //for every point along edge
//           {
                double * coords = axon_point;
            
                float brightness = coords[LOCAL_BRIGHTNESS];
                float std_dev = coords[LOCAL_SIGMA];
                //float avg_bright = coords[AVERAGE_BRIGHTNESS];
                //float avg_rad = coords[AVERAGE_SURFACE];
                float radius = coords[SURFACE];
                
//              //std::cout << "Point at (x, y): ("<< coords[X_COORD] << ", " << coords[Y_COORD] << ") Bright: "<< coords[LOCAL_BRIGHTNESS] << " Rad: " << coords[SURFACE] << std::endl;
                
                
//              if( (brightness > (avg_bright + 10)) && (radius > avg_rad) && (radius >= 3.5) ) 
//                coords[IS_BOUTON] = 1;
//              else
//                coords[IS_BOUTON] = 0;
                
                
                if( coords[IS_BOUTON] && radius > 0.5){ 
                  coords[IS_BOUTON] = 1;
                inBouton = true;
                }
                else
                  coords[IS_BOUTON] = 0;

//          }
          
          return inBouton;
          
        
  
};


void ProximityFinder::markLocalMaximums()
{
    
//    //std::cout<<"in MarkLocalMaximums"<<std::endl;
    double buffer = MAX_DIST;
//    SimpleAxon currentAxon = axon_list.at(axon_count);
//    Edge * currentEdge = currentAxon.getEdge();
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    unsigned int numOfEdges = edges->size();
//    double * currentPoint = currentAxon.getEdgePointCoordinates();
    std::list< double * >::iterator axon_it;
//    unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
    float bright_offset = 5;
    float radius_offset = 0;
    
    for(long pos = numOfEdges -1; pos >= 0; pos--)      //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
    
    
    
    for(axon_it = currentEdge->edgePointCoordinates.begin() ; axon_it != currentEdge->edgePointCoordinates.end(); axon_it++){
        
                        double * current_coords = *axon_it;
                        float diameter = current_coords[SURFACE];
                        float sum = 0;
                        float average = 0;
                        long counter = 1;
                        bool isRadiusMax = false;
                        
//                         int x_pos, y_pos, z_pos;
//                                         
//                                         x_pos = rint(current_coords[X_COORD]/XYSAMPLING);
//                                         y_pos = rint(current_coords[Y_COORD]/XYSAMPLING);
//                                         z_pos = rint(current_coords[Z_COORD]/ZSAMPLING);
//                                         
//                         //std::cout<<"x: "<<x_pos<<"  y: "<<y_pos<<"  z: "<<z_pos<<std::endl;
                        
                        std::list< double * >::iterator next_it = axon_it;              
                        std::list< double * >::iterator prev_it = axon_it;
                        ++next_it;
                        --prev_it;
                        
                        std::list< double * >::iterator next_next_it = next_it;         
                        std::list< double * >::iterator prev_prev_it = prev_it;
                        ++next_next_it;
                        --prev_prev_it; 
         
//                         double *checktheiterator = &axon_it;
//                         //std::cout<<checktheiterator<<std::endl;
    
                        if(axon_it == currentEdge->edgePointCoordinates.begin())     //if at beginning of edge
                        {
                                double * next_coords = *next_it;
                                
//                                //std::cout<<"at the beginning!"<<std::endl;
                                
                                
                        }
                        else if(next_it == currentEdge->edgePointCoordinates.end())       //if at end of edge
                        {
                                double * prev_coords = *prev_it;
                                
//                                //std::cout<<"at the end!"<<std::endl;
                          
                                
                        }
                        else if(prev_it == currentEdge->edgePointCoordinates.begin() && edgeSize >= 4)  //second bouton from beginning
                        {
                                
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                                double * next_next_coords = *next_next_it;
                                
//                                //std::cout<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
//                                //std::cout<<current_coords[SURFACE]<<std::endl;
                                
                          
                                if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
                                  current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
                                {
                                    
                                    
                                    average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
                                    
//                                    //std::cout<<"average:  "<<average<<std::endl;
//                                    //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                    
                                    if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
                                        current_coords[SURFACE] >= next_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(next_coords[SURFACE] >= next_next_coords[SURFACE] && 
                                        next_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    
//                                  if(isRadiusMax)
                                    if(current_coords[SURFACE] > average+radius_offset)
                                    {
                                        
                                        average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
                                        
//                                        //std::cout<<"Surface:   "<<current_coords[SURFACE]<<std::endl;
                                        
//                                        //std::cout<<"average:  "<<average<<std::endl;
//                                        //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                        
                                        if(current_coords[LOCAL_BRIGHTNESS] > (average + bright_offset)){
                                            
                                          current_coords[IS_BOUTON] = 1;
//                                        //std::cout<<"just marked a bouton!"<<std::endl;
                                        }
                                        else{
                                          current_coords[IS_BOUTON] = 0;
//                                          //std::cout<<"just unmarked a bouton!"<<std::endl;
                                        }
                                    }
                                }
                                else
                                  current_coords[IS_BOUTON] = 0;
//                                //std::cout<<"no Bouton!"<<std::endl;
                                
                        }

                        else if(next_next_it == currentEdge->edgePointCoordinates.end() && edgeSize >= 4)  //second bouton from end
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                                double * prev_prev_coords = *prev_prev_it;
                                
                                
                          
                                if(current_coords[LOCAL_BRIGHTNESS] > prev_coords[LOCAL_BRIGHTNESS] && 
                                  current_coords[LOCAL_BRIGHTNESS] > next_coords[LOCAL_BRIGHTNESS] )
                                {
                                    
                                    average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//                                    //std::cout<<"Brightness:    "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                    
//                                    //std::cout<<"average:  "<<average<<std::endl;
//                                   //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                    
                                    
                                    if(current_coords[SURFACE] >= prev_coords[SURFACE] && 
                                        current_coords[SURFACE] >= next_coords[SURFACE] )
                                      isRadiusMax = true;
                                    else if(prev_coords[SURFACE] >= prev_prev_coords[SURFACE] && 
                                        prev_coords[SURFACE] >= current_coords[SURFACE] )
                                      isRadiusMax = true;
                                    
//                                  if(isRadiusMax)
                                    if(current_coords[SURFACE] > average+radius_offset)
                                    {

//                                        //std::cout<<"Surface:   "<<current_coords[SURFACE]<<std::endl;
                                        average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
                                        
//                                        //std::cout<<"average:  "<<average<<std::endl;
//                                        //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                      
                                        if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset)){
                                          current_coords[IS_BOUTON] = 1;
//                                            //std::cout<<"just marked a bouton!"<<std::endl;
                                        }
                                        else{
                                          current_coords[IS_BOUTON] = 0;
//                                          //std::cout<<"just unmarked a bouton!"<<std::endl;
                                            
                                        }  
                                    }
                                }
                                else{
                                  current_coords[IS_BOUTON] = 0;
//                                  //std::cout<<"just unmarked a bouton!"<<std::endl;
                                   
                                }
                                
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

                                    average = averageLocalValues(currentEdge, axon_it, SURFACE, 3);
//                                    //std::cout<<"Brightness:     "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                    
//                                    //std::cout<<"average:  "<<average<<std::endl;
//                                    //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                    
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

                                          average = averageLocalValues(currentEdge, axon_it, LOCAL_BRIGHTNESS, 1);
//                                          //std::cout<<"Surface:   "<<current_coords[SURFACE]<<std::endl;
//                                          //std::cout<<"average:                          "<<average<<std::endl;
//                                          //std::cout<<"current_coors local brightness:   "<<current_coords[LOCAL_BRIGHTNESS]<<std::endl;
                                          
                                          if(current_coords[LOCAL_BRIGHTNESS] > (average+bright_offset)){
                                            current_coords[IS_BOUTON] = 1;
 //                                           //std::cout<<"just marked a bouton!"<<std::endl;
                                          }
                                          else{
                                            current_coords[IS_BOUTON] = 0;
//                                            //std::cout<<"just unmarked a bouton!"<<std::endl;
                                              
                                        }
                                    }
                                }
                                else
                                  current_coords[IS_BOUTON] = 0;{
//                                  //std::cout<<"just unmarked a bouton!"<<std::endl;
                                  }
                                
                        }
        
                }

        
        }
    
            
    //std::vector< Edge * > * edges = amira_graph->edgesPointer();

};


float ProximityFinder::averageLocalValues(Edge * currentEdge, std::list< double * >::iterator centerPoint, int property, int numToAverageInEachDirection)
{
        
#ifdef DEBUG    
        //std::cout<<"in averageLocalValues"<<std::endl;
#endif        
        
        float sum = 0;
        int counter = 0, i=0;
                
        unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
        std::list< double * >::iterator edge_it;

        for(edge_it = centerPoint, i=0; edge_it != currentEdge->edgePointCoordinates.end() && i<=numToAverageInEachDirection; edge_it++, i++) 
        {
                double * coords = *edge_it;
                sum += coords[property];
                counter ++;
                
//                //std::cout<<"sum  :"<<sum<<std::endl;
        }
        
        if(centerPoint != currentEdge->edgePointCoordinates.begin())
        {
          for(edge_it = (--centerPoint), i=0; edge_it != currentEdge->edgePointCoordinates.begin() && i<numToAverageInEachDirection; edge_it--, i++) 
          {
                  double * coords = *edge_it;
                  sum += coords[property];
                  counter++;
                  
//                  //std::cout<<"sum  :"<<sum<<std::endl;
          }
        }

  return sum/(counter);
};






float ProximityFinder::calculateLocalStandardDeviation(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        //std::cout<< "In Calculate Local Standard Deviation" << std::endl;
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
            if(local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0))
                local_size[0] = 0;
            else
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0)-1-local_index[0];
                local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0);
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {
            if(local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1))
                local_size[1] = 0;
            else
                local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1)-1-local_index[1];
                local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1);
        }

#ifdef DEBUG
        //std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
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


float ProximityFinder::calculateLocalBrightness(Image2DType::Pointer image_plane, float x0, float y0, int radius)
{
#ifdef DEBUG    
        //std::cout<< "In Calculate Local Brightness" << std::endl;
#endif
        Image2DType::RegionType local_region;
        Image2DType::IndexType local_index;
        Image2DType::SizeType local_size;

        local_index[0] = x0-radius;
        local_index[1] = y0-radius;

        local_size[0] = 2*radius+1;
        local_size[1] = 2*radius+1;
        
//        //std::cout<<local_size[0]<<std::endl;
//       //std::cout<<local_size[1]<<std::endl;
        
//        //std::cout<<local_index[0]<<std::endl;
        
        
        

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
            if(local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0))
                local_size[0] = 0;
            else
                local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0) - local_index[0] - 1;
 //               local_index[0] = image_plane->GetLargestPossibleRegion().GetSize(0);
                
//             else    
//             local_size[0] = image_plane->GetLargestPossibleRegion().GetSize(0) - local_index[0] - 1;
//             //local_size[0] = local_index[0] - image_plane->GetLargestPossibleRegion().GetSize(0)-1;
//                      
                            
        }
        if(local_index[1] + local_size[1] >= image_plane->GetLargestPossibleRegion().GetSize(1))        
        {   
            if(local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1))
                local_size[1] = 0;
            else
            local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1) - local_index[1] - 1;
//            local_index[1] = image_plane->GetLargestPossibleRegion().GetSize(1);
//             float diff = local_index[1] + 1;
//             if(image_plane->GetLargestPossibleRegion().GetSize(0) <= diff)
//                 local_size[1] = 0;
//            
//             else
//             local_size[1] = image_plane->GetLargestPossibleRegion().GetSize(1) - local_index[1] - 1;
//             //local_size[1] = local_index[1] - image_plane->GetLargestPossibleRegion().GetSize(1)-1;
        }
        
        
        
        

#ifdef DEBUG
        //std::cout<< "x= " << x0<< " y= " << y0 << " index= " << local_index << " size= " << local_size << std::endl << std::flush;
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



float ProximityFinder::calculateThreshold(ImageType::Pointer image, float x0, float y0, float z0)
{
#ifdef DEBUG    
        //std::cout<< "In Calculate Threshold" << std::endl;
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

        if(threshold_index[0]<0)
        {       
                threshold_size[0] = threshold_size[0] + threshold_index[0];
                threshold_index[0]=0;
        }
        if(threshold_index[1]<0)
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
        //std::cout<< "x= " << x0<< " y= " << y0 << " index= " << threshold_index << " size= " << threshold_size << std::endl << std::flush;
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
        float variance = *(covariance.operator[](0));
        float standard_deviation = std::sqrt(variance);

#ifdef DEBUG    
        //std::cout<< "In Calculate Threshold" << mean << std::endl << std::flush;
#endif
        return (mean );
};

float ProximityFinder::getAverageStandardDeviation()
{
//        //std::cout<< "In getAverageStandardDeviation" << std::endl;
  
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();        
        
        float sum = 0;
        int counter = 0;
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();

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

int ProximityFinder::setRadiusAndLocalBrightness(ImageType::Pointer image)
{

        //std::cout<< "In SetRadiusAndLocalBrightness " << std::endl;
        
//        SimpleAxon currentAxon = axon_list.at(axon_count);
//        Edge * currentEdge = currentAxon.getEdge();

        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
//        //std::cout<<"got edges"<<std::endl;
        


        
        unsigned int numOfEdges = edges->size();
        unsigned int end_x = image->GetLargestPossibleRegion().GetSize(0);
//        //std::cout<<end_x<<std::endl;
        unsigned int end_y = image->GetLargestPossibleRegion().GetSize(1);
//        //std::cout<<end_y<<std::endl;
        unsigned int end_z = image->GetLargestPossibleRegion().GetSize(2);
//        //std::cout<<end_z<<std::endl;
        
//        image->Print(//std::cout);
        
        int x_pos, y_pos, z_pos;

        Image2DType::Pointer image_plane_5;
        Image2DType::Pointer image_plane_3;
        
        int count = 1;

        for(int z = 0; z < end_z; z++)
        {
                image_plane_5 = getImagePlane(z, 5, image);
                image_plane_3 = getImagePlane(z,3,image);

//                //std::cout<<"got images"<<std::endl;
//                image_plane_3->Print(//std::cout);
//                image_plane_5->Print(//std::cout);
        
                
                if(numOfEdges == 0)
                {
                        //std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                        return -1;
                }
                else
                {
// #pragma omp parallel for schedule(dynamic,1)
                        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                                
                        {      
                                Edge * currentEdge = edges->at(pos);
//                                //std::cout<<"got currentEdge"<<std::endl;
                                std::list< double * >::iterator edge_it;        
                                std::list< double * >::iterator next_it;
                                
                                int counter = 0;


                                //for every point along edge
                                for(edge_it = currentEdge->edgePointCoordinates.begin();edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                                {
                                        double * coords = *edge_it;
//                                        counter ++;
//                                        //std::cout<<"edge_it:     "<<coords[0]<<std::endl;
                                        
                                        x_pos = rint(coords[X_COORD]/XYSAMPLING);
                                        y_pos = rint(coords[Y_COORD]/XYSAMPLING);
                                        z_pos = rint(coords[Z_COORD]/ZSAMPLING);
                                        
//                                        //std::cout<<"x: "<<x_pos<<"  y: "<<y_pos<<"  z: "<<z_pos<<std::endl;
                                        
                                        if(z_pos == z)
                                        {
                                                /*****************************************************************/
                                                /*set radius info                                                */
                                                /*****************************************************************/
                                          
                                                float x0 = x_pos;
                                                float y0 = y_pos;
                                                float threshold = 0;

//                                              threshold = CalculateThreshold(image_plane_3, x0, y0);
//                                                //std::cout<<"try to set threshold"<<std::endl;
                                                threshold = coords[THRESHOLD];
//                                                //std::cout<<"threshold"<<coords[THRESHOLD]<<std::endl;
                                                
                                                unsigned int nr_of_rays = 20;
                                                float ray_length = 1.5;
                                                int index = 0;

                                                std::vector<VECTOR *> vectors;
                                                VECTOR * adjustment_vect;
                                                float distance = 0;

                                                x0 = x0 + 0.5;
                                                y0 = y0 + 0.5;
        
                                                for(unsigned n=1; n<=nr_of_rays; n++)
                                                {
                                                        VECTOR * tmp;
//                                                        //std::cout<<"try to sendRay"<<std::endl;
                                                        tmp = sendRay(x0, y0, ray_length, nr_of_rays, threshold, n, image_plane_3);
//                                                        //std::cout<<"sent Ray"<<std::endl;
                                                        vectors.push_back(tmp);
                                                        ////std::cout<< "Distance: " << tmp_distance << std::endl;
                                                }
                                                
//                                                //std::cout<<"try to addUpOpposingRaysAndFindShortest"<<std::endl;
                                                index = addUpOpposingRaysAndFindShortest(vectors, &adjustment_vect, &distance);
//                                                //std::cout<<"addedUpOpposingRaysAndFoundShortest"<<std::endl;
//                                              //std::cout<< "Distance: " << distance <<std::endl;
                                                //index = 1;

                                                coords[SURFACE] = distance;
//                                                //std::cout<<"got coords surface"<<std::endl;
                                                
                                                /*****************************************************************/
                                                /*adjust midline based on radius info                            */
                                                /*****************************************************************/
                                                next_it = edge_it;
                                                next_it++;
                                                
//                                                 if(edge_it != currentEdge->edgePointCoordinates.begin() && next_it != currentEdge->edgePointCoordinates.end())
//                                                 {
// //                                                  //std::cout<<"old       "<<coords[X_COORD]<<"         "<<coords[Y_COORD]<<std::endl;  
//                                                   float new_x = coords[X_COORD] + adjustment_vect->coords[X_COORD] * XYSAMPLING;
//                                                   float new_y = coords[Y_COORD] + adjustment_vect->coords[Y_COORD] * XYSAMPLING;
//                                                   
//                                                   if((new_x/XYSAMPLING) < end_x-1 && (new_y/XYSAMPLING) < end_y-1)
//                                                   {
//                                                           coords[X_COORD] = new_x;
//                                                           coords[Y_COORD] = new_y;
//                                                           //std::cout<<"new       "<<coords[X_COORD]<<"         "<<coords[Y_COORD]<<std::endl;  
//                                                   }
//                                                 }
                                                
                                                /*****************************************************************/
                                                /*set brightness info                                            */
                                                /*****************************************************************/
                                                x0 = rint(coords[X_COORD]/XYSAMPLING);
                                                y0 = rint(coords[Y_COORD]/XYSAMPLING);
                                                
//                                                //std::cout<<x0<<"    "<<y0<<std::endl;
                                        
                                                coords[LOCAL_BRIGHTNESS] = calculateLocalBrightness(image_plane_5, x0, y0, 1);
//                                                //std::cout<<"got coords local bright"<<std::endl;
                                                coords[LOCAL_SIGMA] = calculateLocalStandardDeviation(image_plane_5, x0, y0, 2);
//                                                //std::cout<<"got coords local stdv"<<std::endl;
//                                                    coords[LOCAL_SIGMA] = 1;
                                        }
                                        
                                        
                                }
                        }
                }
        }       

        smoothRadii();
        return 0;
};

VECTOR * ProximityFinder::sendRay(float x0, float y0, float ray_length, unsigned int nr_of_rays, float threshold, unsigned int n, Image2DType::Pointer image_plane)
{
#ifdef DEGUB
        //std::cout<< "In SendRay" << std::endl;
#endif
        VECTOR * tmp_vect = new VECTOR;

        float phi = n*(2*PI/nr_of_rays);
        float x_r = 0, y_r = 0, x_f = x0, y_f = y0;
        
//      //std::cout << "Phi: " << n*360/nr_of_rays << std::endl;

        float grey_value = bilinearInterpolation(x_f, y_f, image_plane);
        float old_grey_value = 0;
        unsigned int counter = 0;
        
        bool shouldPrint = false;
        
        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());
        Image2DType::IndexType corner_index;


        while(grey_value >= threshold)
        {
                  
                x_f = x_f + (ray_length * cos(phi));
                y_f = y_f + (ray_length * sin(phi));
                if(x_f <= 1 || y_f <= 1 || x_f >= image_plane->GetLargestPossibleRegion().GetSize(0)-1 || y_f >= image_plane->GetLargestPossibleRegion().GetSize(1)-1)
                        break;
                else
                {
                        corner_index[0] = x_f;
                        corner_index[1] = y_f;
                        it.SetIndex(corner_index);
                        
                        old_grey_value = grey_value;
                        grey_value = it.Get();
                  
                  
                        //grey_value = BilinearInterpolation(x_f, y_f, image_plane);
                        counter++;
                }
                  
        }

        x_f = x_f - x0;
        y_f = y_f - y0;
        
        tmp_vect->coords[X_COORD] = x_f;
        tmp_vect->coords[Y_COORD] = y_f;
        tmp_vect->magnitude = sqrt((x_f*x_f)+(y_f*y_f));
//      tmp_vect->magnitude = ray_length * counter;
        
        return tmp_vect;
};


float ProximityFinder::bilinearInterpolation(float x, float y, Image2DType::Pointer image_plane)
{
#ifdef DEBUG
        //std::cout<< "In BilinearInterpolation" << std::endl;
#endif
        /*unsigned*/ int x1 = (/*unsigned*/ int) (x-1), y1 = (/*unsigned*/ int)(y-1);         //lower corners
        /*unsigned */int x2 = (/*unsigned*/ int)(x+1), y2 = (/*unsigned*/ int)(y+1);  //upper corners

        /*unsigned*/ char Q11 = 0, Q21 = 0, Q22 = 0, Q12 = 0;// grey values of corners  
        
//          //std::cout<<x<<"   "<<y<<std::endl;
//          //std::cout<<x1<<"   "<<y1<<std::endl;
//          //std::cout<<x2<<"   "<<y2<<std::endl;
        
        
//        image_plane->Print(//std::cout);

//        //std::cout<<"GetLargestPossibleRegion"<<std::endl;
        
        int y_size, x_size;
        
        x_size = image_plane->GetLargestPossibleRegion().GetSize(0);
        y_size = image_plane->GetLargestPossibleRegion().GetSize(1);
        
        if(x < x_size && x >= 0 && y < y_size && y >=0){
        
        
        Iterator2DType it(image_plane, image_plane->GetLargestPossibleRegion());
        
//        //std::cout<<"Got the largest possible region"<<std::endl;

        Image2DType::IndexType corner_index;
        corner_index[0] = x1;
//        //std::cout<<corner_index[0]<<std::endl;
        corner_index[1] = y1;
//        //std::cout<<corner_index[1]<<std::endl;
        it.SetIndex(corner_index);
//        //std::cout<<"set index"<<std::endl;
        Q11 = it.Get();
//        //std::cout<<"char gets index"<<std::endl;
//       //std::cout<<" corner 1 1"<<std::endl;
        ;

        corner_index[0] = x2;
        corner_index[1] = y1;
        it.SetIndex(corner_index);
        Q21 = it.Get();
        
//        //std::cout<<" corner 2 1"<<std::endl;

        corner_index[0] = x2;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q22 = it.Get();
        
//        //std::cout<<" corner 2 2"<<std::endl;

        corner_index[0] = x1;
        corner_index[1] = y2;
        it.SetIndex(corner_index);
        Q12 = it.Get();
        
//        //std::cout<<" corner 1 2"<<std::endl;

        float P = 0;
        float R1 = 0, R2 = 0;

        R1=((x2-x)/(x2-x1)*Q11)+((x-x1)/(x2-x1)*Q21);
        R2=((x2-x)/(x2-x1)*Q12)+((x-x1)/(x2-x1)*Q22);
        P =(((y2-y)/(y2-y1))*R1)+(((y-y1)/(y2-y1))*R2);

        return P;
        
        }
    else return 0;    
};

Image2DType::Pointer ProximityFinder::getImagePlane(int z, int depth, ImageType::Pointer input_image)
{
#ifdef DEBUG    
        //std::cout<< "In GetImagePlane" << std::endl;
#endif
  
  
        Image2DType::Pointer image_plane = Image2DType::New();

        Image2DType::IndexType target_index;
        Image2DType::SizeType target_size;

        target_index[0] = 0;
        target_index[1] = 0;

        target_size[0] = input_image->GetLargestPossibleRegion().GetSize(0);
        target_size[1] = input_image->GetLargestPossibleRegion().GetSize(1);
        
        #ifdef DEBUG    
        //std::cout<< target_size << std::endl;
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

        return image_plane;
};


int ProximityFinder::addUpOpposingRaysAndFindShortest(std::vector<VECTOR *> vectors, VECTOR ** adjustment_vect, float* distance )
{
#ifdef DEBUG
        //std::cout<< "AddUpOpposingRaysAndFindShortest" << std::endl;
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
        //std::cout<< "min diam: " << min_distance << std::endl;
#endif
        
        return index;

};

float ProximityFinder::setAverageThreshold(int numOfPointsToAverage)
{

        //std::cout<< "In Set Average threshold " << std::endl;
        
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        
        unsigned int numOfEdges = edges->size();                
        

        if(numOfEdges == 0)
        {
                //std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
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
                                
//                              //std::cout<< "Position "<< primaryPoint[X_COORD] << "  " << primaryPoint[Y_COORD]<< "  "<< primaryPoint[Z_COORD] << std::endl;
//                              //std::cout<< "Position "<< x_pos << "  " << y_pos<< "  "<< z_pos << std::endl;
                                

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


// ImageType::Pointer ProximityFinder::getOriginalImage(char* file_name, int start_index, int end_index)
// {
// 
//         char input_file[1024];
//         strcpy(input_file, file_name);
// 
// 
//         NameGeneratorType::Pointer name_gen = NameGeneratorType::New();         //DECLARE AND INITIALIZE NAME_GENERATOR
//         name_gen->SetSeriesFormat( input_file );
//         name_gen->SetStartIndex( start_index );
//         name_gen->SetEndIndex( end_index );
//         name_gen->SetIncrementIndex( 1 );
//         
//                         
//         SeriesReaderType::Pointer input_reader = SeriesReaderType::New();       //DECLARE AND INITIALIZE INPUT_READER
//         input_reader->SetImageIO( itk::PNGImageIO::New() );
//         input_reader->SetFileNames( name_gen->GetFileNames() );
// 
//         try
//         {
//                 input_reader->Update();
//         }
//         catch( itk::ExceptionObject & err )
//         {
//                 std::cerr << "ImageReaderExceptionObject caught !" << std::endl;
//                 std::cerr << err << std::endl;
//         }
// 
// 
//         ImageType::Pointer original_image;
// 
//         original_image = input_reader->GetOutput();
//         original_image->Update();
//         
//         return original_image;
// };


void ProximityFinder::smoothRadii()
{

//        //std::cout<< "SmoothRadii!! " << std::endl;

                
        std::vector< Edge * > * edges = amira_graph->edgesPointer();
        unsigned int numOfEdges = edges->size();
        double length=0;
        
        
// #pragma omp parallel for schedule(dynamic,1)
        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
        {                       
                Edge * currentEdge = edges->at(pos);
                std::list< double * >::iterator edge_it;        
                unsigned int edgeSize = currentEdge->edgePointCoordinates.size();
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
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                }
                                
                        }
                        else if(next_it == currentEdge->edgePointCoordinates.end())
                        {
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5))
                                {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                }
                        }
                        else
                        {
                                double * next_coords = *next_it;
                                double * prev_coords = *prev_it;
                          
                                if(isWithinEuclideanRange(prev_coords, current_coords, 5))
                                {
                                        diameter = diameter + prev_coords[SURFACE];
                                        counter++;
                                }
                                
                                if(isWithinEuclideanRange(next_coords, current_coords, 5))
                                {
                                        diameter = diameter + next_coords[SURFACE];
                                        counter++;
                                }       
                        }
                                
                        diameter = 0.6666667*(diameter/counter);
                        current_coords[SURFACE] = diameter;
//                      //std::cout<< "Final Diam: " << diameter << std::endl;
        
                }

        }
//        total_axon_length = length; //set the physical length
};

bool ProximityFinder::isWithinEuclideanRange(double* tmp, double* neighbor, unsigned int ldistance)
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


/********************************************************************************************************************
Method: writeProximityLandmarks()

Gets the centers of mass for each proximity and then prints them into a landmark file
*********************************************************************************************************************/
void ProximityFinder::writeProximityLandmarks(char * outputFilename)
{
    std::vector<double *> * center_list = new std::vector<double *>;
    
    
    
    for(int i=0; i<this->proximity_list.size(); i++)
    {
        Proximity * prox = proximity_list.at(i);
        center_list->push_back(prox->getCenterCoords());
    }
    
    if(center_list->size() > 0)
        center_list = transformCoordinates(center_list);
        writeLandmarkFile(center_list, outputFilename);
        writeListFile(center_list, outputFilename);
};


/*******************************************************************************/
/*                                                                             */
/*writeLandmarkFile()                                                          */  
/*      writes a list into a valid landmarkAscii file                          */
/*                                                                             */
/*******************************************************************************/
void ProximityFinder::writeLandmarkFile(std::vector<double*> * list, char* outputFilename)
{
  //std::cout<< "in writeLandmarkFile" << std::endl;
  
// transformToWorldCoordinates(list);
  
  std::string format = outputFilename;
  format += "_locations.landmarkAscii";
  std::ofstream LandMarkData( format.c_str() );
  
  LandMarkData << "# AmiraMesh 3D ASCII 2.0"            << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << ""                                    << std::endl;
  LandMarkData << "define Markers " << list->size()     << std::endl;
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
  
  for(int i=0; i < list->size(); i++)
  {
    LandMarkData << list->at(i)[X_COORD] << " " << list->at(i)[Y_COORD] << " " << list->at(i)[Z_COORD] << std::endl;
  } 
  
  LandMarkData.close();
};

/********************************************************************************/
/*writeListFile()                                                               */
/*      writes a list into a text file                                          */
/********************************************************************************/


void ProximityFinder::writeListFile(std::vector<double*> * list, char* outputFilename)
{
//    transformToWorldCoordinates(list);
    std::string format = outputFilename;
    format += "_locations_list.txt";
    std::ofstream ListData( format.c_str() );
    ImageType::SizeType region_size;
    
    
    
    
    
    ListData << "Proximity location's coordinates"<< std::endl;
    ListData << "Amira coordinates               center coordinates"<< std::endl;
    
    for(int i=0; i < list->size(); i++)
    {
        Proximity * prox = this->proximity_list.at(i);
        
        
        ImageType::RegionType output_region = prox->getOutputRegion();
        region_size = output_region.GetSize();
   
        
        ListData << i << "    " << list->at(i)[X_COORD] - region_size[0]/2*XYSAMPLING << " " << list->at(i)[Y_COORD] - region_size[1]/2*XYSAMPLING << " " << list->at(i)[Z_COORD] - region_size[2]/2*ZSAMPLING << "            "<< list->at(i)[X_COORD] << " " << list->at(i)[Y_COORD] << " " << list->at(i)[Z_COORD]<<std::endl;
    }
    
    ListData.close();
};

void ProximityFinder::removeDoubledProximities(){
    
    //std::cout<<"in removeDoubledProximities"<<std::endl;
    
    
    std::vector<Proximity *> cleaned_prox_list;
    bool merged_proximity = false;
    

    
    for(int i = 0; i < this->proximity_list.size(); i++)
        {
            std::vector<double *> point_list;
            Proximity * prox = this->proximity_list.at(i);         
         
            double * center_coords = prox->getCenterCoords();
/*            //std::cout<<"center_coords:           "<< center_coords<<std::endl;
            //std::cout<<"center_coords:           "<< center_coords[0] << " " << center_coords[1] << " " << center_coords[2] <<std::endl;
  */          
            for(int j = i; j < this->proximity_list.size(); j++)
                {                                     
                    if(j != i)   
                    { 
                        
                        Proximity * other_prox = this->proximity_list.at(j);

                        double * other_center_coords = other_prox->getCenterCoords();
/*                        //std::cout<<"Prev_center_coords:           "<< prev_center_coords<<std::endl;
                        //std::cout<<"Prev_center_coords:           "<< prev_center_coords[0] << " " << prev_center_coords[1] << " " << prev_center_coords[2] <<std::endl;
 */                       
                        double distance_between_points = Utility::euclidean_distance(center_coords, other_center_coords, 3, 1);
                        
                            
                        
//                         //std::cout<<"distance    :  "<<distance_between_points<<std::endl;
/*                        ImageType::Pointer image = prox.getProximityImage();
                        image->Print(//std::cout);*/                       
                        
                        if(distance_between_points < 5 )
                            {                                
                            point_list.push_back(center_coords);
                            point_list.push_back(other_center_coords);
                                                     
                            }
                        else;
                            
                                            
                    }
                    
                    
                
                
                        
                }
                
                if(point_list.size()){
                    Proximity * mergedProximity = new Proximity(&point_list, this->original_image);
                    cleaned_prox_list.push_back(mergedProximity);
                    merged_proximity = true;
                    }
                    
                    if(merged_proximity == false){
                    cleaned_prox_list.push_back(prox); 
                    }
                
            }     

          proximity_list.clear();

          for(int i = 0; i < cleaned_prox_list.size(); i++){
          this->proximity_list.push_back(cleaned_prox_list.at(i));              
          }
}

bool ProximityFinder::readAmiraTransformations()
{
    
       
//     //std::cout<<transformFilename<<std::endl;
//     //std::cout<<filename<<std::endl;
    
//std::cout<<"in readAmiraTransformations"<<std::endl;
        std::ifstream inputStream(inputfilename);
     
        if(!inputStream.fail())
        {
                const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
                const char * numbers = "0123456789";
                const char * signs = "+-";
                const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
                const char * whitespace = "\t ";
                
                std::string currentLine;
                unsigned int line = 0;
                
               
                
                bool parameters = 1;
                bool transform = 0;
                bool correctSection = 1;
                bool correctPrevSection = 0;
                int sectionID = 0;
                unsigned int brackets = 0, transformBrackets = 0;
                unsigned int currentIndex = 0;
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
                        
                        if(currentLine.size())
                                if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
                                        {
//                                              //std::cout << "found correct section transform parameters!" << std::endl;
                                                unsigned int count = 0;
                                                std::string::size_type loc1, loc2, loc3;
                                                loc1 = currentLine.find_first_of(numbers, 0);
                                                loc2 = currentLine.find_first_of(signs, 0);
                                                if(loc2 != std::string::npos)
                                                        if(loc2 < loc1)
                                                                loc1 = loc2;
                                                loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                                while(loc2 != std::string::npos && count < 16)
                                                {
                                                        char * tmp1 = new char[loc2 - loc1];
                                                        currentLine.copy(tmp1, loc2 - loc1, loc1);
                                                        double ftmp1 = atof(tmp1);
                                                        transformation[count%4][count/4]= ftmp1;        // amira files are columns after each other
                                                        loc3 = loc2;
                                                        loc1 = currentLine.find_first_of(numbers, loc3);
                                                        loc2 = currentLine.find_first_of(signs, loc3);
                                                        if(loc2 != std::string::npos)
                                                                if(loc2 < loc1)
                                                                        loc1 = loc2;
                                                        loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                                                        ++count;
                                                        delete [] tmp1;
                                                }
                                             //std::cout << "transformation matrix:" << std::endl;
                                             for(int ii = 0; ii < 4; ++ii)
                                             {
                                                     //std::cout << "[";
                                                     for(int jj = 0; jj < 4; ++jj)
                                                     {
                                                             if(jj < 3)
                                                                     //std::cout << transformation[ii][jj] << ",\t";
                                                             else
                                                                     //std::cout << transformation[ii][jj];
                                                     }
                                                     //std::cout << "]" << std::endl;
                                             }
                                                //remove numeric artifacts from z-axis:
                                                for(int ii = 0; ii < 2; ++ii)
                                                {
                                                        transformation[2][ii] = 0;
                                                        transformation[ii][2] = 0;
                                                }
                                                transformation[2][2] = 1;
                                        }
                       
                }
        }
        
        inputStream.close();
        return 0;
}

std::vector<double* > * ProximityFinder::transformCoordinates(std::vector<double *> * list)
{
    //std::cout<<"in transformCoordinates"<<std::endl;    
    
    std::vector<double *> * new_list = list;
    
    
    
    double ** transform_matrix = this->transformation;

        
        for(int i =0;  i < list->size(); i++) 
                {

                                double tmpPoint[4];
                                double transPoint[4];
                                

                                        tmpPoint[0] = list->at(i)[X_COORD];
                                        tmpPoint[1] = list->at(i)[Y_COORD];
                                        tmpPoint[2] = list->at(i)[Z_COORD];
                                        tmpPoint[3] = 1;
                                        for(int jj = 0; jj < 4; ++jj)
                                        {
                                                transPoint[jj] = 0;
                                                for(int kk = 0; kk < 4; ++kk){
                                                        transPoint[jj] += transform_matrix[jj][kk]*tmpPoint[kk];
//                                                //std::cout<<"Transform      "<<worldTransform[jj][kk]<<std::endl;
                                                }
                                        }
                                        
                                        
                                        new_list->at(i)[X_COORD] = transPoint[0];
                                        new_list->at(i)[Y_COORD] = transPoint[1];
                                        new_list->at(i)[Z_COORD] = transPoint[2];
//                                }
//                        }
//                }     
        }
        
        return new_list;
};
