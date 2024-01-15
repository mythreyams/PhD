

#include "bouton_finder.h"
#include "bouton_params.h"
#include "bouton_cluster.h"
#include "math.h"


int main( int argc , char * argv[])
{
  //#pragma omp parallel
   
  Reader * mergedAmiraReader;
  
  AmiraSpatialGraph * amira_graph;
  AmiraSpatialGraph * resampled_graph;
 
  mergedAmiraReader = new Reader(argv[1],argv[2]);  
  mergedAmiraReader->readSpatialGraphFile(true);
  amira_graph = mergedAmiraReader->getSpatialGraph();

   
  std::vector< Edge * >  destedges;
    //std::list<double *> tempEdgePoints;
    //////std::cout<<"before"<<std::endl;
    // Read the whole cell graph
//     Reader * amiraInputGraphReader = new Reader(input_graph, "/home/mythreya/ResampledSpatialgraph");//"amira_file_withRad_centerline");  
//     amiraInputGraphReader->readSpatialGraphFile(false);
    
     resampled_graph = new AmiraSpatialGraph();
    
    //////std::cout<<"reading done"<<std::endl;
    //Original_SpatialGraph = amiraInputGraphReader->getSpatialGraph();
    
    unsigned int numOfEdges = amira_graph->getNumberOfEdges();
    std::vector< Edge * > * edges = amira_graph->edgesPointer();
    std::vector< Vertex * > * vertices = amira_graph->verticesPointer();
   
    unsigned int numofVertices = vertices->size();
    
    
    for(int i=0; i<numOfEdges; i++)  //for each edge
    {
        bool need_to_exit = 0;
        //////std::cout<<"inside edge loop"<<std::endl;
        
        Edge * currentEdge = edges->at(i);
        std::list< double * >::iterator curr_pt_it;
        std::list< double * >::iterator next_pt_it;
        
        double prev_offset = 0;
        std::list< double * > tempEdgePoints;
        curr_pt_it = currentEdge->edgePointCoordinates.begin();
        next_pt_it = curr_pt_it;
        
        if(currentEdge->edgePointCoordinates.size() == 0)
        {
          // nothing
        }
        else if(currentEdge->edgePointCoordinates.size() == 1)
        {
          
          // Just copy this point as it is 
          tempEdgePoints.push_back(*curr_pt_it);
        }
        else
        {
          //////std::cout<<"size"<<currentEdge->edgePointCoordinates.size()<<std::endl;
            
          next_pt_it++;
          
          // Copy the first point
          tempEdgePoints.push_back(*curr_pt_it);
          
          for( ;next_pt_it != currentEdge->edgePointCoordinates.end(); ) //for every point along edge
          {
              //////std::cout<<"inside point loop"<<std::endl;
              //////std::cout<<"for edge "<<i<<std::endl;
              
              //////std::cout<<(*curr_pt_it)[0]<<" "<<(*curr_pt_it)[1]<<" "<<(*curr_pt_it)[2]<<std::endl;
              //////std::cout<<(*next_pt_it)[0]<<" "<<(*next_pt_it)[1]<<" "<<(*next_pt_it)[2]<<std::endl;
  //             
              
              // if the distance between the pts is less than the delta then ignore those pts
              double dist = Utility::euclidean_distance(*next_pt_it, *curr_pt_it, 3,1);
              
              while (((dist /*+ prev_offset*/)  < RESAMPLING_DISTANCE))
              {
                //////std::cout<<"inside while"<<std::endl;
                //////std::cout<<(*next_pt_it)[0]<<" "<<(*next_pt_it)[1]<<" "<<(*next_pt_it)[2]<<std::endl;
             ////////std::cout<<(*currentEdge->edgePointCoordinates.end())[0]<<" "<<(*currentEdge->edgePointCoordinates.end())[1]<<" "<<(*currentEdge->edgePointCoordinates.end())[2]<<std::endl; 
             
                ++next_pt_it;
                if(next_pt_it != currentEdge->edgePointCoordinates.end())
                {
                    dist = Utility::euclidean_distance(*next_pt_it, *curr_pt_it, 3,1);
                    
                }
                else
                {
                  // we have reached the last edge already.. so write it and exit
                  --next_pt_it;
                  tempEdgePoints.push_back(*(next_pt_it));
                  next_pt_it++;
                  need_to_exit = 1;
                  break;
                }
                
              }
              
              if(!need_to_exit)
              {
                  // Noiw that we are sure that the dist is greater than sampling distance we can sample new points
                  
                  // Scheme is :
                  // Take a unit vector along the edge and fidn the point at 0.2um from the curr point
                  
                  ////////std::cout<<"after while"<<std::endl;
                
                  //////std::cout<<"dist  "<<dist<<std::endl;
                  
                  unsigned int num_of_pts = ((dist / RESAMPLING_DISTANCE) );
                  //////std::cout<<"#pts  "<<num_of_pts<<std::endl;
                  
                  // Find the vector between curr and next point
                  double vec_x = (*next_pt_it)[0] -  (*curr_pt_it)[0];
                  double vec_y = (*next_pt_it)[1] -  (*curr_pt_it)[1];
                  double vec_z = (*next_pt_it)[2] -  (*curr_pt_it)[2];
                  
                  
                    
                  
                  // find how many points are needed in between
                  //dist = Utility::euclidean_distance(*next_pt_it, *curr_pt_it, 3,1);
      //             //////std::cout<<"dist "<<dist<<std::endl;
                  
                  
                  
                  //std::list< double * >::iterator insert_it = curr_pt_it;
                  // Find out the coordinates of these points
                  
                  ////////std::cout<<num_of_pts<<std::endl;
                  for(double j =0; j < num_of_pts; j++)
                  {
                          double * new_pt = new double[3];
                  
                          ////////std::cout<<"in insert loop"<<std::endl;
                          
                          //////std::cout<<(*curr_pt_it)[0]<<" "<<(*curr_pt_it)[1]<<" "<<(*curr_pt_it)[2]<<std::endl;
                          ////////std::cout<<j<<std::endl;
                          new_pt[0] =  (*curr_pt_it)[0] + (((j+1)* RESAMPLING_DISTANCE)/*-prev_offset*/) * (vec_x / dist) ;
                          new_pt[1] =  (*curr_pt_it)[1] + (((j+1)* RESAMPLING_DISTANCE)/*-prev_offset*/) * (vec_y / dist) ;
                          new_pt[2] =  (*curr_pt_it)[2] + (((j+1)* RESAMPLING_DISTANCE)/*-prev_offset*/) * (vec_z / dist) ;
                      
                          //////std::cout<<new_pt[0]<<" "<<new_pt[1]<<" "<<new_pt[2]<<std::endl;
  //             
                          tempEdgePoints.push_back(new_pt);
                          
                          if(j == (num_of_pts-1))
                          {
                            // offset is the distance between the resampled last point and the last point of the original edge
                            prev_offset = Utility::euclidean_distance(*next_pt_it, new_pt, 3,1);
                            
                            // we also need to add the original edge's first and last points
                            
                          }
                      //insert_it++;
                      // insert takes the iterator position where the new entry is to be added and the new entry as parameters
                      //currentEdge->edgePointCoordinates.insert(next_pt_it, new_pt);
                      
                  }
                  
                  
                  
                  // insert to the list
                  curr_pt_it = next_pt_it;
                  next_pt_it++;
                  
                  // if this is the last point of the list
                  if(next_pt_it == currentEdge->edgePointCoordinates.end())
                  {
                      // add the fist and last points to the temp edge points
                      tempEdgePoints.push_back(*curr_pt_it);
                    
                  }
                  
               }
          }
          
        }
        
        // create the new edge
        if(tempEdgePoints.size())
        {
            //std::list< double * >::iterator tempEdgePointIterator;
            //////std::cout<<"in add edge"<<std::endl;
                    
            Edge * tmpEdge = new Edge(currentEdge->edgeConnectivity, tempEdgePoints.size(), currentEdge->label, tempEdgePoints);
            destedges.push_back(tmpEdge);
        }
        
    }
    
   
        
    // Write the new graph to a file
    if(destedges.size())
    {
      
        std::vector< Edge * >::iterator edgeIter;
        
        for(edgeIter = destedges.begin(); edgeIter != destedges.end(); ++edgeIter)
        resampled_graph->addEdge(*edgeIter);
        
        
        //////std::cout<<"New Edges"<<std::endl;
         
    }
    
    for(int i=0; i<numofVertices; i++)  //for each edge
    {
        
        Vertex * currentVertex = vertices->at(i);
        resampled_graph->addVertex(currentVertex);
        double * coords = (currentVertex)->coordinates;
        
    }
    
    //////std::cout<<"before write"<<std::endl;
    // write out the graph
    
    //mergedAmiraReader->writeSpatialGraphFile2(resampled_graph);
    mergedAmiraReader->writeSpatialGraphFile2(argv[1], resampled_graph);
    
    
    return 0;
    

    
//#endif
};
