

#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"
#include "utility.h"
#include <execinfo.h>
#include <signal.h>
#include "string.h"

using namespace std;
#define PROX_DIST 4
#define SPATIAL_GRAPH_SAMPLING_DISTANCE 0.1
#define SOMA_RADIUS_FOR_STARTING_VERTICES 50 // in microns

#define PROX_COL_NUM_ID 0
#define PROX_COL_NUM_COORDS 3 
#define PROX_COL_NUM_DIST 4
#define PROX_COL_NUM_FROM_CELL_ID 5
#define PROX_COL_NUM_TO_CELL_ID 6
#define PROX_COL_NUM_DIST_FROM_AXON 7
#define PROX_COL_NUM_DIST_FROM_DEND 8
#define PROX_COL_NUM_SEC_NUM_1 9
#define PROX_COL_NUM_SEC_NUM_2 10
#define PROX_COL_NUM_DIST_FROM_PRESYNAPTIC_SMOOTHED 11
#define PROX_COL_NUM_DIST_FROM_POSTSYNAPTIC_SMOOTHED 12
#define PROX_COL_NUM_DIST_FROM_PRESYNAPTIC_EUCLIDEAN 13
#define PROX_COL_NUM_DIST_FROM_POSTSYNAPTIC_EUCLIDEAN 14
#define PROX_COL_NUM_DIST_FROM_PRESYNAPTIC_WITHIN_EDGE 15
#define PROX_COL_NUM_DIST_FROM_POSTSYNAPTIC_WITHIN_EDGE 16
#define PROX_COL_NUM_DIST_FROM_PRESYNAPTIC_WITHIN_EDGE_EUCLIDEAN 17
#define PROX_COL_NUM_DIST_FROM_POSTSYNAPTIC_WITHIN_EDGE_EUCLIDEAN 18
#define PROX_COL_NUM_PRESYNAPTIC_LEVEL 19
#define PROX_COL_NUM_POSTSYNAPTIC_LEVEL 20
#define PROX_COL_NUM_PRESYNAPTIC_TYPE 21
#define PROX_COL_NUM_POSTSYNAPTIC_TYPE 22

// globals
//vars to read from csv file
std::vector<string >prox_lines;
std::vector<int >prox_ids;
std::vector<double *>prox_coords;
std::vector<double >prox_dists;
std::vector<int >prox_from_cell_index;
std::vector<int >prox_to_cell_index;
std::vector<double >prox_presynaptic_pathlen;
std::vector<double >prox_postsynaptic_pathlen;
std::vector<int >prox_sec1;
std::vector<int >prox_sec2;
std::vector<double >prox_presynaptic_pathlen_smoothed;
std::vector<double >prox_postsynaptic_pathlen_smoothed;
std::vector<double >prox_presynaptic_dist_euclidean;
std::vector<double >prox_postsynaptic_dist_euclidean;
std::vector<double >prox_presynaptic_pathlen_within_edge;
std::vector<double >prox_postsynaptic_pathlen_within_edge;
std::vector<double >prox_presynaptic_dist_euclidean_within_edge;
std::vector<double >prox_postsynaptic_dist_euclidean_within_edge;
std::vector<double >prox_presynaptic_edge_level;
std::vector<double >prox_postsynaptic_edge_level;
std::vector<string >prox_presynaptic_type;
std::vector<string >prox_postsynaptic_type;


std::vector<double *>soma_landmarks;
std::vector<double *>prox_landmarks;

std::vector<std::vector< double * >*> starting_vertices;

Reader * mergedAmiraReader;
AmiraSpatialGraph * amira_graph;
    
std::vector<AmiraSpatialGraph*>  cell_specific_spatial_graps_vector;
std::vector<std::vector<Edge*> * > cell_specific_axon_edges_list;
std::vector< std::vector<Edge*> * > cell_specific_dendrite_edges_list;
std::vector< std::vector<Edge*> * > cell_specific_soma_edges_list;
std::vector<std::vector<Vertex*> * > cell_specific_axon_vertices_list;
std::vector< std::vector<Vertex*> * > cell_specific_dendrite_vertices_list;
std::vector< std::vector<Vertex*> * > cell_specific_soma_vertices_list;

std::vector<double*>*  soma_locations_vector;
std::vector<double*>*  prox_locations_vector;

std::queue<Vertex*> unprocessed_vertex_que;

int pop_count = 0;
int main_push_cnt = 0;
int push_cnt = 0;

std::vector<Vertex*> pushed_stack;




void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}


bool foundInPushedStack(double * coordinate)
{
  for (int i = 0; i < pushed_stack.size(); i++)
  {
     if( (pushed_stack.at(i)->coordinates[0] == coordinate[0]) &&
	 (pushed_stack.at(i)->coordinates[1] == coordinate[1]) &&
	 (pushed_stack.at(i)->coordinates[2] == coordinate[2]))
	 {
	   
	   return true;
	 }
  }
  
  return false;
  
}




bool isAxon(int ID)
{
    std::vector<int >* id_array = amira_graph->getAxonIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
        }
        
    }
    
    return false;
    
    
};

bool isBasalDendrite(int ID)
{
    std::vector<int >* id_array = amira_graph->getBasalDendriteIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
        }
        
    }
    
    return false;
    
    
};

bool isApicalDendrite(int ID)
{
    std::vector<int >* id_array = amira_graph->getApicalDendriteIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
        }
        
    }
    
    return false;
    
    
};

bool isSoma(int ID)
{
    
    std::vector<int >* id_array = amira_graph->getSomaIDArray();
    
    for(int i =0; i < id_array->size(); i++)
    {
        if(id_array->at(i) == ID)
        {
            return true;
            //std::cout << "isSoma Found Soma with id "<< ID << std::endl;
        }
        
    }
    //std::cout << "isSoma not Found Soma with id "<< ID << std::endl;
    return false;
    
    
};

void readCsvfile(const char* filename)
{
    // read the csv file
    std::ifstream xlsreader(filename);//xlsreader.open(inputxlsfile, std::fstream::in | std::fstream::out);
    std::string currentLine = "";
    int line_number = 0;
    
    //xlsWriter.open (xlsfilename,  std::fstream::in |  std::fstream::app );
    
    if(!xlsreader.fail())
    {
        ////std::cout<<"opened fine "<<inputxlsfile<<std::endl;
	
        while(!std::getline(xlsreader, currentLine).eof())
        {
            line_number++;
            
            if(line_number>1)
            {
                //prox_line.push_back(currentLine);
                char *currLinetmp =  new char[currentLine.length()+1];
                int col_num = 0;
                char * bla;
                strcpy(currLinetmp,currentLine.c_str());
                bla = strtok(currLinetmp, ",");
                prox_ids.push_back( atoi( bla) );
                //std::cout<<"prox_id = "<<bla<<std::endl;
                double * tmp_prox_coords = new double[3];
                int tmp_prox_coords_ind = 0;
                while (bla != NULL)
                {
                    
                    //std::cout<<"currentLine split "<<bla<<std::endl;
                    col_num++;
                    bla = strtok(NULL, ",");
                    if(col_num <= PROX_COL_NUM_COORDS && col_num > PROX_COL_NUM_ID )
                    {
                        //cout << tmp_prox_coords_ind << endl;
                        tmp_prox_coords[tmp_prox_coords_ind] = atof(bla) ;
                        //cout << tmp_prox_coords[tmp_prox_coords_ind]<<endl;
                        tmp_prox_coords_ind++;
                    }
                    else if(col_num == PROX_COL_NUM_DIST)
                    {
                        prox_dists.push_back(atof(bla));
                        //cout<<atof(bla)<<endl;
                        
                    }
                    else if(col_num == PROX_COL_NUM_FROM_CELL_ID)
                    {
                        prox_from_cell_index.push_back(atoi(bla));
                    }
                    else if(col_num == PROX_COL_NUM_TO_CELL_ID) 
                    {
                        prox_to_cell_index.push_back(atoi(bla));
                    }
                    else if(col_num == PROX_COL_NUM_DIST_FROM_AXON) 
                    {
                        prox_presynaptic_pathlen.push_back(atof(bla));
                    }
                    else if(col_num == PROX_COL_NUM_DIST_FROM_DEND) 
                    {
                        prox_postsynaptic_pathlen.push_back(atof(bla));
                    }
                    else if(col_num == PROX_COL_NUM_SEC_NUM_1) 
                    {
                        //prox_sec1.push_back(atoi(bla));
                    }
                    else if(col_num == PROX_COL_NUM_SEC_NUM_2) 
                    {
                        //prox_sec2.push_back(atoi(bla));
                    }
                    
                }
                prox_coords.push_back(tmp_prox_coords);
                
                if(col_num < PROX_COL_NUM_SEC_NUM_1+1)
                {
                    // there is no second sec for this prox
                    // make it zero so that we have same size for all proxes
                    currentLine.append(",0,0");
                    
                }
                else if(col_num < PROX_COL_NUM_SEC_NUM_2+1)
                {
                    // there is no second sec for this prox
                    // make it zero so that we have same size for all proxes
                    currentLine.append(",0");
                    
                }
                prox_presynaptic_pathlen_smoothed.push_back(0.0);
                prox_postsynaptic_pathlen_smoothed.push_back(0.0);
                
                prox_presynaptic_dist_euclidean.push_back(0.0);
                prox_postsynaptic_dist_euclidean.push_back(0.0);
                
                prox_presynaptic_pathlen_within_edge.push_back(0.0);
                prox_postsynaptic_pathlen_within_edge.push_back(0.0);
                
                prox_presynaptic_dist_euclidean_within_edge.push_back(0.0);
                prox_postsynaptic_dist_euclidean_within_edge.push_back(0.0);
                
                prox_presynaptic_edge_level.push_back(0.0);
                prox_postsynaptic_edge_level.push_back(0.0);
                
                prox_presynaptic_type.push_back("Axon");
                prox_postsynaptic_type.push_back("Soma");
                
                // Add this to the global list of lines
                prox_lines.push_back(currentLine);
                //cout << currentLine<<endl;
            }   
        }
    }
    
}


void readLandmarksfile(const char* filename, std::vector<double *>*landmarks)
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
    ////std::cout<<"valid path"<<std::endl;
    while(!std::getline(inputStream, currentLine).eof() && currentLine.find("@1") != 0 ){}
    while(!std::getline(inputStream, currentLine).eof() && currentLine.size() > 0 )
    { 
      double * tmp = new double[ARRAY_LENGTH];
      
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[X_COORD], &tmp[Y_COORD], &tmp[Z_COORD]);
      //std::cout<<" "<<tmp[X_COORD]<<" "<<tmp[Y_COORD]<<" "<<tmp[Z_COORD]<<std::endl;
      
      landmarks->push_back(tmp);
      
    }
    ////std::cout<<"done reading prox landmarks"<<std::endl;
  }

};

void getCellSpecificSpatialGraphs()
{
   int num_of_cells = amira_graph->getNumberofCells();
   
   std::cout<< "Number of Cells "<< num_of_cells<< std::endl;
  
   // init the cell specific spatial graphs
   for(int i = 0; i < num_of_cells; i++)
   {
      cell_specific_spatial_graps_vector.push_back(new AmiraSpatialGraph);
      
      cell_specific_axon_edges_list.push_back( new std::vector<Edge*> );
      cell_specific_dendrite_edges_list.push_back( new std::vector<Edge*> );
      cell_specific_soma_edges_list.push_back( new std::vector<Edge*> );
      cell_specific_axon_vertices_list.push_back( new std::vector<Vertex*> );
      cell_specific_dendrite_vertices_list.push_back( new std::vector<Vertex*> );
      cell_specific_soma_vertices_list.push_back( new std::vector<Vertex*> );
      
   }
   
   
   std::cout<< "intialized ind spa graphs "<<std::endl;
   
   std::vector< Edge * > * edges = amira_graph->edgesPointer();
   std::vector< Vertex * > * vertices = amira_graph->verticesPointer();
   
   int vertexNumberMapping[vertices->size()];
   int cellSpecificVertexNumber[num_of_cells];
   int cellSpecificAxonVertexNumber[num_of_cells];
   int cellSpecificDendriteVertexNumber[num_of_cells];
   int cellSpecificSomaVertexNumber[num_of_cells];
   
   
   // init vertex counts to 0
   for(int i = 0; i < num_of_cells; i++)
   {
      cellSpecificVertexNumber[i] = 0;
      cellSpecificAxonVertexNumber[i] = 0;
      cellSpecificDendriteVertexNumber[i] = 0;
      cellSpecificSomaVertexNumber[i] = 0;
   }
   std::cout<<"edge count"<<edges->size()<<std::endl; 
   std::cout<<"vertex count"<<vertices->size()<<std::endl; 
   
   for(int vertexcount = 0; vertexcount < vertices->size(); vertexcount++)    
   { 
      // get cell index 
      //std::cout<<"edge count"<<edges->size()<<std::endl; 
      //std::cout<<"vertex count"<<vertices->size()<<std::endl;
      //std::cout<< "vertexind" << vertexcount<< std::endl;
      int cell_index = 0;
      bool isItAxon = true;
      Vertex * currVert = vertices->at(vertexcount);
       //std::cout<< " got vertex  "<<  std::endl;
       //std::cout<< "curr vert valid "<< currVert->connectedEdges<<std::endl;
      std::vector<Edge *> * connectedEdges = vertices->at(vertexcount)->connectedEdges;
      //std::cout<< "connectedEdges size " << connectedEdges->size() << std::endl;
      //std::cout<< " Connected edge label "<< vertices->at(vertexcount)->connectedEdges->at(0)->label << std::endl;
      
      // Do stuff only if there are connected edges
      if(connectedEdges->size() == 0)
      {
          break;
      }
      
      
      Edge * connectedFirstEdge = connectedEdges->at(0);
      int connetedEdgeLabel = connectedFirstEdge->label;
      //std::cout<< " Connected edge label "<< connectedFirstEdge->label << std::endl;
      
      if(isAxon(connetedEdgeLabel))
      {
        cell_index = amira_graph->getAxonCellIndex(connetedEdgeLabel);
      }
      else if(isBasalDendrite(connetedEdgeLabel))
      {
        cell_index = amira_graph->getDendriteCellIndex(connetedEdgeLabel);
      }
      else if(isApicalDendrite(connetedEdgeLabel))
      {
        cell_index = amira_graph->getDendriteCellIndex(connetedEdgeLabel);
      }
      else if(isSoma(connetedEdgeLabel))
      {
        cell_index = amira_graph->getSomaCellIndex(connetedEdgeLabel);
      }
      //std::cout<<"cell index "<< cell_index<< std::endl;
      
      // update the vertex number mapping array
      vertexNumberMapping[vertexcount] = cellSpecificVertexNumber[cell_index];
      cellSpecificVertexNumber[cell_index]++;
     
      
      Vertex * tempVertex = new Vertex(vertices->at(vertexcount), false);
      //std::cout<< "got new vertex with label "<< vertices->at(vertexcount)->label<<std::endl;
      cell_specific_spatial_graps_vector.at(cell_index)->verticesPointer()->push_back(tempVertex);
      
   }
   
   
   for(int edgecount = 0; edgecount < edges->size(); edgecount++)    
   {
     
      // get cell index 
      int cell_index = 0;
      bool isItAxon = true;
      
      if(isAxon(edges->at(edgecount)->label))
      {
        cell_index = amira_graph->getAxonCellIndex(edges->at(edgecount)->label);
      }
      else if(isBasalDendrite(edges->at(edgecount)->label))
      {
        cell_index = amira_graph->getDendriteCellIndex(edges->at(edgecount)->label);
      }
      else if(isApicalDendrite(edges->at(edgecount)->label))
      {
        cell_index = amira_graph->getDendriteCellIndex(edges->at(edgecount)->label);
      }
      else if(isSoma(edges->at(edgecount)->label))
      {
        cell_index = amira_graph->getSomaCellIndex(edges->at(edgecount)->label);
      }
      
      //std::cout<< "got new edge with cell ind "<< cell_index<<std::endl;
      
      int new_edge_connectivity_ind_1,new_edge_connectivity_ind_2;
      new_edge_connectivity_ind_1 = vertexNumberMapping[ edges->at(edgecount)->edgeConnectivity[0] ];
      new_edge_connectivity_ind_2 = vertexNumberMapping[ edges->at(edgecount)->edgeConnectivity[1] ];
      
      //std::cout<< "new edge connect 1 = " << new_edge_connectivity_ind_1 <<  " new edge connect 2 = " << new_edge_connectivity_ind_2<<std::endl;
      Edge * tempEdge = new Edge(edges->at(edgecount), new_edge_connectivity_ind_1,new_edge_connectivity_ind_2 );
      
      cell_specific_spatial_graps_vector.at(cell_index)->edgesPointer()->push_back(tempEdge);
      
      // update the connected edges field for the from to to vertex
      int from_vertex = tempEdge->edgeConnectivity[0];
      int to_vertex = tempEdge->edgeConnectivity[1];
      
     
      cell_specific_spatial_graps_vector.at(cell_index)->verticesPointer()->at(from_vertex)->connectedEdges->push_back(tempEdge);
      cell_specific_spatial_graps_vector.at(cell_index)->verticesPointer()->at(to_vertex)->connectedEdges->push_back(tempEdge);
       
      
   }
   
}

void populateProxInfo(Edge * curr_edge, Vertex * parent_vertex, Vertex * child_vertex, double * soma_coords, int cell_index, int basal_neurite_tree_number)
{
    for(int i = 0; i < prox_lines.size(); i++)
    {
        if(((isAxon(curr_edge->label) == true) && (prox_from_cell_index.at(i) == cell_index+1))
            ||((isBasalDendrite(curr_edge->label) == true) && (prox_to_cell_index.at(i) == cell_index+1))
            ||((isApicalDendrite(curr_edge->label) == true) && (prox_to_cell_index.at(i) == cell_index+1))
            ||((isSoma(curr_edge->label) == true) && (prox_to_cell_index.at(i) == cell_index+1)))
        {
            // find the point closest to this prox in this edge
            double min_dist_prox = 1000000;
            double * min_dist_prox_coords = new double[3];
            std::list< double * >::iterator edge_it;
            for( edge_it = curr_edge->edgePointCoordinates.begin();edge_it != curr_edge->edgePointCoordinates.end(); edge_it++) 
            {
                double * coords = *edge_it;
                //cout << "prox_landmark = "<< prox_landmarks.at(i)[0] << " "<< prox_landmarks.at(i)[1]<<" "<<prox_landmarks.at(i)[2]<<endl;
                double dist = Utility::euclidean_distance(prox_landmarks.at(i), coords, 3, 1);
                if(dist < min_dist_prox)
                {
                    min_dist_prox = dist;
                    min_dist_prox_coords[0] = coords[0];
                    min_dist_prox_coords[1] = coords[1];
                    min_dist_prox_coords[2] = coords[2];
                }
            }
            if(min_dist_prox < PROX_DIST)
            {
                cout << "found prox coords " << min_dist_prox_coords[0] << " " << min_dist_prox_coords[1] << " " << min_dist_prox_coords[2] << " " << endl;
                if((isAxon(curr_edge->label) == true) && (prox_from_cell_index.at(i) == cell_index+1))
                {
                    if( (curr_edge->edgePointCoordinates.front()[0] == parent_vertex->coordinates[0])
                        &&(curr_edge->edgePointCoordinates.front()[1] == parent_vertex->coordinates[1])
                        &&(curr_edge->edgePointCoordinates.front()[2] == parent_vertex->coordinates[2]))
                    {
                        double curr_edge_pathlen_from_parent_node = curr_edge->segmentLength(min_dist_prox_coords,true);
                    }
                    else
                    {
                        double curr_edge_pathlen_from_parent_node = curr_edge->segmentLength(min_dist_prox_coords,false);
                    }
                    
                    double curr_edge_pathlen_from_parent_node = curr_edge->segmentLength(min_dist_prox_coords,curr_edge->forward);
                
                    prox_presynaptic_pathlen_smoothed.at(i) = parent_vertex->distanceFromSoma + curr_edge_pathlen_from_parent_node ;
                    
                    prox_presynaptic_dist_euclidean.at(i) = Utility::euclidean_distance(soma_coords,min_dist_prox_coords,3,1);
                    
                    prox_presynaptic_pathlen_within_edge.at(i) = curr_edge_pathlen_from_parent_node;
                    
                    prox_presynaptic_dist_euclidean_within_edge.at(i) = Utility::euclidean_distance(parent_vertex->coordinates,min_dist_prox_coords,3,1);
                    
                    prox_presynaptic_edge_level.at(i) = curr_edge->level;
                    
                    prox_presynaptic_type.at(i) = ("Axon");
                    
                    cout << "prox_presynaptic_pathlen_smoothed for # " << i << " " << prox_presynaptic_pathlen_smoothed.at(i) <<endl;
                    cout << "prox_presynaptic_dist_euclidean for # " << i << " " << prox_presynaptic_dist_euclidean.at(i) <<endl;
                    cout << "prox_presynaptic_pathlen_within_edge for # " << i << " " << prox_presynaptic_pathlen_within_edge.at(i) <<endl;
                    cout << "prox_presynaptic_dist_euclidean_within_edge for # " << i << " " << prox_presynaptic_dist_euclidean_within_edge.at(i) <<endl;
                    cout << "prox_presynaptic_edge_level for # " << i << " " << prox_presynaptic_edge_level.at(i) <<endl;
                    cout << "prox_presynaptic_type for # " << i << " " << prox_presynaptic_type.at(i) <<endl;
                    
                
                }
                else
                {
                    double curr_edge_pathlen_from_parent_node = curr_edge->segmentLength(min_dist_prox_coords,curr_edge->forward);
                
                    prox_postsynaptic_pathlen_smoothed.at(i) = parent_vertex->distanceFromSoma + curr_edge_pathlen_from_parent_node ;
                    
                    prox_postsynaptic_dist_euclidean.at(i) = Utility::euclidean_distance(soma_coords,min_dist_prox_coords,3,1);
                    
                    prox_postsynaptic_pathlen_within_edge.at(i) = curr_edge_pathlen_from_parent_node;
                    
                    prox_postsynaptic_dist_euclidean_within_edge.at(i) = Utility::euclidean_distance(parent_vertex->coordinates,min_dist_prox_coords,3,1);
                    
                    prox_postsynaptic_edge_level.at(i) = curr_edge->level;
                    
                    if(isBasalDendrite(curr_edge->label))
                    {
                        char str[100];
                        strcpy(str, "Dendrite_Basal_");
                        strcat(str, Utility::inttochar(basal_neurite_tree_number));
                        prox_postsynaptic_type.at(i) = str;
                    }
                    else if (isApicalDendrite(curr_edge->label))
                        prox_postsynaptic_type.at(i) = ("Dendrite_Apical");
                    else if(isSoma(curr_edge->label))
                        prox_postsynaptic_type.at(i) = ("Soma");
                    
                    cout << "prox_postsynaptic_pathlen_smoothed for # " << i << " " << prox_postsynaptic_pathlen_smoothed.at(i) <<endl;
                    cout << "prox_postsynaptic_dist_euclidean for # " << i << " " << prox_postsynaptic_dist_euclidean.at(i) <<endl;
                    cout << "prox_postsynaptic_pathlen_within_edge for # " << i << " " << prox_postsynaptic_pathlen_within_edge.at(i) <<endl;
                    cout << "prox_postsynaptic_dist_euclidean_within_edge for # " << i << " " << prox_postsynaptic_dist_euclidean_within_edge.at(i) <<endl;
                    cout << "prox_postsynaptic_edge_level for # " << i << " " << prox_postsynaptic_edge_level.at(i) <<endl;
                    cout << "prox_postsynaptic_type for # " << i << " " << prox_postsynaptic_type.at(i) <<endl;
                    
                }
                 
            }
        }
    
        
    }
}


void setDistanceToVertex( int cell_index, int edgelabel, double * somaCoords, int basal_neurite_tree_number)
{
    std::vector< Vertex * > * vertices = cell_specific_spatial_graps_vector.at(cell_index)->verticesPointer();
    bool is_axon = false;
    while( !unprocessed_vertex_que.empty())
    {
        
        Vertex * curr_vertex = unprocessed_vertex_que.front();
        
        //std::cout<< " popping vertex " << std::endl;
        //std::cout<< " processing vertex "<<curr_vertex->coordinates[0]<< ", " << curr_vertex->coordinates[1]<< ", "<< curr_vertex->coordinates[2]<< " cell id "<<cell_index << " axon "<<is_axon <<std::endl;
        
        unprocessed_vertex_que.pop();
        pop_count++;
        // set distance from soma to this node
        //curr_vertex->distanceFromSoma = 5;
        
        // find this node's distance from the parent node
        //curr_vertex->distanceFromSomaEuclidean = Utility::euclidean_distance(curr_vertex->coordinates,vertices->at(curr_vertex->parentVertexID)->coordinates,3,1)
        
        // check if there are connected unprocessed nodes to this node
        std::vector<Edge*> *connectedEdgesPtr = curr_vertex->getConnectedEdges();
        Edge* min_dist_edge = NULL;
        double min_dist = 0xffffffff;
        //std::cout<< " connectedEdgesPtr->size() " << connectedEdgesPtr->size()<< std::endl;
        
        // Go through all the edges connected this node to grow the tree
        for( int i = 0; i < connectedEdgesPtr->size(); i++)
        {
            Edge * curr_edge = connectedEdgesPtr->at(i);
            
            // process this edge only if the father and child ids are not set
            if((curr_edge->label == edgelabel) && (curr_edge->fatherID == -1) && (curr_edge->childID == -1))
            {
                // set and verify parent and child vertices for this edge
                // amira spatial graph does not guarentee the forward direction
                // from parent to child
                int childVertexID, parentVertexID;
                // assume the forward diretion
                parentVertexID = curr_edge->edgeConnectivity[0];
                childVertexID = curr_edge->edgeConnectivity[1];
                
                // verify whether this assumtion is correct
                // avoid loops
                if( (vertices->at(parentVertexID)->coordinates[0] == vertices->at(childVertexID)->coordinates[0]) &&
                    (vertices->at(parentVertexID)->coordinates[1] == vertices->at(childVertexID)->coordinates[1]) &&
                    (vertices->at(parentVertexID)->coordinates[2] == vertices->at(childVertexID)->coordinates[2]) )
                {
                    // check to see if it is a loop
                    continue;
                }
                else if ((vertices->at(parentVertexID)->coordinates[0] == curr_vertex->coordinates[0]) &&
                        (vertices->at(parentVertexID)->coordinates[1] == curr_vertex->coordinates[1]) &&
                        (vertices->at(parentVertexID)->coordinates[2] == curr_vertex->coordinates[2]) )
                {
                    curr_edge->fatherID = parentVertexID;
                    curr_edge->childID = childVertexID;
                    curr_edge->forward = true;
                    connectedEdgesPtr->at(i)->forward = true;
                }
                else
                {
                    curr_edge->childID = parentVertexID;
                    curr_edge->fatherID = childVertexID;
                    curr_edge->forward = false;
                    connectedEdgesPtr->at(i)->forward = false;
                }
                
                Vertex * parentVertex = vertices->at(curr_edge->fatherID);
                Vertex * childVertex = vertices->at(curr_edge->childID);
                
                if(curr_vertex->vertexLevel == 0)
                {
                    // this is the root node. set father vertex distance
                    parentVertex->distanceFromSoma = Utility::euclidean_distance(parentVertex->coordinates,somaCoords,3,1);
                }
                
                
                childVertex->vertexLevel = parentVertex->vertexLevel + 1;
                // set child node distance
                double curr_edge_length = curr_edge->segmentLength() ;//curr_edge->numEdgePoints * SPATIAL_GRAPH_SAMPLING_DISTANCE;
                
                childVertex->distanceFromSoma = parentVertex->distanceFromSoma + curr_edge_length;
                
                curr_edge->level = parentVertex->vertexLevel;
                
                // Now push the child node for growing the tree
                unprocessed_vertex_que.push(childVertex);
                
                cout << "direction " << curr_edge->forward << endl;
                cout << " Parent Vertex ID "<< curr_edge->fatherID << endl;
                cout << " Parent Vertex coords "<< parentVertex->coordinates[0] << " " << parentVertex->coordinates[1] << " "<<parentVertex->coordinates[2] << endl;
                cout << " Parent Vertex distance "<< parentVertex->distanceFromSoma << endl;
                cout << " child Vertex ID "<< curr_edge->childID << endl;
                cout << " child Vertex coords "<< childVertex->coordinates[0] << " " << childVertex->coordinates[1] << " "<<childVertex->coordinates[2] << endl;
                cout << " Current edge level "<< curr_edge->level;
                cout << " Current edge seg length "<< curr_edge_length;
                cout << " Current edge Euclidean length "<< Utility::euclidean_distance(parentVertex->coordinates,childVertex->coordinates,3,1);
                cout << endl;
                
                populateProxInfo(curr_edge,parentVertex,childVertex,somaCoords,cell_index,basal_neurite_tree_number);
            
            }
            
        }
        
        //edge_level++;
        
        //setDistanceToVertex(cell_index, edgelabel, somaCoords, ++edge_level);
  
    
    }
    
    
  
}







std::list<Vertex *> sortStartingVertices(double * somaLandmark, std::vector<Vertex *> verticestobesorted)
{
    double dist = 0xFFFF;
    double temp_dist = 0;
    std::list<Vertex *> sortedVertices;
    std::vector<double> distancesFromSoma;
    
    std::cout << "given soma landmark "<<somaLandmark[0] << " "<<somaLandmark[1] <<" "<< somaLandmark[2] << std::endl;
    
    std::cout << "vertices to be sorted "<< std::endl;
    for (int i = 0; i <verticestobesorted.size(); i++ )
    {
        //sortedVertices.push_back();
        std::cout << verticestobesorted.at(i)->coordinates[0] << " "<<verticestobesorted.at(i)->coordinates[1] <<" "<< verticestobesorted.at(i)->coordinates[2] << std::endl;
    }
    
    // put the first guy
    Vertex * tmpvertex = new Vertex(verticestobesorted.at(0));
   
    tmpvertex->distanceFromSomaEuclidean = Utility::euclidean_distance(somaLandmark, tmpvertex->coordinates, 3, 1);
    sortedVertices.push_back(tmpvertex);
    
    for(int i = 1; i < verticestobesorted.size(); i++)
    {
        Vertex * tmpvertex = new Vertex(verticestobesorted.at(i));
        tmpvertex->distanceFromSomaEuclidean =  Utility::euclidean_distance(somaLandmark, tmpvertex->coordinates, 3, 1);
        
        /*tmpvertex[0] = verticestobesorted.at(i)->coordinates[0];
        tmpvertex[1] = verticestobesorted.at(i)->coordinates[1];
        tmpvertex[2] = verticestobesorted.at(i)->coordinates[2];
        tmpvertex[3] = temp_dist;*/
        std::list<Vertex *>::iterator it = sortedVertices.begin();
        std::list<Vertex *>::iterator nextit = sortedVertices.begin();
        nextit++;
        for(; it != sortedVertices.end(); it++,nextit++)
        {
            if(nextit == sortedVertices.end())
            {
                //cout << "last guy "<<endl;
                if((*it)->distanceFromSomaEuclidean < tmpvertex->distanceFromSomaEuclidean )
                {
                    sortedVertices.push_back(tmpvertex);
                }
                else
                {
                    sortedVertices.push_front(tmpvertex);
                }
            
            }
            else
            {
                //cout << " not last guy " << endl;
                
                 if(((*it)->distanceFromSomaEuclidean < tmpvertex->distanceFromSomaEuclidean) && ((*nextit)->distanceFromSomaEuclidean >= tmpvertex->distanceFromSomaEuclidean) )
                {
                    sortedVertices.insert(nextit,tmpvertex);
                    break;
                    
                }
                
            }
        }
        
        
    }
    
    std::cout << "sorted vertices "<< std::endl;
    std::list<Vertex *>::iterator listit = sortedVertices.begin();
    for (;listit != sortedVertices.end(); listit++)
    {
        //sortedVertices.push_back();
        std::cout << (*listit)->coordinates[0] << " "<<(*listit)->coordinates[1] <<" "<< (*listit)->coordinates[2] << " "<< (*listit)->distanceFromSomaEuclidean<< std::endl;
    }
    
    return sortedVertices;
    
}

void getStartingVertices()
{
    double min_dist_around_soma = SOMA_RADIUS_FOR_STARTING_VERTICES;
   int pre_ind = 0;
   int post_ind = 0;
   bool is_axon = false;
   
   for(int i = 0; i < soma_landmarks.size(); i++)
   {
        std::cout << "soma # "<< i << std::endl;
        AmiraSpatialGraph * present_graph = cell_specific_spatial_graps_vector.at(i);
        std::vector< Edge * > * edges = present_graph->edgesPointer();
        std::vector< Vertex * > * vertices = present_graph->verticesPointer();
        std::vector<Vertex *> startingverticesforthiscell;
        std::vector<int> basal_neurite_numbers;
        
        
        for (int j = 0; j < present_graph->getNumberOfVertices(); j++)
        {  
            //std::cout<< "init dist from soma "<<vertices->at(j)->distanceFromSoma << std::endl;
            //std::vector<double *> startingverticesforthiscell;
            double dist = Utility::euclidean_distance(soma_landmarks.at(i), vertices->at(j)->coordinates, 3, 1);
            if( ( vertices->at(j)->getConnectedEdges()->size() == 1 ) && ( dist < min_dist_around_soma) )
            {
                if((isSoma(vertices->at(j)->connectedEdges->at(0)->label )==false) && (vertices->at(j)->distanceFromSoma == 0)
                    && (amira_graph->getAxonDendriteOrSomaCellIndex(vertices->at(j)->connectedEdges->at(0)->label) == i))
                {
                    //std::cout<< " node number " << j <<std::endl;
                    //std::cout<< " label " <<vertices->at(j)->connectedEdges->at(0)->label<<std::endl;
                    //std::cout<< " connected edges " << vertices->at(j)->getConnectedEdges()->size()<<std::endl;
                    //std::cout << vertices->at(j)->coordinates[0] << " "<<vertices->at(j)->coordinates[1] << " "<<vertices->at(j)->coordinates[2] << " "<<std::endl;
                    
                    //startingverticesforthiscell.push_back(&vertices->at(j)->coordinates[0]);
                    startingverticesforthiscell.push_back(vertices->at(j));
                        
                }
                        
            }
        }
        
        // start setting up pathlengths for vertices with sorted verices near soma
        if(startingverticesforthiscell.size() > 0)
        {
            int basal_neurite_tree_number = 0;
            
            list<Vertex * > sortedVerticesList = sortStartingVertices(soma_landmarks.at(i),startingverticesforthiscell);
            
            for(list<Vertex * >::iterator it=sortedVerticesList.begin(); it != sortedVerticesList.end(); it++ )
            {
                //if((*it)->distanceFromSoma == 0)
                {
            
                    (*it)->vertexLevel = 0;
                    //(*it)->distanceFromSomaEuclidean = Utility::euclidean_distance((*it)->coordinates,soma_landmarks.at(i),3,1);
                    //(*it)->distanceFromSoma = (*it)->distanceFromSomaEuclidean;
                    unprocessed_vertex_que.push((*it));
                    //std::cout<< "from main pushed vertex number "<< j << std::endl;
                    std::cout<< " pushed from main vertex : "<< (*it)->coordinates[0]<< ", " << (*it)->coordinates[1]<< ", "<< (*it)->coordinates[2]<<std::endl;
                    std::cout<< " pushed from main vertex label : "<< (*it)->connectedEdges->at(0)->label<<std::endl;
                    push_cnt++;
                    main_push_cnt++;
                    pushed_stack.push_back((*it));
                    
                    // start the engine
                    // set the basal dendrite tree number first
                    if(isBasalDendrite((*it)->connectedEdges->at(0)->label))
                    {
                        basal_neurite_tree_number = basal_neurite_tree_number + 1;
                        setDistanceToVertex(i,(*it)->connectedEdges->at(0)->label,soma_landmarks.at(i), basal_neurite_tree_number ); 
                    }
                    else
                    {
                        setDistanceToVertex(i,(*it)->connectedEdges->at(0)->label,soma_landmarks.at(i), 0 );
                    }
                    
                    
                    //break;
                }
                    
            }
                    
        }
        
        //starting_vertices.push_back(&startingverticesforthiscell);
            
   }
    
}



void writeCsvFile(char * outputfile)
{
    std::fstream xlsParamWriter;
    // header
    xlsParamWriter.open (outputfile,   std::fstream::out /*| std::fstream::app*/ );
    

    xlsParamWriter<<"ID"<<","<<"Proximity_Location_X"<<","<<"Proximity_Location_Y"<<","<<"Proximity_Location_Z"<<","
    <<"Distance"<<","<<"From_Axon_of_Cell_#"<<","<<"To_Dendrite_of_Cell_#"<<","<<"Axon_Distance_From_Soma"<<","<<"Dendrite_Distance_From_Soma"<<","<<"Lower_Section_#"<<","<<"Upper_Section_#"<<","
    <<"Presynaptic_Pathlength_From_Soma"<<","<<"Postsynaptic_Pathlength_From_Soma"<<","
    <<"Presynaptic_Euclidean_Distance_From_Soma"<<","<<"Postsynaptic_Euclidean_Distance_From_Soma"<<","
    <<"Presynaptic_Pathlength_From_Parent_Node"<<","<<"Postsynaptic_Pathlength_From_Parent_Node"<<","
    <<"Presynaptic_Euclidean_Distance_From_Parent_Node"<<","<<"Postsynaptic_Euclidean_Distance_From_Parent_Node"<<","
    <<"Presynaptic_Edge_Level_From_Soma"<<","<<"Postsynaptic_Edge_Level_From_Soma"<<","
    <<"Presynaptic_Type"<<","<<"Postsynaptic_Type"<<'\n';
    
    
    for(int i = 0; i<prox_ids.size(); i++)
    {
        //xlsParamWriter<<prox_lines.at(i)<<'\n';
        xlsParamWriter<<prox_lines.at(i)<<","<<prox_presynaptic_pathlen_smoothed.at(i)<<","<<prox_postsynaptic_pathlen_smoothed.at(i)<<","<<
        prox_presynaptic_dist_euclidean.at(i)<<","<<prox_postsynaptic_dist_euclidean.at(i)<<","<<
        prox_presynaptic_pathlen_within_edge.at(i)<<","<<prox_postsynaptic_pathlen_within_edge.at(i)<<","<<
        prox_presynaptic_dist_euclidean_within_edge.at(i)<<","<<prox_postsynaptic_dist_euclidean_within_edge.at(i)<<","<<
        prox_presynaptic_edge_level.at(i)<<","<<prox_postsynaptic_edge_level.at(i)<<","<<
        prox_presynaptic_type.at(i)<<","<<prox_postsynaptic_type.at(i)<<","<<'\n';
        
        
    }
    
    
    xlsParamWriter.close();
    
}

int main(int argc , char * argv[])
{
    signal(SIGSEGV, handler);
    // vars
    //std::vector<double *>prox_landmarks;
    //std::vector<double *>soma_landmarks;
    
    // save input arguments
    char * spatial_graph = argv[1];
    char * prox_landmark_file = argv[2];
    char * soma_landmark_file = argv[3];
    char * prox_csv_file = argv[4];
    char * output_file = argv[5];
    
    // read inputs: spatial graphs and landmarks (scaled ?)
    readCsvfile(prox_csv_file);
    
    mergedAmiraReader = new Reader(argv[1],argv[2]);  
    mergedAmiraReader->readSpatialGraphFile(true);
    amira_graph = mergedAmiraReader->getSpatialGraph();

    getCellSpecificSpatialGraphs();
    
    // read landmarks
    readLandmarksfile(prox_landmark_file, &prox_landmarks );
    readLandmarksfile(soma_landmark_file, &soma_landmarks); 
    cout<<"Read soma landmarks"<<endl;
    
    getStartingVertices();
    
    // convert spatial graph into a tree with lengths, node numbers, segment lengths etc.
    //convertGraphIntoATree();
    cout << "Converted graph into a distance tree"<<endl;
    
    //write csv file with pathlenghts appended
    writeCsvFile(output_file);
    cout << "Wrote new csv file "<<endl;
    
}






