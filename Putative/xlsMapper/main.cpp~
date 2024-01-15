#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace std;

// this small script takes the xls sheet of proximities and adds manual landmarks to 
// the respective proximities. manual landmarks could be overlap, bouton or spine location 
// along with their confidence.
// if a prox has more than 1 matching manual landmark
// then we add the landmarks above or below this line
// based on the distance

// defines
#define MAX_NUM_OF_ARGUMENTS	12	
#define SECOND_MANUAL_LANDMARK  TRUE

enum MasterXlsHeaders
{
  ProxNum = 0,
  ProxCoordX,
  ProxCoordY,
  ProxCoordZ
  
};

// function protos
void printList(list<string *> * listPtr);
void printList(vector<double *> * coordPtr);
void printList(list<vector<double *> *> * listPtr);
void readLandmarkfile(char* filename, int weight, vector<double *>*prox_landmarks);
double euclidean_distance(double * a, double * b, int dim, float z_scalar);
void readCsvAsListofVectors(const char* masterXlsFile, list<vector<double*> *> * masterXlsList );
void writeListofVectorsAsCsv(const char* masterXlsFile, list<vector<double*> *> * masterXlsList );


// global declaration


int main(int argc, char **argv) 
{
    list<vector<double*> *> * masterXlsList = new list<vector<double*> *>;
    vector<double*> *landmarksVector = new vector<double *>;
  
    if( argc < MAX_NUM_OF_ARGUMENTS )
    {
        std::cout << "need "<<MAX_NUM_OF_ARGUMENTS<<"arguments" << std::endl;
	
        return 0;
    
    }
  
    //load the arguments into variables
    const char * masterXlsFile = argv[1]; // input xls file
    const char * masterXlsOutFile = argv[2]; // output xls file
    double landmarkConfidence = atoi(argv[3]); // ranking for the landmark ; now not used
    int ProxRadius = atoi(argv[4]); // size of prox image
    int init_num_cols = atoi(argv[5]); // where does the input xls end
    int final_num_cols = atoi(argv[6]); // end colmn num after adding new data
    int needZeroFilling = atoi(argv[7]); // if its the end of colms (like all boutons have been added); always true now
    int CorrdColNum = atoi(argv[8]); // colm to compare against
    int SecondManualLandmark = atoi(argv[9]); // is there secondary column to compare against
    int SecondaryCoordColNum = atoi(argv[10]); // colmn to compare against after comparing with primary coordinates
    int NoNewLine = atoi(argv[11]); // enabled for bouton and spine; so that no new overlap will be added falsely
    int numOfLandmarkFiles = atoi(argv[12]); // how many manual landmarks to compare against
    //const char * landmarkFile = argv[13]; // list of landmark files as input
    
    // read the master xls 
    readCsvAsListofVectors(masterXlsFile, masterXlsList);
    
cout << "==============================================="<<endl;
cout<< "start printing xls"<<endl; 	
cout << "==============================================="<<endl;
    printList(masterXlsList);
cout << "==============================================="<<endl;
cout<< "end printing xls"<<endl; 	
cout << "==============================================="<<endl;
    
    for(int i = 0; i < numOfLandmarkFiles; i++)
    {
       //char * landmarkfile = new char(argv[12+i]);
       readLandmarkfile(argv[13+i], i+1, landmarksVector);
    }

cout << "==============================================="<<endl;
cout<< "start printing landmarks"<<endl; 	
cout << "==============================================="<<endl;
    
    printList(landmarksVector);
    
cout << "==============================================="<<endl;
cout<< "end printing landmarks"<<endl; 	
cout << "==============================================="<<endl;

    // find the row size before beginning, as the first rows number of cols
    //int init_num_cols = (*(masterXlsList->begin()))->size();
    //int final_num_cols = (*(masterXlsList->begin()))->size();
    
    // run the whole thing twice if we need to compare against two colms; as in the case of second manual landmarks
    for(int j = 0; j < 2; j++)
    {
      int CorrdColNumtoCompare;
      if(j == 0)
	CorrdColNumtoCompare = CorrdColNum;
      else if( SecondaryCoordColNum > 0 )
	CorrdColNumtoCompare = SecondaryCoordColNum;
      else
	break;
      // loop through the landmarks and find the nearest prox
      for(vector<double*>::iterator landmark_it = landmarksVector->begin(); landmark_it != landmarksVector->end(); )
      {
	double min_dist = 0xffffff;
	list<vector<double*> *>::iterator min_iterator_pos;
	
	cout << (*landmark_it)[0]<< " "<< (*landmark_it)[1] << " "<< (*landmark_it)[2]<<" "<< (*landmark_it)[3]<<endl;
	// find the nearest  prox
	// go through the master list 
	for(list<vector<double*> *>::iterator list_it= masterXlsList->begin(); list_it != masterXlsList->end(); list_it++)
	{
	  
	    //double dist = euclidean_distance( &((*landmark_it)[0]), &((*list_it)->at(1)), 3 ,1 );
	    
	    double dist =  sqrt ( ((*landmark_it)[1] - *(*list_it)->at(CorrdColNumtoCompare-1)) * ((*landmark_it)[1] - *(*list_it)->at(CorrdColNumtoCompare-1))
				+ ((*landmark_it)[2] - *(*list_it)->at(CorrdColNumtoCompare)) * ((*landmark_it)[2] - *(*list_it)->at(CorrdColNumtoCompare))  
				+ ((*landmark_it)[3] - *(*list_it)->at(CorrdColNumtoCompare+1)) * ((*landmark_it)[3] - *(*list_it)->at(CorrdColNumtoCompare+1)) );
	    if( dist < min_dist)
	    {
	      min_dist = dist;
	      min_iterator_pos = list_it;
	    }
	    
	}
	
	if( min_dist < ProxRadius )
	{
	    cout << *(*min_iterator_pos)->at(CorrdColNumtoCompare-1) << " " << *(*min_iterator_pos)->at(CorrdColNumtoCompare) << " " << *(*min_iterator_pos)->at(CorrdColNumtoCompare+1) << " " << endl;
	    // see if this prox has been assigned something already, length would have changed
	    if((*min_iterator_pos)->size() > init_num_cols )
	    {
	      // see if we are allowed to add a new line or not
	      // for bouton and spine there shoudl not a new line added
	      // instead need to find the next closest prox and attach this landmark to that
	      if(NoNewLine)
	      {
		  //double * coords_to_be_assigned = new double[3]();
		  // here we need to first see if the current landmark is closer to the 
	 	  // already assigned prox or not 
		  // this is the dist between exsisting and the prox
		  double dist1 = sqrt ( (*(*min_iterator_pos)->at(init_num_cols+1) - *(*min_iterator_pos)->at(CorrdColNumtoCompare-1)) * (*(*min_iterator_pos)->at(init_num_cols+1) - *(*min_iterator_pos)->at(CorrdColNumtoCompare-1))
				+ (*(*min_iterator_pos)->at(init_num_cols+2) - *(*min_iterator_pos)->at(CorrdColNumtoCompare)) * (*(*min_iterator_pos)->at(init_num_cols+2) - *(*min_iterator_pos)->at(CorrdColNumtoCompare))  
				+ ((*(*min_iterator_pos)->at(init_num_cols+3) - *(*min_iterator_pos)->at(CorrdColNumtoCompare+1)) * ((*(*min_iterator_pos)->at(init_num_cols+3) - *(*min_iterator_pos)->at(CorrdColNumtoCompare+1)) )) );
		  // this is the new one
		  double dist2 = sqrt ( ((*landmark_it)[1] - *(*min_iterator_pos)->at(CorrdColNumtoCompare-1)) * ((*landmark_it)[1] - *(*min_iterator_pos)->at(CorrdColNumtoCompare-1))
				+ ((*landmark_it)[2] - *(*min_iterator_pos)->at(CorrdColNumtoCompare)) * ((*landmark_it)[2] - *(*min_iterator_pos)->at(CorrdColNumtoCompare))  
				+ ((*landmark_it)[3] - *(*min_iterator_pos)->at(CorrdColNumtoCompare+1)) * ((*landmark_it)[3] - *(*min_iterator_pos)->at(CorrdColNumtoCompare+1)) );
		  
		  if(dist1 > dist2)
		  {
		     // we have a new winner.. assign the new comer to the current position and take the already assigned guy to be assigned a new overlap coord
		      vector<double * > * copyvector = new vector<double*>;  
		      for(int i = 0; i < init_num_cols; i++)
		      {
			copyvector->push_back( (*min_iterator_pos)->at(i) );
			cout<< *(*min_iterator_pos)->at(i)<<endl; 
		      } 
		      
		      //copyvector->push_back(&(landmarkConfidence));
		      copyvector->push_back(&(*landmark_it)[0]);
		      copyvector->push_back(&(*landmark_it)[1]);
		      copyvector->push_back(&(*landmark_it)[2]);
		      copyvector->push_back(&(*landmark_it)[3]);
		      // insert new vector to the list
		      masterXlsList->insert(min_iterator_pos, copyvector);
		    
		  }
		  
		  
		  
	      }
	      else
	      {
		// lets copy the whole line and add a new line to the list below this one
		// copy the current vector
		vector<double * > * copyvector = new vector<double*>;
		//copyvector->assign(init_num_cols, 0);
		if(SecondManualLandmark)
		{
		    // here we only want to copy till prox details and then fill with zeros
		    // so that we dont duplicate the first manual landmarks
		    for(int i = 0; i < 10; i++)
		    {
		      copyvector->push_back( (*min_iterator_pos)->at(i) );
		      cout<< *(*min_iterator_pos)->at(i)<<endl; 
		    }
		    for(int i = 10; i < init_num_cols; i++)
		    {
		      double *tmp = new double;
		      *tmp = 0;
		      copyvector->push_back( tmp);
		      cout<< *(*min_iterator_pos)->at(i)<<endl; 
		    }
		    
		}
		else
		{
		    for(int i = 0; i < init_num_cols; i++)
		    {
		      copyvector->push_back( (*min_iterator_pos)->at(i) );
		      cout<< *(*min_iterator_pos)->at(i)<<endl; 
		    } 
		}
		
		//copyvector->push_back(&(landmarkConfidence));
		copyvector->push_back(&(*landmark_it)[0]);
		copyvector->push_back(&(*landmark_it)[1]);
		copyvector->push_back(&(*landmark_it)[2]);
		copyvector->push_back(&(*landmark_it)[3]);
		// insert new vector to the list
		masterXlsList->insert(min_iterator_pos, copyvector);
		
		if(final_num_cols < copyvector->size())
		{
		    final_num_cols = copyvector->size();
		}
		
	      }
	    }
	    else
	    {
	      // add this to the end of current vector
	      //(*min_iterator_pos)->push_back(&landmarkConfidence);
	      (*min_iterator_pos)->push_back(&(*landmark_it)[0]);
	      (*min_iterator_pos)->push_back(&(*landmark_it)[1]);
	      (*min_iterator_pos)->push_back(&(*landmark_it)[2]);
	      (*min_iterator_pos)->push_back(&(*landmark_it)[3]);
	      if(final_num_cols < (*min_iterator_pos)->size())
	      final_num_cols = (*min_iterator_pos)->size();
	    }
	}
	
	
	if(landmark_it != landmarksVector->end())
	{
	   landmark_it = landmarksVector->erase(landmark_it);
	}
	
      }
    }
    
    printList(masterXlsList);

    if(needZeroFilling)    
    {
      // fill the rest of the things with zeros
      for(list<vector<double*> *>::iterator list_it= masterXlsList->begin(); list_it != masterXlsList->end(); list_it++)
      {
	if((*list_it)->size() < final_num_cols )
	{
	    for(int i = init_num_cols; i < final_num_cols; i++)
	    {
		double *tmp = new double;
		*tmp = 0;
		(*list_it)->push_back(tmp);
	    }
	}
      }
    }   
    printList(masterXlsList);
    
    // write the updated xls list into a csv file into
    writeListofVectorsAsCsv( masterXlsOutFile, masterXlsList );
    
}

void printList(list<string *> * listPtr)
{
  list<string*>::iterator listit;
  for(listit = listPtr->begin(); listit!= listPtr->end(); listit++)
  {
     cout << *listit<< endl;
  }
}

void printList(vector<double *> * coordPtr)
{
  vector<double *>::iterator listit;
  for(listit = coordPtr->begin(); listit!= coordPtr->end(); listit++)
  {
     cout << (*listit)[0]<< " "<< (*listit)[1] << " "<< (*listit)[2]<<" "<< (*listit)[3]<<endl;
  }
}

void printList(list<vector<double *> *> * listPtr)
{
  list<vector<double *> *>::iterator listit;
  for(listit = listPtr->begin(); listit!= listPtr->end(); listit++)
  {
    vector<double *> *tmpvec =  *listit;
    vector<double*>::iterator vecit;
    
    for(vecit = tmpvec->begin(); vecit!= tmpvec->end(); vecit++)
    {
       cout << (*vecit)[0]<<",";
    }
    cout<<endl;
  }
}


void readLandmarkfile( char* filename, int weight, std::vector<double *>*prox_landmarks)
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
      double * tmp = new double[4];
      
      tmp[0] = weight;
      sscanf(currentLine.c_str(), "%lf %lf %lf ", &tmp[1], &tmp[2], &tmp[3]);
      ////std::cout<<" "<<tmp[X_COORD]<<" "<<tmp[Y_COORD]<<" "<<tmp[Z_COORD]<<std::endl;
      //Proximity * tmpprox = new Proximity(tmp,this->original_image);
      prox_landmarks->push_back(tmp);
      
    }
    ////std::cout<<"done reading prox landmarks"<<std::endl;
  }

}

double euclidean_distance(double * a, double * b, int dim, float z_scalar)
{
    double dist = 0;
    double diff = 0;
    
    for(int i=0; i< dim; i++)
    {
	diff = a[i]-b[i];
	if(i == 2)
	  dist += z_scalar*diff*diff;
	else 
	  dist += diff*diff;
    }
    dist = sqrt(dist);
    return dist;
}

void readCsvAsListofVectors(const char* masterXlsFile, list<vector<double*> *> * masterXlsList )
{
    ifstream ifs(masterXlsFile, ifstream::in);
    
    while(1)
    {
      string currline;
      vector<double *> *currlinevector = new vector<double*>;
      if(getline(ifs,currline,'\n'))
      {
	  stringstream linestream(currline);
	  
	  while(1)
	  {
	    string cell;
	    if(getline(linestream, cell, ','))
	    {
	       double * currCellVal = new double; 
	       *currCellVal = atof(cell.c_str());
	       currlinevector->push_back(currCellVal);
	    }
	    else
	    {
	      break;
	    }
	  }
	  
      }
      else
      {
	break;
      }
       //cout<< currline <<endl;
      masterXlsList->push_back(currlinevector);
    }
}


void writeListofVectorsAsCsv(const char* masterXlsFile, list<vector<double*> *> * listPtr )
{
    ofstream ofs(masterXlsFile, ofstream::out);
      
    list<vector<double *> *>::iterator listit;
    for(listit = listPtr->begin(); listit!= listPtr->end(); listit++)
    {
      vector<double *> *tmpvec =  *listit;
      vector<double*>::iterator vecit;
      
      for(vecit = tmpvec->begin(); vecit!= tmpvec->end(); vecit++)
      {
	ofs << (*vecit)[0]<<",";
      }
      ofs<<endl;
    }

}

