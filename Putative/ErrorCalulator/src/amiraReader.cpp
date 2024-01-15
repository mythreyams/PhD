/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      amiraReader.cpp                                               */
/*                                                                          */
/* Purpose:   interface class between Amira SpatialGraph, Surface and       */
/*            other files and the internal spatial graph data structure     */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   12.5.2013                                                 */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "typedefs.h"
#include "basics.h"
#include "amiraReader.h"


Reader::Reader(const char * filename)
{
        inputSpatialGraph = NULL;
        this->inputFilename = filename;
        letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        numbers = "0123456789";
        signs = "+-";
        otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
        whitespace = "\t ";
        //initializeConstants();
};

Reader::Reader(const char * filename, const char * outputFilename)
{
        inputSpatialGraph = NULL;
        this->inputFilename = filename;
        this->outputFilename = outputFilename;
        letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        numbers = "0123456789";
        signs = "+-";
        otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
        whitespace = "\t ";
        //initializeConstants();
};

Reader::~Reader()
{
        //tbd
};

void Reader::readSpatialGraphFile(bool applyTransform)
{
        std::ifstream inputStream(inputFilename);
        
        if(!inputStream.fail())
        {
//              //////std::cout << "Reading SpatialGraph file " << inputFilename << std::endl;
                std::string currentLine;
                unsigned int line = 0;
                static int dendrtite_count = 0;
                static int axon_count = 0;
                static int reading_Axon_ID = false;
                static int reading_dendrite_ID = false;
                int axon_index = 0;
                int dend_index = 0;
                bool parameters = 0;
                bool transform = 0;
                unsigned int brackets = 0;
                unsigned int vertexTransformIndex = 1000000, edgeTransformIndex = 1000000;
                unsigned int vertexCoordIndex = 1000000, vertexLabelIndex = 1000000, edgeConnectivityIndex = 1000000, edgePointIndex = 1000000, edgeLabelIndex = 1000000, edgePointCoordIndex = 1000000, edgeRadiusIndex = 1000000;
                unsigned int currentIndex = 0;
                unsigned int vertex = 0, edge = 0, point = 0;
                
                std::vector< Vertex * > inputVertices;
                std::vector< Edge * > inputEdges;
                std::list< double * > tmpVertices;
                std::list< int > tmpVertexLabels;
                std::list< int * > tmpEdgeConnections;
                std::list< int > tmpNoEdgePoints;
                std::list< int > tmpEdgeLabels;
                std::list< double * > edgePoints;
                std::vector<int> tmpAxonIDs;
                std::vector<int> tmpDendriteIDs;
                
                
                double ** transformation = new double *[4];
                for(int ii = 0; ii < 4; ++ii)
                {
                        transformation[ii] = new double[4];
                        for(int jj = 0; jj < 4; ++jj)
                        {
                                if(ii != jj)
                                        transformation[ii][jj] = 0;
                                else
                                        transformation[ii][jj] = 1;
                        }
                }
                
                
                
                
                while(!std::getline(inputStream, currentLine).eof() /*&& line < 100*/)
                {
//                      if(!parameters)
//                      {
//                              //////std::cout << currentLine << std::endl;
//                              ++line;
//                      }
                        //////std::cout<<currentLine<<std::endl;
                        if(currentLine.size())
                        {
                                if(currentLine.find("@", 0) == 0)
                                {
                                        char * tmp = new char[currentLine.size() - 1];
                                        currentLine.copy(tmp, currentLine.size() - 1, 1);
                                        currentIndex = atoi(tmp);
//                                      //////std::cout << "Reading data section " << currentIndex << std::endl;
                                        delete [] tmp;
                                        continue;
                                }
                                
                                if(currentIndex == 0)
                                {
                                        std::string::size_type loc = currentLine.find("define", 0);
                                        if(loc != std::string::npos)
                                        {
                                                if(currentLine.find("VERTEX", 7) != std::string::npos)
                                                {
                                                        char * tmp = new char[currentLine.size() - 14];
                                                        currentLine.copy(tmp, currentLine.size() - 14, 14);
                                                        vertex = atoi(tmp);
//                                                      //////std::cout << "vertex = " << vertex << std::endl;
//                                                      inputVertices.resize(vertex);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("EDGE", 7) != std::string::npos)
                                                {
                                                        char * tmp = new char[currentLine.size() - 12];
                                                        currentLine.copy(tmp, currentLine.size() - 12, 12);
                                                        edge = atoi(tmp);
//                                                      //////std::cout << "edges = " << edge << std::endl;
//                                                      inputEdges.resize(edge);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("POINT", 7) != std::string::npos)
                                                {
                                                        char * tmp = new char[currentLine.size() - 13];
                                                        currentLine.copy(tmp, currentLine.size() - 13, 13);
                                                        point = atoi(tmp);
//                                                      //////std::cout << "points = " << point << std::endl;
                                                        delete [] tmp;
                                                        continue;
                                                }
                                        }
                                        
                                        loc = currentLine.find("Parameters", 0);
                                        if(loc == 0)
                                        {
                                                parameters = 1;
                                                if(currentLine.find("{", 0) != std::string::npos)
                                                    {
                                                        //////std::cout<<"parameter true"<<std::endl;
                                                        brackets = 1;
                                                    }
                                                continue;
                                        }
//                                      if(parameters && currentLine.find("{", 0) == std::string::npos && currentLine.find("}", 0) == std::string::npos)
//                                              continue;
                                        if(parameters)
                                        {
                                            //////std::cout << currentLine << std::endl;
                                            
                                                std::string::size_type startPos = 0;
                                                for(std::string::size_type bPos = currentLine.find("{", startPos); bPos != std::string::npos; )
                                                {
                                                        ++brackets;
                                                        if(bPos == currentLine.size() - 1)
                                                                break;
                                                        bPos = currentLine.find("{", bPos+1);
                                                }
                                                for(std::string::size_type bPos = currentLine.find("}", startPos); bPos != std::string::npos; )
                                                {
                                                        --brackets;
                                                        if(bPos == currentLine.size() - 1)
                                                                break;
                                                        bPos = currentLine.find("}", bPos+1);
                                                }
                                                if(!brackets)
                                                {
                                                    //////std::cout<<"parameter true"<<std::endl;
                                                        parameters = 0;
                                                }
                                        }
//                                      if(parameters && currentLine.find("{", 0) != std::string::npos)
//                                      {
//                                              ++brackets;
//                                              continue;
//                                      }
//                                      if(parameters && currentLine.find("}", 0) != std::string::npos)
//                                      {
//                                              --brackets;
//                                              if(!brackets)
//                                                      parameters = 0;
//                                              continue;
//                                      }
                                        
                                        if(parameters && currentLine.find("TransformationMatrix ", 0) != std::string::npos)
                                        {
//                                              //////std::cout << "found correct section transform parameters!" << std::endl;
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
                                                        char * tmp1 = new char[20];
                                                
                                                        for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                        {
                                                            tmp1[i] = 'f';
                                                        }
                                                        //char * tmp1 = new char[loc2 - loc1];
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
//                                              //////std::cout << "transformation matrix:" << std::endl;
//                                              for(int ii = 0; ii < 4; ++ii)
//                                              {
//                                                      //////std::cout << "[";
//                                                      for(int jj = 0; jj < 4; ++jj)
//                                                      {
//                                                              if(jj < 3)
//                                                                      //////std::cout << transformation[ii][jj] << ",\t";
//                                                              else
//                                                                      //////std::cout << transformation[ii][jj];
//                                                      }
//                                                      //////std::cout << "]" << std::endl;
//                                              }
                                                //remove numeric artifacts from z-axis:
                                                for(int ii = 0; ii < 2; ++ii)
                                                {
                                                        transformation[2][ii] = 0;
                                                        transformation[ii][2] = 0;
                                                }
                                                transformation[2][2] = 1;
                                        }
                                        
                                        if(parameters && currentLine.find("Dendrite_", 0) != std::string::npos)
                                        {
                                            //////std::cout<<"reading dendrite true"<<std::endl;
                                            //////std::cout << currentLine << std::endl;
                                            reading_dendrite_ID = true;
                                            std::string::size_type loc1, loc2;
                                            loc1 = currentLine.find_first_of(numbers, 0);
                                            loc2 = currentLine.find_last_of(numbers, currentLine.size());
                                            char * tmp1 = new char[20];
                                                
                                                        for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                        {
                                                            tmp1[i] = 'f';
                                                        }
                                            //char * tmp1 = new char[loc2 - loc1+1];
                                            currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                            //currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                                //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                //////std::cout<<"DendriteID "<<tmp1<<std::endl;
                                                
                                            dend_index = atoi(tmp1);
                                             ////std::cout<<"dend_index "<<dend_index<<std::endl;
                                            //////std::cout<<"DendriteID int "<<ftmp1<<std::endl;
                                            //tmpDendriteIDs.push_back(ftmp1); 
                                            
                                        }
                                        if(parameters && reading_dendrite_ID)
                                        {
                                            if(currentLine.find("Id", 0) != std::string::npos)
                                            {
                                                //////std::cout<<"found id dendrite"<<std::endl;
                                                //////std::cout << currentLine << std::endl;
                                                unsigned int count = 0;
                                                std::string::size_type loc1, loc2, loc3;
                                                loc1 = currentLine.find_first_of(numbers, 0);
                                                loc2 = currentLine.find_last_of(numbers, currentLine.size());
//                                                 if(loc2 != std::string::npos)
//                                                         if(loc2 < loc1)
//                                                                 loc1 = loc2;
//                                                 loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                                //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                         //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                //while(loc2 != std::string::npos && count < 16)
                                                {
                                                        char * tmp1 = new char[20];
                                                
                                                        for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                        {
                                                            tmp1[i] = 'f';
                                                        }
                                                        //char * tmp1 = new char[loc2 - loc1+1];
                                                        currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                                        //currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                                         //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                         //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                         //////std::cout<<"DendriteID "<<tmp1<<std::endl;
                                                         
                                                        int ftmp1 = atoi(tmp1);
                                                        ////std::cout<<"DendriteID int "<<ftmp1<<std::endl;
                                                        tmpDendriteIDs.insert(tmpDendriteIDs.begin()+dend_index-1, ftmp1);         // amira files are columns after each other
                                                       ////std::cout<<"DendriteID int "<<ftmp1<<std::endl;
                                                        delete [] tmp1;
                                                }
                                                
                                                reading_dendrite_ID = false;
                                            }
                                        }
                                        
                                        if(parameters && currentLine.find("Axon_", 0) != std::string::npos)
                                        {
                                            //////std::cout<<"reading axon true"<<std::endl;
                                            //////std::cout << currentLine << std::endl;
                                            reading_Axon_ID = true;
                                            
                                            std::string::size_type loc1, loc2;
                                            loc1 = currentLine.find_first_of(numbers, 0);
                                            loc2 = currentLine.find_last_of(numbers, currentLine.size());
                                            char * tmp1 = new char[20];
                                                
                                            for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                            {
                                                tmp1[i] = 'f';
                                            }
                                            //char * tmp1 = new char[loc2 - loc1+1];
                                            currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                            //currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                                //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                //////std::cout<<"DendriteID "<<tmp1<<std::endl;
                                                
                                            axon_index = atoi(tmp1);
                                            ////std::cout<<"axon_index "<<axon_index<<std::endl;
                                            
                                            
                                        }
                                        if(parameters && reading_Axon_ID)
                                        {
                                            if(currentLine.find("Id", 0) != std::string::npos)
                                            {
                                                //////std::cout<<"found id Axon"<<std::endl;
                                                //////std::cout << currentLine << std::endl;
                                                unsigned int count = 0;
                                                std::string::size_type loc1, loc2, loc3;
                                                loc1 = currentLine.find_first_of(numbers, 0);
                                                loc2 = currentLine.find_last_of(numbers, currentLine.size());
//                                                 if(loc2 != std::string::npos)
//                                                         if(loc2 < loc1)
//                                                                 loc1 = loc2;
//                                                 loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                                                 //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                         //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                          //////std::cout<<"currentLine.size() "<<currentLine.size()<<std::endl;
                                                        
                                                //while(loc2 != std::string::npos && count < 16)
                                                {       
                                                        char * tmp1 = new char[20];
                                                
                                                        for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                                                        {
                                                            tmp1[i] = 'f';
                                                        }
                                                        //char * tmp1 = new char[loc2 - loc1+1];
                                                        currentLine.copy(tmp1, loc2 - loc1+1, loc1);
                                                         //////std::cout<<"loc1 "<<loc1<<std::endl;
                                                         //////std::cout<<"loc2 "<<loc2<<std::endl;
                                                         //////std::cout<<"axonID "<<tmp1<<std::endl;
                                                        int ftmp1 = atoi(tmp1);
                                                        ////std::cout<<"axonID int "<<ftmp1<<std::endl;
                                                        tmpAxonIDs.insert(tmpAxonIDs.begin()+axon_index-1, ftmp1); 
                                                        ////std::cout<<"axonID int "<<ftmp1<<std::endl;
                                                        //transformation[count%4][count/4]= ftmp1;        // amira files are columns after each other
//                                                         loc3 = loc2;
//                                                         loc1 = currentLine.find_first_of(numbers, loc3);
//                                                         loc2 = currentLine.find_first_of(signs, loc3);
//                                                         if(loc2 != std::string::npos)
//                                                                 if(loc2 < loc1)
//                                                                         loc1 = loc2;
//                                                         loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
//                                                         ++count;
//                                                         delete [] tmp1;
                                                }
                                                
                                                reading_Axon_ID = false;
                                            }
                                        }
                                        
                                        loc = currentLine.find("VERTEX", 0);
                                        if(loc == 0)
                                        {
                                                if(currentLine.find("VertexCoordinates", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        vertexCoordIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("VertexLabels", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        vertexLabelIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("TransformInfo", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        vertexTransformIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                        }
                                        
                                        loc = currentLine.find("EDGE", 0);
                                        if(loc == 0)
                                        {
                                                if(currentLine.find("EdgeConnectivity", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgeConnectivityIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("NumEdgePoints", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgePointIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("EdgeLabels", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgeLabelIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("TransformInfo", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgeTransformIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                        }
                                        
                                        loc = currentLine.find("POINT", 0);
                                        if(loc == 0)
                                        {
                                                if(currentLine.find("EdgePointCoordinates", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgePointCoordIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                                if(currentLine.find("Radius", 0) != std::string::npos)
                                                {
                                                        loc = currentLine.find("@", 0);
                                                        char * tmp = new char[currentLine.size() - loc - 1];
                                                        currentLine.copy(tmp, currentLine.size() - loc - 1, loc + 1);
                                                        edgeRadiusIndex = atoi(tmp);
                                                        delete [] tmp;
                                                        continue;
                                                }
                                        }
                                }
                                
                                if(currentIndex == vertexTransformIndex || currentIndex == edgeTransformIndex || currentIndex == edgeRadiusIndex)
                                        continue;
                                
                                if(currentIndex == vertexCoordIndex)
                                {
                                        const char * thisLine = currentLine.c_str();
                                        double * tmpCoords = new double[ARRAY_LENGTH];
                                        char ** endptr = new char*;
                                        tmpCoords[0] = strtod(thisLine, endptr);
                                        tmpCoords[1] = strtod(*endptr, endptr);
                                        tmpCoords[2] = strtod(*endptr, endptr);
                                        
//                                      std::string::size_type loc1 = currentLine.find(" ", 0);
//                                      std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
//                                      char * tmp1 = new char[loc1];
//                                      char * tmp2 = new char[loc2 - loc1];
//                                      char * tmp3 = new char[currentLine.size() - loc2 - 1];
//                                      currentLine.copy(tmp1, loc1, 0);
//                                      currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
//                                      currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
//                                      double * tmpCoords = new double[3];
// //                                   tmpCoords[0] = atof(tmp1);
// //                                   tmpCoords[1] = atof(tmp2);
// //                                   tmpCoords[2] = atof(tmp3);
//                                      char ** endptr = new char*;
//                                      tmpCoords[0] = strtod(tmp1, endptr);
//                                      tmpCoords[1] = strtod(tmp2, endptr);
//                                      tmpCoords[2] = strtod(tmp3, endptr);
                                        
                                        tmpVertices.push_back(tmpCoords);
//                                      delete [] tmp1, delete [] tmp2, delete [] tmp3;
                                        delete endptr;
                                }
                                
                                if(currentIndex == vertexLabelIndex)
                                {
                                        char * tmp = new char[currentLine.size()];
                                        currentLine.copy(tmp, currentLine.size(), 0);
                                        int tmplabel = atoi(tmp);
                                        
                                        tmpVertexLabels.push_back(tmplabel);
                                        delete [] tmp;
                                }
                                
                                if(currentIndex == edgeConnectivityIndex)
                                {
                                        const char * thisLine = currentLine.c_str();
                                        int * tmpCoords = new int[2];
                                        char ** endptr = new char*;
                                        tmpCoords[0] = static_cast< int >(strtol(thisLine, endptr, 10));
                                        tmpCoords[1] = static_cast< int >(strtol(*endptr, endptr, 10));
                                        
//                                      std::string::size_type loc = currentLine.find(" ", 0);
//                                      char * tmp1 = new char[loc];
//                                      char * tmp2 = new char[currentLine.size() - loc - 1];
//                                      currentLine.copy(tmp1, loc, 0);
//                                      currentLine.copy(tmp2, currentLine.size() - loc - 1, loc + 1);
//                                      int * tmpCoords = new int[2];
//                                      tmpCoords[0] = atoi(tmp1);
//                                      tmpCoords[1] = atoi(tmp2);
                                        
                                        tmpEdgeConnections.push_back(tmpCoords);
//                                      delete [] tmp1, delete [] tmp2;
                                        delete endptr;
                                }
                                
                                if(currentIndex == edgePointIndex)
                                {
                                        char * tmp = new char[currentLine.size()];
                                        currentLine.copy(tmp, currentLine.size(), 0);
                                        int tmplabel = atoi(tmp);
                                        
                                        tmpNoEdgePoints.push_back(tmplabel);
                                        delete [] tmp;
                                }
                                
                                if(currentIndex == edgeLabelIndex)
                                {
                                        char * tmp = new char[currentLine.size()];
                                        currentLine.copy(tmp, currentLine.size(), 0);
                                        int tmplabel = atoi(tmp);
                                        
                                        tmpEdgeLabels.push_back(tmplabel);
                                        delete [] tmp;
                                }
                                
                                if(currentIndex == edgePointCoordIndex)
                                {
                                        const char * thisLine = currentLine.c_str();
                                        double * tmpCoords = new double[ARRAY_LENGTH];
                                        char ** endptr = new char*;
                                        tmpCoords[0] = strtod(thisLine, endptr);
                                        tmpCoords[1] = strtod(*endptr, endptr);
                                        tmpCoords[2] = strtod(*endptr, endptr);
                                        
//                                      std::string::size_type loc1 = currentLine.find(" ", 0);
//                                      std::string::size_type loc2 = currentLine.find(" ", loc1 + 1);
//                                      char * tmp1 = new char[loc1];
//                                      char * tmp2 = new char[loc2 - loc1];
//                                      char * tmp3 = new char[currentLine.size() - loc2 - 1];
//                                      currentLine.copy(tmp1, loc1, 0);
//                                      currentLine.copy(tmp2, loc2 - loc1, loc1 + 1);
//                                      currentLine.copy(tmp3, currentLine.size() - loc2 - 1, loc2 + 1);
//                                      double * tmpCoords = new double[3];
// //                                   tmpCoords[0] = atof(tmp1);
// //                                   tmpCoords[1] = atof(tmp2);
// //                                   tmpCoords[2] = atof(tmp3);
//                                      char ** endptr = new char*;
//                                      tmpCoords[0] = strtod(tmp1, endptr);
//                                      tmpCoords[1] = strtod(tmp2, endptr);
//                                      tmpCoords[2] = strtod(tmp3, endptr);
                                        
                                        edgePoints.push_back(tmpCoords);
//                                      delete [] tmp1, delete [] tmp2, delete [] tmp3;
                                        delete endptr;
                                }
                        }
                }
                
                std::list< double * >::iterator tmpvertexit;
                std::list< int >::iterator tmpvertexlabelit;
                std::list< int * >::iterator tmpedgeconnectivityit;
                std::list< int >::iterator tmpnumberedgepointsit;
                std::list< int >::iterator tmpedgelabelit;
                std::list< double * >::iterator edgepointit;
                
                
                edgepointit = edgePoints.begin();
                for(tmpedgeconnectivityit = tmpEdgeConnections.begin(), tmpnumberedgepointsit = tmpNoEdgePoints.begin(), tmpedgelabelit = tmpEdgeLabels.begin();
                tmpedgeconnectivityit != tmpEdgeConnections.end() && tmpnumberedgepointsit != tmpNoEdgePoints.end() && tmpedgelabelit != tmpEdgeLabels.end();
                ++tmpedgeconnectivityit, ++tmpnumberedgepointsit, ++tmpedgelabelit)
                {
                        std::list< double * >::iterator tmpit = edgepointit;
                        for(int ii = 0; ii < *tmpnumberedgepointsit; ++ii)
                                ++tmpit;
                        std::list< double * > tmpPoints(edgepointit, tmpit);
                        Edge * tmpEdge = new Edge(*tmpedgeconnectivityit, *tmpnumberedgepointsit, *tmpedgelabelit, tmpPoints);
                        inputEdges.push_back(tmpEdge);
                        edgepointit = tmpit;
                }
                
                inputSpatialGraph = new AmiraSpatialGraph;
                std::vector< Vertex * >::iterator vertexIter;
                std::vector< Edge * >::iterator edgeIter;
                std::vector< int  >::iterator dendriteIDIter;
                std::vector< int >::iterator axonIDIter;
                
                for(axonIDIter = tmpAxonIDs.begin(); axonIDIter != tmpAxonIDs.end(); ++axonIDIter)
                        inputSpatialGraph->addAxonID(*axonIDIter);
                
                for(dendriteIDIter = tmpDendriteIDs.begin(); dendriteIDIter != tmpDendriteIDs.end(); ++dendriteIDIter)
                        inputSpatialGraph->addDendriteID(*dendriteIDIter);
                
                for(edgeIter = inputEdges.begin(); edgeIter != inputEdges.end(); ++edgeIter)
                        inputSpatialGraph->addEdge(*edgeIter);
                
                
                
                for(tmpvertexit = tmpVertices.begin(), tmpvertexlabelit = tmpVertexLabels.begin(); tmpvertexit != tmpVertices.end() && tmpvertexlabelit != tmpVertexLabels.end(); ++tmpvertexit, ++tmpvertexlabelit)
                {
                        Vertex * tmpVertex = new Vertex(*tmpvertexit, *tmpvertexlabelit);
                        inputVertices.push_back(tmpVertex);
                }
                
                for(vertexIter = inputVertices.begin(); vertexIter != inputVertices.end(); ++ vertexIter)
                        inputSpatialGraph->addVertex(*vertexIter);
                
                
                std::vector< Edge * > * all_edges = inputSpatialGraph->edgesPointer();
                std::vector< Vertex * > * all_vertices = inputSpatialGraph->verticesPointer();
                unsigned int totalNumOfEdges = all_edges->size();       
        
                
                for(long pos = totalNumOfEdges -1; pos >= 0; pos--)     //for each edge in list
                {                       
                        Edge * currentEdge = all_edges->at(pos);
                        int from_vertex = currentEdge->edgeConnectivity[0];
                        int to_vertex = currentEdge->edgeConnectivity[1];
                        
                        all_vertices->at(from_vertex)->connectedEdges->push_back(currentEdge);
                        all_vertices->at(from_vertex)->edgeDirections->push_back(FORWARD);
                        
                        all_vertices->at(to_vertex)->connectedEdges->push_back(currentEdge);
                        all_vertices->at(to_vertex)->edgeDirections->push_back(BACKWARD);
                }
                
                
//              bool isId = 1;
//              for(int ii = 0; ii < 4; ++ii)
//                      for(int jj = 0; jj < 4; ++jj)
//                      {
//                              if(ii != jj)
//                                      if(transformation[ii][jj] != 0)
//                                              isId = 0;
//                              
//                              else
//                                      if(transformation[ii][jj] != 1)
//                                              isId = 0;
//                      }
//              //////std::cout << "transformation matrix:" << std::endl;
//              for(int ii = 0; ii < 4; ++ii)
//              {
//                      //////std::cout << "[";
//                      for(int jj = 0; jj < 4; ++jj)
//                      {
//                              if(jj < 3)
//                                      //////std::cout << transformation[ii][jj] << ",\t";
//                              else
//                                      //////std::cout << transformation[ii][jj];
//                      }
//                      //////std::cout << "]" << std::endl;
//              }
                if(applyTransform)
                {
//                      inputSpatialGraph->printTransformation();
                        //inputSpatialGraph->setTransformation(transformation);
                        //inputSpatialGraph->applyTransformation();
                }
                for(int ii = 0; ii < 4; ++ii)
                        delete [] transformation[ii];
                delete [] transformation;
                
//              //////std::cout << "Vertex number = " << inputSpatialGraph->getNumberOfVertices() << std::endl;
//              //////std::cout << "Edge number = " << inputSpatialGraph->getNumberOfEdges() << std::endl;
//              //////std::cout << "Point number = " << inputSpatialGraph->getNumberOfPoints() << std::endl;
                
//              //////std::cout << "VertexCoordinates @" << vertexCoordIndex << std::endl;
//              //////std::cout << "Vertex GraphLabels @" << vertexLabelIndex << std::endl;
//              //////std::cout << "Vertex TransformInfo @" << vertexTransformIndex << std::endl;
//              //////std::cout << "EdgeConnectivity @" << edgeConnectivityIndex << std::endl;
//              //////std::cout << "NumEdgePoints @" << edgePointIndex << std::endl;
//              //////std::cout << "Edge GraphLabels @" << edgeLabelIndex << std::endl;
//              //////std::cout << "Edge TransformInfo @" << edgeTransformIndex << std::endl;
//              //////std::cout << "EdgePointCoordinates @" << edgePointCoordIndex << std::endl;
//              //////std::cout << "Edge Radius @" << edgeRadiusIndex << std::endl;
                
        
                std::vector< Edge * > * edges = inputSpatialGraph->edgesPointer();
                unsigned int numOfEdges = edges->size();                
                double x_diff, y_diff, z_diff;
                double norm=0;
                bool firstIteration;

                if(numOfEdges == 0)
                {
                        //////std::cout<< "CAUTION!!! Edge List empty!" << std::endl;
                }
                else
                {
                        for(long pos = numOfEdges -1; pos >= 0; pos--)  //for each edge in list
                        {                       
                                Edge * currentEdge = edges->at(pos);
                                std::list< double * >::iterator edge_it;
                                std::list< double * >::iterator prev_it;

                                //for every point along edge
                                for(edge_it = currentEdge->edgePointCoordinates.begin(), firstIteration=true;edge_it != currentEdge->edgePointCoordinates.end(); edge_it++) 
                                {
                                        double * coords = *edge_it;
                                                                        
                                        if(firstIteration)
                                        {
                                            coords[LENGTH_UNTIL_HERE] = 0;
                                            firstIteration = false;
                                        }  
                                        else
                                        {
                                            double * prev_coords = *prev_it;
                                          
                                            x_diff = fabs(coords[X_COORD] - prev_coords[X_COORD]);
                                            y_diff = fabs(coords[Y_COORD] - prev_coords[Y_COORD]);
                                            z_diff = fabs(coords[Z_COORD] - prev_coords[Z_COORD]);
                                            
                                            norm = sqrt( (x_diff*x_diff) + (y_diff*y_diff) + (z_diff*z_diff) );
                                            
                                            currentEdge->physicalLength += norm; //increase length
                                            coords[LENGTH_UNTIL_HERE] = currentEdge->physicalLength;
                                        }
                                        
                                        prev_it = edge_it;
                                }
                                
                        }
                }       

                
                
        }
        
        inputStream.close();
};

void Reader::writeSpatialGraphFile()
{
//      std::list<std::list<Compartment * > >::iterator edge_list_it;
//      std::list<Compartment * >::iterator edge_it;
        std::vector< Vertex * >::iterator vertexIt;
        std::vector< Edge * >::iterator edgeIt;
        
        int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
        
//      for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
//      {
//              number_of_edge_points += (*edge_list_it).size();
//      }
//      for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
//      for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
        
        std::string format = outputFilename;
        format += ".am";
        
        #ifdef DEBUG
        //////std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
        ////////std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
        #endif
        
        std::ofstream NeuroMorphData( format.c_str() );
        
        NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
        NeuroMorphData << "# This SpatialGraph file was created by the automated Bouton and Spine detector BoutonFinder " << std::endl;
        NeuroMorphData << "# BoutonFinder was programmed by Christopher Tull" << std::endl;
        NeuroMorphData << "# Max-Planck-Institute for Biological Cybernetics, Tuebingen, Germany " << std::endl;
        
        NeuroMorphData << "define VERTEX " << inputSpatialGraph->getNumberOfVertices() << std::endl;
        NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
//      NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
        NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
        
        NeuroMorphData << "Parameters {GraphLabels {"                           <<std::endl;
        NeuroMorphData << "        Neuron { "                                   <<std::endl;
        NeuroMorphData << "            Dendrite {"                              <<std::endl;
        NeuroMorphData << "                ApicalDendrite {"                    <<std::endl;
        NeuroMorphData << "                    Color 1 0.5 0.5,"                <<std::endl;
        NeuroMorphData << "                    Id 4 }"                          <<std::endl;
        NeuroMorphData << "                BasalDendrite {"                     <<std::endl;
        NeuroMorphData << "                    Color 0.8 0.4 0.4,"              <<std::endl;
        NeuroMorphData << "                    Id 5 }"                          <<std::endl;
        NeuroMorphData << "                Color 1 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 3 }"                              <<std::endl;
        NeuroMorphData << "            Axon {"                                  <<std::endl;
        NeuroMorphData << "                Color 0 0 1,"                        <<std::endl;
        NeuroMorphData << "                Id 6 }"                              <<std::endl;
        NeuroMorphData << "            Soma {"                                  <<std::endl;
        NeuroMorphData << "                Color 1 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 7 }"                              <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2 }"                                  <<std::endl;
        NeuroMorphData << "        Landmark {"                                  <<std::endl;
        NeuroMorphData << "            Pia {"                                   <<std::endl;
        NeuroMorphData << "                Color 0 1 0.5,"                      <<std::endl;
        NeuroMorphData << "                Id 9 }"                              <<std::endl;
        NeuroMorphData << "            Vessel {"                                <<std::endl;
        NeuroMorphData << "                Color 1 0.5 0,"                      <<std::endl;
        NeuroMorphData << "                Id 10 }"                             <<std::endl;
        NeuroMorphData << "            Barrel {"                                <<std::endl;
        NeuroMorphData << "                aRow {"                              <<std::endl;
        NeuroMorphData << "                    A1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 13 }"                     <<std::endl;
        NeuroMorphData << "                    A2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 14 }"                     <<std::endl;
        NeuroMorphData << "                    A3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 15 }"                     <<std::endl;
        NeuroMorphData << "                    A4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 16 }"                     <<std::endl;
        NeuroMorphData << "                Color 1 0.2 0.2,"                    <<std::endl;
        NeuroMorphData << "                Id 12 }"                             <<std::endl;
        NeuroMorphData << "                bRow {"                              <<std::endl;
        NeuroMorphData << "                    B1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 18 }"                     <<std::endl;
        NeuroMorphData << "                    B2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 19 }"                     <<std::endl;
        NeuroMorphData << "                    B3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 20 }"                     <<std::endl;
        NeuroMorphData << "                    B4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 21 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.25 0.25,"              <<std::endl;
        NeuroMorphData << "                    Id 17 }"                         <<std::endl;
        NeuroMorphData << "                cRow {"                              <<std::endl;
        NeuroMorphData << "                    C1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 23 }"                     <<std::endl;
        NeuroMorphData << "                    C2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 24 }"                     <<std::endl;
        NeuroMorphData << "                    C3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 25 }"                     <<std::endl;
        NeuroMorphData << "                    C4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 26 }"                     <<std::endl;
        NeuroMorphData << "                    C5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 27 }"                     <<std::endl;
        NeuroMorphData << "                    C6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 28 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.3 0.3,"                <<std::endl;
        NeuroMorphData << "                    Id 22 }"                         <<std::endl;
        NeuroMorphData << "                dRow {"                              <<std::endl;
        NeuroMorphData << "                    D1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 30 }"                     <<std::endl;
        NeuroMorphData << "                    D2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 31 }"                     <<std::endl;
        NeuroMorphData << "                    D3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 32 }"                     <<std::endl;
        NeuroMorphData << "                    D4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 33 }"                     <<std::endl;
        NeuroMorphData << "                    D5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 34 }"                     <<std::endl;
        NeuroMorphData << "                    D6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 35 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.35 0.35,"              <<std::endl;
        NeuroMorphData << "                    Id 29 }"                         <<std::endl;
        NeuroMorphData << "                eRow {"                              <<std::endl;
        NeuroMorphData << "                    E1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 37 }"                     <<std::endl;
        NeuroMorphData << "                    E2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 38 }"                     <<std::endl;
        NeuroMorphData << "                    E3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 39 }"                     <<std::endl;
        NeuroMorphData << "                    E4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 40 }"                     <<std::endl;
        NeuroMorphData << "                    E5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 41 }"                     <<std::endl;
        NeuroMorphData << "                    E6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 42 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.4 0.4,"                <<std::endl;
        NeuroMorphData << "                    Id 36 }"                         <<std::endl;
        NeuroMorphData << "                greekRow {"                          <<std::endl;
        NeuroMorphData << "                    Alpha {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 44 }"                     <<std::endl;
        NeuroMorphData << "                    Beta {"                          <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 45 }"                     <<std::endl;
        NeuroMorphData << "                    Gamma {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 46 }"                     <<std::endl;
        NeuroMorphData << "                    Delta {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 47 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.1 0.1,"                <<std::endl;
        NeuroMorphData << "                    Id 43 }"                         <<std::endl;
        NeuroMorphData << "                Color 0 1 0,"                        <<std::endl;
        NeuroMorphData << "                Id 11 }"                             <<std::endl;
        NeuroMorphData << "            WhiteMatter {"                           <<std::endl;
        NeuroMorphData << "                Color 0.5 1 0.75,"                   <<std::endl;
        NeuroMorphData << "                Id 48 }"                             <<std::endl;
        NeuroMorphData << "            OtherBarrels {"                          <<std::endl;
        NeuroMorphData << "                Color 1 0 1,"                        <<std::endl;
        NeuroMorphData << "                Id 49 }"                             <<std::endl;
        NeuroMorphData << "            ZAxis {"                                 <<std::endl;
        NeuroMorphData << "                Color 0 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 50 }"                             <<std::endl;
        NeuroMorphData << "            Color 0 1 1,"                            <<std::endl;
        NeuroMorphData << "            Id 8 }"                                  <<std::endl;
        NeuroMorphData << "        Id 0,"                                       <<std::endl;
        NeuroMorphData << "        Color 0 0 0 }"                               <<std::endl;
        NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                    <<std::endl;
  
        
        NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 "           << std::endl;
        NeuroMorphData << "VERTEX {int GraphLabels } @2 "                       << std::endl;
        
        NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 "                << std::endl;
        NeuroMorphData << "EDGE { int NumEdgePoints } @4 "                      << std::endl;
        NeuroMorphData << "EDGE { int GraphLabels } @5 "                        << std::endl;
        
        NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 "         << std::endl;
        NeuroMorphData << "POINT { float Radius } @7 "                          << std::endl;
        
        if(inputSpatialGraph->getNumberOfVertices())
        {
                NeuroMorphData << "\n@1 # Vertices xyz coordinates"                     << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                        NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
                
                NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                {
                        NeuroMorphData << (*vertexIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfEdges())
        {
                NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
                int last_index = 0;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
                }
                
                NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
                }
                
                NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfPoints())
        {
                NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
                        }
                        
                }
                
                NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[SURFACE] << std::endl;
                        }
                        
                }
        }
        
        NeuroMorphData.close();
};


void Reader::writeSpatialGraphFile2()
{
//      std::list<std::list<Compartment * > >::iterator edge_list_it;
//      std::list<Compartment * >::iterator edge_it;
        std::vector< Vertex * >::iterator vertexIt;
        std::vector< Edge * >::iterator edgeIt;
        
        int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
        
//      for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
//      {
//              number_of_edge_points += (*edge_list_it).size();
//      }
//      for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
//      for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
        
        std::string format = outputFilename;
        format += ".am";
        
        #ifdef DEBUG
        //////std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
        ////////std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
        #endif
        
        std::ofstream NeuroMorphData( format.c_str() );
        
        NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
        NeuroMorphData << "# This SpatialGraph file was created by the automated Bouton and Spine detector BoutonFinder " << std::endl;
        NeuroMorphData << "# BoutonFinder was programmed by Christopher Tull" << std::endl;
        NeuroMorphData << "# Max-Planck-Institute for Biological Cybernetics, Tuebingen, Germany " << std::endl;
        
        NeuroMorphData << "define VERTEX " << inputSpatialGraph->getNumberOfVertices() << std::endl;
        NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
//      NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
        NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
        
        NeuroMorphData << "Parameters {"                                        <<std::endl;
        NeuroMorphData << "    VertexLabels {"                          <<std::endl;
        NeuroMorphData << "        LowEndingVertexLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        HighEndingVertexLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntersecVertexLabel {"                       <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        NormalEndingVertexLabel {"                   <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;  
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    EdgeLabels {"                                    <<std::endl;
        NeuroMorphData << "        BottomEdgeLabel {"                           <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopEdgeLabel {"                              <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntermediateEdgeLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopToBottomEdgeLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_1 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.0833333 1 0,"                    <<std::endl;
        NeuroMorphData << "            Id 8"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_2 {"                                <<std::endl;
        NeuroMorphData << "            Color 0 1 0.958333,"                     <<std::endl;
        NeuroMorphData << "            Id 9"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_3 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.840278 0.351575 0.0583526,"      <<std::endl;
        NeuroMorphData << "            Id 10"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_1 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.375 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 15"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_2 {"                                    <<std::endl;
        NeuroMorphData << "            Color 1 1 0.5,"                          <<std::endl;
        NeuroMorphData << "            Id 16"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_3 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.875 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 17"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    ContentType \"HxSpatialGraph\""                  <<std::endl;
        NeuroMorphData << "}"                                                   <<std::endl;
  
        
        NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 "           << std::endl;
        NeuroMorphData << "VERTEX {int VertexLabels } @2 "                       << std::endl;
        
        NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 "                << std::endl;
        NeuroMorphData << "EDGE { int NumEdgePoints } @4 "                      << std::endl;
        NeuroMorphData << "EDGE { int EdgeLabels } @5 "                        << std::endl;
        
        NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 "         << std::endl;
        NeuroMorphData << "POINT { float Radius } @7 "                          << std::endl;
        
        if(inputSpatialGraph->getNumberOfVertices())
        {
                NeuroMorphData << "\n@1 # Vertices xyz coordinates"                     << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                        NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
                
                NeuroMorphData << "\n@2 # Vertex Label" << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                {
                        NeuroMorphData << (*vertexIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfEdges())
        {
                NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
                int last_index = 0;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
                }
                
                NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
                }
                
                NeuroMorphData << "\n@5 # Edge Labels" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfPoints())
        {
                NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
                        }
                        
                }
                
                NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->radiusList.size())
                        {
                                std::list< double >::const_iterator radiusIt;
                                for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
                                        NeuroMorphData << *radiusIt << std::endl;
                        }
                        else
                        {
                                for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
                                {
                                        NeuroMorphData << (*edgeIt)->radius << std::endl;
                                }
                        }
                }
        }
        
        NeuroMorphData.close();
};

void Reader::writeSpatialGraphFile2(AmiraSpatialGraph * given_spatial_graph)
{
    ////std::cout<< " in writing spatial graph " <<std::endl;
//      std::list<std::list<Compartment * > >::iterator edge_list_it;
//      std::list<Compartment * >::iterator edge_it;
        std::vector< Vertex * >::iterator vertexIt;
        std::vector< Edge * >::iterator edgeIt;
        
        int number_of_edge_points = given_spatial_graph->getNumberOfPoints();
        
//      for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
//      {
//              number_of_edge_points += (*edge_list_it).size();
//      }
//      for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
//      for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
        
        std::string format = outputFilename;
        //format += ".am";
        
        #ifdef DEBUG
        //////std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
        ////////std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
        #endif
        
        std::ofstream NeuroMorphData( format.c_str() );
        
        NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
        NeuroMorphData << "# This SpatialGraph file was created by the automated Bouton and Spine detector BoutonFinder " << std::endl;
        NeuroMorphData << "# BoutonFinder was programmed by Christopher Tull" << std::endl;
        NeuroMorphData << "# Max-Planck-Institute for Biological Cybernetics, Tuebingen, Germany " << std::endl;
        
        NeuroMorphData << "define VERTEX " << given_spatial_graph->getNumberOfVertices() << std::endl;
        NeuroMorphData << "define EDGE " << given_spatial_graph->getNumberOfEdges()  << std::endl;
//      NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
        NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
        
        NeuroMorphData << "Parameters {"                                        <<std::endl;
        NeuroMorphData << "    VertexLabels {"                          <<std::endl;
        NeuroMorphData << "        LowEndingVertexLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        HighEndingVertexLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntersecVertexLabel {"                       <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        NormalEndingVertexLabel {"                   <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;  
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    EdgeLabels {"                                    <<std::endl;
        NeuroMorphData << "        BottomEdgeLabel {"                           <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopEdgeLabel {"                              <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntermediateEdgeLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopToBottomEdgeLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_1 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.0833333 1 0,"                    <<std::endl;
        NeuroMorphData << "            Id 8"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_2 {"                                <<std::endl;
        NeuroMorphData << "            Color 0 1 0.958333,"                     <<std::endl;
        NeuroMorphData << "            Id 9"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_3 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.840278 0.351575 0.0583526,"      <<std::endl;
        NeuroMorphData << "            Id 10"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_1 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.375 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 15"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_2 {"                                    <<std::endl;
        NeuroMorphData << "            Color 1 1 0.5,"                          <<std::endl;
        NeuroMorphData << "            Id 16"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_3 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.875 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 17"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    ContentType \"HxSpatialGraph\""                  <<std::endl;
        NeuroMorphData << "}"                                                   <<std::endl;
  
        
        NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 "           << std::endl;
        NeuroMorphData << "VERTEX {int VertexLabels } @2 "                       << std::endl;
        
        NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 "                << std::endl;
        NeuroMorphData << "EDGE { int NumEdgePoints } @4 "                      << std::endl;
        NeuroMorphData << "EDGE { int EdgeLabels } @5 "                        << std::endl;
        
        NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 "         << std::endl;
        NeuroMorphData << "POINT { float Radius } @7 "                          << std::endl;
        
        if(given_spatial_graph->getNumberOfVertices())
        {
            ////std::cout<< given_spatial_graph->getNumberOfVertices() <<std::endl;
                NeuroMorphData << "\n@1 # Vertices xyz coordinates"                     << std::endl;
                for(vertexIt = given_spatial_graph->verticesBegin(); vertexIt != given_spatial_graph->verticesEnd(); ++vertexIt)
                        NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
                
                NeuroMorphData << "\n@2 # Vertex Label" << std::endl;
                for(vertexIt = given_spatial_graph->verticesBegin(); vertexIt != given_spatial_graph->verticesEnd(); ++vertexIt)
                {
                        NeuroMorphData << (*vertexIt)->label << std::endl;
                }
        }
        
        if(given_spatial_graph->getNumberOfEdges())
        {
             ////std::cout<< given_spatial_graph->getNumberOfEdges() <<std::endl;
           
                NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
                int last_index = 0;
                for(edgeIt = given_spatial_graph->edgesBegin(); edgeIt != given_spatial_graph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
                }
                
                NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
                for(edgeIt = given_spatial_graph->edgesBegin(); edgeIt != given_spatial_graph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
                }
                
                NeuroMorphData << "\n@5 # Edge Labels" << std::endl;
                for(edgeIt = given_spatial_graph->edgesBegin(); edgeIt != given_spatial_graph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->label << std::endl;
                }
        }
        
        if(given_spatial_graph->getNumberOfPoints())
        {
            ////std::cout<< given_spatial_graph->getNumberOfPoints() <<std::endl;
           
                NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
                for(edgeIt = given_spatial_graph->edgesBegin(); edgeIt != given_spatial_graph->edgesEnd(); ++edgeIt)
                {
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
                        }
                        
                }
                 ////std::cout<< "done with points" <<std::endl;
                
                NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
                for(edgeIt = given_spatial_graph->edgesBegin(); edgeIt != given_spatial_graph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->radiusList.size())
                        {
                                std::list< double >::const_iterator radiusIt;
                                for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
                                        NeuroMorphData << *radiusIt << std::endl;
                        }
                        else
                        {
                                for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
                                {
                                        NeuroMorphData << (*edgeIt)->radius << std::endl;
                                }
                        }
                }
        }
        
        NeuroMorphData.close();
};

void Reader::writeSpatialGraphFileFromEdges()
{
        std::vector< Edge * >::iterator edgeIt;
        
        int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
        
        std::string format = outputFilename;
        format += ".am";
        
        #ifdef DEBUG
        //////std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
        ////////std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
        #endif
        
        std::ofstream NeuroMorphData( format.c_str() );
        
        NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
        NeuroMorphData << "# This SpatialGraph file was created by the Neuron Reconstruction Tool NeuroMorph " << std::endl;
        NeuroMorphData << "# NeuroMorph was programmed by Marcel Oberlaender and Philip J. Broser," << std::endl;
        NeuroMorphData << "# Max-Planck-Institute for Medical Research Heidelberg, Germany " << std::endl;
        
        int nrOfVertices = 0;
        for(int ii = 0; ii < inputSpatialGraph->getNumberOfEdges(); ++ii)
        {
                if((*(inputSpatialGraph->edgesPointer()))[ii]->edgeConnectivity[0] == (*(inputSpatialGraph->edgesPointer()))[ii]->edgeConnectivity[1])
                        nrOfVertices += 1;
                else
                        nrOfVertices += 2;
        }
        
        NeuroMorphData << "define VERTEX " << nrOfVertices << std::endl;
        NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
        NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
        
        NeuroMorphData << "Parameters {GraphLabels {"                           <<std::endl;
        NeuroMorphData << "        Neuron { "                                   <<std::endl;
        NeuroMorphData << "            Dendrite {"                              <<std::endl;
        NeuroMorphData << "                ApicalDendrite {"                    <<std::endl;
        NeuroMorphData << "                    Color 1 0.5 0.5,"                <<std::endl;
        NeuroMorphData << "                    Id 4 }"                          <<std::endl;
        NeuroMorphData << "                BasalDendrite {"                     <<std::endl;
        NeuroMorphData << "                    Color 0.8 0.4 0.4,"              <<std::endl;
        NeuroMorphData << "                    Id 5 }"                          <<std::endl;
        NeuroMorphData << "                Color 1 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 3 }"                              <<std::endl;
        NeuroMorphData << "            Axon {"                                  <<std::endl;
        NeuroMorphData << "                Color 0 0 1,"                        <<std::endl;
        NeuroMorphData << "                Id 6 }"                              <<std::endl;
        NeuroMorphData << "            Soma {"                                  <<std::endl;
        NeuroMorphData << "                Color 1 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 7 }"                              <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2 }"                                  <<std::endl;
        NeuroMorphData << "        Landmark {"                                  <<std::endl;
        NeuroMorphData << "            Pia {"                                   <<std::endl;
        NeuroMorphData << "                Color 0 1 0.5,"                      <<std::endl;
        NeuroMorphData << "                Id 9 }"                              <<std::endl;
        NeuroMorphData << "            Vessel {"                                <<std::endl;
        NeuroMorphData << "                Color 1 0.5 0,"                      <<std::endl;
        NeuroMorphData << "                Id 10 }"                             <<std::endl;
        NeuroMorphData << "            Barrel {"                                <<std::endl;
        NeuroMorphData << "                aRow {"                              <<std::endl;
        NeuroMorphData << "                    A1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 13 }"                     <<std::endl;
        NeuroMorphData << "                    A2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 14 }"                     <<std::endl;
        NeuroMorphData << "                    A3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 15 }"                     <<std::endl;
        NeuroMorphData << "                    A4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.2 0.2,"            <<std::endl;
        NeuroMorphData << "                        Id 16 }"                     <<std::endl;
        NeuroMorphData << "                Color 1 0.2 0.2,"                    <<std::endl;
        NeuroMorphData << "                Id 12 }"                             <<std::endl;
        NeuroMorphData << "                bRow {"                              <<std::endl;
        NeuroMorphData << "                    B1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 18 }"                     <<std::endl;
        NeuroMorphData << "                    B2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 19 }"                     <<std::endl;
        NeuroMorphData << "                    B3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 20 }"                     <<std::endl;
        NeuroMorphData << "                    B4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.25 0.25,"          <<std::endl;
        NeuroMorphData << "                        Id 21 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.25 0.25,"              <<std::endl;
        NeuroMorphData << "                    Id 17 }"                         <<std::endl;
        NeuroMorphData << "                cRow {"                              <<std::endl;
        NeuroMorphData << "                    C1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 23 }"                     <<std::endl;
        NeuroMorphData << "                    C2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 24 }"                     <<std::endl;
        NeuroMorphData << "                    C3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 25 }"                     <<std::endl;
        NeuroMorphData << "                    C4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 26 }"                     <<std::endl;
        NeuroMorphData << "                    C5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 27 }"                     <<std::endl;
        NeuroMorphData << "                    C6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.3 0.3,"            <<std::endl;
        NeuroMorphData << "                        Id 28 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.3 0.3,"                <<std::endl;
        NeuroMorphData << "                    Id 22 }"                         <<std::endl;
        NeuroMorphData << "                dRow {"                              <<std::endl;
        NeuroMorphData << "                    D1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 30 }"                     <<std::endl;
        NeuroMorphData << "                    D2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 31 }"                     <<std::endl;
        NeuroMorphData << "                    D3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 32 }"                     <<std::endl;
        NeuroMorphData << "                    D4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 33 }"                     <<std::endl;
        NeuroMorphData << "                    D5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 34 }"                     <<std::endl;
        NeuroMorphData << "                    D6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.35 0.35,"          <<std::endl;
        NeuroMorphData << "                        Id 35 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.35 0.35,"              <<std::endl;
        NeuroMorphData << "                    Id 29 }"                         <<std::endl;
        NeuroMorphData << "                eRow {"                              <<std::endl;
        NeuroMorphData << "                    E1 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 37 }"                     <<std::endl;
        NeuroMorphData << "                    E2 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 38 }"                     <<std::endl;
        NeuroMorphData << "                    E3 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 39 }"                     <<std::endl;
        NeuroMorphData << "                    E4 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 40 }"                     <<std::endl;
        NeuroMorphData << "                    E5 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 41 }"                     <<std::endl;
        NeuroMorphData << "                    E6 {"                            <<std::endl;
        NeuroMorphData << "                        Color 1 0.4 0.4,"            <<std::endl;
        NeuroMorphData << "                        Id 42 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.4 0.4,"                <<std::endl;
        NeuroMorphData << "                    Id 36 }"                         <<std::endl;
        NeuroMorphData << "                greekRow {"                          <<std::endl;
        NeuroMorphData << "                    Alpha {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 44 }"                     <<std::endl;
        NeuroMorphData << "                    Beta {"                          <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 45 }"                     <<std::endl;
        NeuroMorphData << "                    Gamma {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 46 }"                     <<std::endl;
        NeuroMorphData << "                    Delta {"                         <<std::endl;
        NeuroMorphData << "                        Color 1 0.1 0.1,"            <<std::endl;
        NeuroMorphData << "                        Id 47 }"                     <<std::endl;
        NeuroMorphData << "                    Color 1 0.1 0.1,"                <<std::endl;
        NeuroMorphData << "                    Id 43 }"                         <<std::endl;
        NeuroMorphData << "                Color 0 1 0,"                        <<std::endl;
        NeuroMorphData << "                Id 11 }"                             <<std::endl;
        NeuroMorphData << "            WhiteMatter {"                           <<std::endl;
        NeuroMorphData << "                Color 0.5 1 0.75,"                   <<std::endl;
        NeuroMorphData << "                Id 48 }"                             <<std::endl;
        NeuroMorphData << "            OtherBarrels {"                          <<std::endl;
        NeuroMorphData << "                Color 1 0 1,"                        <<std::endl;
        NeuroMorphData << "                Id 49 }"                             <<std::endl;
        NeuroMorphData << "            ZAxis {"                                 <<std::endl;
        NeuroMorphData << "                Color 0 0 0,"                        <<std::endl;
        NeuroMorphData << "                Id 50 }"                             <<std::endl;
        NeuroMorphData << "            Color 0 1 1,"                            <<std::endl;
        NeuroMorphData << "            Id 8 }"                                  <<std::endl;
        NeuroMorphData << "        Id 0,"                                       <<std::endl;
        NeuroMorphData << "        Color 0 0 0 }"                               <<std::endl;
        NeuroMorphData << "ContentType \"HxSpatialGraph\" }"                    <<std::endl;
        
        NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 "           << std::endl;
        NeuroMorphData << "VERTEX {int GraphLabels } @2 "                       << std::endl;
        
        NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 "                << std::endl;
        NeuroMorphData << "EDGE { int NumEdgePoints } @4 "                      << std::endl;
        NeuroMorphData << "EDGE { int GraphLabels } @5 "                        << std::endl;
        
        NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 "         << std::endl;
        NeuroMorphData << "POINT { float Radius } @7 "                          << std::endl;
        
        if(inputSpatialGraph->getNumberOfEdges())
        {
                NeuroMorphData << "\n@1 # Vertices xyz coordinates"                     << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
                                NeuroMorphData << (*edgeIt)->edgePointCoordinates.front()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Z_COORD] << std::endl;
                        else
                        {
                                NeuroMorphData << (*edgeIt)->edgePointCoordinates.front()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.front()[Z_COORD] << std::endl;
                                NeuroMorphData << (*edgeIt)->edgePointCoordinates.back()[X_COORD] << " " << (*edgeIt)->edgePointCoordinates.back()[Y_COORD] << " " << (*edgeIt)->edgePointCoordinates.back()[Z_COORD] << std::endl;
                        }
                }
                
                NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
                                NeuroMorphData << (*edgeIt)->label << std::endl;
                        else
                        {
                                NeuroMorphData << (*edgeIt)->label << std::endl;
                                NeuroMorphData << (*edgeIt)->label << std::endl;
                        }
                }
                
                NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
                int lastID = 0;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->edgeConnectivity[0] == (*edgeIt)->edgeConnectivity[1])
                        {
                                NeuroMorphData << lastID << " " << lastID << std::endl;
                                ++lastID;
                        }
                        else
                        {
                                NeuroMorphData << lastID << " " << lastID + 1 << std::endl;
                                lastID += 2;
                        }
                }
                
                NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
                }
                
                NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        NeuroMorphData << (*edgeIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfPoints())
        {
                NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
                        }
                        
                }
                
                NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                        if((*edgeIt)->radiusList.size())
                        {
                                std::list< double >::const_iterator radiusIt;
                                for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
                                        NeuroMorphData << *radiusIt << std::endl;
                        }
                        else
                                for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
                                {
                                        NeuroMorphData << (*edgeIt)->radius << std::endl;
                                }
                                
                }
        }
        
        NeuroMorphData.close();
};

void Reader::writeSpatialGraphFile(int Axon_ID, int Dendrite_ID)
{
//      std::list<std::list<Compartment * > >::iterator edge_list_it;
//      std::list<Compartment * >::iterator edge_it;
        std::vector< Vertex * >::iterator vertexIt;
        std::vector< Edge * >::iterator edgeIt;
        
        int number_of_edge_points = inputSpatialGraph->getNumberOfPoints();
        
//      for(edge_list_it = amira_spatial_graph->edge_list.begin(); edge_list_it != amira_spatial_graph->edge_list.end(); ++edge_list_it)
//      {
//              number_of_edge_points += (*edge_list_it).size();
//      }
//      for(edge_list_contour_it = amira_contour_graph->edge_list.begin(); edge_list_contour_it != amira_contour_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
//      for(edge_list_contour_it = amira_bvp_graph->edge_list.begin(); edge_list_contour_it != amira_bvp_graph->edge_list.end(); ++edge_list_contour_it)
//      {
//              number_of_edge_points += (*edge_list_contour_it).size();
//      }
        
        std::string format = outputFilename;
        format += ".am";
        
        #ifdef DEBUG
        //////std::cout << "WriteSpatialGraphFile: " << format.c_str()  << std::endl;
        ////////std::cout<< "Vertex List Size: " << amira_spatial_graph-> vertice_list.size() << " Edge List Size: "<< amira_spatial_graph->edge_list.size() <<std::endl;
        #endif
        
        std::ofstream NeuroMorphData( format.c_str() );
        
        NeuroMorphData << "# AmiraMesh 3D ASCII 2.0" << std::endl;
        NeuroMorphData << "# This SpatialGraph file was created by the automated Bouton and Spine detector BoutonFinder " << std::endl;
        NeuroMorphData << "# BoutonFinder was programmed by Christopher Tull" << std::endl;
        NeuroMorphData << "# Max-Planck-Institute for Biological Cybernetics, Tuebingen, Germany " << std::endl;
        
        NeuroMorphData << "define VERTEX " << inputSpatialGraph->getNumberOfVertices() << std::endl;
        NeuroMorphData << "define EDGE " << inputSpatialGraph->getNumberOfEdges()  << std::endl;
//      NeuroMorphData << "define GRAPH " << /*amira_spatial_graph->vertice_list.size() +*/ /*amira_contour_graph->vertice_list.size() +*/ inputSpatialGraph->getNumberOfVertices() + /*amira_spatial_graph->edge_list.size() +*/ /*amira_contour_graph->edge_list.size() +*/ inputSpatialGraph->getNumberOfEdges() << std::endl;
        NeuroMorphData << "define POINT " << number_of_edge_points << std::endl;
        
        NeuroMorphData << "Parameters {"                                        <<std::endl;
        NeuroMorphData << "    VertexLabels {"                          <<std::endl;
        NeuroMorphData << "        LowEndingVertexLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        HighEndingVertexLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntersecVertexLabel {"                       <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        NormalEndingVertexLabel {"                   <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;  
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    EdgeLabels {"                                    <<std::endl;
        NeuroMorphData << "        BottomEdgeLabel {"                           <<std::endl;
        NeuroMorphData << "            Color 0 0 1,"                            <<std::endl;
        NeuroMorphData << "            Id 1"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopEdgeLabel {"                              <<std::endl;
        NeuroMorphData << "            Color 0 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 2"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        IntermediateEdgeLabel {"                     <<std::endl;
        NeuroMorphData << "            Color 1 1 0,"                            <<std::endl;
        NeuroMorphData << "            Id 3"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        TopToBottomEdgeLabel {"                      <<std::endl;
        NeuroMorphData << "            Color 1 0 0,"                            <<std::endl;
        NeuroMorphData << "            Id 4"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_1 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.0833333 1 0,"                    <<std::endl;
        NeuroMorphData << "            Id 8"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_2 {"                                <<std::endl;
        NeuroMorphData << "            Color 0 1 0.958333,"                     <<std::endl;
        NeuroMorphData << "            Id 9"                                    <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Dendrite_3 {"                                <<std::endl;
        NeuroMorphData << "            Color 0.840278 0.351575 0.0583526,"      <<std::endl;
        NeuroMorphData << "            Id 10"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_1 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.375 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 15"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_2 {"                                    <<std::endl;
        NeuroMorphData << "            Color 1 1 0.5,"                          <<std::endl;
        NeuroMorphData << "            Id 16"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Axon_3 {"                                    <<std::endl;
        NeuroMorphData << "            Color 0.875 0.5 1,"                      <<std::endl;
        NeuroMorphData << "            Id 17"                                   <<std::endl;
        NeuroMorphData << "        }"                                           <<std::endl;
        NeuroMorphData << "        Color 1 1 1,"                                <<std::endl;
        NeuroMorphData << "        Id 0"                                        <<std::endl;
        NeuroMorphData << "    }"                                               <<std::endl;
        NeuroMorphData << "    ContentType \"HxSpatialGraph\""                  <<std::endl;
        NeuroMorphData << "}"                                                   <<std::endl;
  
        
        NeuroMorphData << "VERTEX { float[3] VertexCoordinates } @1 "           << std::endl;
        NeuroMorphData << "VERTEX {int VertexLabels } @2 "                      << std::endl;
        
        NeuroMorphData << "EDGE { int[2] EdgeConnectivity } @3 "                << std::endl;
        NeuroMorphData << "EDGE { int NumEdgePoints } @4 "                      << std::endl;
        NeuroMorphData << "EDGE { int EdgeLabels } @5 "                         << std::endl;
        
        NeuroMorphData << "POINT { float[3] EdgePointCoordinates } @6 "         << std::endl;
        NeuroMorphData << "POINT { float Radius } @7 "                          << std::endl;
        
        if(inputSpatialGraph->getNumberOfVertices())
        {
                NeuroMorphData << "\n@1 # Vertices xyz coordinates"                     << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                    if((*vertexIt)->label == Axon_ID || (*vertexIt)->label == Dendrite_ID)    
                    NeuroMorphData << (*vertexIt)->coordinates[X_COORD] << " " << (*vertexIt)->coordinates[Y_COORD]  << " " << (*vertexIt)->coordinates[Z_COORD]  << std::endl;
                
                NeuroMorphData << "\n@2 # Vertex Graph Label" << std::endl;
                for(vertexIt = inputSpatialGraph->verticesBegin(); vertexIt != inputSpatialGraph->verticesEnd(); ++vertexIt)
                {
                    if((*vertexIt)->label == Axon_ID || (*vertexIt)->label == Dendrite_ID)
                        NeuroMorphData << (*vertexIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfEdges())
        {
                NeuroMorphData << "\n@3 # Edge Identifiers" << std::endl;
                int last_index = 0;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                    if((*edgeIt)->label == Axon_ID || (*edgeIt)->label == Dendrite_ID)
                        NeuroMorphData << (*edgeIt)->edgeConnectivity[0] << " " << (*edgeIt)->edgeConnectivity[1] << std::endl;
                }
                
                NeuroMorphData << "\n@4 # Number of Points per Edge" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                    if((*edgeIt)->label == Axon_ID || (*edgeIt)->label == Dendrite_ID)
                        NeuroMorphData <<  (*edgeIt)->numEdgePoints <<std::endl;
                }
                
                NeuroMorphData << "\n@5 # Edge Graph Labels" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                    if((*edgeIt)->label == Axon_ID || (*edgeIt)->label == Dendrite_ID)
                        NeuroMorphData << (*edgeIt)->label << std::endl;
                }
        }
        
        if(inputSpatialGraph->getNumberOfPoints())
        {
                NeuroMorphData << "\n@6 # Point xyz coordinates" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                    if((*edgeIt)->label == Axon_ID || (*edgeIt)->label == Dendrite_ID){
                        std::list< double * >::iterator contour_it;
                        for(contour_it = (*edgeIt)->edgePointCoordinates.begin(); contour_it != (*edgeIt)->edgePointCoordinates.end(); ++contour_it)
                        {
                                NeuroMorphData << (*contour_it)[X_COORD] << " " << (*contour_it)[Y_COORD] << " " << (*contour_it)[Z_COORD] << std::endl;
                        }
                    }   
                }
                
                NeuroMorphData << "\n@7 # Radius at Point" << std::endl;
                for(edgeIt = inputSpatialGraph->edgesBegin(); edgeIt != inputSpatialGraph->edgesEnd(); ++edgeIt)
                {
                    if((*edgeIt)->label == Axon_ID || (*edgeIt)->label == Dendrite_ID){    
                    if((*edgeIt)->radiusList.size())
                        {
                                std::list< double >::const_iterator radiusIt;
                                for(radiusIt = (*edgeIt)->radiusList.begin(); radiusIt != (*edgeIt)->radiusList.end(); ++radiusIt)
                                        NeuroMorphData << *radiusIt << std::endl;
                        }
                        else
                                for(int ii = 0; ii < (*edgeIt)->edgePointCoordinates.size(); ++ii)
                                {
                                        NeuroMorphData << (*edgeIt)->radius << std::endl;
                                }
                }             
                }
        }
        
        NeuroMorphData.close();
};









