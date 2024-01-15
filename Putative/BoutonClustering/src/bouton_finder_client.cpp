

#include "bouton_finder.h"
#include "bouton_params.h"
#include "math.h"

// iterator for the params
//std::list< bouton_params * > BoutonParamsList;

// since we are collecting params for all the points in the ROI
// it would be a good idea to change this to list again
// for now putting 500 as the size
bouton_params BoutonParamsArray[50000];// = new bouton_params[20];
int NumOfLandmarks = 0;
int AxonID = 0;


void writeToXls(char * ouptputXlsPath, int section_number, int roi_number)
{
    

    bool file_exsists = 0;
    std::string s1("BoutonParams.xls");
    std::string s2("BoutonParams.csv");
    
    std::cout<<"In writexls"<<std::endl;
    
    if(chdir(ouptputXlsPath) != 0)
            perror("Couldn't open image drirectory!");
    
    // Iterate over the params list and write each param of each bouton one by one
    std::list< bouton_params * >::iterator paramIterator;
    //std::ofstream xlsParamWriter("BoutonParams.xls",+a);
    
    std::ifstream f1(s1.c_str());
    if (f1.good()) {
        f1.close();
        file_exsists = true;
    } else {
        f1.close();
        file_exsists = false;
    } 
    
    std::ifstream f2(s2.c_str());
    if (f2.good()) {
        f2.close();
        file_exsists = true;
    } else {
        f2.close();
        file_exsists = false;
    } 
    
    
    std::fstream xlsParamWriter;
    std::fstream csvParamWriter;
    
    
    xlsParamWriter.open ("BoutonParams.xls", std::fstream::in | std::fstream::out | std::fstream::app);
    csvParamWriter.open ("BoutonParams.csv", std::fstream::in | std::fstream::out | std::fstream::app);
    
    // write header if it does not exists
    if(!file_exsists)
    {
        xlsParamWriter<<"Section Number"<<'\t'<<"ROI Number"<<'\t'
        <<"Actual Landmark X"<<'\t'<<"Actual Landmark Y"<<'\t'<<"Actual Landmark Z"<<'\t'
        <<"Transfomed Landmark X"<<'\t'<<"Transformed Landmark Y"<<'\t'<<"Transformed Landmark Z"<<'\t'
        <<"Min Radius XY 1 plane"<<'\t'<<"Min Radius XY 3 plane"<<'\t'<<"Min Radius XY 5 plane"<<'\t' 
        <<"Avg Radius XY 1 plane"<<'\t'<<"Avg Radius XY 3 plane"<<'\t'<<"Avg Radius XY 5 plane"<<'\t' 
        <<"Brightness Total XY 1 plane"<<'\t'
        <<"Brightness Contour XY 1 plane"<<'\t'<<"Brightness Contour XY 3 plane"<<'\t'<<"Brightness Contour XY 5 plane"<<'\t'
        <<"Brightness Min Contour XY 1 plane"<<'\t'<<"Brightness Min Contour XY 3 plane"<<'\t'<<"Brightness Min Contour XY 5 plane"<<'\t'
        <<"Brightness ratio Contour XY 1 plane"<<'\t'<<"Brightness ratio Contour XY 3 plane"<<'\t'<<"Brightness ratio Contour XY 5 plane"<<'\t'
        <<"Roi Orientation XZ Axon 1 plane"<<'\t'<<"Roi Orientation XZ Axon 2 plane"<<'\t'<<"Roi Orientation XZ Axon 3 plane"<<'\t' 
        <<"1 pt Orientation XY"<<'\t'<<"1 pt Orientation YZ"<<'\t'<<"1 pt Orientation XZ"<<'\t'
        <<"2 pt Orientation XY"<<'\t'<<"2 pt Orientation YZ"<<'\t'<<"2 pt Orientation XZ"<<'\t'
        <<"3 pt Orientation XY"<<'\t'<<"3 pt Orientation YZ"<<'\t'<<"3 pt Orientation XZ"<<'\t'
        <<"4 pt Orientation XY"<<'\t'<<"4 pt Orientation YZ"<<'\t'<<"4 pt Orientation XZ"<<'\t'
        <<"5 pt Orientation XY"<<'\t'<<"5 pt Orientation YZ"<<'\t'<<"5 pt Orientation XZ"<<'\n';
        
         csvParamWriter<<"Section Number"<<","<<"ROI Number"<<","
        <<"Actual Landmark X"<<'\t'<<"Actual Landmark Y"<<'\t'<<"Actual Landmark Z"<<'\t'
	 <<"Transfomed Landmark X"<<","<<"Transformed Landmark Y"<<","<<"Transformed Landmark Z"<<","
        <<"Min Radius XY 1 plane"<<","<<"Min Radius XY 3 plane"<<","<<"Min Radius XY 5 plane"<<"," 
        <<"Avg Radius XY 1 plane"<<'\t'<<"Avg Radius XY 3 plane"<<'\t'<<"Avg Radius XY 5 plane"<<'\t' 
        <<"Brightness Total XY 1 plane"<<","
        <<"Brightness Contour XY 1 plane"<<","<<"Brightness Contour XY 3 plane"<<","<<"Brightness Contour XY 5 plane"<<","
        <<"Brightness Min Contour XY 1 plane"<<","<<"Brightness Min Contour XY 3 plane"<<","<<"Brightness Min Contour XY 5 plane"<<","
        <<"Brightness ratio Contour XY 1 plane"<<","<<"Brightness ratio Contour XY 3 plane"<<","<<"Brightness ratio Contour XY 5 plane"<<","
        <<"Roi Orientation XZ Axon 1 plane"<<'\t'<<"Roi Orientation XZ Axon 2 plane"<<'\t'<<"Roi Orientation XZ Axon 3 plane"<<'\t' 
	<<"1 pt Orientation XY"<<","<<"1 pt Orientation YZ"<<","<<"1 pt Orientation XZ"<<","
        <<"2 pt Orientation XY"<<","<<"2 pt Orientation YZ"<<","<<"2 pt Orientation XZ"<<","
        <<"3 pt Orientation XY"<<","<<"3 pt Orientation YZ"<<","<<"3 pt Orientation XZ"<<","
        <<"4 pt Orientation XY"<<","<<"4 pt Orientation YZ"<<","<<"4 pt Orientation XZ"<<","
        <<"5 pt Orientation XY"<<","<<"5 pt Orientation YZ"<<","<<"5 pt Orientation XZ"<<'\n';
    }
    
    
    for(int i = 0; i<NumOfLandmarks; i++)
    {
        //std::cout<< (*paramIterator)->landmark[0]<<std::endl;
        
        xlsParamWriter<<section_number<<'\t'<<roi_number<<'\t'
        <<BoutonParamsArray[i].landmark[0]<<'\t'<<BoutonParamsArray[i].landmark[1]<<'\t'<<BoutonParamsArray[i].landmark[2]<<'\t'
        <<BoutonParamsArray[i].transformed_landmark[0]<<'\t'<<BoutonParamsArray[i].transformed_landmark[1]<<'\t'<<BoutonParamsArray[i].transformed_landmark[2]<<'\t'
        <<BoutonParamsArray[i].radius_1_plane_XY[0]<<'\t'<<BoutonParamsArray[i].radius_3_plane_XY[0]<<'\t'<<BoutonParamsArray[i].radius_5_plane_XY[0]<<'\t'
        <<BoutonParamsArray[i].radius_1_plane_XY[1]<<'\t'<<BoutonParamsArray[i].radius_3_plane_XY[1]<<'\t'<<BoutonParamsArray[i].radius_5_plane_XY[1]<<'\t'
        <<BoutonParamsArray[i].Brightness_1_plane_Total<<'\t'
        <<BoutonParamsArray[i].Brightness_1_plane_XY_Contour[0]<<'\t'<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour[0]<<'\t'<<BoutonParamsArray[i].Brightness_5_plane_XY_Contour[0]<<'\t'
        <<BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0]<<'\t'<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour_Min[0]<<'\t'<<BoutonParamsArray[i].Brightness_5_plane_XY_Contour_Min[0]<<'\t'
        <<BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0]/BoutonParamsArray[i].Brightness_1_plane_XY_BB[0]<<'\t'<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour_Min[0] / BoutonParamsArray[i].Brightness_3_plane_XY_BB[0]<<'\t'
        <<BoutonParamsArray[i].Brightness_5_plane_XY_Contour_Min[0] / BoutonParamsArray[i].Brightness_5_plane_XY_BB[0]<<'\t'
        <<BoutonParamsArray[i].Axon_1_XZ_Angle<<'\t'<<BoutonParamsArray[i].Axon_2_XZ_Angle<<'\t'<<BoutonParamsArray[i].Axon_3_XZ_Angle<<'\t';
        
       
        
        for(int it = 0; it < 15; it++)
        {
            //std::cout<<(*paramIterator)->orientation[it]<<std::endl;
            xlsParamWriter<< BoutonParamsArray[i].orientation[it]<<'\t';
        }
        
        
        xlsParamWriter<<'\n';
        
        std::fstream xlsParamWriter;
        csvParamWriter<<section_number<<","<<roi_number<<","
        <<BoutonParamsArray[i].landmark[0]<<","<<BoutonParamsArray[i].landmark[1]<<","<<BoutonParamsArray[i].landmark[2]<<","
	<<BoutonParamsArray[i].transformed_landmark[0]<<","<<BoutonParamsArray[i].transformed_landmark[1]<<","<<BoutonParamsArray[i].transformed_landmark[2]<<","
        <<BoutonParamsArray[i].radius_1_plane_XY[0]<<","<<BoutonParamsArray[i].radius_3_plane_XY[0]<<","<<BoutonParamsArray[i].radius_5_plane_XY[0]<<","
        <<BoutonParamsArray[i].radius_1_plane_XY[1]<<","<<BoutonParamsArray[i].radius_3_plane_XY[1]<<","<<BoutonParamsArray[i].radius_5_plane_XY[1]<<","
        <<BoutonParamsArray[i].Brightness_1_plane_Total<<","
	<<BoutonParamsArray[i].Brightness_1_plane_XY_Contour[0]<<","<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour[0]<<","<<BoutonParamsArray[i].Brightness_5_plane_XY_Contour[0]<<","
        <<BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0]<<","<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour_Min[0]<<","<<BoutonParamsArray[i].Brightness_5_plane_XY_Contour_Min[0]<<","
        <<BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0]/BoutonParamsArray[i].Brightness_1_plane_XY_BB[0]<<","<<BoutonParamsArray[i].Brightness_3_plane_XY_Contour_Min[0] / BoutonParamsArray[i].Brightness_3_plane_XY_BB[0]<<","
        <<BoutonParamsArray[i].Brightness_5_plane_XY_Contour_Min[0] / BoutonParamsArray[i].Brightness_5_plane_XY_BB[0]<<","
        <<BoutonParamsArray[i].Axon_1_XZ_Angle<<","<<BoutonParamsArray[i].Axon_2_XZ_Angle<<","<<BoutonParamsArray[i].Axon_3_XZ_Angle<<",";
        
        
        for(int it = 0; it < 15; it++)
        {
            //std::cout<<(*paramIterator)->orientation[it]<<std::endl;
            csvParamWriter<< BoutonParamsArray[i].orientation[it]<<",";
        }
        
        
        csvParamWriter<<'\n';
        
        
    }
    
    
    
    xlsParamWriter.close();
    csvParamWriter.close();
    
#if 0

    std::fstream xlsParamWriter;
    xlsParamWriter.open ("BoutonContours.xls", std::fstream::in | std::fstream::out | std::fstream::app);
    
    std::cout<<"In writexls"<<std::endl;
    
     if(chdir(ouptputXlsPath) != 0)
            perror("Couldn't open image drirectory!");
    
      
        for(int i = 0; i < NumOfLandmarks; i++)
        {
        
            xlsParamWriter<< BoutonParamsArray[i].landmark[0]<<'\t'<<BoutonParamsArray[i].landmark[1]<<'\t'<<BoutonParamsArray[i].landmark[2]<<'\t'
            <<BoutonParamsArray[i].transformed_landmark[0]/XYSAMPLING<<'\t'<<BoutonParamsArray[i].transformed_landmark[1]/XYSAMPLING<<'\t'<<BoutonParamsArray[i].transformed_landmark[2]/ZSAMPLING<<'\n';
            
            for(int it = 0; it < 72; it++)
            {
                //std::cout<<(*paramIterator)->orientation[it]<<std::endl;
                xlsParamWriter<< BoutonParamsArray[i].contour_1_plane_XY[it]<<'\t';
            }
            
            
        
            xlsParamWriter<<'\n';
            
            xlsParamWriter<<BoutonParamsArray[i].Brightness_1_plane_XY_Contour[0];
            xlsParamWriter<<'\n';
            
            std::list<double*>::iterator itr;
            
            for(itr= BoutonParamsArray[i].brightness_points_1.begin(); itr != BoutonParamsArray[i].brightness_points_1.end(); itr++)
            {
                xlsParamWriter<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
                std::cout<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
               
            }
        }
#endif

}

#if 0
// Read the landmark locations
bool getLandmarks(char * landmark_folder)
{
    int pram_array_ind =0;
    // open the bouton landmark file
    if(chdir(landmark_folder) != 0)
            perror("Couldn't open image drirectory!");
    
    
    std::string s("Bouton.landmarkAscii");
    
    
    std::ifstream boutonStream(s.c_str());
    
    
    if (!boutonStream.good())
    {
        boutonStream.close();
        return false;
    } 
    else 
    {
       
      //std::ifstream boutonStream("Bouton.landmarkAscii");
            
    // read all the landmarks into a list
    
    
    //readLandmarkCoordinates(&boutonStream, 1);
    
    
    
           // const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
            const char * numbers = "0123456789";
            const char * signs = "+-";
            //const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
            const char * whitespace = "\t\n ";
            
            std::string currentLine;
            unsigned int line = 0;
                            
            //bool parameters = 1;
            //bool transform = 0;
            //unsigned int brackets = 0, transformBrackets = 0;
            //unsigned int currentIndex = 0;
            bool inValidSection = 0;
            //int landmarkCount = 0;
            
            while(!std::getline(boutonStream, currentLine).eof() )
            {
                // set the flag for the right section in the file
                //std::cout<<currentLine<<std::endl;
                    
                // in the right section read the coordinates
                if(inValidSection && currentLine.size())
                {
                    //std::cout<<currentLine<<std::endl;
                    double bouton_coord[3];
                    //bouton_params* b;
                    unsigned int count = 0;
                    std::string::size_type loc1, loc2, loc3;
                    loc1 = 0;
                    loc2 = currentLine.find_first_of(signs, 0);
                    if(loc2 != std::string::npos)
                            if(loc2 < loc1)
                                    loc1 = loc2;
                    loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                    while(loc2 != std::string::npos && count<3)
                    {
                            
                            char * tmp1 = new char[20];
                            
                            for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                            {
                                // 
                                tmp1[i] = 'f';
                            }
                            currentLine.copy(tmp1, loc2 - loc1, loc1);
                            double ftmp1 = atof(tmp1);
                            
                            
                            bouton_coord[count]= ftmp1;
                            //std::cout << ftmp1 << std::endl;
                            loc3 = loc2;
                            loc1 = currentLine.find_first_of(numbers, loc3);
                            loc2 = currentLine.find_first_of(signs, loc3);
                            if(loc2 != std::string::npos)
                                    if(loc2 < loc1)
                                            loc1 = loc2;
                            loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                            //std::cout << loc2 <<'\t'<<count<< std::endl;
                            ++count;
                            delete [] tmp1;
                    }
                    
                   
                    
                     
                    /*
                    // initialize the bouton landmark thru constructor
                    else
                    {
                        b = new bouton_params(bouton_coord, 1);
                    }*/
                    //inValidSection = 0;
                    
                    
                    
                    
                    //std::cout<<(*It)->landmark<<'\n';
                    
                    //bouton_params * tmpParam = new bouton_params(bouton_coord, 1);
                    //BoutonParamsArray[pram_array_ind] = *tmpParam;
                    
                    BoutonParamsArray[pram_array_ind].landmark[0] = bouton_coord[0];
                    BoutonParamsArray[pram_array_ind].landmark[1] = bouton_coord[1];
                    BoutonParamsArray[pram_array_ind].landmark[2] = bouton_coord[2];
                    
                    //transformed_landmark = new double [3];
                    BoutonParamsArray[pram_array_ind].transformed_landmark[0] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[1] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[2] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[3] = 0;
                    
                    BoutonParamsArray[pram_array_ind].isBouton = 1;
                    
                    //nearestPoint = new double [3];
                    BoutonParamsArray[pram_array_ind].nearestPoint[0] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[1] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[2] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[3] = 0;
                    
                    BoutonParamsArray[pram_array_ind].distanceFromLandmark = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_1_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_1_plane_XY[1] = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_3_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_3_plane_XY[1] = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_5_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_5_plane_XY[1] = 0;
                    
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                       BoutonParamsArray[pram_array_ind].contour_1_plane_XY[j] = 0;
                    }
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                        BoutonParamsArray[pram_array_ind].contour_3_plane_XY[j] = 0;
                    }
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                        BoutonParamsArray[pram_array_ind].contour_5_plane_XY[j] = 0;
                    }
                     
                    BoutonParamsArray[pram_array_ind].Brightness_1_plane_XY_BB[0] = 0;
                    BoutonParamsArray[pram_array_ind].Brightness_3_plane_XY_BB[0] = 0;
                    BoutonParamsArray[pram_array_ind].Brightness_5_plane_XY_BB[0] = 0;
                    
                    
                    pram_array_ind++;
                        //BoutonParamsList.push_back(tmpParam);
                        
                        
                        //delete [] tmpParam;
                     //new bouton_params(bouton_coord, 1);
                    //std::cout<<tmpParam->landmark[0]<<'\t'<<tmpParam->landmark[1]<<'\t'<<tmpParam->landmark[2]<<'\n';
                    //std::cout<<&tmpParam->landmark[0]<<'\t'<<&tmpParam->landmark[1]<<'\t'<<&tmpParam->landmark[2]<<'\n';
                    
                }
                
                else if(currentLine.find("@1", 0) != std::string::npos)
                {
                    //std::cout<<currentLine<<std::endl;
                    ++line;
                    if(line == 2)
                    {
                        inValidSection = 1;
                    }
                }
                    
                    
            }
            
            boutonStream.close();
            
        NumOfLandmarks = pram_array_ind; 
        
        return true;
            
    }
    
    
    
    
    
    // open non bouton file and read into a the same list
#if 0
    std::ifstream NonboutonStream("NonBouton.landmarkAscii");
                    
    if(!NonboutonStream.fail())
    {
            //const char * letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
            const char * numbers = "0123456789";
            const char * signs = "+-";
            //const char * otherChars = ":;\'\"\\()[]{}!@#$%^&_=|<>?";
            const char * whitespace = "\t\n ";
            
            std::string currentLine;
            unsigned int line = 0;
                            
            //bool parameters = 1;
            //bool transform = 0;
           // unsigned int brackets = 0, transformBrackets = 0;
           // unsigned int currentIndex = 0;
            bool inValidSection = 0;
            
            
            while(!std::getline(NonboutonStream, currentLine).eof() )
            {
                // set the flag for the right section in the file
                
                    //std::cout<<currentLine<<std::endl;
                // in the right section read the coordinates
                if(inValidSection && currentLine.size())
                {
                    //std::cout<<currentLine<<std::endl;
                    double bouton_coord[3];
                    //bouton_params* nb;
                    unsigned int count = 0;
                    std::string::size_type loc1, loc2, loc3;
                    loc1 = 0;
                    loc2 = currentLine.find_first_of(signs, 0);
                    if(loc2 != std::string::npos)
                            if(loc2 < loc1)
                                    loc1 = loc2;
                    loc2 = currentLine.find_first_of(whitespace, loc1 + 1); //ignores last value: is always 1 anyways
                    while(loc2 != std::string::npos && count<3)
                    {
                            
                            char * tmp1 = new char[20];
                            
                            for(int i = 0; i < 20/*loc2-loc1+1*/; i++ )
                            {
                                // 
                                tmp1[i] = 'f';
                            }
                            currentLine.copy(tmp1, loc2 - loc1, loc1);
                            double ftmp1 = atof(tmp1);
                            
                            
                            bouton_coord[count]= ftmp1;
                            //std::cout << ftmp1 << std::endl;
                            loc3 = loc2;
                            loc1 = currentLine.find_first_of(numbers, loc3);
                            loc2 = currentLine.find_first_of(signs, loc3);
                            if(loc2 != std::string::npos)
                                    if(loc2 < loc1)
                                            loc1 = loc2;
                            loc2 = currentLine.find_first_of(whitespace, loc1 + 1);
                            //std::cout << loc2 <<'\t'<<count<< std::endl;
                            ++count;
                            delete [] tmp1;
                    }
                    
                    // initialize the bouton landmark thru constructor
                    
                    //bouton_params * tmpParam = new bouton_params(bouton_coord, 0);
                    //BoutonParamsList.push_back(tmpParam);
                    
                    
                    //BoutonParamsArray[pram_array_ind] = *tmpParam;
                    BoutonParamsArray[pram_array_ind].landmark[0] = bouton_coord[0];
                    BoutonParamsArray[pram_array_ind].landmark[1] = bouton_coord[1];
                    BoutonParamsArray[pram_array_ind].landmark[2] = bouton_coord[2];
                    
                    //transformed_landmark = new double [3];
                    BoutonParamsArray[pram_array_ind].transformed_landmark[0] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[1] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[2] = 0;
                    BoutonParamsArray[pram_array_ind].transformed_landmark[3] = 0;
                    
                    BoutonParamsArray[pram_array_ind].isBouton = 0;
                    
                    //nearestPoint = new double [3];
                    BoutonParamsArray[pram_array_ind].nearestPoint[0] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[1] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[2] = 0;
                    BoutonParamsArray[pram_array_ind].nearestPoint[3] = 0;
                    
                    BoutonParamsArray[pram_array_ind].distanceFromLandmark = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_1_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_1_plane_XY[1] = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_3_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_3_plane_XY[1] = 0;
                    
                    BoutonParamsArray[pram_array_ind].radius_5_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].radius_5_plane_XY[1] = 0;
                    
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                       BoutonParamsArray[pram_array_ind].contour_1_plane_XY[j] = 0;
                    }
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                        BoutonParamsArray[pram_array_ind].contour_3_plane_XY[j] = 0;
                    }
                    for(int j = 0; j< NUM_RAYS * 4 ; j++)
                    {
                        BoutonParamsArray[pram_array_ind].contour_5_plane_XY[j] = 0;
                    }
                     
                    BoutonParamsArray[pram_array_ind].Brightness_1_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].Brightness_3_plane_XY[0] = 0;
                    BoutonParamsArray[pram_array_ind].Brightness_5_plane_XY[0] = 0;
                    
                    
                    pram_array_ind++;
                    
                   //std::cout<<tmpParam->landmark[0]<<'\t'<<tmpParam->landmark[1]<<'\t'<<tmpParam->landmark[2]<<'\n';
                    
                }
                
                else if(currentLine.find("@1", 0) != std::string::npos)
                {
                    //std::cout<<currentLine<<std::endl;
                    ++line;
                    if(line == 2)
                    {
                        inValidSection = 1;
                    }
                }
                    
                    
            }
            
    }
    
    NonboutonStream.close();
    
#endif
    
    
}


#endif

int main( int argc , char * argv[])
{
    
    int start_index = atoi(argv[3]);
    int end_index = atoi(argv[4]);
    // bool time_to_transformgraph = true;
    
    int section_number = atoi(argv[7]);
    int roi_number = atoi(argv[8]);
    
    std::cout << std::endl << argc <<std::endl;
   std::cout << std::endl << "roi num"<<roi_number <<std::endl;
   std::cout << std::endl << "roi num"<<argv[8] <<std::endl;
    
    if(argc != 9)
    {
        std::cout << std::endl << "Proper format:  ./BoutonFinder 'axonfile.am' 'proximity_cropped_image_stack_and_hxfile_path' start_index end_index 'loadhxfile_with_translation' 'outputpath'" << std::endl << std::endl;
        std::cout << std::endl << argc <<std::endl;
        return -1;
    }
  

    BoutonFinder* finder = new BoutonFinder(argv[1], argv[2], start_index, end_index, argv[5], argv[6]);
    
    //finder->getAxonPoints(finder->getOriginalImage());
    // constructor reads t
    //bouton_params * params;
    // Read landmark file and get the location into bouton and non bouton structures
    /*if(!getLandmarks(argv[2]))
    {
        return -1;
    }*/
    
     
//
    
    //std::list< bouton_params * > boutonParamsList;
    // iterate over all the landmarks read and find parameters for each of them
    
    
    //finder->findSpatialGraphPoint(/*&boutonParamsList*/);
    
    //std::cout<<"Axon ID"<<" "<<AxonID<<std::endl;
    
    
    
    //finder->transformGraphandLandmark(/*&boutonParamsList*/);
    
    
    
    /*for(int i = 0; i< NumOfLandmarks; i++)
    {
        
        std::cout<<i<<'\t'<<BoutonParamsArray[i].transformed_landmark[0]<<'\t'<<BoutonParamsArray[i].transformed_landmark[1]<<'\t'<<BoutonParamsArray[i].transformed_landmark[2]<<std::endl;
        
        
    }*/
    
    std::cout<<"Num of Landmarks "<<NumOfLandmarks<<std::endl;
   
    // Distance from landmark only applies to boutons not non boutons
    //for(int i = 0; i< NumOfLandmarks; i++)
    {
       
        
        //BoutonParamsArray[i].distanceFromLandmark = Utility::euclidean_distance(BoutonParamsArray[i].landmark,BoutonParamsArray[i].nearestPoint,3,1);
       
    }
    
    // Populate the non bouton points which fall within this ROI region
    
    
    
    //std::cout<<NumOfLandmarks<<std::endl;
    
    
    
    
  //#if 0  
    for(int i = 0; i< NumOfLandmarks; i++)
    {
       
        std::cout<< i <<std::endl;
        //BoutonParamsArray[i].distanceFromLandmark = Utility::euclidean_distance(BoutonParamsArray[i].landmark,BoutonParamsArray[i].nearestPoint,3,1);
        
        
        //(*paramIterator)->orientation = new std::list<double *>;
        finder->findOrientationOfAxon(BoutonParamsArray[i].transformed_landmark,5,BoutonParamsArray[i].orientation);     
        
        
        // Find the radius at each point in 1,3 and 5 planes
        
        
        int z_pos = rint(BoutonParamsArray[i].transformed_landmark[Z_COORD]/ZSAMPLING);
	std::cout<< BoutonParamsArray[i].transformed_landmark[X_COORD]/XYSAMPLING<<'\t'<<BoutonParamsArray[i].transformed_landmark[Y_COORD]/XYSAMPLING
	<<BoutonParamsArray[i].transformed_landmark[Z_COORD]/ZSAMPLING<<'\n';  
	
	std::cout<< BoutonParamsArray[i].landmark[X_COORD]<<'\t'<<BoutonParamsArray[i].landmark[Y_COORD]
	<<BoutonParamsArray[i].landmark[Z_COORD]<<'\n';  
            
        std::cout<<z_pos<<std::endl;
        //std::cout<<i<<'\t'<<BoutonParamsArray[i].transformed_landmark[Z_COORD]<<std::endl;
        
        
        Image2DType::Pointer image_plane_1 = finder->getImagePlane(z_pos,1,finder->original_image);
        finder->sendRay(BoutonParamsArray[i].transformed_landmark[0],BoutonParamsArray[i].transformed_landmark[1],image_plane_1,BoutonParamsArray[i].radius_1_plane_XY,BoutonParamsArray[i].contour_1_plane_XY, &BoutonParamsArray[i].contour_1_plane_XY_min_rad, BoutonParamsArray[i].orientation[0]);
        
        
        Image2DType::Pointer image_plane_3 = finder->getImagePlane(z_pos,3,finder->original_image);
        finder->sendRay(BoutonParamsArray[i].transformed_landmark[0],BoutonParamsArray[i].transformed_landmark[1],image_plane_3,BoutonParamsArray[i].radius_3_plane_XY,BoutonParamsArray[i].contour_3_plane_XY, &BoutonParamsArray[i].contour_3_plane_XY_min_rad, BoutonParamsArray[i].orientation[0]);
        
        
        Image2DType::Pointer image_plane_5 = finder->getImagePlane(z_pos,5,finder->original_image);
        finder->sendRay(BoutonParamsArray[i].transformed_landmark[0],BoutonParamsArray[i].transformed_landmark[1],image_plane_5,BoutonParamsArray[i].radius_5_plane_XY,BoutonParamsArray[i].contour_5_plane_XY, &BoutonParamsArray[i].contour_5_plane_XY_min_rad, BoutonParamsArray[i].orientation[0]);
                
        
        finder->findBrightness(image_plane_1, BoutonParamsArray[i].contour_1_plane_XY, BoutonParamsArray[i].contour_1_plane_XY_min_rad,BoutonParamsArray[i].Brightness_1_plane_XY_BB, BoutonParamsArray[i].Brightness_1_plane_XY_Contour,BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min,&BoutonParamsArray[i].brightness_points_1);
        //std::cout<<BoutonParamsArray[i].Brightness_1_plane_XY[0]<<std::endl;
        
//         std::list<double*>::iterator itr;
//             
//         for(itr= BoutonParamsArray[i].brightness_points_1.begin(); itr != BoutonParamsArray[i].brightness_points_1.end(); itr++)
//         {
//             //xlsParamWriter<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
//             std::cout<< (*itr)[0]<<'\t'<<(*itr)[1]<<'\t'<<(*itr)[2]<<'\n';
//             
//         }
        
        /*for(int ind = 0; ind < NUM_RAYS*4; ind=ind+2)
        {
            
            std::cout<< BoutonParamsArray[i].contour_1_plane_XY[ind]<<'\t'<<BoutonParamsArray[i].contour_1_plane_XY[ind+1]<<'\n';             
            
            
        }*/
        
        
        finder->findBrightness(image_plane_3, BoutonParamsArray[i].contour_3_plane_XY, BoutonParamsArray[i].contour_3_plane_XY_min_rad,BoutonParamsArray[i].Brightness_3_plane_XY_BB,BoutonParamsArray[i].Brightness_3_plane_XY_Contour,BoutonParamsArray[i].Brightness_3_plane_XY_Contour_Min,&BoutonParamsArray[i].brightness_points_3);
        //std::cout<<BoutonParamsArray[i].Brightness_3_plane_XY[0]<<std::endl;
              
        finder->findBrightness(image_plane_5, BoutonParamsArray[i].contour_5_plane_XY, BoutonParamsArray[i].contour_5_plane_XY_min_rad,BoutonParamsArray[i].Brightness_5_plane_XY_BB,BoutonParamsArray[i].Brightness_5_plane_XY_Contour,BoutonParamsArray[i].Brightness_5_plane_XY_Contour_Min,&BoutonParamsArray[i].brightness_points_5);
        //std::cout<<BoutonParamsArray[i].Brightness_5_plane_XY[0]<<std::endl;
      
        
        
    }
    
    
    
    for (int j = 0; j < NumOfLandmarks; j++)
    {
      
      // Get the avg brightness of the roi along the axonfile
      BoutonParamsArray[0].Brightness_1_plane_Total += BoutonParamsArray[j].Brightness_1_plane_XY_Contour_Min[0];
      
      
    }
    
    
    for (int j = 0; j < NumOfLandmarks; j++)
    {
      
      // Get the avg brightness of the roi along the axonfile
      BoutonParamsArray[j].Brightness_1_plane_Total = BoutonParamsArray[0].Brightness_1_plane_Total;
      
      
    }
    
    
    
    if(NumOfLandmarks)
    {
       // write the parameters into an xls file
       writeToXls(argv[6],section_number,roi_number);
    }
    
    //finder->writeBoutonLandmarkfile();
        
    
    //finder->findBoutons();
    

    
    delete finder;

    //delete params;

    return 0;
};
