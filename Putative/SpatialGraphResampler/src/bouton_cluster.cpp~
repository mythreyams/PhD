/****************************************************************************/
/*                                                                          */
/* File:      bouton_cluster.cpp                                            */
/*                                                                          */
/* Purpose:                                                                 */
/*                                                                          */
/*                                                                          */
/* Author:    Mythreya Seetharama                                           */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail: marcel.oberlaender@tuebingen.mpg.de                               */
/*                                                                          */
/* History:   15-02-2016                                                    */
/*                                                                          */
/* Remarks: This module colsters the putative bouton points based on        */
/*          radius and brightness and finds the centriod                    */
/*                                                                          */
/****************************************************************************/
// #include "bouton_params.h"
#include "bouton_cluster.h"
#include "bouton_finder.h"
#include "string.h"

using namespace std;

#define MAX_INTER_BOUTON_DIST   0.3
#define MIN_BOUTON_STRENGTH     4

extern bouton_params BoutonParamsArray[];
//extern int NumOfLandmarks;


bouton_cluster::bouton_cluster(int NumOfLandmarks, std::list<double *>* bouton_centroids, std::list<double *>* overlap_centroids)
{

  ////std:://cout<<NumOfLandmarks<<std::endl;
    double *centroid1 = new double[2];
    double *centroid2 = new double[2];
    bool done_with_clustering = false;
    
    // set the intial centroids
    // find the brightest and dimmest points as init centroids
    getInitialCentroids( NumOfLandmarks, centroid1, centroid2 );
    
     
    while(done_with_clustering == false)
    {
        
        // Assign each point to the nearest centroid
        assignClusters( NumOfLandmarks, centroid1, centroid2, &done_with_clustering);
        
        // get new centroids
        computeNewCentriods( NumOfLandmarks, centroid1, centroid2);
       
        
    }
    
    createRowBoutonsList(NumOfLandmarks);
    
    groupBoutonPoints();
    
    // find the centre of mass for valid bouton groups
    findBoutonandOverlapCentroids(  NumOfLandmarks, bouton_centroids, overlap_centroids);
    
    
};


// Private member functions

void bouton_cluster:: getInitialCentroids( int NumOfLandmarks, double *centroid1, double *centroid2)
{
    double max_bri = 0;
    double min_bri = 0xffff;
    int centroididx1,centroididx2;
  
    for (int i = 0; i < NumOfLandmarks; i++)
    {
        if(BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0] < min_bri)
        {
            min_bri = BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
            centroididx2 = i;
        }
        
        if(BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0] > max_bri)
        {
            max_bri = BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
            centroididx1 = i;
        }
        
        // init the previous cluster state here
        BoutonParamsArray[i].prev_which_cluster = 0;
    }
    
    
    centroid1[0] = BoutonParamsArray[centroididx1].radius_1_plane_XY_perpendicular_otsu[0];
    centroid1[1] = BoutonParamsArray[centroididx1].Brightness_1_plane_XY_Contour_Min[0];
    
    
    centroid2[0] = BoutonParamsArray[centroididx2].radius_1_plane_XY_perpendicular_otsu[0];
    centroid2[1] = BoutonParamsArray[centroididx2].Brightness_1_plane_XY_Contour_Min[0];
    
     
};

void bouton_cluster:: assignClusters( int NumOfLandmarks, double *centroid1, double *centroid2, bool *done_with_clustering)
{
    double dist1=0;
    double dist2=0;
    
    *done_with_clustering = true;
    
    for (int i = 0; i < NumOfLandmarks; i++)
    {
        double *param = new double[2];
        
        param[1] = BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
        param[0] = BoutonParamsArray[i].radius_1_plane_XY_perpendicular_otsu[0];
        
       
        
        dist1 = Utility::euclidean_distance(param,centroid1,2,1);
        dist2 = Utility::euclidean_distance(param,centroid2,2,1);
                  
        
        if (dist1 > dist2)
        {
            // this point belongs to bouton cluster
            BoutonParamsArray[i].which_cluster = 2;
	    
            
        }
        else
        {
            // this one belongs to non bouton cluster
            BoutonParamsArray[i].which_cluster = 1;
	    ////std:://cout << BoutonParamsArray[i].landmark[0]  << " " << BoutonParamsArray[i].landmark[Y_COORD] << " " << BoutonParamsArray[i].landmark[Z_COORD] << std::endl;
	    
            
            
        }
        
        // wait for 3 consecutive non membership changes to conclude that the clustering is done
        if(BoutonParamsArray[i].which_cluster != BoutonParamsArray[i].prev_which_cluster)
        {
            *done_with_clustering = false;
            
        }
        BoutonParamsArray[i].prev_which_cluster = BoutonParamsArray[i].which_cluster;
    }
    
    
};



void bouton_cluster:: computeNewCentriods( int NumOfLandmarks, double *centroid1, double *centroid2)
{
    double *sumcl1 = new double[2];
    double *sumcl2 = new double[2];
    double cl1cnt = 0;
    double cl2cnt = 0;
   
    
    for (int i = 0; i < NumOfLandmarks; i++)
    {
        
        if(BoutonParamsArray[i].which_cluster == 1)
        {
            sumcl1[0] += BoutonParamsArray[i].radius_1_plane_XY_perpendicular_otsu[0];
            sumcl1[1] += BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
            cl1cnt++;
         
        }
        else
        {
            sumcl2[0] += BoutonParamsArray[i].radius_1_plane_XY_perpendicular_otsu[0];
            sumcl2[1] += BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
            cl2cnt++;
             
             
        }
        
    }
   
    centroid1[0] = sumcl1[0] / cl1cnt;
    centroid1[1] = sumcl1[1] / cl1cnt;
    
    centroid2[0] = sumcl2[0] / cl2cnt;
    centroid2[1] = sumcl2[1] / cl2cnt;
   
};

void bouton_cluster::createRowBoutonsList( int NumOfLandmarks)
{
    

    for (int i = 0; i < NumOfLandmarks; i++)
    {
        
        if(BoutonParamsArray[i].which_cluster == 1)
        {
           ////std:://cout << BoutonParamsArray[i].landmark[0]  << " " << BoutonParamsArray[i].landmark[Y_COORD] << " " << BoutonParamsArray[i].landmark[Z_COORD] << std::endl;
           double *tempcoor = new double[12];
           tempcoor[0] = BoutonParamsArray[i].landmark[0];
           tempcoor[1] = BoutonParamsArray[i].landmark[1];
           tempcoor[2] = BoutonParamsArray[i].landmark[2];
           tempcoor[3] = BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
           tempcoor[4] = BoutonParamsArray[i].overlap[0]; 
           tempcoor[5] = BoutonParamsArray[i].overlap[1];
           tempcoor[6] = BoutonParamsArray[i].perpendicular_otsu_front[0];
           tempcoor[7] = BoutonParamsArray[i].perpendicular_otsu_front[1];
           tempcoor[8] = BoutonParamsArray[i].perpendicular_otsu_front[2];
           tempcoor[9] = BoutonParamsArray[i].perpendicular_otsu_back[0];
           tempcoor[10] = BoutonParamsArray[i].perpendicular_otsu_back[1];
           tempcoor[11] = BoutonParamsArray[i].perpendicular_otsu_back[2];
           
            if(BoutonParamsArray[i].overlap[0])
            {
              //////////std:://cout<<"bouton overlap front"<< tempcoor[6] << " "<< " "<< tempcoor[7] <<" "<<tempcoor[8]<< std::endl;
            }
            
            if(BoutonParamsArray[i].overlap[1])
            {
              //////////std:://cout<<"bouton overlap back"<< tempcoor[9] << " "<< " "<< tempcoor[10] <<" "<<tempcoor[11]<< std::endl;
            }
           
            rawBoutonList.push_back(tempcoor); 
         
        }
        
    }
    
};

void bouton_cluster:: groupBoutonPoints(void)
{
    
    std::list<double *> tempList;
    
    std::list<double *>::iterator rawBoutonsIterator;// = rawBoutons->begin();
    
    std::list<double *>::iterator tempListIterator;// = tempList.begin();
    
    
    // state machine
    state = Init;
    while(state != End)
    {
        switch(state)
        {
            case Init:
            {  
                //cout<< state<<endl;
                if(!rawBoutonList.empty())
                {
                    rawBoutonsIterator = rawBoutonList.begin();
                    tempList.push_back(*rawBoutonsIterator);
                    //tempListIterator = tempList.begin();
                    
                    rawBoutonsIterator = rawBoutonList.erase(rawBoutonsIterator);
                    
                    state = Compare;
                    
                }
                else
                {
                    state = End;
                }
            }
            break;
            
            case Compare:
            {
                //cout<< state<<endl;
                for(tempListIterator = tempList.begin(); tempListIterator != tempList.end(); tempListIterator++)
                {
                    
                    for(rawBoutonsIterator = rawBoutonList.begin(); rawBoutonsIterator != rawBoutonList.end(); rawBoutonsIterator++)
                    {
                        double dist = Utility::euclidean_distance(*tempListIterator,*rawBoutonsIterator, 3,1 );
                        
                        if(dist < MAX_INTER_BOUTON_DIST)
                        {
                            ////std:://cout<<"inner while if"<<std::endl;
                            // push this entry into the temp list 
                            tempList.push_back(*rawBoutonsIterator);
                            
                            // remove the element from the  
                            rawBoutonsIterator = rawBoutonList.erase(rawBoutonsIterator);
                            
                        }
                         
                    }
                    
                }
                
                state = Push;
                
            }
            break;
            
            case Push:
            {
               //cout<< state<<endl;
               if(!tempList.empty())
               {
                   //std::list<double *> insertlist = tempList;
                   //name_list2 = new List<string>(name_list1);
                   list<double *> *insertlist = new list<double *>(tempList);
                    
                   //cout<<"pointer "<<insertlist<<" "<<insertlist<<endl;
                   groupedBoutonLists.push_back(insertlist);
                    
                   tempList.erase(tempList.begin(),tempList.end());
               }
               
               state = Init;
                
            }
            break;
            
            default:
                break;
            
        }
    }
    
    
};

void bouton_cluster:: findBoutonandOverlapCentroids( int NumOfLandmarks, std::list<double *>* bouton_centroids, std::list<double *>* overlap_centroids)
{
    BoutonFinder *bfinder = new BoutonFinder();
    // get the valid bouton lists and find centroid
    list< list< double *> *>::iterator it;
    for(it = groupedBoutonLists.begin(); it != groupedBoutonLists.end(); it++)
    {
        // this is a contigous bouton list
        list<double *> *il= *it;
        
        if (il->size() >= MIN_BOUTON_STRENGTH)
        {
            std::list<double *>::iterator temp_bouton_coords_it;
            double *bouton_centroid_overlap = new double[11];
            double brisum = 0;
            int overlap_front_count = 0;
            int overlap_back_count = 0;
            
            
            for (int i =0; i < 11; i++)
            {
                bouton_centroid_overlap[i] = 0;
            }
            
            //////////std:://cout<< "init "<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
            
            for( temp_bouton_coords_it = il->begin(); temp_bouton_coords_it != il->end() ; temp_bouton_coords_it++)
            {
                bouton_centroid_overlap[0] += (*temp_bouton_coords_it)[0] * (*temp_bouton_coords_it)[3];
                bouton_centroid_overlap[1] += (*temp_bouton_coords_it)[1] * (*temp_bouton_coords_it)[3];
                bouton_centroid_overlap[2] += (*temp_bouton_coords_it)[2] * (*temp_bouton_coords_it)[3]; 
                brisum += (*temp_bouton_coords_it)[3];
                
                if((*temp_bouton_coords_it)[4] == true)
                {
                    
                    bouton_centroid_overlap[3] = true;
                    bouton_centroid_overlap[5] += (*temp_bouton_coords_it)[6];
                    bouton_centroid_overlap[6] += (*temp_bouton_coords_it)[7];
                    bouton_centroid_overlap[7] += (*temp_bouton_coords_it)[8];
                    overlap_front_count++;
                    ////////std:://cout<< "valid "<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
                    
                    double * tempOverlap = new double[3]; 
                    tempOverlap[0] = (*temp_bouton_coords_it)[6];
                    tempOverlap[1] = (*temp_bouton_coords_it)[7];
                    tempOverlap[2] = (*temp_bouton_coords_it)[8];
                    //overlap_centroids->push_back(tempOverlap);
                    
                    
                }
                if((*temp_bouton_coords_it)[5] == true)
                {
                    
                    bouton_centroid_overlap[4] = true;
                    bouton_centroid_overlap[8] += (*temp_bouton_coords_it)[9];
                    bouton_centroid_overlap[9] += (*temp_bouton_coords_it)[10];
                    bouton_centroid_overlap[10] += (*temp_bouton_coords_it)[11];
                    overlap_back_count++;
                    ////////std:://cout<< "valid "<<(*temp_bouton_coords_it)[9] << " "<< (*temp_bouton_coords_it)[10] << " "<< (*temp_bouton_coords_it)[11] << std::endl;
                    
                    double * tempOverlap = new double[3]; 
                    tempOverlap[0] = (*temp_bouton_coords_it)[9];
                    tempOverlap[1] = (*temp_bouton_coords_it)[10];
                    tempOverlap[2] = (*temp_bouton_coords_it)[11];
                    //overlap_centroids->push_back(tempOverlap);
                }
                ////std:://cout<< (*temp_bouton_coords_it)[0] << " "<< (*temp_bouton_coords_it)[1] << " "<< (*temp_bouton_coords_it)[2] <<(*temp_bouton_coords_it)[3]<< std::endl;
                
            }
            
            //////////std:://cout<< "front count "<<overlap_front_count << std::endl;
            //////////std:://cout<< "back cnt "<<overlap_back_count << std::endl;
            
            ////////////std:://cout <<"sum"<< std::endl;
            ////////////std:://cout<< bouton_centroid[0] << " "<< bouton_centroid[1] << " "<< bouton_centroid[2] << std::endl;
            //////////////std:://cout<< "divisor"<< ( ((*temp_bouton_coords_it)[3])) << std::endl;
            
            
            bouton_centroid_overlap[0] = bouton_centroid_overlap[0]/brisum;
            bouton_centroid_overlap[1] = bouton_centroid_overlap[1]/brisum;
            bouton_centroid_overlap[2] = bouton_centroid_overlap[2]/brisum;
            
            if(bouton_centroid_overlap[3] == true)
            {
                double * tempOverlap = new double[3]; 
                //////////std:://cout<< "sum"<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
                bouton_centroid_overlap[5] /= overlap_front_count;
                bouton_centroid_overlap[6] /= overlap_front_count;
                bouton_centroid_overlap[7] /= overlap_front_count;
                //////////std:://cout<< "divide"<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
                
                tempOverlap[0] = bouton_centroid_overlap[5];
                tempOverlap[1] = bouton_centroid_overlap[6];
                tempOverlap[2] = bouton_centroid_overlap[7];
                overlap_centroids->push_back(tempOverlap);
                
                //delete tempOverlap;
                //std:://cout<< "overlap centroid front"<<tempOverlap[0] << " "<< tempOverlap[1] << " "<< tempOverlap[2] << std::endl;
                
                
            }
            
            if(bouton_centroid_overlap[4] == true)
            {
                double * tempOverlap = new double[3];
                bouton_centroid_overlap[8]  /= overlap_back_count;
                bouton_centroid_overlap[9]  /= overlap_back_count;
                bouton_centroid_overlap[10] /= overlap_back_count;
                
                tempOverlap[0] = bouton_centroid_overlap[8];
                tempOverlap[1] = bouton_centroid_overlap[9];
                tempOverlap[2] = bouton_centroid_overlap[10];
                overlap_centroids->push_back(tempOverlap);
                
                
                //delete tempOverlap;
                //std:://cout<< "overlap centroid back "<<tempOverlap[0] << " "<< tempOverlap[1] << " "<< tempOverlap[2] << std::endl;
                
                
            }
            
            ////////////std:://cout << std::endl;
            ////////////std:://cout<< bouton_centroid[0] << " "<< bouton_centroid[1] << " "<< bouton_centroid[2] << std::endl;
//                 double * centroid_only = new double[3];
//              centroid_only[0] = bouton_centroid_overlap[0];
//                 centroid_only[1] = bouton_centroid_overlap[1];
//              centroid_only[2] = bouton_centroid_overlap[2];
            bouton_centroids->push_back(bouton_centroid_overlap);
            
            
        
        }
        
        
    }
    //bfinder->writeListinAmiraCoords( *bouton_centroids, "/home/mythreya/Documents/output/S05_S28/diameter", "BoutonCentroids",0,0,0,1);
    //bfinder->writeListinAmiraCoords( *overlap_centroids, "/home/mythreya/Documents/output/S05_S28/diameter", "OverlapCentroids",0,1,1,1);
    delete bfinder;
    
};
    



    
  /*
   int k = 0;
    ////////////std:://cout<<"in findcentroids"<<std::endl;
    std::list<double *> bouton_coords; 
    std::list<double *> temp_bouton_coords;
    double groupCnt = 0;
    
    
    
    // create a list of bouton coords
    for (int i = 0; i < NumOfLandmarks; i++)
    {
        
        if(BoutonParamsArray[i].which_cluster == 1)
        {
	   ////std:://cout << BoutonParamsArray[i].landmark[0]  << " " << BoutonParamsArray[i].landmark[Y_COORD] << " " << BoutonParamsArray[i].landmark[Z_COORD] << std::endl;
           double *tempcoor = new double[12];
           tempcoor[0] = BoutonParamsArray[i].landmark[0];
           tempcoor[1] = BoutonParamsArray[i].landmark[1];
           tempcoor[2] = BoutonParamsArray[i].landmark[2];
           tempcoor[3] = BoutonParamsArray[i].Brightness_1_plane_XY_Contour_Min[0];
	   tempcoor[4] = BoutonParamsArray[i].overlap[0]; 
	   tempcoor[5] = BoutonParamsArray[i].overlap[1];
	   tempcoor[6] = BoutonParamsArray[i].perpendicular_otsu_front[0];
	   tempcoor[7] = BoutonParamsArray[i].perpendicular_otsu_front[1];
	   tempcoor[8] = BoutonParamsArray[i].perpendicular_otsu_front[2];
	   tempcoor[9] = BoutonParamsArray[i].perpendicular_otsu_back[0];
	   tempcoor[10] = BoutonParamsArray[i].perpendicular_otsu_back[1];
	   tempcoor[11] = BoutonParamsArray[i].perpendicular_otsu_back[2];
	   
	    if(BoutonParamsArray[i].overlap[0])
	    {
	      //////////std:://cout<<"bouton overlap front"<< tempcoor[6] << " "<< " "<< tempcoor[7] <<" "<<tempcoor[8]<< std::endl;
	    }
	    
	    if(BoutonParamsArray[i].overlap[1])
	    {
	      //////////std:://cout<<"bouton overlap back"<< tempcoor[9] << " "<< " "<< tempcoor[10] <<" "<<tempcoor[11]<< std::endl;
	    }
	   
            bouton_coords.push_back(tempcoor); 
         
        }
        
    }
    
    bfinder->writeListinAmiraCoords( bouton_coords, "/home/mythreya/Documents/output/S05_S28/diameter", "before_sort",0,0,0,0);
   
    
    // need to sort raw bouton coords, since nodes usually mess up the sequence of points
    std::list<double *>sorted_bouton_coords;
    
    std::list<double *>::iterator bouton_coords_it;
    std::list<double *>::iterator sorted_bouton_coords_it;
    std::list<double *>::iterator temp_it;
    
    
    bouton_coords_it = bouton_coords.begin();
    
    
    // put the first guy in sorted list
    sorted_bouton_coords.push_back(*bouton_coords_it);
    sorted_bouton_coords_it = sorted_bouton_coords.begin();
    
    bouton_coords.erase(bouton_coords_it);
    
    
    
    
    // printing only block
    int loopi = 0;
    list< list< double *> *>::iterator it;
    for(it = groupedBoutonLists.begin(); it != groupedBoutonLists.end(); it++)
    {
        
        list<double *> *il= *it;
        //cout<<"pointer main"<<*it<<endl;
        list<double *>::iterator ilit;
        
        for (ilit = il->begin();  ilit != il->end() ; ilit++)
        {
            //cout<<"printing back"<<(*ilit)[0]<<(*ilit)[1]<<(*ilit)[2]<<endl;
        }
        
        loopi++;
        std::string format = "valid_boutons_list_";
        if(loopi ==1)
        {
            format += "1";
        }
        if(loopi ==2)
        {
            format += "2";
        }
        if(loopi ==3)
        {
            format += "3";
        }
        if(loopi ==4)
        {
            format += "4";
        }
        if(loopi ==5)
        {
            format += "5";
        }
        if(loopi ==6)
        {
            format += "6";
        }
        if(loopi ==7)
        {
            format += "7";
        }
        //cout<<*it<<endl;
        //cout<<format<<endl;
        //cout<<loopi<<endl;
        bfinder->writeListinAmiraCoords( **it, "/home/mythreya/Documents/output/S05_S28/diameter", format.c_str(),0,0,0,0);
    }
    
    //bfinder->writeListinAmiraCoords( sorted_bouton_coords, "/home/mythreya/Documents/output/S05_S28/diameter", "after_sort",0,0,0,0);*/
    
//     BoutonFinder *bfinder  = new BoutonFinder();
//     std::list<double *>::iterator new_bouton_coords_it;
//     std::list<double *>::iterator new_next_bouton_coords_it;
//     
//     new_bouton_coords_it = sorted_bouton_coords.begin();
//     new_next_bouton_coords_it = sorted_bouton_coords.begin();
//     new_next_bouton_coords_it++;
//     
//     for( ; new_next_bouton_coords_it != sorted_bouton_coords.end() ; new_bouton_coords_it++,new_next_bouton_coords_it++)
//     {
//         
//         ////////////std:://cout<<"in main loop"<<std::endl;
//         
//         temp_bouton_coords.push_back(*new_bouton_coords_it);
//         
//         
//         double dist = Utility::euclidean_distance(*new_bouton_coords_it,*new_next_bouton_coords_it,3,1);
// 	////std:://cout<<" "<<(*new_bouton_coords_it)[0]<<" "<<(*new_bouton_coords_it)[1]<<" "<<(*new_bouton_coords_it)[2]<<std::endl;
// 	////std:://cout<<" "<<(*new_next_bouton_coords_it)[0]<<" "<<(*new_next_bouton_coords_it)[1]<<" "<<(*new_next_bouton_coords_it)[2]<<std::endl;
//         ////std:://cout << " "<< dist<<std::endl;
//         if (dist > MAX_INTER_BOUTON_DIST)
//         {
//              
//             // this is a seperate bouton grp.. lets see if it is big enuf
//             if( temp_bouton_coords.size() > MIN_BOUTON_STRENGTH)
//             {
//                 // this a valid grp of boutons.. lets find the centroid
//                 std::list<double *>::iterator temp_bouton_coords_it;
//                 double *bouton_centroid_overlap = new double[11];
//                 double brisum = 0;
// 		int overlap_front_count = 0;
//                 int overlap_back_count = 0;
// 		
// 		
// 		for (int i =0; i < 11; i++)
// 		{
// 		  bouton_centroid_overlap[i] = 0;
// 		}
// 		
// 		//////////std:://cout<< "init "<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
// 		
//                 for( temp_bouton_coords_it = temp_bouton_coords.begin(); temp_bouton_coords_it != temp_bouton_coords.end() ; temp_bouton_coords_it++)
//                 {
//                     bouton_centroid_overlap[0] += (*temp_bouton_coords_it)[0] * (*temp_bouton_coords_it)[3];
//                     bouton_centroid_overlap[1] += (*temp_bouton_coords_it)[1] * (*temp_bouton_coords_it)[3];
//                     bouton_centroid_overlap[2] += (*temp_bouton_coords_it)[2] * (*temp_bouton_coords_it)[3]; 
//                     brisum += (*temp_bouton_coords_it)[3];
// 		    
// 		    if((*temp_bouton_coords_it)[4] == true)
// 		    {
// 		      
// 		      bouton_centroid_overlap[3] = true;
// 		      bouton_centroid_overlap[5] += (*temp_bouton_coords_it)[6];
// 		      bouton_centroid_overlap[6] += (*temp_bouton_coords_it)[7];
// 		      bouton_centroid_overlap[7] += (*temp_bouton_coords_it)[8];
// 		      overlap_front_count++;
// 		      ////////std:://cout<< "valid "<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
// 		      
// 		      double * tempOverlap = new double[3]; 
// 		      tempOverlap[0] = (*temp_bouton_coords_it)[6];
// 		      tempOverlap[1] = (*temp_bouton_coords_it)[7];
// 		      tempOverlap[2] = (*temp_bouton_coords_it)[8];
// 		      overlap_centroids->push_back(tempOverlap);
// 		      
// 		      
// 		    }
// 		    if((*temp_bouton_coords_it)[5] == true)
// 		    {
// 		      
// 		      bouton_centroid_overlap[4] = true;
// 		      bouton_centroid_overlap[8] += (*temp_bouton_coords_it)[9];
// 		      bouton_centroid_overlap[9] += (*temp_bouton_coords_it)[10];
// 		      bouton_centroid_overlap[10] += (*temp_bouton_coords_it)[11];
// 		      overlap_back_count++;
// 		      ////////std:://cout<< "valid "<<(*temp_bouton_coords_it)[9] << " "<< (*temp_bouton_coords_it)[10] << " "<< (*temp_bouton_coords_it)[11] << std::endl;
// 		      
// 		      double * tempOverlap = new double[3]; 
// 		      tempOverlap[0] = (*temp_bouton_coords_it)[9];
// 		      tempOverlap[1] = (*temp_bouton_coords_it)[10];
// 		      tempOverlap[2] = (*temp_bouton_coords_it)[11];
// 		      overlap_centroids->push_back(tempOverlap);
// 		    }
//                     //////////////std:://cout<< (*temp_bouton_coords_it)[0] << " "<< (*temp_bouton_coords_it)[1] << " "<< (*temp_bouton_coords_it)[2] <<(*temp_bouton_coords_it)[3]<< std::endl;
//                     
//                 }
//                 
//                 //////////std:://cout<< "front count "<<overlap_front_count << std::endl;
// 		//////////std:://cout<< "back cnt "<<overlap_back_count << std::endl;
//                 
//                 ////////////std:://cout <<"sum"<< std::endl;
//                 ////////////std:://cout<< bouton_centroid[0] << " "<< bouton_centroid[1] << " "<< bouton_centroid[2] << std::endl;
//                 //////////////std:://cout<< "divisor"<< ( ((*temp_bouton_coords_it)[3])) << std::endl;
//                 
//                 
//                 bouton_centroid_overlap[0] = bouton_centroid_overlap[0]/brisum;
//                 bouton_centroid_overlap[1] = bouton_centroid_overlap[1]/brisum;
//                 bouton_centroid_overlap[2] = bouton_centroid_overlap[2]/brisum;
// 		
// 		if(bouton_centroid_overlap[3] == true)
// 		{
// 		  double * tempOverlap = new double[3]; 
// 		  //////////std:://cout<< "sum"<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
// 		  bouton_centroid_overlap[5] /= overlap_front_count;
// 		  bouton_centroid_overlap[6] /= overlap_front_count;
// 		  bouton_centroid_overlap[7] /= overlap_front_count;
// 		  //////////std:://cout<< "divide"<<bouton_centroid_overlap[5] << " "<< bouton_centroid_overlap[6] << " "<< bouton_centroid_overlap[7]<< std::endl;
// 		  
// 		  tempOverlap[0] = bouton_centroid_overlap[5];
// 		  tempOverlap[1] = bouton_centroid_overlap[6];
// 		  tempOverlap[2] = bouton_centroid_overlap[7];
// 		  //overlap_centroids->push_back(tempOverlap);
// 		  
// 		  //delete tempOverlap;
// 		  ////////std:://cout<< "overlap centroid front"<<tempOverlap[0] << " "<< tempOverlap[1] << " "<< tempOverlap[2] << std::endl;
// 		  
// 		  
// 		}
// 		
// 		if(bouton_centroid_overlap[4] == true)
// 		{
// 		  double * tempOverlap = new double[3];
// 		  bouton_centroid_overlap[8]  /= overlap_back_count;
// 		  bouton_centroid_overlap[9]  /= overlap_back_count;
// 		  bouton_centroid_overlap[10] /= overlap_back_count;
// 		  
// 		  tempOverlap[0] = bouton_centroid_overlap[8];
// 		  tempOverlap[1] = bouton_centroid_overlap[9];
// 		  tempOverlap[2] = bouton_centroid_overlap[10];
// 		  //overlap_centroids->push_back(tempOverlap);
// 		  
// 		  
// 		  //delete tempOverlap;
// 		  ////////std:://cout<< "overlap centroid back "<<tempOverlap[0] << " "<< tempOverlap[1] << " "<< tempOverlap[2] << std::endl;
// 		  
// 		   
// 		}
// 		
//                 ////////////std:://cout << std::endl;
//                 ////////////std:://cout<< bouton_centroid[0] << " "<< bouton_centroid[1] << " "<< bouton_centroid[2] << std::endl;
// //                 double * centroid_only = new double[3];
// // 		centroid_only[0] = bouton_centroid_overlap[0];
// //                 centroid_only[1] = bouton_centroid_overlap[1];
// // 		centroid_only[2] = bouton_centroid_overlap[2];
// 		bouton_centroids->push_back(bouton_centroid_overlap);
// // 		delete centroid_only;
// 		//get_auto_boutons()->push_back(bouton_centroid);
// 		//auto_boutons->push_back(bouton_centroid);
// 		std::string format = "valid_boutons_";
// 		k++;
// 		if(k ==1)
// 		{
// 		  format += "1";
// 		}
// 		if(k ==2)
// 		{
// 		  format += "2";
// 		}
// 		if(k ==3)
// 		{
// 		  format += "3";
// 		}
// 		if(k ==4)
// 		{
// 		  format += "4";
// 		}
// 		//format += std::to_string(k);
//                 bfinder->writeListinAmiraCoords( temp_bouton_coords, "/home/mythreya/Documents/output/S05_S28/diameter", format.c_str(),0,0,0,0);
//                 // empty the temp list so that we can start afresh
//                 temp_bouton_coords.erase(temp_bouton_coords.begin(),temp_bouton_coords.end());
//                
//                 ////////////std:://cout<< "size after delete"<<temp_bouton_coords.size() << std::endl;
//                 
//             }
//             else
//             {
//                 temp_bouton_coords.erase(temp_bouton_coords.begin(),temp_bouton_coords.end());
//             }
//             
//         }
//        
//     }
//     delete bfinder;
//    
// };



