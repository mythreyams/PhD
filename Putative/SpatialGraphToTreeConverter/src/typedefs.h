/****************************************************************************/
/*                                                                          */
/* File:      typedefs.h                                                    */
/*                                                                          */
/* Purpose:   header file for all necessary inculdes and typedefs for files */
/*            associated with the NeuroMorph project                        */
/*                                                                          */
/* Author:    Marcel Oberlaender                                            */
/*            Max-Planck-Institut fuer medizinische Forschung               */
/*            Jahnstrasse 29                                                */
/*            D-69120 Heidelberg                                            */
/*                                                                          */
/*                                                                          */
/* EMail: marcel.oberlaender@mpimf-heidelberg.mpg.de                        */
/*                                                                          */
/* History:   09.02.2006                                                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/
#ifndef TYPEDEF
#define TYPEDEF

//#define DEBUG


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


//#include "string"
#include <sstream>
#include <string.h>
//#include <direct.h>
#include <stdio.h>
#include "unistd.h"
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <iomanip>
// #include "fann.h"
// #include "floatfann.h"


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <map>


//stats regarding image resolution
#define CONFOCAL

#ifdef CONFOCAL
  #define XYSAMPLING  	0.092  //confocal
#else
  #define XYSAMPLING  	0.184    //brightfield
#endif

#define ZSAMPLING 	 0.5

//length of the double array used to represent points  
#define ARRAY_LENGTH 32
  
//VoxelFeatures of the Measurement Vector Type
#define X_COORD			0
#define Y_COORD			1
#define Z_COORD			2
#define IDENTIFIER		3
#define SURFACE			4
#define AVERAGE_STD_DEV		5
#define LOCAL_BRIGHTNESS	6
#define THRESHOLD		7
#define LENGTH_UNTIL_HERE	8
#define IS_BOUTON		9
#define IS_TRUE_BOUTON		10
#define LOCAL_SIGMA		11
#define IS_PAIRED		12
#define IS_USER1_BOUTON		13
#define IS_USER2_BOUTON		14
#define IS_USER3_BOUTON		15
#define CONFIDENCE		16
  
#define LOCAL_BRIGHTNESS_RAD2	17
#define LOCAL_BRIGHTNESS_RAD3	18
#define LOCAL_SIGMA_RAD3	19
#define AVG_BRIGHT_RANGE3	20
#define AVG_BRIGHT_RANGE50	21
#define AVG_RAD_RANGE3		22
#define AVG_RAD_RANGE4		23
#define AVG_RAD_RANGE5		24
  
#define VALID_FOR_NN		25
#define AVG_SURF                26
#define AVG_BRIGHT              27  
  
#define ZAxis			50
  
//Edge directions
#define FORWARD			true
#define BACKWARD		false
  
  
#define AXON_ID_1       15
#define AXON_ID_2       16
#define AXON_ID_3       17
#define UNKNOWN_ID      18
#define DENDRITE_ID_1   8
#define DENDRITE_ID_2   9
#define DENDRITE_ID_3   10  
  
// #define PRINT_TRAINING_DATA	1
// #define TRAIN_NETWORK		

  
typedef struct{
  float coords[3];
  float magnitude;
  
}VECTOR;
#endif
