/****************************************************************************/
/*                                                                          */
/* File:      utility.h                                                     */
/*                                                                          */
/* Purpose:   Header file for the Utility class                             */
/*                                                                          */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   31.03.2014                                                    */
/*                                                                          */
/* Remarks:   Contains a number of static utility methods defining          */
/*            common computations and needed tasks                          */
/*                                                                          */
/****************************************************************************/


#include "typedefs.h"

#ifndef UTILITY
#define UTILITY

class Utility       
{
 public:
  Utility();
  
  static int calc_mode(long* histogram, int size);
  static float calc_mean(long* histogram, int size, int lower_threshold);
  static int calc_max(long* histogram, int size);
  static float calc_variance(long* histogram, int size, float mean, int lower_threshold);
  double calc_vector_mean(std::vector<double> list);
  double calc_vector_variance(std::vector<double> list, double mean);
  static int calc_median(long* histogram, int lower, int upper);
  static double euclidean_distance(double * a, double * b, int dim, float z_scalar);
  static int calc_dotProduct_int(int a[], int b[], int dim);
  static double calc_dotProduct_double(double * a, double * b, int dim);
  static std::string get_cwd();

  static std::string makeString(int num) 
  {
    std::ostringstream ss;
    ss << num;
    return ss.str();
  }

 private:

};

#endif