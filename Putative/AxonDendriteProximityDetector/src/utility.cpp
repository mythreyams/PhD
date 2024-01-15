
#include "utility.h"
#include <list>


/****************************************************************************/
/*get_cwd() utility function that returns string with current working dir   */
/****************************************************************************/

std::string Utility::get_cwd()
{
    std::size_t buf_size = 1024;
    char* buf = NULL;
    char* r_buf;
  
    buf = static_cast<char*>(realloc(buf, buf_size));
    r_buf = getcwd(buf, buf_size);

    std::string str(buf);
    free(buf);
    return str;
};


int Utility::calc_mode(long* histogram, int size)
{
    int temp = 0, mode = 0;
    
    for(int i=0;i<size;i++)    									//calculate mode of pixels
    {
      ////std::cout << i << " "<< mode_array[i] << std::endl;
      if(histogram[i] > temp)
      {
	temp = histogram[i];
	mode = i;
      }
    }
    
    return mode;
};


float Utility::calc_mean(long* histogram, int size, int lower_threshold)
{
  double sum=0, count=0;
  
  for(int i=lower_threshold+1;i<size;i++)
  {
    sum += histogram[i] * i;
    count += histogram[i];
  }
  
  return (float) (sum/count);
};

double Utility::calc_vector_mean(std::vector<double> list)
{
  double sum=0;
  
  for(int i=0;i<list.size();i++)
  {
    sum += list.at(i);
  }
  
  return (double) (sum/list.size());
};

double Utility::calc_vector_variance(std::vector<double> list, double mean)
{
  double sum=0, temp=0;
  
  
  for(int i=0;i<list.size();i++)
  {
    temp = list.at(i) - mean;
    sum += (temp*temp);
  }
  
  return (double) (sum/list.size());
};



int Utility::calc_max(long* histogram, int size)
{
  int max=0;
  for(int i=0;i<size;i++)
  {
    if(histogram[i] != 0)
      max = i;
  }
  
  return max;
};


float Utility::calc_variance(long* histogram, int size, float mean, int lower_threshold)
{
  double sum=0, count=0, temp=0;
  
  
  for(int i=lower_threshold+1;i<size;i++)
  {
    temp = i - mean;
    sum += histogram[i] * (temp*temp);
    count += histogram[i];
  }
  
  return (float) (sum/count);
};

int Utility::calc_median(long* histogram, int lower, int upper)
{
  long sum=0;
  int i;
  
  
  for(i=lower;i<upper;i++)
  {
      sum += histogram[i];
  }
  sum /= 2;
  
  for(i=lower;i<upper;i++)
  {   
      if(sum <=0)
	break;
	
      sum -= histogram[i];
  }
  
  return i;
};

int Utility::calc_dotProduct_int(int a[], int b[], int dim)
{
    int sum = 0;
    
    for(int i=0; i< dim; i++)
    {
	sum += a[i]*b[i];
    }
    return sum;
};

double Utility::calc_dotProduct_double(double * a, double * b, int dim)
{
    double sum = 0;
    
    for(int i=0; i< dim; i++)
    {
	sum += a[i]*b[i];
    }
    return sum;
};

double Utility::euclidean_distance(double * a, double * b, int dim, float z_scalar)
{
    double dist = 0;
    double diff = 0;
    
    for(int i=0; i< dim; i++)
    {
	diff = a[i]-b[i];
	if(i == Z_COORD)
	  dist += z_scalar*diff*diff;
	else 
	  dist += diff*diff;
    }
    return sqrt(dist);
};


char * Utility::inttochar(int i)
{
  /* Room for INT_DIGITS digits, - and '\0' */
  static char buf[3 + 2];
  char *p = buf + 3 + 1;       /* points to terminating '\0' */
  if (i >= 0) {
    do {
      *--p = '0' + (i % 10);
      i /= 10;
    } while (i != 0);
    return p;
  }
  else {                        /* i < 0 */
    do {
      *--p = '0' - (i % 10);
      i /= 10;
    } while (i != 0);
    *--p = '-';
  }
  return p;
};












