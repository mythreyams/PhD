
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
      ////////std::cout << i << " "<< mode_array[i] << std::endl;
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

// this function appends given data to the end of given row
void Utility::appendToXls(const char * inputxlsfile, const char * outputxlsfile, std::vector<int> *data, int section_number)
{
    ////std::cout<<"inputxlsfile "<<inputxlsfile<<std::endl;
    ////std::cout<<"outputxlsfile "<<outputxlsfile<<std::endl;
    std::ifstream xlsreader(inputxlsfile);//xlsreader.open(inputxlsfile, std::fstream::in | std::fstream::out);
    std::ofstream xlsWriter(outputxlsfile);
    std::string currentLine = "";
    int line_number = 0;
    
    //xlsWriter.open (xlsfilename,  std::fstream::in |  std::fstream::app );
    
    if(!xlsreader.fail() && ! xlsWriter.fail())
    {
        ////std::cout<<"opened fine "<<inputxlsfile<<std::endl;
	
	while(!std::getline(xlsreader, currentLine).eof())
	{
	  //xlsWriter = std::getline(xlsWriter, currentLine);
	   
	  
	  if( foundInVector(data, line_number))
	  {
	    //////std::cout<<"out of while "<<xlsfilename<<std::endl;
	    ////std::cout<<"line_number "<<line_number<<std::endl;
	    ////std::cout<<"currentLine.length() "<<currentLine.length()<<std::endl;
	    ////std::cout<<"currentLine.lasof() "<<currentLine.at(currentLine.length()-1)<<std::endl;
	    //////std::cout<<"data "<<data<<" "<<Utility::inttochar(data)<<std::endl;
	    ////std::cout<<"currentLine "<<currentLine<<std::endl;
         
	    currentLine.insert(currentLine.length(), ",");
	    currentLine.insert(currentLine.length(), Utility::inttochar(section_number));
	    ////std::cout<<"after insert "<<currentLine<<std::endl;
	    //xlsWriter<<currentLine;
	  }
	  
	  xlsWriter<<currentLine<<'\n';
	  
	  line_number++;
	  
	}
	
        
        
    }
    // write header if it does not exists
    
    xlsreader.close();
    xlsWriter.close();
    
    
};

bool Utility::foundInVector(  std::vector<int> *data, int line_number)
{
   for(int i = 0; i < data->size(); i++)
   {
     
     if(data->at(i) == line_number)
     {
        return true;
     }
     
     
   }
   
   return false;
  
  
};

// Applies given tranformation to the given point
void Utility::applyTransformationToPoint(double * point, double **given_transformation)
{
    
    double oldCoords[4], newCoords[4];
    for(int ii = 0; ii < 3; ++ii)
    {
            oldCoords[ii] = point[ii];
            newCoords[ii] = 0;
    }
    oldCoords[3] = 1;
    newCoords[3] = 1;
    for(int ii = 0; ii < 3; ++ii)
            for(int jj = 0; jj < 4; ++jj)
                    newCoords[ii] += given_transformation[ii][jj]*oldCoords[jj];
    
    for(int ii = 0; ii < 3; ++ii)
            point[ii] = newCoords[ii];
    
};

bool Utility::Match3DPoints(double * point1, double * point2, int precision)
{
  // accuracy tells how many decimal places to match
  
  for(int i = 0; i < 3; i++)
  {
      point1[i] *= pow(10, precision);
      point1[i] = round(point1[i]);
      
      point2[i] *= pow(10, precision);
      point2[i] = round(point2[i]);
      
      if (point1[i] != point2[i])
      {
	return false;
      }
      
  }
  
  return true;
      

};








