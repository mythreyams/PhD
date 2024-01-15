/****************************************************************************/
/*                                                                          */
/* File:      bounding_box.h                                                */
/*                                                                          */
/* Purpose:   Header file containing declarations and implementation of     */
/*            the methods for the BoundingBox class                         */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   26.03.2014                                                    */
/*                                                                          */
/* Remarks:    Defines a 3d voxel box and contains methods for updating the */
/*             coordinates of this box. Each new point that lies ouside the */
/*             box expands the coordinates. Points that already lie inside  */
/*             the box have no affect.                                      */
/*                                                                          */
/****************************************************************************/


#ifndef BOUNDINGBOX
#define BOUNDINGBOX

#include "typedefs.h"

#define INFTY 10000000


//****************************


class BoundingBox
{
    
public:
    BoundingBox()
    {
        this->max_coords[X_COORD] = 0;
        this->max_coords[Y_COORD] = 0;
        this->max_coords[Z_COORD] = 0;
        
        this->min_coords[X_COORD] = INFTY;
        this->min_coords[Y_COORD] = INFTY;
        this->min_coords[Z_COORD] = INFTY;
    };
    
    //scales the coordinates into voxel coordinates before expanding
    void expandBox_amiraCoords(double * coords)
    {
         int x = rint(coords[X_COORD]/XYSAMPLING);
         int y = rint(coords[Y_COORD]/XYSAMPLING);
        int z = rint(coords[Z_COORD]/ZSAMPLING);
        
        expandBox(x,y,z);
    };
    
    void expandBox_voxelCoords( int coords[])
    {
        expandBox(coords[X_COORD], coords[Y_COORD], coords[Z_COORD]);
    };   
    
    unsigned int X_dim()                        //returns the width of the box in the x dimension
    {
            int x_dim = max_coords[X_COORD] - min_coords[X_COORD];
            
            if(x_dim > 0)
                return x_dim;
            else
                return 0;
    }
    
    unsigned int Y_dim()                        //returns the width of the box in the y dimension
    {
            int y_dim = max_coords[Y_COORD] - min_coords[Y_COORD];
            
            if(y_dim > 0)
                return y_dim;
            else
                return 0;
    }
    
    unsigned int Z_dim()                        //returns the width of the box in the z dimension
    {
            int z_dim = max_coords[Z_COORD] - min_coords[Z_COORD];
            
            if(z_dim > 1)
                return z_dim;
            else
                return 1;
    }
    
    unsigned int X_dim_Amira()                        //returns the width of the box in the x dimension
    {
            int x_dim = (max_coords[X_COORD] - min_coords[X_COORD])*XYSAMPLING;
            
            if(x_dim > 0)
                return x_dim;
            else
                return 0;
    }
    
    unsigned int Y_dim_Amira()                        //returns the width of the box in the y dimension
    {
            int y_dim = (max_coords[Y_COORD] - min_coords[Y_COORD])*XYSAMPLING;
            
            if(y_dim > 0)
                return y_dim;
            else
                return 0;
    }
    
    unsigned int Z_dim_Amira()                        //returns the width of the box in the z dimension
    {
            int z_dim = (max_coords[Z_COORD] - min_coords[Z_COORD])*ZSAMPLING;
            
            if(z_dim > 1)
                return z_dim;
            else
                return 1;
    }
    
    ImageType::IndexType getStartIndex()        //used externally to define a region based on the bounding box
    {
        ImageType::IndexType start;
        start[0] =   min_coords[0];  // first index on X
        start[1] =   min_coords[1];  // first index on Y
        start[2] =   min_coords[2];  // first index on Z
        
        return start;
    };
    
    double * getCenter_amiraCoords()
    {
        double * center_coords = (double *)calloc(3, sizeof(double));
        center_coords[X_COORD] = ( ( (double)(max_coords[X_COORD] + min_coords[X_COORD])) / 2 )*XYSAMPLING;
        center_coords[Y_COORD] = ( ( (double)(max_coords[Y_COORD] + min_coords[Y_COORD])) / 2 )*XYSAMPLING;
        center_coords[Z_COORD] = ( ( (double)(max_coords[Z_COORD] + min_coords[Z_COORD])) / 2 )*ZSAMPLING;
        
        return center_coords;
    }
    
      int * getCenter_voxelCoords()
    {
          int * center_coords = ( int *)calloc(3, sizeof( int));
        center_coords[X_COORD] = ( (max_coords[X_COORD] + min_coords[X_COORD]) / 2 );
        center_coords[Y_COORD] = ( (max_coords[Y_COORD] + min_coords[Y_COORD]) / 2 );
        center_coords[Z_COORD] = ( (max_coords[Z_COORD] + min_coords[Z_COORD]) / 2 );
        
        return center_coords;
    }
    
    void getBoundingBox(double bounds[6]){
        
        bounds[0] = this->max_coords[X_COORD]*XYSAMPLING;
        bounds[1] = this->max_coords[Y_COORD]*XYSAMPLING;
        bounds[2] = this->max_coords[Z_COORD]*ZSAMPLING;
        bounds[3] = this->min_coords[X_COORD]*XYSAMPLING;
        bounds[4] = this->min_coords[Y_COORD]*XYSAMPLING;
        bounds[5] = this->min_coords[Z_COORD]*ZSAMPLING;
        
        
    }
    
    
private:
    int max_coords[3];
    int min_coords[3];
    
    void expandBox(int x, int y, int z)
    {
        if(x > this->max_coords[X_COORD])
        {
            this->max_coords[X_COORD] = x;
        }
        if(y > this->max_coords[Y_COORD])
        {
            this->max_coords[Y_COORD] = y;
        }
        if(z > this->max_coords[Z_COORD])
        {
            this->max_coords[Z_COORD] = z;
        }
        
        if(x < this->min_coords[X_COORD] /*&& x >= 0*/ )
        {
            this->min_coords[X_COORD] = x;
        }
        if(y < this->min_coords[Y_COORD] /*&& y >= 0*/  )
        {
            this->min_coords[Y_COORD] = y;
        }
        if(z < this->min_coords[Z_COORD] /*&& z >= 0*/ )
        {
            this->min_coords[Z_COORD] = z;
        }
    };
};
  


#endif

