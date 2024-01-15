/****************************************************************************/
/*                                                                          */
/* File:      edge_wrapper.h                                                */
/*                                                                          */
/* Purpose:   Defines and implements a the methods of the EdgeWrapper class */
/*                                                                          */
/*                                                                          */
/* Author:    Christopher Tull                                              */
/*            Max-Planck-Institut fuer biologische Kybernetik               */
/*                                                                          */
/*                                                                          */
/* EMail:     christopher.tull@tuebingen.mpg.de                             */
/*                                                                          */
/* History:   25.03.2014                                                    */
/*                                                                          */
/* Remarks:   Serves as a wrapper for the previously implemented Edge class.*/
/*            This greatly simplifies the interaction with the Edge class   */
/*            by using iterator-style hasNext() and getNext() methods       */
/*            This also serves to improve readability.                      */
/*                                                                          */
/****************************************************************************/


#ifndef EDGEWRAPPER
#define EDGEWRAPPER

#include "typedefs.h"
#include "typedefs_two.h"
#include "basics.h"

//#define DEBUG

//****************************

class EdgeWrapper
{
 public:
    int edgeID;
     
    EdgeWrapper(Edge * edge, int ID)
    {
        this->edgeID = ID;
        this->edge = edge; 
        setToBegin();
    };
    

    void setToBegin()
    {
        edge_it = edge->edgePointCoordinates.begin();
    };
     
    bool hasNext()
    {
        if(edge_it == edge->edgePointCoordinates.end())
            return false;
        else
            return true;
    };
    
    double * getNext()
    {
        double * currentPoint = *edge_it;
        edge_it++;
        
        return currentPoint;  
    };
    
    bool hasPrevious()
    {
        if(edge_it == edge->edgePointCoordinates.begin())
            return false;
        else
            return true;
    };
    
    double * getPrevious()
    {
        edge_it--;
        double * currentPoint = *edge_it;
        
        return currentPoint;  
    };
    
    double * getEdgePointCoordinates(){
        
        double * currentPoint = *edge_it;
        
        return currentPoint;
    }
    
    Edge * getEdge(){
        return this->edge;        
    }
    
    int  getID(){
        return this->edge->label;
    }

  
 private:
    Edge * edge;
    std::list< double * >::iterator edge_it;
     
};

  


#endif

