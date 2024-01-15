/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      basics.cpp                                                    */
/*                                                                          */
/* Purpose:   basic classes for the internal data structure                 */
/*            SpatialGraph, Edge, Vertex(deprecated)                        */
/*                                                                          */
/* Author:    Robert Egger                                                  */
/*            Max-Planck-Florida Institut                                   */
/*            5353 Parkside Drive                                           */
/*            Jupiter, Florida 33458                                        */
/*            USA                                                           */
/*                                                                          */
/* EMail:     Robert.Egger@maxplanckflorida.org                             */
/*                                                                          */
/* History:   22.12.2010                                                    */
/*                                                                          */
/* Remarks:   All rights are reserved by the Max-Planck-Society             */
/*                                                                          */
/****************************************************************************/

#include "basics.h"

//extern std::list< bouton_params * > BoutonParamsList;
extern bouton_params  BoutonParamsArray[20];
extern int NumOfLandmarks;
Vertex::Vertex(float * _coordinates, int _label)
{
	coordinates[0] = (double)_coordinates[0];
	coordinates[1] = (double)_coordinates[1];
	coordinates[2] = (double)_coordinates[2];
	label = _label;
	connectedEdges = new std::vector<Edge *>;
	edgeDirections = new std::vector<bool>;
};

Vertex::Vertex(double * _coordinates, int _label)
{
	coordinates[0] = _coordinates[0];
	coordinates[1] = _coordinates[1];
	coordinates[2] = _coordinates[2];
	label = _label;
	connectedEdges = new std::vector<Edge *>;
	edgeDirections = new std::vector<bool>;
};

Vertex::Vertex(Vertex * otherVertex)
{
	this->coordinates[0] = otherVertex->coordinates[0];
	this->coordinates[1] = otherVertex->coordinates[1];
	this->coordinates[2] = otherVertex->coordinates[2];
	this->label = otherVertex->label;
	this->connectedEdges = otherVertex->connectedEdges;
	this->edgeDirections = otherVertex->edgeDirections;
}

Vertex::~Vertex()
{
	//tbd
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = 0.0;
	fatherID = -1;
	physicalLength = 0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, double _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radius = _radius;
	fatherID = -1;
	physicalLength = 0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, std::list< double > radList )
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	edgePointCoordinates = _edgePointCoordinates;
	radiusList = radList;
	fatherID = -1;
	physicalLength = 0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = 0.0;
	fatherID = -1;
	physicalLength = 0;
};

Edge::Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius)
{
	edgeConnectivity[0] = _edgeConnectivity[0];
	edgeConnectivity[1] = _edgeConnectivity[1];
	numEdgePoints = _numEdgePoints;
	label = _label;
	std::list< float * >::iterator iter;
	for(iter = _edgePointCoordinates.begin(); iter != _edgePointCoordinates.end(); ++iter)
		edgePointCoordinates.push_back((double *)(*iter));
	radius = _radius;
	fatherID = -1;
	physicalLength = 0;
};

Edge::Edge(Edge * otherEdge)
{
	this->edgeConnectivity[0] = otherEdge->edgeConnectivity[0];
	this->edgeConnectivity[1] = otherEdge->edgeConnectivity[1];
	this->label = otherEdge->label;
	this->numEdgePoints = otherEdge->numEdgePoints;
	this->radius = otherEdge->radius;
	std::list< double >::const_iterator radiusListIt;
	for(radiusListIt = otherEdge->radiusList.begin(); radiusListIt != otherEdge->radiusList.end(); ++radiusListIt)
		this->radiusList.push_back(*radiusListIt);
	std::list< double *  >::const_iterator edgePtListIt;
	for(edgePtListIt = otherEdge->edgePointCoordinates.begin(); edgePtListIt != otherEdge->edgePointCoordinates.end(); ++edgePtListIt)
	{
		double * pt = new double [3];
		pt[0] = (*edgePtListIt)[0], pt[1] = (*edgePtListIt)[1], pt[2] = (*edgePtListIt)[2];// pt[3] = (*edgePtListIt)[3], pt[4] = (*edgePtListIt)[4];
		this->edgePointCoordinates.push_back(pt);
	}
	this->fatherID = otherEdge->fatherID;
	this->physicalLength = otherEdge->physicalLength;
}

Edge::~Edge()
{
	edgePointCoordinates.clear();
	radiusList.clear();
}

double Edge::segmentLength()
{
	if(!edgePointCoordinates.size())
		return 0;
	
	double length = 0, * lastPt = edgePointCoordinates.front();
	std::list< double * >::const_iterator ptIt = edgePointCoordinates.begin();
	++ptIt;
	while(ptIt != edgePointCoordinates.end())
	{
		double * nextPt = *ptIt;
		length += sqrt((nextPt[0] - lastPt[0])*(nextPt[0] - lastPt[0]) + (nextPt[1] - lastPt[1])*(nextPt[1] - lastPt[1]) + (nextPt[2] - lastPt[2])*(nextPt[2] - lastPt[2]));
		lastPt = nextPt;
		++ptIt;
	}
	return length;
};

AmiraSpatialGraph::AmiraSpatialGraph()
{
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
		{
			if(ii != jj)
				transformation[ii][jj] = 0;
			else
				transformation[ii][jj] = 1;
		}
	isIdentity = 1;
	transformationApplied = 0;
        inverseTransformationApplied = 0;
	homeBarrel = 0;
};

AmiraSpatialGraph::~AmiraSpatialGraph()
{
	std::vector< Edge * >::iterator it1;
	std::vector< Vertex * >::iterator it2;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
		delete *it1;
	for(it2 = vertices.begin(); it2 != vertices.end(); ++it2)
		delete *it2;
	
	vertices.clear();
	edges.clear();
};

unsigned int AmiraSpatialGraph::getNumberOfPoints()
{
	unsigned int number = 0;
	std::vector< Edge * >::iterator it1;
	for(it1 = edges.begin(); it1 != edges.end(); ++it1)
	{
		number += (*it1)->numEdgePoints;
	}
	
	return number;
}

void AmiraSpatialGraph::addVertex(Vertex * newVertex)
{
	vertices.push_back(newVertex);
};

void AmiraSpatialGraph::addEdge(Edge * newEdge)
{
	edges.push_back(newEdge);
};

void AmiraSpatialGraph::getBoundingBox(double bounds[6])
{
        double xMin = 1E09, xMax = -1E09, yMin = 1E09, yMax = -1E09, zMin = 1E09, zMax = -1E09;
        bool hasPoints = 0;
        for(int ii = 0; ii < edges.size(); ++ii)
        {
                std::list< double * >::const_iterator edgePtListIt;
                for(edgePtListIt = edges[ii]->edgePointCoordinates.begin(); edgePtListIt != edges[ii]->edgePointCoordinates.end(); ++edgePtListIt)
                {
                        hasPoints = 1;
                        double * pt = *edgePtListIt;
                        if(pt[0] < xMin)
                                xMin = pt[0];
                        if(pt[0] > xMax)
                                xMax = pt[0];
                        if(pt[1] < yMin)
                                yMin = pt[1];
                        if(pt[1] > yMax)
                                yMax = pt[1];
                        if(pt[2] < zMin)
                                zMin = pt[2];
                        if(pt[2] > zMax)
                                zMax = pt[2];
                }
        }
        
        if(hasPoints)
        {
                bounds[0] = xMin;
                bounds[1] = xMax;
                bounds[2] = yMin;
                bounds[3] = yMax;
                bounds[4] = zMin;
                bounds[5] = zMax;
        }
        else
        {
                bounds[0] = 0;
                bounds[1] = 0;
                bounds[2] = 0;
                bounds[3] = 0;
                bounds[4] = 0;
                bounds[5] = 0;
        }
}


/*
void AmiraSpatialGraph::setTransformation(double ** newTransform)
{
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			transformation[ii][jj] = newTransform[ii][jj];
	
	bool isId = 1;
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
		{
			if(ii != jj)
				if(fabs(transformation[ii][jj]) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
			
			else
				if(fabs(transformation[ii][jj] - 1) > 1E-6)
				{
					isId = 0;
// 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
				}
		}
	
	isIdentity = isId;
	transformationApplied = 0;
};*/

// void AmiraSpatialGraph::setTransformation(TransformPointerType newTransform)
// {
// 	HomogeneousMatrixPointerType mat = newTransform->GetMatrix();
// 	for(int ii = 0; ii < 4; ++ii)
// 		for(int jj = 0; jj < 4; ++jj)
// 			transformation[ii][jj] = mat->GetElement(ii, jj);
// 	
// 	bool isId = 1;
// 	for(int ii = 0; ii < 4; ++ii)
// 		for(int jj = 0; jj < 4; ++jj)
// 		{
// 			if(ii != jj)
// 				if(fabs(transformation[ii][jj]) > 1E-6)
// 				{
// 					isId = 0;
// // 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
// 				}
// 			
// 			else
// 				if(fabs(transformation[ii][jj] - 1) > 1E-6)
// 				{
// 					isId = 0;
// // 					std::cout << "T[" << ii << "][" << jj << "] = " << transformation[ii][jj] << std::endl;
// 				}
// 		}
// 	
// 	isIdentity = isId;
// 	transformationApplied = 0;
// };

void AmiraSpatialGraph::applyTransformation()
{
	if(!isIdentity && !transformationApplied)
	{
// 		printTransformation();
		transformationApplied = 1;
// 		std::flush(std::cout << "Transforming " << vertices.size() << " vertices..." << std::endl);
		std::vector< Vertex * >::iterator vertexIt;
// 		int vertexCount = 1;
		for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
		{
// 			std::cout << "Transforming Vertex " << vertexCount << " of " << vertices.size() << std::endl;
// 			++vertexCount;
			double oldCoords[4], newCoords[4];
			for(int ii = 0; ii < 3; ++ii)
			{
				oldCoords[ii] = (*vertexIt)->coordinates[ii];
				newCoords[ii] = 0;
			}
			oldCoords[3] = 1;
			newCoords[3] = 1;
			for(int ii = 0; ii < 3; ++ii)
				for(int jj = 0; jj < 4; ++jj)
					newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
			
			for(int ii = 0; ii < 3; ++ii)
				(*vertexIt)->coordinates[ii] = newCoords[ii];
		}
		
// 		std::flush(std::cout << "Transforming " << edges.size() << " edges..." << std::endl);
		std::vector< Edge * >::iterator edgeIt;
// 		int edgeCount = 1;
		for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
		{
// 			std::cout << "Transforming Edge " << edgeCount << " of " << edges.size() << std::endl;
// 			std::cout << (*edgeIt)->edgePointCoordinates.size() << " points" << std::endl;
// 			++edgeCount;
			std::list< double * >::iterator edgeListIt;
			for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
			{
				double oldCoords[4], newCoords[4];
				for(int ii = 0; ii < 3; ++ii)
				{
					oldCoords[ii] = (*edgeListIt)[ii];
					newCoords[ii] = 0;
				}
				oldCoords[3] = 1;
				newCoords[3] = 1;
				for(int ii = 0; ii < 3; ++ii)
					for(int jj = 0; jj < 4; ++jj)
						newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
				
				for(int ii = 0; ii < 3; ++ii)
					(*edgeListIt)[ii] = newCoords[ii];
			}
		}
// 		std::flush(std::cout << "done!" << std::endl);
	}
};

//  Mythreya: Modified to keep track of a point of interest after transformation 
void AmiraSpatialGraph::applyInverseTransformation(double inverse_transformation[4][4] /*std::list< bouton_params * > *boutonParamsList*/)
{
    //std::list<bouton_params*>::iterator paramIterator;
    /*for(paramIterator = BoutonParamsList.begin(); paramIterator != BoutonParamsList.end();  paramIterator++)
    {
            std::cout<<'\t'<<'\t'<<(*paramIterator)->nearestPoint[0]<<'\t'<<(*paramIterator)->nearestPoint[1]<<'\t'<<(*paramIterator)->nearestPoint[2]<<std::endl;
    }*/
    

    //double * tracked_point = new double[3];
        if(true/*!isIdentity && !inverseTransformationApplied*/)
        {
//              printTransformation();
                inverseTransformationApplied = 1;
//              std::flush(std::cout << "Transforming " << vertices.size() << " vertices..." << std::endl);
                std::vector< Vertex * >::iterator vertexIt;
//              int vertexCount = 1;
                for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
                {
//                      std::cout << "Transforming Vertex " << vertexCount << " of " << vertices.size() << std::endl;
//                      ++vertexCount;
                    
                        double oldCoords[4], newCoords[4];
                        for(int ii = 0; ii < 3; ++ii)
                        {
                                oldCoords[ii] = (*vertexIt)->coordinates[ii];
                                newCoords[ii] = 0;
                        }
                        oldCoords[3] = 1;
                        newCoords[3] = 1;
                        for(int ii = 0; ii < 3; ++ii)
                                for(int jj = 0; jj < 4; ++jj)
                                        newCoords[ii] += inverse_transformation[ii][jj]*oldCoords[jj];
                        
                        for(int ii = 0; ii < 3; ++ii)
                                (*vertexIt)->coordinates[ii] = newCoords[ii];
                }
                
                
//              std::flush(std::cout << "Transforming " << edges.size() << " edges..." << std::endl);
                std::vector< Edge * >::iterator edgeIt;
                
//                 std::cout<<"wtf begin"<<std::endl;
//              int edgeCount = 1;
                for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
                {
//                      std::cout << "Transforming Edge " << edgeCount << " of " << edges.size() << std::endl;
//                      std::cout << (*edgeIt)->edgePointCoordinates.size() << " points" << std::endl;
//                      ++edgeCount;
                        std::list< double * >::iterator edgeListIt;
                        for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
                        {
                                double oldCoords[4], newCoords[4];
                                for(int ii = 0; ii < 3; ++ii)
                                {
                                        oldCoords[ii] = (*edgeListIt)[ii];
                                        newCoords[ii] = 0;
                                }
                                oldCoords[3] = 1;
                                newCoords[3] = 1;
                                for(int ii = 0; ii < 3; ++ii)
                                        for(int jj = 0; jj < 4; ++jj)
                                                newCoords[ii] += inverse_transformation[ii][jj]*oldCoords[jj];
                                
                                for(int ii = 0; ii < 3; ++ii)
                                {
                                        (*edgeListIt)[ii] = newCoords[ii];
                                }
                                
                                // iterate thr the landmark list in the params list
                                // see if it matches the current point
                                // if so update the corresponding transformed landmark
                                
                                //std::cout<<'\t'<<'\t'<<'\t'<<"begin"<<std::endl;
                                
                                
                                
                                //std::cout<<'\t'<<'\t'<<'\t'<<"begin"<<std::endl;
                                
                                
                                /*
                                // if the this point is the point being tracked
                                // set its new coordinates as well
                                if((oldCoords[0] == point_tracked[0]) && (oldCoords[0] == point_tracked[0]) && (oldCoords[0] == point_tracked[0]))
                                {
                                    tracked_point[0] =  newCoords[0];
                                    tracked_point[1] =  newCoords[1];
                                    tracked_point[2] =  newCoords[2];
                                    
                                }*/
                        }
                }
//                 std::cout<<"wtf end"<<std::endl;
                //std::cout << "ok here 3" << std::endl;
//              std::flush(std::cou< "done!" << std::endl);
        }
        
        //return tracked_point;

};
// 
// // apply transformation only to one specific label
// void AmiraSpatialGraph::applyTransformation(int label)
// {
// 	if(!isIdentity && !transformationApplied)
// 	{
// // 		printTransformation();
// 		transformationApplied = 1;
// 		std::vector< Vertex * >::iterator vertexIt;
// 		int vertexCount = 0;
// 		for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
// 		{
// 			if((*vertexIt)->label == label)
// 			{
// 				++vertexCount;
// 				double oldCoords[4], newCoords[4];
// 				for(int ii = 0; ii < 3; ++ii)
// 				{
// 					oldCoords[ii] = (*vertexIt)->coordinates[ii];
// 					newCoords[ii] = 0;
// 				}
// 				oldCoords[3] = 1;
// 				newCoords[3] = 1;
// 				for(int ii = 0; ii < 3; ++ii)
// 					for(int jj = 0; jj < 4; ++jj)
// 						newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
// 				
// 				for(int ii = 0; ii < 3; ++ii)
// 					(*vertexIt)->coordinates[ii] = newCoords[ii];
// 			}
// 		}
// // 		std::cout << "Transforming " << vertexCount << " vertices..." << std::endl;
// 		
// 		std::vector< Edge * >::iterator edgeIt;
// 		int edgeCount = 0;
// 		for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
// 		{
// 			if((*edgeIt)->label == label)
// 			{
// 				++edgeCount;
// 				std::list< double * >::iterator edgeListIt;
// 				for(edgeListIt = (*edgeIt)->edgePointCoordinates.begin(); edgeListIt != (*edgeIt)->edgePointCoordinates.end(); ++edgeListIt)
// 				{
// 					double oldCoords[4], newCoords[4];
// 					for(int ii = 0; ii < 3; ++ii)
// 					{
// 						oldCoords[ii] = (*edgeListIt)[ii];
// 						newCoords[ii] = 0;
// 					}
// 					oldCoords[3] = 1;
// 					newCoords[3] = 1;
// 					for(int ii = 0; ii < 3; ++ii)
// 						for(int jj = 0; jj < 4; ++jj)
// 							newCoords[ii] += transformation[ii][jj]*oldCoords[jj];
// 					
// 					for(int ii = 0; ii < 3; ++ii)
// 						(*edgeListIt)[ii] = newCoords[ii];
// 				}
// 			}
// 		}
// // 		std::cout << "Transforming " << edgeCount << " edges..." << std::endl;
// 	}
// };
// 
// void AmiraSpatialGraph::printTransformation()
// {
// 	std::cout << "SpatialGraph transformation matrix:" << std::endl << std::endl;
// 	for(int ii = 0; ii < 4; ++ii)
// 	{
// 		std:: cout << "[";
// 		for(int jj = 0; jj < 3; ++jj)
// 			std::cout << transformation[ii][jj] << ", ";
// 		std::cout << transformation[ii][3] << "]" << std::endl;
// 	}
// 	std::cout << std::endl;
// 	std::cout << "isIdentity = " << isIdentity << std::endl;
// 	std::cout << "transformationApplied = " << transformationApplied << std::endl;
// };

std::vector< Vertex * >::iterator AmiraSpatialGraph::verticesBegin()
{
	std::vector< Vertex * >::iterator it = vertices.begin();
	return it;
};

std::vector< Edge * >::iterator AmiraSpatialGraph::edgesBegin()
{
	std::vector< Edge * >::iterator it = edges.begin();
	return it;
};

std::vector< Vertex * >::iterator AmiraSpatialGraph::verticesEnd()
{
	std::vector< Vertex * >::iterator it = vertices.end();
	return it;
};

std::vector< Edge * >::iterator AmiraSpatialGraph::edgesEnd()
{
	std::vector< Edge * >::iterator it = edges.end();
	return it;
};

void AmiraSpatialGraph::mergeSpatialGraph(AmiraSpatialGraph * otherSpatialGraph)
{
	unsigned int oldVertexNr = vertices.size();
	std::vector< Vertex * >::iterator vertexIt;
	for(vertexIt = otherSpatialGraph->verticesBegin(); vertexIt != otherSpatialGraph->verticesEnd(); ++vertexIt)
	{
		Vertex * newVertex = new Vertex(*vertexIt);
		this->vertices.push_back(newVertex);
	}
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = otherSpatialGraph->edgesBegin(); edgeIt != otherSpatialGraph->edgesEnd(); ++edgeIt)
	{
		Edge * newEdge = new Edge(*edgeIt);
		newEdge->edgeConnectivity[0] += oldVertexNr;
		newEdge->edgeConnectivity[1] += oldVertexNr;
		this->edges.push_back(newEdge);
	}
};

void AmiraSpatialGraph::copyEdgesWithLabel(AmiraSpatialGraph * otherSpatialGraph, int edgeLabel)
{
//         std::flush(std::cout << "Start copying edge with label " << edgeLabel << std::endl);
        unsigned int oldVertexNr = vertices.size();
//         std::flush(std::cout << "oldVertexNr: " << oldVertexNr << std::endl);
        unsigned int vertexCount = 0;
//         std::flush(std::cout << "otherSpatialGraph->getNumberOfEdges(): " << otherSpatialGraph->getNumberOfEdges() << std::endl);
        std::vector< Edge * >::iterator edgeIt;
        for(edgeIt = otherSpatialGraph->edgesBegin(); edgeIt != otherSpatialGraph->edgesEnd(); ++edgeIt)
        {
            if((*edgeIt)->label == edgeLabel)
            {
//                 std::flush(std::cout << "copying edge with label " << edgeLabel << std::endl);
//                 std::flush(std::cout << "oldVertexNr: " << oldVertexNr << std::endl);
//                 std::flush(std::cout << "vertexCount: " << vertexCount << std::endl);
                Edge * newEdge = new Edge(*edgeIt);
                int vertex1 = newEdge->edgeConnectivity[0];
                int vertex2 = newEdge->edgeConnectivity[1];
//                 std::flush(std::cout << "vertex1: " << vertex1 << std::endl);
//                 std::flush(std::cout << "vertex2: " << vertex2 << std::endl);
                if(vertex1 != vertex2)
                {
                    Vertex * newVertex1 = new Vertex(otherSpatialGraph->verticesPointer()->at(vertex1));
                    Vertex * newVertex2 = new Vertex(otherSpatialGraph->verticesPointer()->at(vertex2));
                    newEdge->edgeConnectivity[0] = oldVertexNr + vertexCount;
                    newEdge->edgeConnectivity[1] = oldVertexNr + vertexCount + 1;
//                     std::flush(std::cout << "newEdge->edgeConnectivity[0]: " << newEdge->edgeConnectivity[0] << std::endl);
//                     std::flush(std::cout << "newEdge->edgeConnectivity[1]: " << newEdge->edgeConnectivity[1] << std::endl);
                    this->vertices.push_back(newVertex1);
                    this->vertices.push_back(newVertex2);
                    vertexCount += 2;
                }
                else
                {
                    Vertex * newVertex1 = new Vertex(otherSpatialGraph->verticesPointer()->at(vertex1));
                    newEdge->edgeConnectivity[0] = oldVertexNr + vertexCount;
                    newEdge->edgeConnectivity[1] = oldVertexNr + vertexCount;
//                     std::flush(std::cout << "newEdge->edgeConnectivity[0]: " << newEdge->edgeConnectivity[0] << std::endl);
//                     std::flush(std::cout << "newEdge->edgeConnectivity[1]: " << newEdge->edgeConnectivity[1] << std::endl);
                    this->vertices.push_back(newVertex1);
                    vertexCount += 1;
                }
                this->edges.push_back(newEdge);
                
                // Copy radius information
                
            }
        }
};

// resets SpatialGraph. use with extreme caution!!! all data will be gone
void AmiraSpatialGraph::clear()
{
	std::vector< Edge * >::iterator edgeIt;
	std::vector< Vertex * >::iterator vertexIt;
	for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
		delete *edgeIt;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
		delete *vertexIt;
	edges.clear(), vertices.clear();
	homeBarrel = 0;
	isIdentity = 1;
	transformationApplied = 0;
	for(int ii = 0; ii < 4; ++ii)
		for(int jj = 0; jj < 4; ++jj)
			transformation[ii][jj] = (ii == jj) ? 1 : 0;
};
/*
void AmiraSpatialGraph::vesselsToPoints()
{
	std::vector< Edge * >::iterator edgeIt;
	std::vector< Vertex * >::iterator vertexIt;
	std::vector< Vertex * > verticesVec;
	
// 	std::flush(std::cout << "converting blood vessels to points..." << std::endl);
	
	for(vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt)
	{
		verticesVec.push_back(*vertexIt);
	}
	
	for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
	{
		if((*edgeIt)->label == Vessel)
		{
			double circ = 0;
			double * center = new double[3];
			center[0] = 0, center[1] = 0, center[2] = 0;
			double * curr;
			double * last;
			std::list< double * >::iterator it = (*edgeIt)->edgePointCoordinates.begin();
			last = *it;
			++it;
			while(it != (*edgeIt)->edgePointCoordinates.end())
			{
				curr = *it;
				circ += std::sqrt((curr[0] - last[0])*(curr[0] - last[0]) + (curr[1] - last[1])*(curr[1] - last[1]) + (curr[2] - last[2])*(curr[2] - last[2]));
				center[0] += curr[0];
				center[1] += curr[1];
				center[2] += curr[2];
				last = curr;
				++it;
			}
			center[0] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[1] /= (double)((*edgeIt)->numEdgePoints - 1);
			center[2] /= (double)((*edgeIt)->numEdgePoints - 1);
			(*edgeIt)->numEdgePoints = 1;
			(*edgeIt)->edgePointCoordinates.clear();
			(*edgeIt)->edgePointCoordinates.push_back(center);
			(*edgeIt)->radius = circ/(2*PI);
			
			for(int ii = 0; ii < 3; ++ii)
			{
				verticesVec[(*edgeIt)->edgeConnectivity[0]]->coordinates[ii] = center[ii];
				verticesVec[(*edgeIt)->edgeConnectivity[1]]->coordinates[ii] = center[ii];
			}
		}
	}
	
	int ii = 0;
	for(vertexIt = vertices.begin(); vertexIt != vertices.end() && ii < verticesVec.size(); ++vertexIt, ++ii)
	{
		for(int jj = 0; jj < 3; ++jj)
			(*vertexIt)->coordinates[jj] = verticesVec[ii]->coordinates[jj];
	}
	
	removeDuplicatePoints();
};

void AmiraSpatialGraph::removeDuplicatePoints()
{
	std::vector< Edge * >::iterator edgeIt = edges.begin();
	std::vector< Vertex * >::iterator vertexIt = vertices.begin();
	bool duplicate = 0;
	while(edgeIt != edges.end() && vertexIt != vertices.end())
	{
		bool duplicateVertex = (bool)((*edgeIt)->edgeConnectivity[0] != (*edgeIt)->edgeConnectivity[1]);
		duplicate = 0;
		
		if((*edgeIt)->label == Vessel)
		{
			std::vector< Edge * >::iterator checkIt = edgeIt;
			std::vector< Vertex * >::iterator checkIt2 = vertexIt;
			++checkIt;
			++checkIt2;
			if(duplicateVertex)
				++checkIt2;
			
			while(checkIt != edges.end() && checkIt2 != vertices.end())
			{
				if((*edgeIt)->label == Vessel)
				{
					double * tmp1 = (*edgeIt)->edgePointCoordinates.front();
					double * tmp2 = (*checkIt)->edgePointCoordinates.front();
					
					if(std::abs(tmp1[2] - tmp2[2]) > 1.0)
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
						continue;
					}
					
					float dist = std::sqrt((tmp1[0] - tmp2[0])*(tmp1[0] - tmp2[0]) + (tmp1[1] - tmp2[1])*(tmp1[1] - tmp2[1]) + (tmp1[2] - tmp2[2])*(tmp1[2] - tmp2[2]));
					if(dist < (*edgeIt)->radius || dist < (*checkIt)->radius)
					{
						if((*edgeIt)->radius < (*checkIt)->radius)
						{
							duplicate = 1;
							edgeIt = edges.erase(edgeIt);
							vertexIt = vertices.erase(vertexIt);
							if(duplicateVertex)
								vertexIt = vertices.erase(vertexIt);
							
							std::vector< Edge * >::iterator updateIt;
							for(updateIt = edgeIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
							++checkIt;
							++checkIt2;
							if(duplicateVertex)
								++checkIt2;
						}
						else
						{
							duplicate = 0;
							checkIt = edges.erase(checkIt);
							checkIt2 = vertices.erase(checkIt2);
							if(duplicateVertex)
								checkIt2 = vertices.erase(checkIt2);
							
							std::vector< Edge * >::iterator updateIt;
							for(updateIt = checkIt; updateIt != edges.end(); ++updateIt)
							{
								if(duplicateVertex)
								{
									(*updateIt)->edgeConnectivity[0] -= 2;
									(*updateIt)->edgeConnectivity[1] -= 2;
								}
								else
								{
									(*updateIt)->edgeConnectivity[0] -= 1;
									(*updateIt)->edgeConnectivity[1] -= 1;
								}
							}
						}
					}
					else
					{
						++checkIt;
						++checkIt2;
						if(duplicateVertex)
							++checkIt2;
					}
				}
				else
				{
					++checkIt;
					++checkIt2;
					if(duplicateVertex)
						++checkIt2;
				}
			}
		}
		
		if(!duplicate)
		{
			++edgeIt;
			++vertexIt;
			if(duplicateVertex)
				++vertexIt;
		}
	}
};*/

// calculate total segment length from edge/* startID to soma
// double AmiraSpatialGraph::totalSegmentLength(int startID)
// {
// 	if(startID >= edges.size() || startID < 0)
// 		return 0;
// // 	double totalLength = edges[startID]->segmentLength();
// 	double totalLength = 0;
// 	double thisVec[3], thisNorm = 0;
// 	for(int ii = 0; ii < 3; ++ii)
// 	{
// 		thisVec[ii] = edges[startID]->edgePointCoordinates.back()[ii] - edges[startID]->edgePointCoordinates.front()[ii];
// 		thisNorm += thisVec[ii]*thisVec[ii];
// 	}
// 	totalLength = sqrt(thisNorm);
// 	int fatherID = edges[startID]->fatherID;
// // 	int segmentCnt = 1;
// 	while(fatherID != -1)
// 	{
// 		int currID = fatherID;
// // 		totalLength += edges[currID]->segmentLength();
// 		thisNorm = 0;
// 		for(int ii = 0; ii < 3; ++ii)
// 		{
// 			thisVec[ii] = edges[currID]->edgePointCoordinates.back()[ii] - edges[currID]->edgePointCoordinates.front()[ii];
// 			thisNorm += thisVec[ii]*thisVec[ii];
// 		}
// 		totalLength += sqrt(thisNorm);
// 		fatherID = edges[currID]->fatherID;
// // 		++segmentCnt;
// 	}
// // 	std::cout << "nr of segments to soma = " << segmentCnt << std::endl;
// 	return totalLength;
// };
/*
double AmiraSpatialGraph::cumulatedSegmentAngle(int startID)
{
	if(startID >= edges.size() || startID < 0)
		return 0;
	double totalAngle = 0;
	double thisVec[3], thisNorm = 0;
	for(int ii = 0; ii < 3; ++ii)
	{
		thisVec[ii] = edges[startID]->edgePointCoordinates.back()[ii] - edges[startID]->edgePointCoordinates.front()[ii];
		thisNorm += thisVec[ii]*thisVec[ii];
	}
	thisNorm = sqrt(thisNorm);
	if(thisNorm)
		for(int ii = 0; ii < 3; ++ii)
			thisVec[ii] /= thisNorm;
	int fatherID = edges[startID]->fatherID;
	int segmentCnt = 1;
	while(fatherID != -1)
	{
		++segmentCnt;
		int currID = fatherID;
		double nextVec[3], nextNorm = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			nextVec[ii] = edges[currID]->edgePointCoordinates.back()[ii] - edges[currID]->edgePointCoordinates.front()[ii];
			nextNorm += nextVec[ii]*nextVec[ii];
		}
		nextNorm = sqrt(nextNorm);
		if(nextNorm)
			for(int ii = 0; ii < 3; ++ii)
				nextVec[ii] /= nextNorm;
		double angle = 0;
		for(int ii = 0; ii < 3; ++ii)
		{
			angle += thisVec[ii]*nextVec[ii];
			thisVec[ii] = nextVec[ii];
		}
		angle = acos(angle);
		totalAngle += angle;
		fatherID = edges[currID]->fatherID;
	}
	return totalAngle;
// 	return totalAngle/(double)segmentCnt;
};*/

bool AmiraSpatialGraph::isLabelInSpatialGraph(int checkLabel)
{
	for(int ii = 0; ii < edges.size(); ++ii)
		if(edges[ii]->label == checkLabel)
		{
			return 1;
		}
	return 0;
};

// //extract all points of landmark 'label' planewise; return their plane indices in zIndexList; return false if empty
// bool AmiraSpatialGraph::extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList)
// {
// 	std::vector< Edge * >::iterator edgeIt;
// 	edgeIt = this->edges.begin();
// 	if(edgeIt != this->edges.end())
// 	{
// 		int zIndex = 0;
// 		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
// 			++edgeIt;
// 		
// 		zIndex = lround((*edgeIt)->edgePointCoordinates.back()[2]);
// 		zIndexList.push_back(zIndex);
// 		std::list< double * > planeCoordinates;
// 		std::list< double * >::iterator pointListIt;
// 		while(edgeIt != this->edges.end())
// 		{
// 			if((*edgeIt)->label != label)
// 			{
// 				++edgeIt;
// 				continue;
// 			}
// 			int tmpZ = lround((*edgeIt)->edgePointCoordinates.back()[2]);
// 			// 			std::flush(std::cout << "in plane " << tmpZ << std::endl);
// 			if(tmpZ == zIndex)
// 			{
// 				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
// 				{
// 					double * tmpPoint = *pointListIt;
// // 					tmpPoint[2] = (double)(int)(tmpPoint[2] + 0.5);
// 					tmpPoint[2] = round(tmpPoint[2]);
// 					planeCoordinates.push_back(*pointListIt);
// 					planeCoordinates.back()[2] = tmpPoint[2];
// 				}
// 			}
// 			else
// 			{
// 				std::list< double * > tmpList(planeCoordinates);
// 				planewisePointList.push_back(tmpList);
// 				planeCoordinates.clear();
// 				
// 				for(pointListIt = (*edgeIt)->edgePointCoordinates.begin(); pointListIt != (*edgeIt)->edgePointCoordinates.end(); ++pointListIt)
// 				{
// 					double * tmpPoint = *pointListIt;
// // 					tmpPoint[2] = (double)(int)(tmpPoint[2] + 0.5);
// 					tmpPoint[2] = round(tmpPoint[2]);
// 					planeCoordinates.push_back(*pointListIt);
// 					planeCoordinates.back()[2] = tmpPoint[2];
// 				}
// 				zIndex = tmpZ;
// 				zIndexList.push_back(zIndex);
// 			}
// 			++edgeIt;
// 		}
// 		zIndexList.sort();
// 		planewisePointList.push_back(planeCoordinates);
// 		std::list< std::list< double * > >::const_iterator controlit;
// 		for(controlit = planewisePointList.begin(); controlit != planewisePointList.end(); ++controlit)
// 			if(controlit->size())
// 				return 1;
// 		
// 		return 0;
// 	}
// 	
// 	return 0;
// };
/*
//extract all points of landmark 'label' as PolyData; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, PolyDataPointerType polyData)
{
	if(!polyData)
		polyData = PolyDataPointerType::New();
	polyData->Allocate(1);
	PointsPointerType points = PointsPointerType::New();
	points->SetDataTypeToFloat();
	int lastID = 0;
	std::vector< Edge * >::iterator edgeIt;
	for(edgeIt = this->edges.begin(); edgeIt != this->edges.end(); ++edgeIt)
	{
		if((*edgeIt)->label != label)
			continue;
		
		int end = (*edgeIt)->edgePointCoordinates.size();
		if(label >= Landmark && label < ZAxis)	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
			--end;
		PolygonPointerType poly = PolygonPointerType::New();
		poly->GetPointIds()->SetNumberOfIds(end);
		
		std::list< double * >::iterator pointListIt;
		pointListIt = (*edgeIt)->edgePointCoordinates.begin();
		for(int ii = 0; ii < end; ++pointListIt, ++ii)
		{
			double * tmp = new double[3];
			tmp[0] = (*pointListIt)[0], tmp[1] = (*pointListIt)[1], tmp[2] = (*pointListIt)[2];
			points->InsertNextPoint(tmp);
			poly->GetPointIds()->SetId(ii, ii + lastID);
		}
		lastID += end;
		polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
	}
	polyData->SetPoints(points);
	polyData->Update();
	
	if(polyData->GetNumberOfPoints())
		return 1;
	
	return 0;
};

//extract all points of landmark 'label' as PolyData; return their plane indices in zIndexList; return false if empty
bool AmiraSpatialGraph::extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList)
{
	std::vector< Edge * >::iterator edgeIt;
	edgeIt = this->edges.begin();
	if(edgeIt != this->edges.end())
	{
		polyData->Allocate(1);
		PointsPointerType points = PointsPointerType::New();
		points->SetDataTypeToFloat();
		int lastID = 0;
		while(edgeIt != this->edges.end() && (*edgeIt)->label != label)
			++edgeIt;
		
		std::list< double * >::iterator pointListIt;
		while(edgeIt != this->edges.end())
		{
			if((*edgeIt)->label != label)
			{
				++edgeIt;
				continue;
			}
			int end = (*edgeIt)->edgePointCoordinates.size();
			if(label >= Landmark)	// vtkPolygon does NOT use the same point twice on a contour as a SpatialGraph does
				--end;
			PolygonPointerType poly = PolygonPointerType::New();
			poly->GetPointIds()->SetNumberOfIds(end);
			pointListIt = (*edgeIt)->edgePointCoordinates.begin();
			for(int ii = 0; ii < end; ++pointListIt, ++ii)
			{
				double * tmp = *pointListIt;
				points->InsertNextPoint(tmp);
				poly->GetPointIds()->SetId(ii, ii + lastID);
			}
			lastID += end;
			polyData->InsertNextCell(poly->GetCellType(), poly->GetPointIds());
			++edgeIt;
		}
		polyData->SetPoints(points);
		polyData->Update();
		
		for(int ii = 0; ii < polyData->GetNumberOfCells(); ++ii)
		{
			double * bounds = polyData->GetCell(ii)->GetBounds();
			if(lround(bounds[4]) == lround(bounds[5]))
				zIndexList.push_back(lround(bounds[4]));
		}
		zIndexList.sort();
		zIndexList.unique();
		
		if(polyData->GetNumberOfPoints())
			return 1;
			
		return 0;
	}
	
	return 0;
};

//writes all cells in object as closed loops in SpatialGraph
void AmiraSpatialGraph::addPolyDataObject(PolyDataPointerType object, int label)
{
	if(object->GetNumberOfCells())
	{
		for(int ii = 0; ii < object->GetNumberOfCells(); ++ii)
		{
			PointsPointerType cellPoints = object->GetCell(ii)->GetPoints();
			if(cellPoints->GetNumberOfPoints())
			{
				double * start = new double[3];
				double * end = new double[3];
				cellPoints->GetPoint(0, start);
				int endIndex = (label >= Landmark) ? 0 : cellPoints->GetNumberOfPoints()-1;
				cellPoints->GetPoint(0, end);
				std::list< double * > edgePts;
				for(int jj = 0; jj < cellPoints->GetNumberOfPoints(); ++jj)
				{
					double * pt = new double[3];
					cellPoints->GetPoint(jj, pt);
					edgePts.push_back(pt);
				}
				if(!endIndex)
					edgePts.push_back(end);
				int connectivity[2];
				int nrEdgePts = edgePts.size();
				connectivity[0] = getNumberOfVertices();
				connectivity[1] = getNumberOfVertices() + 1;
				Vertex * startVertex = new Vertex(start, label);
				Vertex * endVertex = new Vertex(end, label);
				Edge * newEdge = new Edge(connectivity, nrEdgePts, label, edgePts);
				addVertex(startVertex);
				addVertex(endVertex);
				addEdge(newEdge);
			}
		}
	}
};*/

void AmiraSpatialGraph::addLine(double start[3], double end[3])
{
	double * endPoint = new double[3];
	double * bottomPoint = new double[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		endPoint[ii] = end[ii];
		bottomPoint[ii] = start[ii];
	}
	Vertex * newVert1 = new Vertex(endPoint, ZAxis);
	Vertex * newVert2 = new Vertex(bottomPoint, ZAxis);
	this->addVertex(newVert2);
	this->addVertex(newVert1);
	int connectionIndex[2];
	if(!this->getNumberOfVertices())
	{
		connectionIndex[0] = 0;
		connectionIndex[1] = 1;
	}
	else
	{
		connectionIndex[0] = this->getNumberOfVertices() - 2;
		connectionIndex[1] = this->getNumberOfVertices() - 1;
	}
	std::list< double * > axisCoords;
	axisCoords.push_back(bottomPoint);
	axisCoords.push_back(endPoint);
	Edge * newAxis = new Edge(connectionIndex, 2, ZAxis, axisCoords);
	this->addEdge(newAxis);
};

void AmiraSpatialGraph::addLine(double start[3], double end[3], int ID)
{
	double * endPoint = new double[3];
	double * bottomPoint = new double[3];
	for(int ii = 0; ii < 3; ++ii)
	{
		endPoint[ii] = end[ii];
		bottomPoint[ii] = start[ii];
	}
	Vertex * newVert1 = new Vertex(endPoint, ID);
	Vertex * newVert2 = new Vertex(bottomPoint, ID);
	this->addVertex(newVert2);
	this->addVertex(newVert1);
	int connectionIndex[2];
	if(!this->getNumberOfVertices())
	{
		connectionIndex[0] = 0;
		connectionIndex[1] = 1;
	}
	else
	{
		connectionIndex[0] = this->getNumberOfVertices() - 2;
		connectionIndex[1] = this->getNumberOfVertices() - 1;
	}
	std::list< double * > axisCoords;
	axisCoords.push_back(bottomPoint);
	axisCoords.push_back(endPoint);
	Edge * newAxis = new Edge(connectionIndex, 2, ID, axisCoords);
	this->addEdge(newAxis);
};

//careful! impossible to undo
void AmiraSpatialGraph::removeLabel(int label)
{
	std::vector< Vertex * >::iterator vertexIt;
	vertexIt = this->vertices.begin();
	while(vertexIt != this->vertices.end())
	{
		if((*vertexIt)->label == label)
		{
			delete *vertexIt;
			vertexIt = this->vertices.erase(vertexIt);
		}
		else
			++vertexIt;
	}
	std::vector< Edge * >::iterator edgeIt;
	edgeIt = this->edges.begin();
	while(edgeIt != this->edges.end())
	{
		if((*edgeIt)->label == label)
		{
			delete *edgeIt;
			edgeIt = this->edges.erase(edgeIt);
		}
		else
			++edgeIt;
	}
};
/*
Column::Column()
{
	contours = NULL;
	this->top = new double[3];
	this->bottom =  new double[3];
};

Column::Column ( Column* otherColumn )
{
	this->contours = PolyDataPointerType::New();
	this->contours->DeepCopy(otherColumn->contours);
	this->top = new double[3];
	this->top[0] = otherColumn->top[0], this->top[1] = otherColumn->top[1], this->top[2] = otherColumn->top[2];
	this->bottom =  new double[3];
	this->bottom[0] = otherColumn->bottom[0], this->bottom[1] = otherColumn->bottom[1], this->bottom[2] = otherColumn->bottom[2];
}

Column::Column(PolyDataPointerType contours, double * top, double * bottom)
{
	this->contours = contours;
	this->top = new double[3];
	this->top[0] = top[0], this->top[1] = top[1], this->top[2] = top[2];
	this->bottom =  new double[3];
	this->bottom[0] = bottom[0], this->bottom[1] = bottom[1], this->bottom[2] = bottom[2];
};

Column::~Column()
{
	if(this->top)
		delete top;
	if(this->bottom)
		delete bottom;
};

void Column::getCenter ( double center[3] )
{
	center[0] = 0.5*(top[0] + bottom[0]);
	center[1] = 0.5*(top[1] + bottom[1]);
	center[2] = 0.5*(top[2] + bottom[2]);
}

void Column::translateColumn(const double * shift)
{
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] += shift[ii];
		this->bottom[ii] += shift[ii];
	}
	
	if(this->contours)
	{
		TransformFilterType transform = TransformFilterType::New();
		TransformPointerType translation = TransformPointerType::New();
		translation->Translate(shift);
		transform->SetTransform(translation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
}

void Column::rotateColumn(gsl_matrix * rot)
{
	TransformFilterType transform = TransformFilterType::New();
	TransformPointerType rotation = TransformPointerType::New();
	HomogeneousMatrixPointerType mat = HomogeneousMatrixPointerType::New();
	for(int ii = 0; ii < 3; ++ii)
		for(int jj = 0; jj < 3; ++jj)
			mat->SetElement(ii, jj, gsl_matrix_get(rot, ii, jj));
	for(int ii = 0; ii < 3; ++ii)
	{
		mat->SetElement(ii, 3, 0);
		mat->SetElement(3, ii, 0);
	}
	mat->SetElement(3, 3, 1);
	
	double hTop[4], hBottom[4];
	for(int ii = 0; ii < 3; ++ii)
	{
		hTop[ii] = this->top[ii];
		hBottom[ii] = this->bottom[ii];
	}
	hTop[3] = hBottom[3] = 1;
	mat->MultiplyPoint(hTop, hTop);
	mat->MultiplyPoint(hBottom, hBottom);
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] = hTop[ii];
		this->bottom[ii] = hBottom[ii];
	}
	
	if(this->contours)
	{
		rotation->SetMatrix(mat);
		transform->SetTransform(rotation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
};

void Column::rotateColumn(HomogeneousMatrixPointerType mat)
{
	TransformFilterType transform = TransformFilterType::New();
	TransformPointerType rotation = TransformPointerType::New();
	double hTop[4], hBottom[4];
	for(int ii = 0; ii < 3; ++ii)
	{
		hTop[ii] = this->top[ii];
		hBottom[ii] = this->bottom[ii];
	}
	hTop[3] = hBottom[3] = 1;
	mat->MultiplyPoint(hTop, hTop);
	mat->MultiplyPoint(hBottom, hBottom);
	for(int ii = 0; ii < 3; ++ii)
	{
		this->top[ii] = hTop[ii];
		this->bottom[ii] = hBottom[ii];
	}
	
	if(this->contours)
	{
		rotation->SetMatrix(mat);
		transform->SetTransform(rotation);
		transform->SetInput(this->contours);
		transform->Update();
		this->contours->DeepCopy(transform->GetOutput());
		this->contours->Update();
	}
};*/
/*
Surface::Surface(PolyDataPointerType mesh)
{
	dataValid = 0;
	intersectionFound = 0;
	surfaceMesh = PolyDataPointerType::New();
	surfaceMesh->Allocate(1);
	surfaceMesh->DeepCopy(mesh);
	locator = CellLocatorPointerType::New();
	locator->AutomaticOn();
	locator->CacheCellBoundsOn();
	locator->SetDataSet(surfaceMesh);
	locator->BuildLocator();
};

void Surface::intersectLine(double * axis, double * center)
{
	double * pt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*4000;
		a1[jj] -= axis[jj]*4000;
	}
	int intersection = locator->IntersectWithLine(a0, a1, tol, t, pt, pcoords, subId, cellID, intersectCell);
	dataValid = 1;
	if(intersection)
	{
		intersectPt = pt;
		intersectID = cellID;
		intersectionFound = 1;
	}
	else
	{
		delete [] pt;
		intersectPt = NULL;
		intersectID = -1;
		intersectionFound = 0;
	}
};

void Surface::intersectLineInDirection ( double* axis, double* center )
{
	double * pt = new double[3];
	double a0[3], a1[3], tol = 1E-03, t, pcoords[3];
	int subId;
	vtkIdType cellID;
	GenericCellPointerType intersectCell = GenericCellPointerType::New();
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] = center[jj];
		a1[jj] = center[jj];
	}
	for(int jj = 0; jj < 3; ++jj)
	{
		a0[jj] += axis[jj]*4000;
	}
	int intersection = locator->IntersectWithLine(a0, a1, tol, t, pt, pcoords, subId, cellID, intersectCell);
	dataValid = 1;
	if(intersection)
	{
		intersectPt = pt;
		intersectID = cellID;
		intersectionFound = 1;
	}
	else
	{
		delete [] pt;
		intersectPt = NULL;
		intersectID = -1;
		intersectionFound = 0;
	}
}

double * Surface::getLastIntersectPoint()
{
	if(dataValid)
	{
		if(intersectPt)
		{
			double * pt = new double[3];
			pt[0] = intersectPt[0], pt[1] = intersectPt[1], pt[2] = intersectPt[2];
			return pt;
		}
		return NULL;
	}
	else
		return NULL;
};

void Surface::getLastIntersectPoint(double pt[3])
{
	if(dataValid && intersectPt)
		pt[0] = intersectPt[0], pt[1] = intersectPt[1], pt[2] = intersectPt[2];
};

vtkIdType Surface::getLastIntersectCellID()
{
	if(dataValid)
	{
		return intersectID;
	}
	return 0;
};

ClosedSurface::ClosedSurface ( PolyDataPointerType mesh ) : Surface ( mesh )
{
	insideSurfaceFilter = SelectEnclosedPointsFilterType::New();
	insideSurfaceFilter->Initialize(mesh);
	insideSurfaceFilter->CheckSurfaceOn();
// 	insideSurfaceFilter->Print(std::cout);
}

ClosedSurface::~ClosedSurface()
{
	if(insideSurfaceFilter) insideSurfaceFilter->Complete();
}

// use vtkSelectEnclosedPoints filter
bool ClosedSurface::isPointInsideSurface ( double pt[3] )
{
	return insideSurfaceFilter->IsInsideSurface(pt);
}


*/








