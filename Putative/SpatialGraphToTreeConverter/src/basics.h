/****************************************************************************/
/*                                                                          */
/* Program:   CortexCoordinates                                             */
/*                                                                          */
/* File:      basics.h                                                      */
/*                                                                          */
/* Purpose:   basic classes for the internal data structure                 */
/*            SpatialGraph, Edge, Vertex(deprecated) etc.                   */
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

#pragma once
#include "typedefs.h"
//#include "typedefs_two.h"

#ifndef BASICS
#define BASICS



class Edge
{
	public:
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< float * > _edgePointCoordinates, float _radius);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, double _radius);
		Edge(int * _edgeConnectivity, int _numEdgePoints, int _label, std::list< double * > _edgePointCoordinates, std::list< double > radList );
		Edge(Edge * otherEdge);
		Edge(Edge * otherEdge, int edgeConnectivity1, int edgeConnectivity2 );
		
		~Edge();
// 	private:
		int edgeConnectivity[2];
		int numEdgePoints;
		double physicalLength;
		int label;
		std::list< double * > edgePointCoordinates;
		double radius;
		std::list< double > radiusList;
		// for hoc file structure
		int fatherID;
        int childID;
        int level;
        bool forward;
		double segmentLength();
        double segmentLength(double *,bool );
        
};


class Vertex
{
	public:
		Vertex(float * _coordinates, int _label);
		Vertex(double * _coordinates, int _label);
		Vertex(Vertex * otherVertex);
		Vertex(Vertex * otherVertex, double _distanceFromSoma);
		Vertex(Vertex * otherVertex, bool copyConnectedEdge);
		~Vertex();
		
		std::vector<Edge *> * getConnectedEdges(){return connectedEdges;}
		
// 	private:
		double coordinates[3];
		int label;
		std::vector<Edge *> * connectedEdges;
		std::vector<bool> * edgeDirections;
		double distanceFromSoma;
        double distanceFromSomaEuclidean;
        int parentVertexID;
        int vertexLevel;
};

class AmiraSpatialGraph
{
	public:
		AmiraSpatialGraph();
		~AmiraSpatialGraph();
		
		// appends otherSpatialGraph to this SpatialGraph
		void mergeSpatialGraph(AmiraSpatialGraph * otherSpatialGraph);
                // appends edges from otherSpatialGraph with label edgeLabel to this SpatialGraph
                void copyEdgesWithLabel(AmiraSpatialGraph * otherSpatialGraph, int edgeLabel);
		// clears all data from internal data storage
		// use with caution!
		void clear();
		
		void addVertex( Vertex * newVertex );
		void addEdge( Edge * newEdge );
        void addAxonID( int axonID);
        void addBasalDendriteID( int dendriteID);
        void addApicalDendriteID( int dendriteID);
		void addSomaID( int somaID);
        std::vector<int>* getAxonIDArray( void);
        std::vector<int>* getBasalDendriteIDArray( void);
        std::vector<int>* getApicalDendriteIDArray( void);
        std::vector<int>* getSomaIDArray( void);
        int getNumberofCells( void);
        int getDendriteCellIndex( int ID);
        int getAxonCellIndex( int ID);
		int getSomaCellIndex( int ID);
        int getDendriteOrSomaCellIndex( int ID);
        int getAxonDendriteOrSomaCellIndex( int ID);
                
		//void addPolyDataObject(PolyDataPointerType object, int label);
		void addLine(double start[3], double end[3]);
		void addLine(double start[3], double end[3], int ID);
		// clears all structures of type 'label' from internal data storage
		// use with caution!
		void removeLabel(int label);
		
		// sets internal transformation matrix
		// from double[4][4]
// 		void setTransformation(double ** newTransform);
// 		// sets internal transformation matrix
// 		// from vtkTransform
// 		void setTransformation(TransformPointerType newTransform);
// 		// applies transformation matrix to all structures.
// 		// can only be applied once; to apply more than once,
// 		// set transformation before each use of this method
 		void applyTransformation();
        void applyInverseTransformation(double ** inverse_transformation);
// 		// applies transformation matrix to all structures of type 'label'.
// 		// can only be applied once; to apply more than once,
// 		// set transformation before each use of this method
// 		void applyTransformation(int label);
// 		
// 		void printTransformation();
		
		std::vector< Vertex * >::iterator verticesBegin();
		std::vector< Edge * >::iterator edgesBegin();
		std::vector< Vertex * >::iterator verticesEnd();
		std::vector< Edge * >::iterator edgesEnd();
		std::vector< Vertex * > * verticesPointer(){return &vertices;}
		std::vector< Edge * > * edgesPointer(){return &edges;}
		
		bool isLabelInSpatialGraph(int checkLabel);
		// DEPRECATED: extract all points of landmark 'label' planewise;
		// return their plane indices in zIndexList;
		// sanity check: returns false if no points found
// 		bool extractLandmark(int label, std::list< std::list< double * > >& planewisePointList, std::list< int >& zIndexList);
// 		// stores structures of type 'label' in polyData and their plane indices in zIndexList
// 		// sanity check: returns false if no points found
// 		bool extractLandmark(int label, PolyDataPointerType polyData, std::list< int >& zIndexList);
// 		// stores structures of type 'label' in polyData
// 		// preferred method to interface SpatialGraph structures with VTK
// 		// sanity check: returns false if no points found
// 		bool extractLandmark(int label, PolyDataPointerType polyData);
		
		unsigned int getNumberOfVertices() { return vertices.size(); }
		unsigned int getNumberOfEdges() { return edges.size(); }
		unsigned int getNumberOfPoints();
        void getBoundingBox(double bounds[6]);
		double getTotalGraphLength() {return total_graph_length;}
		void setTotalGraphLength(double length) {total_graph_length = length;}
		
// 		void setHomeBarrel(int ID){this->homeBarrel = ID;}
// 		int getHomeBarrel(){return this->homeBarrel;}
// 		// calculate total segment length from edge startID to soma
// 		// calculate total segment angle from edge startID to soma
// 		double totalSegmentLength(int startID);
// 		double cumulatedSegmentAngle(int startID);
// 		
// 		void vesselsToPoints();
	
	private:
// 		unsigned int numberOfVertices, numberOfEdges, numberOfPoints;
		
		std::vector< Vertex * > vertices;
		std::vector< Edge * > edges;
		std::vector< int > AxonIDs;
        std::vector< int > BasalDendriteIDs;
        std::vector< int > ApicalDendriteIDs;
		std::vector< int > SomaIDs;
                
		int homeBarrel;
		double total_graph_length;
		
		bool isIdentity;	//avoid going through all points if transformation == 1
		bool transformationApplied, inverseTransformationApplied;	//flag to ensure transformation is at most applied once
		double transformation[4][4];
                double inverse_transformation[4][4];
                
		
		void removeDuplicatePoints();
};

// contains top and bottom center pts as double *
// // top and bottom contours in vtkPolyData in order top, bottom
// class Column
// {
// 	public:
// 		Column();
// 		Column(Column * otherColumn);
// 		Column(PolyDataPointerType contours, double * top, double * bottom);
// 		~Column();
// 		
// 		double getHeight(){return sqrt((top[0] - bottom[0])*(top[0] - bottom[0]) + (top[1] - bottom[1])*(top[1] - bottom[1]) + (top[2] - bottom[2])*(top[2] - bottom[2]));}
// 		void getCenter ( double center[3] );
// 		void translateColumn(const double * shift);
// 		void rotateColumn(gsl_matrix * rot);
// 		void rotateColumn(HomogeneousMatrixPointerType mat);
// 		
// 		PolyDataPointerType contours;
// 		double * top;
// 		double * bottom;
// };

// class Surface
// {
// 	public:
// 		Surface(PolyDataPointerType mesh);
// 		virtual ~Surface() {};
// 		
// 		void intersectLine(double * axis, double * center);
// 		void intersectLineInDirection(double * axis, double * center);
// 		
// 		double * getLastIntersectPoint();
// 		void getLastIntersectPoint(double pt[3]);
// 		
// 		vtkIdType getLastIntersectCellID();
// 		bool isValid(){return dataValid;}
// 		bool isIntersectionFound(){return intersectionFound;}
// 		
// 		PolyDataPointerType ptr(){return surfaceMesh;}
// // 		void setSurfaceMesh(PolyDataPointerType mesh){surfaceMesh = mesh;}
// 		
// 	protected:
// 	       PolyDataPointerType surfaceMesh;
// 	       CellLocatorPointerType locator;
// 	       
// 	       bool dataValid, intersectionFound;
// 	       double * intersectPt;
// 	       vtkIdType intersectID;
// };
/*
class ClosedSurface : public Surface
{
	public:
		ClosedSurface ( PolyDataPointerType mesh );
		~ClosedSurface();
		
		bool isPointInsideSurface ( double pt[3] );
		
	private:
		SelectEnclosedPointsFilterType insideSurfaceFilter;
};
*/

#endif
