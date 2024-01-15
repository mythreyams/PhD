# hocclasses.py

import math
from decimal import Decimal

############################################################################################
############################################################################################
############################################################################################
class Segment:
    points        = None
    segmentName   = None
    father        = None
    sons          = None
    label         = None
    # numpoints     = 0    # Not used?
    # branchorder   = -1   # Not used?

    def __init__(self):
        
        # initializes the data members
        self.points      = None
        self.segmentName = None
        self.father      = None
        self.sons        = None
        self.label       = None
        # self.numpoints   = 0
        # self.branchorder = -1

    def length(self):
        
        # Compute segment length using Euclidean distance between each point
        # and adding up all the sections
        
        x = [row[1] for row in self.points]
        y = [row[2] for row in self.points]
        z = [row[3] for row in self.points]
        r = [row[4] for row in self.points]

        length       = 0
        surface_area = 0
        polarR       = 0
        
##         print "x = ", x
##         print "y = ", y
##         print "z = ", z

##         print
##         print "max(x), min(x) = ", max(x), min(x)
##         print "max(y), min(y) = ", max(y), min(y)
##         print "max(z), min(z) = ", max(z), min(z)
##         print

        # Euclidean length of whole segment

        Euclidlength = math.sqrt(
            (max(x)-min(x))**2 +
            (max(y)-min(y))**2 +
            (max(z)-min(z))**2)

        # print "Euclidean Length     = ", Euclidlength

##         print
        
##         print "i, i+1, section length, range", range(len(self.points)-1)

        for i in range(len(self.points)-1):

            length = length + math.sqrt(
                (x[i]-x[i+1])**2 +
                (y[i]-y[i+1])**2 +
                (z[i]-z[i+1])**2
                )
            
            surface_area = surface_area + \
                           2*math.pi*(abs(r[i+1]+r[i])/2)*math.sqrt(
                (x[i]-x[i+1])**2 +
                (y[i]-y[i+1])**2 +
                (z[i]-z[i+1])**2
                )

##             polarR(i) = math.sqrt(
##                 x[i]**2 +
##                 y[i]**2 +
##                 z[i]**2
##                 )
        
            
##             print i, i+1, x[i], x[i+1], "   =", math.sqrt(
##                 (x[i]-x[i+1])**2 +
##                 (y[i]-y[i+1])**2 +
##                 (z[i]-z[i+1])**2
##                 )
            
        # r_ave = sum(r)/len(r)
        # surface_area = 2*math.pi*r_ave*length
        # print "Trace Length         = ", length
        
        #return round(float(length), 3)
        #decimal.getcontext().prec = 4
        #return decimal.Decimal(float(length))
        return [float(length), float(surface_area)]
        # return [float(length), float(polarR)]

    def printSegment(self):
        
        print
        print "*********************************"
        print "Segment Characteristics:"
        print "SegmentName          = ", self.segmentName
        print "Number of points     = ", len(self.points)
        # self.numpoints = len(self.points)
        # print "Number of points (numpoints)     = ", self.numpoints
        # print "Range of points     = ", range(1,(len(self.points)+1))
        print "father               = ", self.father
        print "sons                 = ", self.sons
        print "length               = ", self.length()
        # print self.length()
        print "label                = ", self.label
        print "Points:                 "
        for point in self.points:
            print "  ", point
        # print "*********************************"
        print
        

############################################################################################
############################################################################################
############################################################################################
class Neuron(Segment):

    name            = None
    segments        = None
    segmentLengths  = None
    segment_counter = 0
    transMatrix     = 0


    def __init__(self):
        
        self.name            = None
        self.segments        = None
        self.segmentLengths  = None
        self.segment_counter = None
    
    def printCell(self):
        
        print
        print "*********************************"
        print "Cell Characteristics:"
        print
        print "Name                                        = ", self.name
        print "Number of Segments                          = ", self.segment_counter
        print
        print "Fathers                                     = ", self.father
        print
        print "Number of Fathers                           = ", len(self.father)
        print
        print "Sons                                        = ", self.sons
        print
        print "Number of Sons                              = ", len(self.sons)
        print
        print "Segments                                    = ", self.segments
        print
        keysofpoints = self.points.keys()
        # print "Segments as keys of points  = ", keysofpoints
        print
        # sortedkeys = keysofpoints.sort()
        # print "Segments as keys of points (sorted) = ", sortedkeys
        # print "Points                = ", self.points
        print "Segment Lengths and surface areas            = ", self.segmentLengths
        print
        print "Number of Points per Segment                 = ", self.numpoints
        print
        
        totalcelllength = 0
        totalcellsurfarea = 0
        # Had to dink around with for axon files.  Ended up setting
        # self.segment_counter to self.segment_counter-1
        #
        # print range(self.segment_counter)
        # print "self.segmentLength = ", self.segmentLengths[1117][1]
        for i in range(self.segment_counter-1):
                       totalcelllength = totalcelllength + self.segmentLengths[i][0]
                       totalcellsurfarea = totalcellsurfarea + self.segmentLengths[i][1]
                       
        print "Total Number of points                       = ", sum(self.numpoints)
        print
        print "Total Cell Length                            = ", totalcelllength
        print
        print "Total Cell Surface Area                      = ", totalcellsurfarea
        print
        print "*********************************"
        print

