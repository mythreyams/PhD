#import pandas as pd
import os
from enum import Enum
#from anytree import Node, RenderTree
import numpy as np
from scipy import spatial
#import vtk
import glob
import sys
import pathlib
import pandas as pd
#import networkx as nx
import json
from scipy.spatial import distance
from spatial_graphs.Vectors import Vectors

class Edge:
    Id = 0
    color = [0,0,0]
    nodes = []
    pts = []
    def __init__(self,filename):
        self.Id = []
        self.nodes = []
        self.pts = []


class EdgeList(Edge):
    edge_list = []
    def __init__(self,filename='',):
        self.edge_list = [Edge.__init__(self)]
        

class Contours(EdgeList):
    contours = []
    Id = 0
    color = [0,0,0]
    '''Read the contours from the spatial graph which have one of the given labels'''
    def __init__(self,filename,label=[]):
        self.contours = [EdgeList.__init__(self)]
        

class GraphDecorators():
    '''
    This is a data structure for string the '@' decorators of the amira graph
    '''
    def __init__(self):
        self.decorator_invalid = 0
        self.decorator_vertices = 1
        self.decorator_vertex_labels = 2
        self.decorator_edge_connectivity = 3
        self.decorator_num_pts_per_edge = 4
        self.decorator_edge_labels = 5
        self.decorator_edge_pt_coords = 6
        self.decorator_edge_pt_radii = 7


class SpatialGraphData():
    '''
    This is a data structure class, used for labelled subgraph of a spatial graph(leaf nodes of the label tree)
    such as : axon, apical_dend, pia, wm
    '''
    def __init__(self,vertices=[],vertex_labels=[],edge_connectivity=[],num_pts_per_edge=[],edge_labels=[],edge_pt_coords=[],edge_pt_radii=[],edge_list=[],transformation_matrix = [],label_name='',label_id=0,label_color=[0,0,0]):
        self.vertices = []#vertices
        self.vertex_labels = []#vertex_labels
        self.edge_connectivity = []#edge_connectivity
        self.num_pts_per_edge = []#num_pts_per_edge
        self.edge_labels = []#edge_labels
        self.edge_pt_coords = []#edge_pt_coords
        self.edge_pt_radii = []#edge_pt_radii
        self.transformation_matrix = []#transformation_matrix
        self.edge_list = []#edge_list
        self.length = 0
        self.vert_id = 0 # this should be used for creating graph from scatch, with continous vert numbers
        self.label_name = ''#label_name
        self.label_id = 0#label_id
        self.label_color = [0,0,0]#label_color
        
    def set_label(self,new_label):
        for i in range(len(self.vertices)):
            self.vertex_labels[i] = new_label
        for i in range(len(self.edge_list)):
            self.edge_labels[i] = new_label
            
    def apply_scaling(self,scaling):
        for i in range(len(self.vertices)):
            pt = self.vertices[i]
            self.vertices[i] = [pt[0]*scaling[0],pt[1]*scaling[1],pt[2]*scaling[2]]
        for i in range(len(self.edge_pt_coords)):
            pt = self.edge_pt_coords[i]
            self.edge_pt_coords[i] = [pt[0]*scaling[0],pt[1]*scaling[1],pt[2]*scaling[2]]
            self.edge_pt_radii[i] = 0
         
    def set_z_coord(self,z):
        for i in range(len(self.vertices)):
            pt = self.vertices[i]
            if len(pt)==3:
                self.vertices[i] = [pt[0],pt[1],z]
            elif len(pt)==4:
                self.vertices[i] = [pt[0],pt[1],z,pt[3]]
        for i in range(len(self.edge_pt_coords)):
            pt = self.edge_pt_coords[i]
            if len(pt)==3:
                self.edge_pt_coords[i] = [pt[0],pt[1],z]
            elif len(pt)==4:
                self.edge_pt_coords[i] = [pt[0],pt[1],z,pt[3]]

    def copy_spatial_graph_data(self,graph_data_to_be_copied):
        self.vertices = []#vertices
        self.vertex_labels = []#vertex_labels
        self.edge_connectivity = []#edge_connectivity
        self.num_pts_per_edge = []#num_pts_per_edge
        self.edge_labels = []#edge_labels
        self.edge_pt_coords = []#edge_pt_coords
        self.edge_pt_radii = []#edge_pt_radii
        self.transformation_matrix = []#transformation_matrix
        self.edge_list = []#edge_list
        self.length = 0
        
        for i in range(len(graph_data_to_be_copied.vertices)):
            self.vertices.append(graph_data_to_be_copied.vertices[i])
        for i in range(len(graph_data_to_be_copied.vertex_labels)):
            self.vertex_labels.append(graph_data_to_be_copied.vertex_labels[i])
        for i in range(len(graph_data_to_be_copied.edge_connectivity)):
            self.edge_connectivity.append(graph_data_to_be_copied.edge_connectivity[i])
        for i in range(len(graph_data_to_be_copied.num_pts_per_edge)):
            self.num_pts_per_edge.append(graph_data_to_be_copied.num_pts_per_edge[i])
        for i in range(len(graph_data_to_be_copied.edge_labels)):
            self.edge_labels.append(graph_data_to_be_copied.edge_labels[i])
        for i in range(len(graph_data_to_be_copied.edge_pt_coords)):
            self.edge_pt_coords.append(graph_data_to_be_copied.edge_pt_coords[i])
        for i in range(len(graph_data_to_be_copied.edge_pt_radii)):
            self.edge_pt_radii.append(graph_data_to_be_copied.edge_pt_radii[i])
        for i in range(len(graph_data_to_be_copied.edge_list)):
            self.edge_list.append(graph_data_to_be_copied.edge_list[i])
        self.transformation_matrix = graph_data_to_be_copied.transformation_matrix
        self.length = graph_data_to_be_copied.length
        
    def add_edge(self,start_vert,end_vert,label_id=0):
        self.vertices.append(start_vert)
        self.vertices.append(end_vert)
        self.vertex_labels.append(label_id)
        self.vertex_labels.append(label_id)
        
        self.edge_connectivity.append([len(self.vertices)-2,len(self.vertices)-1])
        
        self.num_pts_per_edge.append(2)
        self.edge_labels.append(label_id)
        
        self.edge_pt_coords.append(start_vert)
        self.edge_pt_coords.append(end_vert)
        
        self.edge_list.append([start_vert,end_vert])
        
        self.edge_pt_radii.append(0)
        self.edge_pt_radii.append(0)

    def add_contour(self,contour,label_id=0,sort_pts=True):
        # make sure the contour is closed
        # sort the contour pts
        if sort_pts:
            sorted_pts = []
            sorted_pts.append(contour[0])

            def is_in_list(ele,given_list):
                for thing in given_list:
                    if thing[0]==ele[0] and thing[1]==ele[1] and thing[2]==ele[2]:
                        return True
                return False
            for i in range(len(contour)):
                print(contour[i])
                min_pt = []
                min_dist = 999999
                for j in range(len(contour)):

                    dist = distance.euclidean(sorted_pts[len(sorted_pts)-1],contour[j])
                    #print(dist)
                    if dist<min_dist and dist>0 and not is_in_list(contour[j],sorted_pts):
                        min_dist=dist
                        min_pt = contour[j]
                if len(min_pt) > 0:
                    sorted_pts.append(min_pt)

            sorted_pts.append(contour[0])
            contour = sorted_pts

        self.vertices.append(contour[0])
        self.vertices.append(contour[len(contour)-1])
        self.vertex_labels.append(label_id)
        self.vertex_labels.append(label_id)

        self.edge_connectivity.append([len(self.vertices)-2,len(self.vertices)-1])

        self.num_pts_per_edge.append(len(contour))
        self.edge_labels.append(label_id)

        for pt in contour:
            self.edge_pt_coords.append(pt)
            self.edge_pt_radii.append(0)

        self.edge_list.append(contour)

    def apply_transformation(self,tr_mat,):
        if len(self.vertices)>0:
            tmp = np.matmul((self.add_ones_column(self.vertices)),tr_mat)
            self.vertices = (tmp)
        if len(self.edge_pt_coords)>0:
            tmp = np.matmul((self.add_ones_column(self.edge_pt_coords)),tr_mat)
            self.edge_pt_coords = (tmp)
        if len(self.edge_list) > 0:
            for i in range(len(self.edge_list)):
                tmp = np.matmul((self.add_ones_column(self.edge_list[i])),tr_mat)
                self.edge_list[i] = (tmp)

    def add_ones_column(self,data):
            a = np.array(data).transpose()
            b = [a[0],a[1],a[2],np.ones(len(data))]
            c = np.array(b).transpose()
            return c


class Dendrite(SpatialGraphData):
    def __init__(self,):
        self.apical_dendrite = SpatialGraphData()
        self.basal_dendrite = SpatialGraphData()

    def apply_transformation(self,tr_mat,):
        self.apical_dendrite.apply_transformation(tr_mat)
        self.basal_dendrite.apply_transformation(tr_mat)


class Neuron(Dendrite,SpatialGraphData):
    def __init__(self,):
        self.soma = SpatialGraphData()
        self.axon = SpatialGraphData()
        self.dendrite = Dendrite()
        self.all_neurites_subgraphdata = SpatialGraphData()

    def apply_transformation(self,tr_mat,):
        self.soma.apply_transformation(tr_mat)
        self.axon.apply_transformation(tr_mat)
        self.dendrite.apply_transformation(tr_mat)
        self.all_neurites_subgraphdata.apply_transformation(tr_mat)


class Barrel(SpatialGraphData):

    def __init__(self,data, validate=False, barrel_projections_present = False, reverse_direction = False, row_name='', column_name='',):
        # validate barrels if needed 
        # validation involves:
        # checking barrel height : shoudld be betweeen 350+-100
        self.barrel_name = row_name+column_name
        self.contours = data
        self.projections_present = barrel_projections_present
        self.sorted_vertices = []
        self.sorted_vertex_labels = []
        self.sorted_contours_list = []
        self.sorted_edge_num_pts = []
        self.sorted_edge_connectivity = []
        self.sorted_edge_labels = []
        self.sorted_edge_pt_coords = []
        self.sorted_edge_pt_radii =[]
        
        self.pia_projection_barrel = []
        self.pia_projection_barrel_centroid = []
        self.wm_projection_barrel = []
        self.wm_projection_barrel_centroid = []
        self.top_barrel = []
        self.top_barrel_centroid = []
        self.bottom_barrel = []
        self.bottom_barrel_centroid = []
        
        self.barrel_centroid = []
        self.barrel_column_centroid = []
        self.barrel_height = 0
        self.column_height = 0

        self.reverse_direction = reverse_direction

        if len(data.edge_pt_coords)>0:
            self.get_sorted_contours(validate)
        
    def get_sorted_contours(self,validate):
        inds = np.argsort((np.array(self.contours.vertices)[:,2]))
        for i in range(len(inds)):
            
            #self.sorted_vertex_labels.append(self.contours.vertex_labels[inds[i]])
            if i%2==0:
                self.sorted_contours_list.append(self.contours.edge_list[int(inds[i]/2)])
                self.sorted_vertices.append(self.contours.vertices[inds[i]])
                #self.sorted_edge_num_pts.append(self.contours.)

                ##print(self.sorted_contours_list)
        # sorted list has contours from z -ve to z +ve.. for now assume wm to pia
        if self.projections_present and len(self.sorted_contours_list) >=3:
            if not self.reverse_direction:
                self.wm_projection_barrel = self.sorted_contours_list[0]
                self.wm_projection_barrel_centroid = np.array(self.wm_projection_barrel).mean(axis=0)

                self.pia_projection_barrel = self.sorted_contours_list[len(self.sorted_contours_list)-1]
                self.pia_projection_barrel_centroid = np.array(self.pia_projection_barrel).mean(axis=0)
            else:
                self.wm_projection_barrel = self.sorted_contours_list[len(self.sorted_contours_list)-1]
                self.wm_projection_barrel_centroid = np.array(self.wm_projection_barrel).mean(axis=0)

                self.pia_projection_barrel = self.sorted_contours_list[0]
                self.pia_projection_barrel_centroid = np.array(self.pia_projection_barrel).mean(axis=0)

            self.column_height = self.pia_projection_barrel_centroid[2] - self.wm_projection_barrel_centroid[2]

            self.bottom_barrel = self.sorted_contours_list[1]
            self.bottom_barrel_centroid = np.array(self.bottom_barrel).mean(axis=0)

            self.top_barrel = self.sorted_contours_list[len(self.sorted_contours_list)-2]
            self.top_barrel_centroid = np.array(self.top_barrel).mean(axis=0)

            self.barrel_height = self.top_barrel_centroid[2] - self.bottom_barrel_centroid[2]

            self.barrel_centroid = (self.top_barrel_centroid + self.bottom_barrel_centroid) / 2.0
            
            self.barrel_column_centroid = (self.pia_projection_barrel_centroid + self.top_barrel_centroid +\
                                           self.bottom_barrel_centroid + self.wm_projection_barrel_centroid) / 4.0
            
            ##print(self.barrel_height,self.column_height)

            if validate:

                if (self.barrel_height < 450 and self.barrel_height > 250) \
                        and (self.column_height < 2500 and self.column_height > 1300):
                    # barrel valid 
                    # write only top, bottom barrels leaving the mid ones in graph data
                    new_graphdata = SpatialGraphData()

                    new_graphdata.vertices.append(self.sorted_vertices[0])
                    new_graphdata.vertices.append(self.sorted_vertices[0])
                    new_graphdata.vertices.append(self.sorted_vertices[1])
                    new_graphdata.vertices.append(self.sorted_vertices[1])
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-2])
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-2])
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-1])
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-1])

                    for i in range(8):
                        new_graphdata.vertex_labels.append(self.contours.label_id)

                    new_graphdata.edge_list.append(self.sorted_contours_list[0])
                    new_graphdata.edge_connectivity.append([0,1])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[0]))
                    for pt in (self.sorted_contours_list[0]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    new_graphdata.edge_list.append(self.sorted_contours_list[1])
                    new_graphdata.edge_connectivity.append([2,3])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[1]))
                    for pt in (self.sorted_contours_list[1]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    new_graphdata.edge_list.append(self.sorted_contours_list[len(self.sorted_contours_list)-2])
                    new_graphdata.edge_connectivity.append([4,5])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[len(self.sorted_contours_list)-2]))
                    for pt in (self.sorted_contours_list[len(self.sorted_contours_list)-2]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    new_graphdata.edge_list.append(self.sorted_contours_list[len(self.sorted_contours_list)-1])
                    new_graphdata.edge_connectivity.append([6,7])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[len(self.sorted_contours_list)-1]))
                    for pt in (self.sorted_contours_list[len(self.sorted_contours_list)-1]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    self.contours = new_graphdata
                    
                else:
                    #print('invalid barrel {}'.format(self.contours.label_name))
                    self.invalidate_barrel()
        
        else:
            if not self.reverse_direction:
                self.bottom_barrel = self.sorted_contours_list[0]
                self.bottom_barrel_centroid = np.array(self.bottom_barrel).mean(axis=0)

                self.top_barrel = self.sorted_contours_list[len(self.sorted_contours_list)-1]
                self.top_barrel_centroid = np.array(self.top_barrel).mean(axis=0)
            else:
                self.bottom_barrel = self.sorted_contours_list[0]
                self.bottom_barrel_centroid = np.array(self.bottom_barrel).mean(axis=0)

                self.top_barrel = self.sorted_contours_list[len(self.sorted_contours_list)-1]
                self.top_barrel_centroid = np.array(self.top_barrel).mean(axis=0)


            self.barrel_height = self.top_barrel_centroid[2] - self.bottom_barrel_centroid[2]

            #print('{} Barrel height {}'.format(self.contours.label_name,self.barrel_height))
            if validate:

                if (self.barrel_height <= 450 and self.barrel_height >= 150):
                    # barrel valid 
                    # write only top, bottom barrels leaving the mid ones in graph data
                    #print('valid barrel')
                    #print(self.contours.label_name, self.barrel_height,self.column_height)
                    ##print(self.top_barrel_centroid[2],self.bottom_barrel_centroid[2])
                    
                    new_graphdata = SpatialGraphData()

                    new_graphdata.vertices.append(self.sorted_vertices[0])
                    new_graphdata.vertices.append(self.sorted_vertices[0])
                    
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-1])
                    new_graphdata.vertices.append(self.sorted_vertices[len(self.sorted_contours_list)-1])

                    for i in range(4):
                        new_graphdata.vertex_labels.append(self.contours.label_id)

                    new_graphdata.edge_list.append(self.sorted_contours_list[0])
                    new_graphdata.edge_connectivity.append([0,1])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[0]))
                    for pt in (self.sorted_contours_list[0]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    new_graphdata.edge_list.append(self.sorted_contours_list[len(self.sorted_contours_list)-1])
                    new_graphdata.edge_connectivity.append([2,3])
                    new_graphdata.edge_labels.append(self.contours.label_id)
                    new_graphdata.num_pts_per_edge.append(len(self.sorted_contours_list[len(self.sorted_contours_list)-1]))
                    for pt in (self.sorted_contours_list[len(self.sorted_contours_list)-1]):
                        new_graphdata.edge_pt_coords.append(pt)
                        new_graphdata.edge_pt_radii.append(0)

                    self.contours = new_graphdata
                    
                else:
                    ##print('invalid barrel')
                    ##print(self.contours.label_name, self.barrel_height,self.column_height)
                    ##print(self.top_barrel_centroid[2],self.bottom_barrel_centroid[2])
                    self.invalidate_barrel()

    def invalidate_barrel(self,):
        #print('invalidating barrel')
        self.contours = SpatialGraphData()
        self.sorted_contours_list = []
        self.pia_projection_barrel = []
        self.pia_projection_barrel_centroid = []
        self.wm_projection_barrel = []
        self.wm_projection_barrel_centroid = []
        self.top_barrel = []
        self.top_barrel_centroid = []
        self.bottom_barrel = []
        self.bottom_barrel_centroid = []

        self.barrel_centroid = []
        self.barrel_column_centroid = []
        self.barrel_height = 0
        self.column_height = 0

    def apply_transformation(self,tr_mat):
        self.contours.apply_transformation(tr_mat)
        for i in range(len(self.sorted_contours_list)):
            self.sorted_contours_list[i] = self.apply_transformation_to_graph_data(tr_mat,self.sorted_contours_list[i])
        self.sorted_vertices = self.apply_transformation_to_graph_data(tr_mat,self.sorted_vertices)
        self.sorted_edge_pt_coords = self.apply_transformation_to_graph_data(tr_mat,self.sorted_edge_pt_coords)

        self.pia_projection_barrel = self.apply_transformation_to_graph_data(tr_mat,self.pia_projection_barrel)
        self.pia_projection_barrel_centroid = self.apply_transformation_to_pt(tr_mat,self.pia_projection_barrel_centroid)
        self.wm_projection_barrel = self.apply_transformation_to_graph_data(tr_mat,self.wm_projection_barrel)
        self.wm_projection_barrel_centroid = self.apply_transformation_to_pt(tr_mat,self.wm_projection_barrel_centroid)
        self.top_barrel = self.apply_transformation_to_graph_data(tr_mat,self.top_barrel)
        self.top_barrel_centroid = self.apply_transformation_to_pt(tr_mat,self.top_barrel_centroid)
        self.bottom_barrel = self.apply_transformation_to_graph_data(tr_mat,self.bottom_barrel)
        self.bottom_barrel_centroid = self.apply_transformation_to_pt(tr_mat,self.bottom_barrel_centroid)

        self.barrel_centroid = self.apply_transformation_to_pt(tr_mat,self.barrel_centroid)
        self.barrel_column_centroid = self.apply_transformation_to_pt(tr_mat,self.barrel_column_centroid)

    def apply_transformation_to_graph_data(self,tr_mat,graph_data):
        if len(graph_data) > 0:
            tmp = np.matmul((self.add_ones_column(graph_data)),tr_mat)
            graph_data = (tmp)
        return graph_data

    def add_ones_column(self,data):
        #if len(data) > 1:
        #reshaped = np.array(data).reshape(len(data),3)
        a = np.array(data).transpose()
        b = [a[0],a[1],a[2],np.ones(len(data))]
        c = np.array(b).transpose()
        return c

    def apply_transformation_to_pt(self,tr_mat,pt):
        if len(pt) > 0:
            tmp = np.matmul((self.add_one(pt)),tr_mat)
            pt = (tmp)
        return pt

    def add_one(self,pt):
        return [pt[0],pt[1],pt[2],1]


class BarrelRow(Barrel):
     def __init__(self,row_name=''):
        self.row_name = row_name
        self.one = Barrel(SpatialGraphData(),row_name=row_name,column_name='1')
        self.two = Barrel(SpatialGraphData(),row_name=row_name,column_name='2')
        self.three = Barrel(SpatialGraphData(),row_name=row_name,column_name='3')
        self.four = Barrel(SpatialGraphData(),row_name=row_name,column_name='4')
        self.five = Barrel(SpatialGraphData(),row_name=row_name,column_name='5')
        self.six = Barrel(SpatialGraphData(),row_name=row_name,column_name='6')
        self.seven = Barrel(SpatialGraphData(),row_name=row_name,column_name='7')
        self.all_barrels_graphdata = SpatialGraphData()
        self.single_barrels_list = []

     def apply_transformation(self,tr_mat):
        self.one.apply_transformation(tr_mat)
        self.two.apply_transformation(tr_mat)
        self.three.apply_transformation(tr_mat)
        self.four.apply_transformation(tr_mat)
        self.five.apply_transformation(tr_mat)
        self.six.apply_transformation(tr_mat)
        self.seven.apply_transformation(tr_mat)
        self.all_barrels_graphdata.apply_transformation(tr_mat)


class Barrels(BarrelRow,):
     def __init__(self,):
        self.a_row = BarrelRow(row_name='A')
        self.b_row = BarrelRow(row_name='B')
        self.c_row = BarrelRow(row_name='C')
        self.d_row = BarrelRow(row_name='D')
        self.e_row = BarrelRow(row_name='E')
        self.greek_row = BarrelRow(row_name='Greek')
        self.all_rows_graphdata = SpatialGraphData()
        self.rows_list = []

     def apply_transformation(self,tr_mat):
        self.a_row.apply_transformation(tr_mat)
        self.b_row.apply_transformation(tr_mat)
        self.c_row.apply_transformation(tr_mat)
        self.d_row.apply_transformation(tr_mat)
        self.e_row.apply_transformation(tr_mat)
        self.greek_row.apply_transformation(tr_mat)
        self.all_rows_graphdata.apply_transformation(tr_mat)


class MatchBarrels(Barrels):
    def __init__(self,ref_barrels,barrels,use_barrels_alone=False,use_projections_alone=False):
        self.ref_barrels = ref_barrels
        self.barrels = barrels
        self.ref_matching_barrels = []
        self.matching_barrels = []
        self.ref_barrel_centroids = []
        self.barrel_centroids = []
        self.ref_barrel_pia_centroids = []
        self.barrel_pia_centroids = []
        self.ref_barrel_wm_centroids = []
        self.barrel_wm_centroids = []
        self.ref_barrel_top_centroids = []
        self.barrel_top_centroids = []
        self.ref_barrel_bottom_centroids = []
        self.barrel_bottom_centroids = []
        for ref_barrel_row in ref_barrels.rows_list:
            ##print(barrel_row.one.contours.edge_list)
            for ref_barrel in ref_barrel_row.single_barrels_list:

                for barrel_row in barrels.rows_list:
                    ##print(barrel_row.one.contours.edge_list)
                    for barrel in barrel_row.single_barrels_list:
                        if len(ref_barrel.contours.vertices)>0 and len(barrel.contours.vertices)>0 and\
                                (ref_barrel.contours.label_id == barrel.contours.label_id) and \
                                barrel.contours.label_id!=0:

                            self.ref_matching_barrels.append(ref_barrel)
                            self.matching_barrels.append(barrel)
                            ##print(ref_barrel.top_barrel_centroid,barrel.top_barrel_centroid)

        if use_barrels_alone:
            for i in range(len(self.ref_matching_barrels)):
                if len(self.ref_matching_barrels[i].top_barrel_centroid)>0 and len(self.matching_barrels[i].top_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].top_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].top_barrel_centroid)
                if len(self.ref_matching_barrels[i].bottom_barrel_centroid)>0 and len(self.matching_barrels[i].bottom_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].bottom_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].bottom_barrel_centroid)

        elif use_projections_alone:
            for i in range(len(self.ref_matching_barrels)):
                if len(self.ref_matching_barrels[i].pia_projection_barrel_centroid)>0 and len(self.matching_barrels[i].pia_projection_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].pia_projection_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].pia_projection_barrel_centroid)
                    self.ref_barrel_pia_centroids.append(self.ref_matching_barrels[i].pia_projection_barrel_centroid)
                    self.barrel_pia_centroids.append(self.matching_barrels[i].pia_projection_barrel_centroid)
                if len(self.ref_matching_barrels[i].wm_projection_barrel_centroid)>0 and len(self.matching_barrels[i].wm_projection_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].wm_projection_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].wm_projection_barrel_centroid)
                    self.ref_barrel_wm_centroids.append(self.ref_matching_barrels[i].wm_projection_barrel_centroid)
                    self.barrel_wm_centroids.append(self.matching_barrels[i].wm_projection_barrel_centroid)

        else:# use both barrels and column(projections)
            for i in range(len(self.ref_matching_barrels)):
                if len(self.ref_matching_barrels[i].top_barrel_centroid)>0 and len(self.matching_barrels[i].top_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].top_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].top_barrel_centroid)
                    self.ref_barrel_top_centroids.append(self.ref_matching_barrels[i].top_barrel_centroid)
                    self.barrel_top_centroids.append(self.matching_barrels[i].top_barrel_centroid)
                if len(self.ref_matching_barrels[i].bottom_barrel_centroid)>0 and len(self.matching_barrels[i].bottom_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].bottom_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].bottom_barrel_centroid)
                    self.ref_barrel_bottom_centroids.append(self.ref_matching_barrels[i].bottom_barrel_centroid)
                    self.barrel_bottom_centroids.append(self.matching_barrels[i].bottom_barrel_centroid)
                if len(self.ref_matching_barrels[i].pia_projection_barrel_centroid)>0 and len(self.matching_barrels[i].pia_projection_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].pia_projection_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].pia_projection_barrel_centroid)
                    self.ref_barrel_pia_centroids.append(self.ref_matching_barrels[i].pia_projection_barrel_centroid)
                    self.barrel_pia_centroids.append(self.matching_barrels[i].pia_projection_barrel_centroid)
                if len(self.ref_matching_barrels[i].wm_projection_barrel_centroid)>0 and len(self.matching_barrels[i].wm_projection_barrel_centroid)>0:
                    self.ref_barrel_centroids.append(self.ref_matching_barrels[i].wm_projection_barrel_centroid)
                    self.barrel_centroids.append(self.matching_barrels[i].wm_projection_barrel_centroid)
                    self.ref_barrel_wm_centroids.append(self.ref_matching_barrels[i].wm_projection_barrel_centroid)
                    self.barrel_wm_centroids.append(self.matching_barrels[i].wm_projection_barrel_centroid)

    def get_matching_barrels(self):
        return self.ref_matching_barrels,self.matching_barrels

    def get_matching_barrel_centroids(self):
        return self.ref_barrel_centroids,self.barrel_centroids

    def get_matching_barrel_pia_centroids(self):
        return self.ref_barrel_pia_centroids,self.barrel_pia_centroids

    def get_matching_barrel_wm_centroids(self):
        return self.ref_barrel_wm_centroids,self.barrel_wm_centroids

    def get_matching_barrel_top_centroids(self):
        return self.ref_barrel_top_centroids,self.barrel_top_centroids

    def get_matching_barrel_bottom_centroids(self):
        return self.ref_barrel_bottom_centroids,self.barrel_bottom_centroids

    def get_matching_standard_barrel_for_given(self,barrel):
        for ref_barrel_row in self.ref_barrels.rows_list:
            ##print(barrel_row.one.contours.edge_list)
            for ref_barrel in ref_barrel_row.single_barrels_list:

                if ref_barrel.contours.label_id == barrel.contours.label_id:
                    return ref_barrel

    def get_matching_axes(self,):
        ref_barrel_axes = []
        barrel_axes = []
        # return barrel coumns
        for i in range(len(self.matching_barrels)):
            ref_barrel_axes.append([self.ref_barrel_wm_centroids[i],self.ref_barrel_pia_centroids[i]])
            barrel_axes.append([self.barrel_wm_centroids[i],self.barrel_pia_centroids[i]])

        return ref_barrel_axes,barrel_axes

    def get_errors(self,exp_name='',alignment_type=''):
        barrel_names = []
        pia_center_errors = []
        barrel_center_errors = []
        wm_center_errors = []
        column_center_errors = []
        pia_centers_means = []
        pia_centers_stds = []
        barrel_centers_means = []
        barrel_centers_stds = []
        wm_centers_means = []
        wm_centers_stds = []
        column_centers_means = []
        column_centers_stds = []

        angles_means = []
        angles_stds = []

        for i in range(len(self.matching_barrels)):
            barrel_names.append(self.matching_barrels[i].barrel_name)
            #print(self.matching_barrels[i].barrel_name)
            pia_center_errors.append(spatial.distance.euclidean(self.ref_barrel_pia_centroids[i],self.barrel_pia_centroids[i]))
            wm_center_errors.append(spatial.distance.euclidean(self.ref_barrel_wm_centroids[i],self.barrel_wm_centroids[i]))

            ref_barrel_centroid = (np.array(self.ref_matching_barrels[i].top_barrel_centroid) + np.array(self.ref_matching_barrels[i].bottom_barrel_centroid)) / 2.0
            barrel_centroid = (np.array(self.matching_barrels[i].top_barrel_centroid) + np.array(self.matching_barrels[i].bottom_barrel_centroid)) / 2.0
            barrel_center_errors.append(spatial.distance.euclidean(ref_barrel_centroid,barrel_centroid))

            ref_column_centroid = (np.array(self.ref_barrel_pia_centroids[i]) + np.array(self.ref_barrel_wm_centroids[i])) / 2.0
            column_centroid = (np.array(self.barrel_pia_centroids[i]) + np.array(self.barrel_wm_centroids[i])) / 2.0
            column_center_errors.append(spatial.distance.euclidean(ref_column_centroid,column_centroid))

            # Find angluar precision of barrel column axis
            #angles,andgles_mean,angles_std,mean_uv = self.get_mean_angle_between_vectors\
            #    ([[self.ref_barrel_wm_centroids[i],self.ref_barrel_pia_centroids[i]],[self.barrel_wm_centroids[i],self.barrel_pia_centroids[i]]])
            #angles_means.append(angles_mean)
            #angles_stds.append(angles_std)

            #print(np.array([np.array(self.ref_barrel_pia_centroids[i]),np.array(self.barrel_pia_centroids[i])]).mean(axis=0))
            pia_centers_means.append(np.array([np.array(self.ref_barrel_pia_centroids[i]),np.array(self.barrel_pia_centroids[i])]).mean(axis=0))
            pia_centers_stds.append(np.array([np.array(self.ref_barrel_pia_centroids[i]),np.array(self.barrel_pia_centroids[i])]).std(axis=0))
            barrel_centers_means.append(np.array([np.array(ref_barrel_centroid),np.array(barrel_centroid)]).mean(axis=0))
            barrel_centers_stds.append(np.array([np.array(ref_barrel_centroid),np.array(barrel_centroid)]).std(axis=0))
            wm_centers_means.append(np.array([np.array(self.ref_barrel_wm_centroids[i]),np.array(self.barrel_wm_centroids[i])]).mean(axis=0))
            wm_centers_stds.append(np.array([np.array(self.ref_barrel_wm_centroids[i]),np.array(self.barrel_wm_centroids[i])]).std(axis=0))
            column_centers_means.append(np.array([np.array(ref_column_centroid),np.array(column_centroid)]).mean(axis=0))
            column_centers_stds.append(np.array([np.array(ref_column_centroid),np.array(column_centroid)]).std(axis=0))


        exp_name_list = []
        alignment_type_list = []
        for i in range(len(barrel_names)):
            exp_name_list.append(exp_name)
        for i in range(len(barrel_names)):
            alignment_type_list.append(alignment_type)
        df = pd.DataFrame()
        df['Exp_Name'] = exp_name_list
        df['Alignment_Type'] = alignment_type_list
        df['Barrel_Names'] = barrel_names
        df['Pia_Errors'] = pia_center_errors
        df['Barrel_Center_Errors'] = barrel_center_errors
        df['WM_Errors'] = wm_center_errors
        df['Column_Center_Errors'] = column_center_errors

        df['Pia_Centers_Means'] = pia_centers_means
        df['Pia_Centers_Stds'] = pia_centers_stds
        df['WM_Centers_Means'] = wm_centers_means
        df['WM_Centers_Stds'] = wm_centers_stds
        df['Barrel_Centers_Means'] = barrel_centers_means
        df['Barrel_Centers_Stds'] = barrel_centers_stds
        df['Column_Centers_Means'] = column_centers_means
        df['Column_Centers_Stds'] = column_centers_stds

        return df,pia_centers_stds,wm_centers_stds,barrel_centers_stds,column_centers_stds

    def get_matching_centroids_for_given_barrel_name(self,barrel_name='D2'):
        for i in range(len(self.matching_barrels)):
            if self.matching_barrels[i].barrel_name == barrel_name:
                return self.barrel_pia_centroids[i],self.barrel_wm_centroids[i],\
                        (np.array(self.barrel_pia_centroids[i]) + np.array(self.barrel_wm_centroids[i])) / 2.0



class AmiraSpatialGraph(Neuron,Barrels,SpatialGraphData,GraphDecorators) :
    ''' 
    This class has the following functionalities
        1) It can read / write Amira Spatial Graphs
        2) It can split / join spatial graphs based on 
            - Labelled objects like pia,wm,neuron etc
            - z coords: used to split 3D graph into section graphs and vice-versa
        3) Apply transformation to the spatial graph
    '''
    
    def __init__(self, filename = '',labels=[], read_header_only =False, axis_directions=[1,1,1],\
                 scaling=[1,1,1],resample_resolution = None,compute_edge_length=False,create_section_graphs=False,\
                 validate_barrels=False,read_transformation_matrix=False,barrel_projections_present=False,generic_graph=False,\
                 labelled_graphs_to_be_read=[],file_format = 'am',reverse_direction=False,reverse_barrel_direction=False):
        ##print('init called')
        
        if not filename or (read_header_only and filename is not None ):
            if read_header_only and filename is not None:
                self.input_filename = filename
                self.read_transformation_matrix = read_transformation_matrix 
                self.label_tree_text = self.read_label_tree_text()#'/nas1/Data_Mythreya/MotorCortexProject/V0/0_Inputs/SpatialGraphs/Default_Spatial_Graph_Label_Header.txt'
            else:
                self.input_filename = ''
                #self.default_header_file = '/nas1/Data_Mythreya/MotorCortexProject/V0/0_Inputs/SpatialGraphs/Default_Spatial_Graph_Label_Header.txt'
                self.default_header_file = os.path.join(pathlib.Path(__file__).parent.absolute() , pathlib.Path('./Default_Spatial_Graph_Label_Header.txt'))
                #print(self.default_header_file)
                self.read_transformation_matrix = read_transformation_matrix 
                self.label_tree_text = self.read_text(self.default_header_file)
            
            self.label_tree = []
            self.label_tree_depth_level = 0
            self.label_tree.append(Node(['Root',0,[0,0,0]]))
            self.barrel_projections_present = barrel_projections_present
            self.axis_directions = axis_directions
            self.scaling = scaling
            self.validate_barrels = validate_barrels
            self.decorators = GraphDecorators()
            self.graph_data = SpatialGraphData()
            self.labelMap = None
            self.pia = SpatialGraphData()
            self.labelled_subgraphs_list = []
            self.wm = SpatialGraphData()
            self.bvs = SpatialGraphData()
            self.barrels = Barrels()
            self.neuron = Neuron()
            
            self.number_z_sections = 0
            self.section_subgraphs = []#SpatialGraphData()[]#self.split_spatial_graph_into_z_sections()
            
            self.transformation_matrix = []
            self.transformation_applied = False
            self.resample_resolution = None
            self.compute_edge_length = False
            self.create_section_graphs = False
            self.labelled_graphs_to_be_read = []
            self.generic_graph = generic_graph
            self.reverse_direction = reverse_direction
            self.reverse_barrel_direction = reverse_barrel_direction
            
        else:
            # Read from the given spatial graph file
            self.input_filename = filename
            self.file_format = file_format
            #default_header_file = '/nas1/Data_Mythreya/MotorCortexProject/V0/0_Inputs/SpatialGraphs/Default_Spatial_Graph_Label_Header.txt'
            self.default_header_file = os.path.join(pathlib.Path(__file__).parent.absolute() , pathlib.Path('./Default_Spatial_Graph_Label_Header.txt'))
            self.axis_directions = axis_directions
            self.scaling = scaling
            self.validate_barrels = validate_barrels
            self.compute_edge_length = compute_edge_length
            self.create_section_graphs = create_section_graphs
            self.transformation_matrix = []
            self.transformation_applied = False
            # first read the label tree from the file header
            self.read_transformation_matrix = read_transformation_matrix 
            self.label_tree_text = self.read_label_tree_text()
            self.label_tree = []
            self.label_tree_depth_level = 0
            self.label_tree.append(Node(['Root',0,[0,0,0]]))
            self.labelMap = None
            # parse the label tree text to find label id and color for the leaf nodes
            self.parse_label_tree(0)
            self.barrel_projections_present = barrel_projections_present
            # read the '@' decorators
            self.decorators = self.read_decorators()
            self.labelled_subgraphs_list = []
            self.labelled_graphs_to_be_read = labelled_graphs_to_be_read
            self.generic_graph = generic_graph
            self.reverse_direction = reverse_direction
            self.reverse_barrel_direction = reverse_barrel_direction

            if self.file_format == 'ASC':
                # read the amira default header
                #self.default_header_file = os.curdir+'/Default_Spatial_Graph_Label_Header.txt'
                self.default_header_file = os.path.join(pathlib.Path(__file__).parent.absolute() , pathlib.Path('./Default_Spatial_Graph_Label_Header.txt'))
                self.label_tree_text = self.read_text(self.default_header_file)
                
            else:
                # now read the raw spatial graph
                self.graph_data = self.read_spatial_graph()
            ##print(self.graph_data.transformation_matrix)
            
            #self.remove_isolated_nodes()
            self.resample_resolution = resample_resolution
            
            # get only those labelled sugraphs that the user needs
            # other wise takes too long
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['Pia','alpha','pia']):
                self.pia = self.get_labelled_subgraph(['Pia','alpha','pia'],label_id=9,label_color=[0,1,0.5])
                self.labelled_subgraphs_list.append(self.pia)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['WM','wm','WhiteMatter','Label0','cc']):
                self.wm = self.get_labelled_subgraph(['WM','wm','WhiteMatter','Label0','cc'],label_id=10,label_color=[1,0.5,0] )
                self.labelled_subgraphs_list.append(self.wm)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['Vessels','vessels','Vessel','vessel',]):
                self.bvs = self.get_labelled_subgraph(['Vessels','vessels','Vessel','vessel',])
                self.labelled_subgraphs_list.append(self.bvs)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['Neuron','neuron']):
                self.neuron = self.get_neuron_subgraph(['Neuron','neuron'])
                self.labelled_subgraphs_list.append(self.neuron.all_neurites_subgraphdata)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['Barrels','barrels','Barrel','barrel']):
                self.barrels = self.get_barrels_subgraph(['Barrels','barrels','Barrel','barrel'])
                self.labelled_subgraphs_list.append(self.barrels.all_rows_graphdata)
            
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_landmarks','Rabies_Landmarks',]):
                self.rabies_landmarks = self.get_labelled_subgraph(['rabies_landmarks','Rabies_Landmarks',])
                self.labelled_subgraphs_list.append(self.rabies_landmarks)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vS1_landmarks','Rabies_vS1_Landmarks',]):
                self.rabies_S1_landmarks = self.get_labelled_subgraph(['rabies_vS1_landmarks','Rabies_vS1_Landmarks',])
                self.labelled_subgraphs_list.append(self.rabies_S1_landmarks)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vM1_landmarks','Rabies_vM1_Landmarks',]):
                self.rabies_M1_landmarks = self.get_labelled_subgraph(['rabies_vM1_landmarks','Rabies_vM1_Landmarks',])
                self.labelled_subgraphs_list.append(self.rabies_M1_landmarks)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vM1_l5_landmarks','Rabies_vM1_L5_Landmarks',]):
                self.rabies_M1_l5_landmarks = self.get_labelled_subgraph(['rabies_vM1_l5_landmarks','Rabies_vM1_L5_Landmarks',])
                self.labelled_subgraphs_list.append(self.rabies_M1_l5_landmarks)
                
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_landmarks_centroid','Rabies_Landmarks_Centroid',]):
                self.rabies_landmarks_centroid = self.get_labelled_subgraph(['rabies_landmarks_centroid','Rabies_Landmarks_Centroid',])
                self.labelled_subgraphs_list.append(self.rabies_landmarks_centroid)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vS1_landmarks_centroid','Rabies_vS1_Landmarks_Centroid',]):
                self.rabies_S1_landmarks_centroid = self.get_labelled_subgraph(['rabies_vS1_landmarks_centroid','Rabies_vS1_Landmarks_Centroid',])
                self.labelled_subgraphs_list.append(self.rabies_S1_landmarks_centroid)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vM1_landmarks_centroid','Rabies_vM1_Landmarks_Centroid',]):
                self.rabies_M1_landmarks_centroid = self.get_labelled_subgraph(['rabies_vM1_landmarks_centroid','Rabies_vM1_Landmarks_Centroid',])
                self.labelled_subgraphs_list.append(self.rabies_M1_landmarks_centroid)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['rabies_vM1_L5_landmarks_centroid','Rabies_vM1_L5_Landmarks_Centroid',]):
                self.rabies_M1_L5_landmarks_centroid = self.get_labelled_subgraph(['rabies_vM1_L5_landmarks_centroid','Rabies_vM1_L5_Landmarks_Centroid',])
                self.labelled_subgraphs_list.append(self.rabies_M1_L5_landmarks_centroid)
                
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['pia_to_wm_axis_field','Pia_To_WM_Axis_Field',]):
                self.pia_wm_axis_field = self.get_labelled_subgraph(['pia_to_wm_axis_field','Pia_To_WM_Axis_Field',])
                self.labelled_subgraphs_list.append(self.pia_wm_axis_field)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['wm_to_pia_axis_field','WM_To_Pia_Axis_Field',]):
                self.wm_pia_axis_field = self.get_labelled_subgraph(['wm_to_pia_axis_field','WM_To_Pia_Axis_Field',])
                self.labelled_subgraphs_list.append(self.wm_pia_axis_field)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['s1_axis_field','S1_Axis_Field',]):
                self.S1_axis_field = self.get_labelled_subgraph(['s1_axis_field','S1_Axis_Field',])
                self.labelled_subgraphs_list.append(self.S1_axis_field)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['m1_axis_field','M1_Axis_Field',]):
                self.M1_axis_field = self.get_labelled_subgraph(['m1_axis_field','M1_Axis_Field',])
                self.labelled_subgraphs_list.append(self.M1_axis_field)
        
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['M1_PCA_0','m1_pca_0']):
                self.m1_pca_1 = self.get_labelled_subgraph(['M1_PCA_0','m1_pca_0',])
                self.labelled_subgraphs_list.append(self.m1_pca_1)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['M1_PCA_1','m1_pca_1']):
                self.m1_pca_2 = self.get_labelled_subgraph(['M1_PCA_1','m1_pca_1',])
                self.labelled_subgraphs_list.append(self.m1_pca_2)
            if len(labelled_graphs_to_be_read) == 0 or self.label_present_in_to_read_list(['M1_PCA_2','m1_pca_2']):
                self.m1_pca_3 = self.get_labelled_subgraph(['M1_PCA_2','m1_pca_2',])
                self.labelled_subgraphs_list.append(self.m1_pca_3)
            
            # now stitch back together the whole graph from labelled subgraphs
            if not generic_graph:
                self.update_graph_data()
            #self.graph_data = self.combine_subgraphs([self.pia,self.wm,\
            #                                           self.bvs,self.barrels.all_rows_graphdata,\
            #                                           self.neuron.all_neurites_subgraphdata])
            
            # split the graph into z sections
            if len(self.pia.edge_pt_coords) > 0 and self.create_section_graphs:
                z_sorted_coords = sorted( self.pia.edge_pt_coords, key=lambda coord: coord[2])
                unique_z_list = np.unique(np.transpose(np.array(z_sorted_coords))[2])
                if len(unique_z_list) > 1:
                    self.number_z_sections = len(unique_z_list)
                    #print(self.number_z_sections)
                    self.section_subgraphs = self.split_spatial_graph_into_z_sections()
            else:
                self.number_z_sections = 1
                self.section_subgraphs = [self.graph_data]#self.split_spatial_graph_into_z_sections()
      
    def label_present_in_to_read_list(self,alternate_labels):
        for label in alternate_labels:
            for label_ in self.labelled_graphs_to_be_read:
                if label == label_:
                    return True
        return False
            
    def apply_transformation(self,transformation_matrix=[],z_desired=None,trans_2d_only=False):
        
        if len(transformation_matrix) > 0:
            self.transformation_matrix = transformation_matrix
        
        # transformation matrix 4x4 as a row first list
        if  len(self.transformation_matrix)> 0:
            self.transformation_applied = True
            tr_mat = np.reshape(self.transformation_matrix,[4,4])
            if z_desired is not None:
                # to translate in z to the desired position for section graphs
                ##print(self.pia.edge_pt_coords[0][2])
                #tr_mat[2,3] = z_desired - self.pia.edge_pt_coords[0][2]
                tr_mat[2,3] = z_desired - self.graph_data.edge_pt_coords[0][2]
            if trans_2d_only:
                tr_mat[2,2] = 1
                tr_mat[0,2] = 0
                tr_mat[1,2] = 0
                tr_mat[2,0] = 0
                tr_mat[2,1] = 0
            
            if self.generic_graph:
                self.graph_data = self.apply_transformation_to_labelled_subgraph(self.graph_data,tr_mat)
            else:
                #self.pia = self.apply_transformation_to_labelled_subgraph(self.pia,tr_mat)
                #self.wm = self.apply_transformation_to_labelled_subgraph(self.wm,tr_mat)
                #self.bvs = self.apply_transformation_to_labelled_subgraph(self.bvs,tr_mat)
                #self.barrels.all_rows_graphdata = self.apply_transformation_to_labelled_subgraph(self.barrels.all_rows_graphdata,tr_mat)
                #self.neuron.all_neurites_subgraphdata = self.apply_transformation_to_labelled_subgraph(self.neuron.all_neurites_subgraphdata,tr_mat)
                #self.update_graph_data()
                self.pia.apply_transformation(tr_mat)
                self.wm.apply_transformation(tr_mat)
                self.bvs.apply_transformation(tr_mat)
                self.barrels.apply_transformation(tr_mat)
                self.neuron.apply_transformation(tr_mat)
                self.update_graph_data()

    def apply_transformation_to_labelled_subgraph(self,graph_data,tr_mat):
        #if len(graph_data.vertices)>0:
        #    np.array(graph_data.vertices).transpose()
        ##print(graph_data.vertices)
        if len(graph_data.vertices)>0:
            #print(graph_data.vertices)
            tmp = np.matmul((self.add_ones_column(graph_data.vertices)),tr_mat)
            #print(tmp)
            graph_data.vertices = (tmp)
        if len(graph_data.edge_pt_coords)>0:
            tmp = np.matmul((self.add_ones_column(graph_data.edge_pt_coords)),tr_mat)
            graph_data.edge_pt_coords = (tmp)
        #for i in range(len(graph_data.vertices)):
        #    pt = graph_data.vertices[i]
        #    graph_data.vertices[i] = np.matmul(tr_mat,[pt[0],pt[1],pt[2],1])
        #for i in range(len(graph_data.edge_pt_coords)):
        #    pt = graph_data.edge_pt_coords[i]
        #    graph_data.edge_pt_coords[i] = np.matmul(tr_mat,[pt[0],pt[1],pt[2],1])
            
        return graph_data
    
    def apply_rotation(self,rot_mat=[],translation=None,trans_2d_only=False):
        if self.generic_graph == True:
            # usually the case for axis fields
            self.graph_data = self.apply_rotation_to_labelled_subgraph(self.graph_data,rot_mat,translation)
        else:
            transformed_subgraphs_list = []
            for subgraph in self.labelled_subgraphs_list:
                transformed_subgraphs_list.append(self.apply_rotation_to_labelled_subgraph(subgraph,rot_mat,translation))
            #self.pia = self.apply_rotation_to_labelled_subgraph(self.pia,rot_mat,translation)
            #self.wm = self.apply_rotation_to_labelled_subgraph(self.wm,rot_mat,translation)
            #self.bvs = self.apply_rotation_to_labelled_subgraph(self.bvs,rot_mat,translation)
            #self.barrels.all_rows_graphdata = self.apply_rotation_to_labelled_subgraph\
            #                                    (self.barrels.all_rows_graphdata,rot_mat,translation)
            #self.neuron.all_neurites_subgraphdata = self.apply_rotation_to_labelled_subgraph\
            #                                            (self.neuron.all_neurites_subgraphdata,rot_mat,translation)
            self.labelled_subgraphs_list = transformed_subgraphs_list
            self.update_graph_data()    
        
    def apply_rotation_to_labelled_subgraph(self,graph_data,rot_mat,translation=None):
        if len(graph_data.vertices)>0:
            if translation is not None:
                graph_data.vertices = np.array(graph_data.vertices) - np.array(translation)
            tmp = np.matmul(rot_mat,np.transpose((graph_data.vertices)))
            graph_data.vertices = np.transpose(tmp)

        if len(graph_data.edge_pt_coords)>0:
            if translation is not None:
                graph_data.edge_pt_coords = np.array(graph_data.edge_pt_coords) - np.array(translation)
            tmp = np.matmul(rot_mat,np.transpose((graph_data.edge_pt_coords)))
            graph_data.edge_pt_coords = np.transpose(tmp)
        return graph_data
            
    def add_ones_column(self,data):
        a = np.array(data).transpose()
        b = [a[0],a[1],a[2],np.ones(len(data))]
        c = np.array(b).transpose()
        return c

    def read_label_tree_text(self):
        label_tree_text = []
        read_flag = False
        lines = []
        with open(self.input_filename,'r') as sp:
            lines = sp.readlines()
        for line in lines:
            if line.startswith('Parameters {'):
                read_flag = True

            if read_flag:
                ##print(line)
                if self.read_transformation_matrix == False and line.rfind("TransformationMatrix ") > -1: 
                    ##print('found tx line')
                    continue
                elif line.rfind("@") == 0 or line.isspace() or line.startswith('VERTEX {') : 
                        return label_tree_text
                else:
                      
                    label_tree_text.append(line)
        
        
        return label_tree_text
    
    def parse_label_tree(self,line_number):
        ''' Read in the whole label tree hierarchy structure of the spatial graph 
            This defines label ids. which will be used to fetch edge pts
            belonging to a certain labelled object such as pia, wm or neuron
        '''
        # look for a specific label (pia,wm etc)
        # if found put them in the respective object
        # read only the labels under "GraphLabel"
        valid_tree = False
        lines = []
        #with open(self.input_filename,'r') as sp:
        #    lines = sp.readlines()
        lines = self.label_tree_text
        for i in range(len(lines)):
            line = lines[i]
            ##print(line)
            if (line.rfind("{")>-1) and i > line_number :
                self.label_tree_depth_level = self.label_tree_depth_level + 1
                if  self.label_tree_depth_level >= 0:

                    # fetch node name 
                    # this ends in eiter '{' or a line feed
                    #line.split("\n")
                    label_name = line[0:line.rfind('{')].rstrip().lstrip()
                    #print(label_name)
                    Id = 0
                    color = [0,0,0]
                    # check the next lines to see if this is a leaf node
                    # fetch the labels if thats tha case
                    # next line should be 'color' next to that Id
                    if (lines[i+1].rfind('Color') > -1) and (lines[i+2].rfind('Id') > -1) :
                        ##print(lines[i+1],lines[i+2])
                        for sub_string in lines[i+2].split():
                            if sub_string.isdigit():
                                Id = int(sub_string)
                                break
                        #Id = int(lines[i+2][lines[i+2].rfind("Id")+3:])
                        color_str = lines[i+1][lines[i+1].rfind("Color")+6:len(lines[i+1])-2]
                        color = list(map(float,color_str.split(" ")))

                    # find parent node
                    # parent node is the guy who is one level above and while doing
                    # a back search in the linear list
                    node_id = len(self.label_tree)
                    while node_id > 0 :
                        node_id = node_id -1
                        node = self.label_tree[node_id]
                        if node.depth == self.label_tree_depth_level-1:
                            break

                    self.label_tree.append(Node([label_name,Id,color],parent=node))
                    self.parse_label_tree(i)

            if line.rfind("}")>-1 and i > line_number:
                ##print(line)
                ##print('exiting block')
                ##print(self.label_tree_depth_level)
                self.label_tree_depth_level = self.label_tree_depth_level - 1
    
    def resample_contour(self,cnt,res,downscale=1):
        resptlist=[]
        ##print(cnt)
        for i in range(len(cnt)-1):
       
            pt1 = np.reshape(np.array(cnt)[i,:],[1,3])
            pt2 = np.reshape(np.array(cnt)[i+1,:],[1,3])

            dist = (spatial.distance.cdist(pt1,pt2))
            num_pts_inserted = int(dist / res)
            ##print('num_pts_to be inserted: {}'.format(num_pts_inserted))
            #m = (pt2[:,1]-pt2[:,0]) / (pt1[:,1]-pt1[:,0])
            #if num_pts_inserted > 0:
                ##print('needed_{}'.format(num_pts_inserted))
            ptlist = []
            ptlist.append(cnt[i])
            ptlist.append(cnt[i+1])
            ##print(ptlist)
            
            ptlist = self.get_mid_line(ptlist,num_pts_inserted)
            ##print(ptlist)
            for i in range(len(ptlist)-1):
                resptlist.append(ptlist[i])
        # add the last point now
        resptlist.append(cnt[len(cnt)-1])
        ##print(resptlist)
        return resptlist

    def get_mid_line(self,ptlist,limit):
        newptlist = []
        ##print('in getmidline')
        ##print(ptlist)
        if len(ptlist) < (limit/2):
            ##print('inside check')
            for i in range(len(ptlist)-1):
                newptlist.append(ptlist[i])
                ##print(ptlist[i])
                newptlist.append([(ptlist[i][0]+ptlist[i+1][0])/2,                                  (ptlist[i][1]+ptlist[i+1][1])/2,                                  (ptlist[i][2]+ptlist[i+1][2])/2])
            newptlist.append(ptlist[i+1])
            ptlist = newptlist
            
            return self.get_mid_line(ptlist,limit)
            #return newptlist
        else:
            
            return ptlist
    
    def get_edge_length(self,edge):
        length = 0
        for i in range(len(edge)-1):
            length = length + spatial.distance.cdist(np.array(edge[i]).reshape(1,3),np.array(edge[i+1]).reshape(1,3))
        return length
            
    def read_decorators(self):
        decorators = GraphDecorators()
        lines = []
        with open(self.input_filename, 'r') as sp_gr:
            lines = sp_gr.readlines()
            # for each line
        for line in lines:
            #read all the tags
            if line.rfind("VERTEX { float[3] VertexCoordinates } @")>-1:
                decorators.decorator_vertices = int(line[line.rfind("VERTEX { float[3] VertexCoordinates } @")+len("VERTEX { float[3] VertexCoordinates } @"):]) 
                ##print(self.verts_tag)
            if line.rfind("VERTEX { int GraphLabels } @")>-1:
                decorators.decorator_vertex_labels = int(line[line.rfind("VERTEX { int GraphLabels } @")+len("VERTEX { int GraphLabels } @"):]) 
                ##print(verts_labels_tag)
            if line.rfind("EDGE { int[2] EdgeConnectivity } @")>-1:
                decorators.decorator_edge_connectivity = int(line[line.rfind("EDGE { int[2] EdgeConnectivity } @")+len("EDGE { int[2] EdgeConnectivity } @"):]) 
                ##print(edge_connectivity_tag)
            if line.rfind("EDGE { int NumEdgePoints } @")>-1:
                decorators.decorator_num_pts_per_edge = int(line[line.rfind("EDGE { int NumEdgePoints } @")+len("EDGE { int NumEdgePoints } @"):]) 
                ##print(pts_per_edge_tag)
                #print edge_num_points_id
            if line.rfind("EDGE { int EdgeLabels } @")>-1:
                decorators.decorator_edge_labels = int(line[line.rfind("EDGE { int EdgeLabels } @")+len("EDGE { int EdgeLabels } @"):])
                ##print(edge_labels_tag)
            if line.rfind("EDGE { int GraphLabels } @")>-1:
                decorators.decorator_edge_labels = int(line[line.rfind("EDGE { int GraphLabels } @")+len("EDGE { int GraphLabels } @"):]) 
                ##print(edge_labels_tag)
                #print edge_ids_id
            if line.rfind("POINT { float[3] EdgePointCoordinates } @")>-1:
                decorators.decorator_edge_pt_coords = int(line[line.rfind("POINT { float[3] EdgePointCoordinates } @")+len("POINT { float[3] EdgePointCoordinates } @"):]) 
                ##print(edge_pts_tag)
                #print edge_points_id
            if line.rfind("POINT { float Radius } @")>-1:
                decorators.edge_pts_radii_tag = int(line[line.rfind("POINT { float Radius } @")+len("POINT { float Radius } @"):])
                ##print(edge_pts_radii_tag)
        
        return decorators
        
    def read_spatial_graph(self):
        ''' Just read all the vert and edge coords in the sp graph '''
        
        class GraphReadingStates(Enum):
            reading_invalid = 0
            reading_vertices = 1
            reading_vertex_labels = 2
            reading_edge_connectivity = 3
            reading_num_pts_per_edge = 4
            reading_edge_labels = 5
            reading_edge_coords = 6
            reading_pt_radii = 7
        #reading_state = GraphReadingStates
        
        # instantiate the graph datastructure
        graph_data = SpatialGraphData()
        ##print(len(graph_data.vertices))
        lines = []
        with open(self.input_filename, 'r') as sp_gr:
            lines = sp_gr.readlines()
            
        # for each line
        for line in lines:
            ##print(line)
            # read transformation matrix if at all
            if not line.find("TransformationMatrix ") == -1:
                graph_data.transformation_matrix = (line[line.rfind("TransformationMatrix ")+len("TransformationMatrix "):-1])
                self.transformation_matrix = list(map(float, graph_data.transformation_matrix.split(" ")[:]))
                #print (self.transformation_matrix)

            # set the states as per tags encountered
            if line.startswith('@{}'.format(self.decorators.decorator_vertices)):
                ##print('found')
                GraphReadingStates = GraphReadingStates.reading_vertices
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_vertex_labels)) == 0:
                GraphReadingStates = GraphReadingStates.reading_vertex_labels
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_edge_connectivity)) == 0:
                GraphReadingStates = GraphReadingStates.reading_edge_connectivity
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_num_pts_per_edge)) == 0:
                GraphReadingStates = GraphReadingStates.reading_num_pts_per_edge
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_edge_labels)) == 0:
                GraphReadingStates = GraphReadingStates.reading_edge_labels
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_edge_pt_coords)) == 0:
                GraphReadingStates = GraphReadingStates.reading_edge_coords
                #print(GraphReadingStates)
                continue
            if line.find("@{}".format(self.decorators.decorator_edge_pt_radii)) == 0:
                GraphReadingStates = GraphReadingStates.reading_pt_radii
                #print(GraphReadingStates)
                continue


            # now read the actual elements 
            # skip this iteration if encountered a black line or next decorator
            if GraphReadingStates != GraphReadingStates.reading_invalid:
                if (line.rfind("@") == 0 or line.isspace()): 
                    GraphReadingStates = GraphReadingStates.reading_invalid
                    continue

            if GraphReadingStates == GraphReadingStates.reading_vertices :
                pt = list(map(float,line.split()))
                pt = [pt[0]*self.axis_directions[0]*self.scaling[0],pt[1]*self.axis_directions[1]*self.scaling[1],pt[2]*self.axis_directions[2]*self.scaling[2]]
                graph_data.vertices.append(pt)
                ##print(line)

            if GraphReadingStates == GraphReadingStates.reading_vertex_labels :
                graph_data.vertex_labels.append(int(line))

            if GraphReadingStates == GraphReadingStates.reading_edge_connectivity :
                graph_data.edge_connectivity.append(list(map(int,line.split())))

            if GraphReadingStates == GraphReadingStates.reading_num_pts_per_edge :
                graph_data.num_pts_per_edge.append( int(line))

            if GraphReadingStates == GraphReadingStates.reading_edge_labels :
                graph_data.edge_labels.append( int(line))

            if GraphReadingStates == GraphReadingStates.reading_edge_coords :
                pt = list(map(float,line.split()))
                pt = [pt[0]*self.axis_directions[0]*self.scaling[0],pt[1]*self.axis_directions[1]*self.scaling[1],pt[2]*self.axis_directions[2]*self.scaling[2]]
                graph_data.edge_pt_coords.append(pt)

            if GraphReadingStates == GraphReadingStates.reading_pt_radii :
                graph_data.edge_pt_radii.append( float(line))
                    
        self.create_edges(graph_data)
                
        return graph_data
    
    def create_edges(self,graph_data):
        # populate edges from the raw edge coords list 
        edge_pt_start_ind = 0
        for edge_id in range(len(graph_data.num_pts_per_edge)):
            edge = graph_data.edge_pt_coords[edge_pt_start_ind: edge_pt_start_ind + graph_data.num_pts_per_edge[edge_id]]
            graph_data.edge_list.append(graph_data.edge_pt_coords[edge_pt_start_ind: edge_pt_start_ind + graph_data.num_pts_per_edge[edge_id]])
            ##print(edge_point_coords[edge_pt_start_ind:edge_pt_start_ind+edge_num_points[edge_id]])
            edge_pt_start_ind = edge_pt_start_ind + graph_data.num_pts_per_edge[edge_id]
            
    def get_labelled_subgraph(self,alternate_labels,label_id=0,label_color=[0,0,0]):
        '''
        This function will return the labelled subgraph
        ''' 
        if self.file_format == 'ASC':
            return self.get_labelled_subgraph_from_asc_file(alternate_labels,label_id,label_color)
            
        else:
            # get the label id
            #for alternate_labels in label_list:
            label_tree_node = self.get_label_tree_node(alternate_labels)
            if label_tree_node != None:
                #label_name,label_id,label_color = label_tree_node 
                return self.populate_labelled_subgraph(label_tree_node)            
            return SpatialGraphData()
    
    def get_label_tree_node(self,alternate_labels):
        for label_name in alternate_labels:
            label_tree_node = self.search_label_tree_node(label_name) 
            #print(label_tree_node)
            if label_tree_node != None:
                #label_name,label_id,label_color = label_tree_node 
                return label_tree_node            
        return None
        
    def get_labelled_subgraph_from_asc_file(self,label_names,label_id=0,label_color=[0,0,0]):
        class HocReadingStates(Enum):
            hoc_reading_idle = 0
            hoc_reading_root_edge = 1
            hoc_reading_child_edge = 2
            hoc_reading_contour = 3

        tree_list  = []
        master_graphdata = SpatialGraphData()
        with open(self.input_filename,'r') as fp:
            HocReadingStates = HocReadingStates.hoc_reading_idle
            edge =[]
            vert_from = []
            first_edge = False
            lines = fp.readlines()
            graph_data = SpatialGraphData()
            edge_names_list = []
            for i in range(len(lines)):
                line = lines[i]
                #print(line)
                #print(HocReadingStates)
                if HocReadingStates == HocReadingStates.hoc_reading_idle:
                    if line.startswith('  ('+(label_names[0])) or line.startswith('  ('+(label_names[1])):
                        graph_data = SpatialGraphData()
                        edge_names_list = []
                        edge = []
                        #print('found tree')
                        if lines[i+1].endswith('Root\n'):
                            # this is the root edge
                            #print('first_edge')
                            HocReadingStates = HocReadingStates.hoc_reading_root_edge
                            edge_names_list.append('R')

                    elif line.startswith('("'+label_names[0]) or line.startswith('("'+label_names[1]):
                            # reading the contour
                            #print('reading contour')
                            HocReadingStates = HocReadingStates.hoc_reading_contour
                            graph_data = SpatialGraphData()
                            edge = []
                            edge_names_list = []
                elif HocReadingStates == HocReadingStates.hoc_reading_contour:
                    if line.find('End of contour')>0:
                        HocReadingStates = HocReadingStates.hoc_reading_idle
                        # process the edge that we just read
                        vert_from = edge[0]
                        #print(len(edge))
                        self.process_edge(graph_data,edge,label_id,vert_from=vert_from,vert_to=edge[len(edge)-1],                                     first_edge=True)
                        edge =[]
                        vert_from = []

                        if len(graph_data.vertices) > 0:
                            tree_list.append(graph_data)
                            graph_data = SpatialGraphData()
                    elif len(line.split(';')) > 1:
                        #print(line)
                        x = line.split('(')[1].lstrip().split(' ')[0]
                        y = (line[line.find(x)+len(x):].lstrip()).split(' ')[0]
                        z = (line[line.find(y)+len(y):].lstrip()).split(' ')[0]
                        pt = [float(x),float(y),float(z)]
                        edge.append(pt)

                elif HocReadingStates == HocReadingStates.hoc_reading_root_edge:
                    if line.find('(\n')>0 :
                        HocReadingStates = HocReadingStates.hoc_reading_child_edge
                        # process the edge that we just read
                        vert_from = edge[0]
                        #print(len(edge))
                        self.process_edge(graph_data,edge,label_id,vert_from=vert_from,vert_to=edge[len(edge)-1],                                     first_edge=True)
                        edge =[]
                        vert_from = []

                    elif line.startswith('\n') or line.find('End of tree')>0:
                        HocReadingStates = HocReadingStates.hoc_reading_idle
                        # process the edge that we just read
                        if len(edge)>0:
                            vert_from = edge[0]
                            #print(len(edge))
                            self.process_edge(graph_data,edge,label_id,vert_from=vert_from,vert_to=edge[len(edge)-1],                                         first_edge=True)
                            edge =[]
                            vert_from = []
                            tree_list.append(graph_data)
                            graph_data = SpatialGraphData()
                            edge_names_list = []
                    elif len(line.split(';')) > 1:
                        #print(line)
                        x = line.split('(')[1].lstrip().split(' ')[0]
                        y = (line[line.find(x)+len(x):].lstrip()).split(' ')[0]
                        z = (line[line.find(y)+len(y):].lstrip()).split(' ')[0]
                        pt = [float(x),float(y),float(z)]
                        edge.append(pt)

                elif HocReadingStates == HocReadingStates.hoc_reading_child_edge:
                    #print(len(line.split(' ')))
                    if line.startswith('\n') or line.endswith('|\n') or line.endswith('Normal\n') or                        line.endswith('(\n')  or line.find('End of split')>0 or line.find('Generated')>0:
                        #print(line)
                        if len(edge)>0:
                            #print(edge_names_list)
                            current_edge_name = edge_names_list[len(edge_names_list)-1]
                            parent_edge_name = current_edge_name[0:-2]
                            #print(current_edge_name,parent_edge_name)
                            parent_edge_ind = self.find_ind_for_element(edge_names_list,parent_edge_name)
                            #print(parent_edge_ind)
                            parent_edge = graph_data.edge_list[parent_edge_ind]
                            vert_from = parent_edge[len(parent_edge)-1]
                            self.process_edge(graph_data,edge,label_id,vert_from=vert_from,vert_to=edge[len(edge)-1])
                            edge =[]
                            vert_from = []
                    elif line.find('End of tree')>0:
                        HocReadingStates = HocReadingStates.hoc_reading_idle

                        # all this tree to the sg 
                        if len(graph_data.vertices) > 0:
                            tree_list.append(graph_data)

                    elif len(line.split(' ')) > 9:
                        if line.find('R') > 0:
                            # this is the first pt of edge.. there store the prev edge connection info
                            ind = line.find('R-')
                            edge_names_list.append(str(line[ind:-1]))

                        #print(line)
                        x = line.split('(')[1].lstrip().split(' ')[0]
                        y = (line[line.find(x)+len(x):].lstrip()).split(' ')[0]
                        z = (line[line.find(y)+len(y):].lstrip()).split(' ')[0]
                        pt = [float(x),float(y),float(z)]
                        edge.append(pt)

        for i in range(len(tree_list)):
            master_graphdata = self.combine_subgraphs([master_graphdata,tree_list[i]])
        
        return master_graphdata             
   
    def process_edge(self,graph_data,edge,label_id,vert_from_ind=0,vert_to_ind=0,vert_from=[],vert_to=[],first_edge=False,contour=False):
        #if len(edge)==1:
        #    edge.append(edge[0])#(edge,label_id,vert_from,vert_from_ind,vert_to,vert_to_ind)
        if not first_edge:
            edge.insert(0,vert_from)
        graph_data.edge_list.append(edge)
        graph_data.edge_labels.append(label_id)
        #print(vert_from,vert_to)
        graph_data.num_pts_per_edge.append(len(edge))

    #     if contour:
    #         # close the loop as it is a contour
    #         graph_data.vertices.append(edge[0])
    #         graph_data.vertices.append(edge[len(edge)-1])
    #         graph_data.edge_connectivity.append([len(graph_data.vertices)-2,len(graph_data.vertices)-1])
    #         graph_data.vertex_labels.append(label_id)
    #         graph_data.vertex_labels.append(label_id)

    #     else:
        if first_edge:
            graph_data.vertices.append(vert_from)
            graph_data.vertex_labels.append(label_id)

        graph_data.vertices.append(vert_to)
        graph_data.vertex_labels.append(label_id)

        vert_from_ind = self.find_ind_for_verts(graph_data.vertices,vert_from)
        vert_to_ind = self.find_ind_for_verts(graph_data.vertices,vert_to)
        #print(vert_from_ind,vert_to_ind)
        graph_data.edge_connectivity.append([vert_from_ind,vert_to_ind])
            #graph_data.vertex_labels.append(label_id)

        for pt in edge:
            graph_data.edge_pt_coords.append(pt)
            graph_data.edge_pt_radii.append(0)
        
    def find_ind_for_verts(self,vertices,vert):
        for i in range(len(vertices)):
            if vertices[i][0] == vert[0] and vertices[i][1] == vert[1] and vertices[i][2] == vert[2]:
                return i
        return -1

    def find_ind_for_element(self,element_list,element):
        for i in range(len(element_list)):
            if element_list[i] == element :
                return i
        return -1

    def get_neuron_subgraph(self,alternate_labels):
        neuron = Neuron()
        subgraphs_list = []
        neuron.axon = self.get_labelled_subgraph(['Axon','Axon'],label_id=6,label_color=[0, 0, 1])
        subgraphs_list.append(neuron.axon)
        neuron.soma = self.get_labelled_subgraph(['Soma','CellBody'],label_id=7,label_color=[1, 0, 0])
        subgraphs_list.append(neuron.soma)
        neuron.dendrite.apical_dendrite = self.get_labelled_subgraph(['ApicalDendrite','Apical'],label_id=4,label_color=[1, 0.5, 0.5])
        subgraphs_list.append(neuron.dendrite.apical_dendrite)
        neuron.dendrite.basal_dendrite = self.get_labelled_subgraph(['BasalDendrite','Dendrite'],label_id=5,label_color=[0.8, 0.4, 0.4])
        subgraphs_list.append(neuron.dendrite.basal_dendrite)
        neuron.all_neurites_subgraphdata = self.combine_subgraphs(subgraphs_list)
        return neuron
        
    def get_barrels_subgraph(self,alternate_labels):
        barrels = Barrels()
        barrels.a_row = self.get_barrel_row('A',label_id=12,label_color=[1,0.2,0.2])
        barrels.rows_list.append(barrels.a_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,barrels.a_row.all_barrels_graphdata])
        barrels.b_row = self.get_barrel_row('B',label_id=17,label_color=[1,0.25,0.25])
        barrels.rows_list.append(barrels.b_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,barrels.b_row.all_barrels_graphdata])
        barrels.c_row = self.get_barrel_row('C',label_id=22,label_color=[1,0.3,0.3])
        barrels.rows_list.append(barrels.c_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,barrels.c_row.all_barrels_graphdata])
        barrels.d_row = self.get_barrel_row('D',label_id=29,label_color=[1,0.35,0.35])
        barrels.rows_list.append(barrels.d_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,barrels.d_row.all_barrels_graphdata])
        barrels.e_row = self.get_barrel_row('E',label_id=36,label_color=[1,0.4,0.4])
        barrels.rows_list.append(barrels.e_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,barrels.e_row.all_barrels_graphdata])
        barrels.greek_row = self.get_barrel_greek_row(['Alpha','Beta','Gamma','Delta'],label_id=43,label_color=[1,0.1,0.1])
        barrels.rows_list.append(barrels.greek_row)
        barrels.all_rows_graphdata = self.combine_subgraphs([barrels.all_rows_graphdata,\
                                                             barrels.greek_row.all_barrels_graphdata])
        barrels.all_rows_graphdata
        
        return barrels
        
    def get_barrel_row(self,row_label,label_id=0,label_color=[0,0,0]):
        barrel_row = BarrelRow()
        barrel_row.one = Barrel(self.get_labelled_subgraph([row_label+'1',row_label+'1'],label_id+1,label_color),\
                                 self.validate_barrels, self.barrel_projections_present, self.reverse_barrel_direction,row_name=row_label, column_name='1')
        barrel_row.single_barrels_list.append(barrel_row.one)
        barrel_row.two = Barrel(self.get_labelled_subgraph([row_label+'2',row_label+'2'],label_id+2,label_color), \
                                self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='2')
        barrel_row.single_barrels_list.append(barrel_row.two)
        barrel_row.three = Barrel(self.get_labelled_subgraph([row_label+'3',row_label+'3'],label_id+3,label_color), \
                                  self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='3')
        barrel_row.single_barrels_list.append(barrel_row.three)
        barrel_row.four = Barrel(self.get_labelled_subgraph([row_label+'4',row_label+'4'],label_id+4,label_color), \
                                 self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='4')
        barrel_row.single_barrels_list.append(barrel_row.four)
        barrel_row.five = Barrel(self.get_labelled_subgraph([row_label+'5',row_label+'5'],label_id+5,label_color), \
                                 self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='5')
        barrel_row.single_barrels_list.append(barrel_row.five)
        barrel_row.six = Barrel(self.get_labelled_subgraph([row_label+'6',row_label+'6'],label_id+6,label_color), \
                                self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='6')
        barrel_row.single_barrels_list.append(barrel_row.six)
        barrel_row.seven = Barrel(self.get_labelled_subgraph([row_label+'7',row_label+'7'],label_id+7,label_color), \
                                  self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name=row_label, column_name='7')
        barrel_row.single_barrels_list.append(barrel_row.seven)
        #append_list = []
        #for i in range(1,8):
        #    append_list.append(self.get_labelled_subgraph([row_label+'{}'.format(i)]))
            #barrel_row.single_barrels_list.append(self.get_labelled_subgraph([row_label+'{}'.format(i)]))
        barrel_row.all_barrels_graphdata  = self.combine_subgraphs([barrel_row.one.contours,barrel_row.two.contours,\
                                                                    barrel_row.three.contours,barrel_row.four.contours,\
                                                                    barrel_row.five.contours,barrel_row.six.contours,
                                                                   barrel_row.seven.contours]) 
        return barrel_row
    
    def get_barrel_greek_row(self,greekrow_labels,label_id=0,label_color=[0,0,0]):
        barrel_row = BarrelRow()
        barrel_row.one = Barrel(self.get_labelled_subgraph([greekrow_labels[0],'closed_Alpha'],label_id+1,label_color),\
                                self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name='Greek_Row', column_name='Alpha')
        barrel_row.single_barrels_list.append(barrel_row.one)
        barrel_row.two = Barrel(self.get_labelled_subgraph([greekrow_labels[1],greekrow_labels[1]],label_id+2,label_color),\
                                self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name='Greek_Row', column_name='Beta',)
        barrel_row.single_barrels_list.append(barrel_row.two)
        barrel_row.three = Barrel(self.get_labelled_subgraph([greekrow_labels[2],greekrow_labels[2]],label_id+3,label_color),\
                                  self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name='Greek_Row', column_name='Gamma')
        barrel_row.single_barrels_list.append(barrel_row.three)
        barrel_row.four = Barrel(self.get_labelled_subgraph([greekrow_labels[3],greekrow_labels[3]],label_id+4,label_color),\
                                 self.validate_barrels,self.barrel_projections_present,self.reverse_barrel_direction,row_name='Greek_Row', column_name='Delta')
        barrel_row.single_barrels_list.append(barrel_row.four)
        barrel_row.single_barrels_list.append(barrel_row.five)
        barrel_row.single_barrels_list.append(barrel_row.six)
        barrel_row.single_barrels_list.append(barrel_row.seven)
        barrel_row.all_barrels_graphdata  = self.combine_subgraphs([barrel_row.one.contours,barrel_row.two.contours,\
                                                                    barrel_row.three.contours,barrel_row.four.contours])
        return barrel_row
                                 
    def combine_subgraphs(self,subgraphlist):
        combined_subgraph = SpatialGraphData()
        num_vertices = 0
        for subgraph in subgraphlist:
            ##print('combining subgraph ')
            ##print(subgraph.label_id)
            if subgraph is not None and len(subgraph.vertices)>0:
                for vert in subgraph.vertices:
                    ##print('combining vert ')
                    ##print(vert)
                    combined_subgraph.vertices.append(vert)
                for vert_label in subgraph.vertex_labels:
                    combined_subgraph.vertex_labels.append(vert_label)
                for edge_connectivity in subgraph.edge_connectivity:
                    combined_subgraph.edge_connectivity.append([edge_connectivity[0]+num_vertices,edge_connectivity[1]+num_vertices])
                for num_pts_per_edge in subgraph.num_pts_per_edge:
                    combined_subgraph.num_pts_per_edge.append(num_pts_per_edge)
                for edge_label in subgraph.edge_labels:
                    combined_subgraph.edge_labels.append(edge_label)
                for pt in subgraph.edge_pt_coords:
                    combined_subgraph.edge_pt_coords.append(pt)
                for radius in subgraph.edge_pt_radii:
                    combined_subgraph.edge_pt_radii.append(radius)
                for edge in subgraph.edge_list:
                    combined_subgraph.edge_list.append(edge)
                num_vertices = num_vertices + len(subgraph.vertices)
        return combined_subgraph
                
    def split_spatial_graph_into_z_sections(self,):
        section_subgraphs = []
        neuron_section_not_found = True
        # find the number of z sections in pia and set number and z_step of sections
        z_sorted_coords = sorted( self.pia.edge_pt_coords, key=lambda coord: round(coord[2]))
        unique_z_list = np.unique(np.transpose(np.array(z_sorted_coords))[2])
        z_step_raw = np.abs(unique_z_list[0]-unique_z_list[1])
        #print(z_step_raw)
        if np.abs(z_step_raw-50) < np.abs(z_step_raw-100):
            self.z_step = 50
        else:
            self.z_step = 100
        self.number_z_sections = len(np.unique(np.transpose(np.array(z_sorted_coords))[2]))
        #print(self.z_step,self.number_z_sections)
        for z in unique_z_list:
            
            section_subgraph_pia = self.populate_section_subgraph(z,self.pia)
            #print(section_subgraph_pia.length)
            ##print('Pia subgraph points {}'.format(len(section_subgraph_pia.edge_pt_coords)))
            section_subgraph_wm = self.populate_section_subgraph(z,self.wm)
            ##print('WM subgraph points {}'.format(len(section_subgraph_wm.edge_pt_coords)))
            section_subgraph_bvs = self.populate_section_subgraph(z,self.bvs)
            ##print(len(self.barrels.all_rows_graphdata.edge_pt_coords))
            section_subgraph_barrels = self.populate_section_subgraph(z,self.barrels.all_rows_graphdata)
            ##print(len(section_subgraph_barrels.edge_pt_coords))
            if len(self.neuron.soma.edge_pt_coords) > 0:
                soma_center = np.mean(np.transpose(np.array(self.neuron.soma.edge_pt_coords))[2])
                ##print(soma_center,z,self.z_step)
                if np.abs(soma_center - z) <=  self.z_step and neuron_section_not_found:
                    neuron_section_not_found = False
                    ##print('soma section is {}'.format(z))
                    # this is the section to which neuron soma belongs
                    # so append this neuron as part of this spatial graph
                    subgraph = self.combine_subgraphs([section_subgraph_pia,section_subgraph_wm,section_subgraph_barrels,\
                                                       section_subgraph_bvs,self.neuron.all_neurites_subgraphdata])
                else:
                    subgraph = self.combine_subgraphs([section_subgraph_pia,section_subgraph_wm,section_subgraph_barrels,\
                                                       section_subgraph_bvs])
            else:
                subgraph = self.combine_subgraphs([section_subgraph_pia,section_subgraph_wm,section_subgraph_barrels,\
                                                   section_subgraph_bvs])
            subgraph.length = section_subgraph_pia.length
            #print(subgraph.length)
            section_subgraphs.append(subgraph)
        return section_subgraphs
        
    def search_label_tree_node(self,label_name):
        for node in self.label_tree:
            if node.name[0] == label_name:
                ##print(node.name[0])
                return node.name
            
    def populate_section_subgraph(self,z_coord,graph_data):
        ##print('graph data '+str(len(graph_data.edge_pt_coords)))
        sub_graph = SpatialGraphData() 
        if graph_data is not None:
            ##print(sub_graph.vertices)
            tmp_vertices = []
            #sub_graph.label_name = label_node[0]#label_name
            #sub_graph.label_id = label_node[1]#label_id
            #sub_graph.label_color = label_node[2]#label_color

            # loop through vertices for this id
            #for i in range(len(graph_data.vertex_labels)):
            #    if round(graph_data.vertices[i][2]) == z_coord:
            #        #print(graph_data.vertices[i])
            #        sub_graph.vertices.append(graph_data.vertices[i])
            #        sub_graph.vertex_labels.append(graph_data.vertex_labels[i])

            start_ind = 0
            end_ind = 0
            sub_graph_vert_ind = 0
            for i in range(len(graph_data.num_pts_per_edge)):
                if len(graph_data.edge_list) > 0:
                    edge_pt_z = graph_data.edge_list[i][0][2] 
                    if round(edge_pt_z) <= z_coord+self.z_step/2 and round(edge_pt_z) >= z_coord-self.z_step/2:
                        start_ind = end_ind
                        end_ind = start_ind + graph_data.num_pts_per_edge[i]
                        #print start_ind
                        #sub_graph.edge_connectivity.append(self.graph_data.edge_connectivity[i])
                        sub_graph.num_pts_per_edge.append(graph_data.num_pts_per_edge[i])
                        sub_graph.edge_labels.append(graph_data.edge_labels[i])
                        sub_graph.edge_list.append(graph_data.edge_list[i])
                        for pt in graph_data.edge_list[i]:
                            sub_graph.edge_pt_coords.append(pt)
                        ##print(len(sub_graph.edge_pt_coords))
                        for radius in range(len(graph_data.edge_list[i])):
                            sub_graph.edge_pt_radii.append(0)

                        # map from original graphs verts to the new one
                        #vert_from_ind_original = graph_data.edge_connectivity[i][0]
                        #vert_to_ind_original = graph_data.edge_connectivity[i][1]
                        #from_inds = self.get_vert_index(graph_data.vertices[vert_from_ind_original],sub_graph.vertices)
                        #to_inds = self.get_vert_index(graph_data.vertices[vert_to_ind_original],sub_graph.vertices)

                        #sub_graph.vertices,sub_graph.vertex_labels,from_index,to_index = self.validate_vertex_index(from_inds,to_inds,vert_from_ind_original,vert_to_ind_original,\
                        #                                           sub_graph.vertices,sub_graph.vertex_labels,sub_graph.label_id)

                         # map from original graphs verts to the new one
                        vert_from_ind_original = graph_data.edge_connectivity[i][0]
                        vert_to_ind_original = graph_data.edge_connectivity[i][1]

                        pt1 = graph_data.vertices[vert_from_ind_original]
                        pt2 = graph_data.vertices[vert_to_ind_original]

                        from_ind = self.get_vert_index_new(graph_data.vertices[vert_from_ind_original],sub_graph.vertices)
                        if from_ind < 0 or self.compare_pts(pt1,pt2):# not found add this to vert list or its a loop
                            sub_graph.vertices.append(graph_data.vertices[vert_from_ind_original])
                            sub_graph.vertex_labels.append(graph_data.label_id)
                            from_ind = len(sub_graph.vertices)-1

                        to_ind = self.get_vert_index_new(graph_data.vertices[vert_to_ind_original],sub_graph.vertices)
                        if to_ind < 0 or self.compare_pts(pt1,pt2):# not found add this to vert list
                            sub_graph.vertices.append(graph_data.vertices[vert_to_ind_original])
                            sub_graph.vertex_labels.append(graph_data.label_id)
                            to_ind = len(sub_graph.vertices)-1
                        # now add the edge connectivity with the new vertex indices
                        sub_graph.edge_connectivity.append([from_ind,to_ind])
                        if self.compute_edge_length:
                            sub_graph.length = sub_graph.length + self.get_edge_length(graph_data.edge_list[i])

            sub_graph.transformation_matrix = graph_data.transformation_matrix
        
        return sub_graph
    
    def populate_labelled_subgraph(self,label_node):
        sub_graph = SpatialGraphData() 
        tmp_vertices = []
        sub_graph.label_name = label_node[0]#label_name
        sub_graph.label_id = label_node[1]#label_id
        sub_graph.label_color = label_node[2]#label_color
        
        ##print(label_node[0],label_node[1],label_node[2])
        # loop through vertices for this id
        #for i in range(len(self.graph_data.vertex_labels)):
        #    if sub_graph.label_id == self.graph_data.vertex_labels[i]:
        #        sub_graph.vertices.append(self.graph_data.vertices[i])
        #        sub_graph.vertex_labels.append(self.graph_data.vertex_labels[i])
        
        start_ind = 0
        end_ind = 0
        sub_graph_vert_ind = 0
        for i in range(len(self.graph_data.num_pts_per_edge)):
            start_ind = end_ind
            end_ind = start_ind + self.graph_data.num_pts_per_edge[i]
            if self.graph_data.edge_labels[i] == sub_graph.label_id:
                
                #sub_graph.edge_connectivity.append(self.graph_data.edge_connectivity[i])
                edge = self.graph_data.edge_pt_coords[start_ind:end_ind]
                ##print('before resampling {}'.format(len(edge)))
                if self.resample_resolution is not None:
                    edge = self.resample_contour(edge,res=self.resample_resolution)
                sub_graph.edge_list.append(edge)
                ##print('after resampling {}'.format(len(edge)))
                ##print(edge)
                #Landmarks(pts=edge).write_landmarks('edge_{}')
                
                sub_graph.num_pts_per_edge.append(len(edge))
                sub_graph.edge_labels.append(self.graph_data.edge_labels[i])
                
                if self.compute_edge_length:
                    sub_graph.length = sub_graph.length + self.get_edge_length(edge)
                
                for pt in edge:
                    sub_graph.edge_pt_coords.append(pt)
                    sub_graph.edge_pt_radii.append(0)
                #for radius in self.graph_data.edge_pt_radii[start_ind:end_ind]:
                #    sub_graph.edge_pt_radii.append(radius)
                     
                # map from original graphs verts to the new one
                vert_from_ind_original = self.graph_data.edge_connectivity[i][0]
                vert_to_ind_original = self.graph_data.edge_connectivity[i][1]
                
                pt1 = self.graph_data.vertices[vert_from_ind_original]
                pt2 = self.graph_data.vertices[vert_to_ind_original]
                ##print(pt1,pt2)
                from_ind = self.get_vert_index_new(self.graph_data.vertices[vert_from_ind_original],sub_graph.vertices)
                if from_ind < 0 or self.compare_pts(pt1,pt2):# not found add this to vert list
                    sub_graph.vertices.append(self.graph_data.vertices[vert_from_ind_original])
                    sub_graph.vertex_labels.append(sub_graph.label_id)
                    from_ind = len(sub_graph.vertices)-1
                ##print(sub_graph.vertices)
                to_ind = self.get_vert_index_new(self.graph_data.vertices[vert_to_ind_original],sub_graph.vertices)
                if to_ind < 0 or self.compare_pts(pt1,pt2):# not found add this to vert list
                    sub_graph.vertices.append(self.graph_data.vertices[vert_to_ind_original])
                    sub_graph.vertex_labels.append(sub_graph.label_id)
                    to_ind = len(sub_graph.vertices)-1
                ##print(sub_graph.vertices)
                #sub_graph.vertices,sub_graph.vertex_labels,from_index,to_index = self.validate_vertex_index(from_inds,to_inds,vert_from_ind_original,vert_to_ind_original,\
                #                                           sub_graph.vertices,sub_graph.vertex_labels,sub_graph.label_id)
                
                # now add the edge connectivity with the new vertex indices
                sub_graph.edge_connectivity.append([from_ind,to_ind])
                ##print(from_ind,to_ind)
        sub_graph.transformation_matrix = self.graph_data.transformation_matrix
        
        ##print(sub_graph.edge_pt_coords)
        return sub_graph                      
        
    def compare_pts(self,pt1,pt2):
        if pt1[0] == pt2[0] and  pt1[1] == pt2[1] and pt1[2] == pt2[2]:
            return True
        else:
            return False
        
    def validate_vertex_index(self,from_inds,to_inds,vert_from_ind_original,vert_to_ind_original,vertices_list,vertex_labels_list,label_id):
        '''
        This function validates the new vertex indices for the subgraph:
        - takes cares of missing vertex by adding it to the vertex list
        - takes care of double vertices (which would be the case for loops)
        '''
        from_index = -1
        to_index = -1
        if len(from_inds)==0:
            # did not find the vertex from this label.. may be it is with some other label
            # so add it to this list
            # 
            vertices_list.append(self.graph_data.vertices[vert_from_ind_original])
            vertex_labels_list.append(label_id)
            from_index = len(vertices_list)-1
        if len(to_inds)==0:
            # did not find the vertex from this label.. may be it is with some other label
            # so add it to this list
            vertices_list.append(self.graph_data.vertices[vert_to_ind_original])
            vertex_labels_list.append(label_id)
            to_index = len(vertices_list)-1
        
        if len(from_inds)==1:
            from_index = from_inds[0]
        if len(to_inds)==1:
            to_index = to_inds[0]

        if len(from_inds)>1:
            # this would be the case for loops... let the lowe
            from_index = from_inds[0]
                
        if len(to_inds)>1:
            # this would be the case for loops... let the lowe
            to_index = to_inds[1]
            
        return vertices_list,vertex_labels_list,from_index,to_index
    
    def get_vert_index(self,vert,vert_list):
        inds = []
        #if len(vert)>0 and len(vert_list)>0:
        #    dist = spatial.distance_matrix(np.reshape(np.array(vert),[1,3]),vert_list)
        #    inds = np.argmin(dist,axis=1)
        #    #print(inds)
        for i in range(len(vert_list)):
            if vert_list[i][0] == vert[0] and  vert_list[i][1] == vert[1] and vert_list[i][2] == vert[2]:
                inds.append(i)
        return inds
              
    def get_vert_index_new(self,vert,vert_list):
        #inds = []
        #if len(vert)>0 and len(vert_list)>0:
        #    dist = spatial.distance_matrix(np.reshape(np.array(vert),[1,3]),vert_list)
        #    inds = np.argmin(dist,axis=1)
        #    #print(inds)
        for i in range(len(vert_list)):
            if vert_list[i][0] == vert[0] and  vert_list[i][1] == vert[1] and vert_list[i][2] == vert[2]:
                return(i)
        return -1

    def read_text(self,filename):
        with open(filename,'r') as reader:
            #for line in reader.readlines():
            return reader.readlines()
        
    def write_text(self,lines):
        with open(filename,'w') as writer:
            #for line in reader.readlines():
            writer.writelines(lines)
            
    def write_spatial_graph(self,filename):
        ''' 
        File includes following sections
        1) DEFINE section for the vertex,edge and point numbers
        2) Label Tree
        3) various tags
        4) Actual data containing verts,edges and poits
        '''
        with open(filename,'w') as wr:
            wr.writelines('# Avizo 3D ASCII 2.0\n\n\n')
            wr.writelines('define VERTEX {}\n'.format(len(self.graph_data.vertices)))
            wr.writelines('define EDGE {}\n'.format(len(self.graph_data.num_pts_per_edge)))
            wr.writelines('define POINT {}\n'.format(len(self.graph_data.edge_pt_coords)))
            wr.writelines('\n')
            
            wr.writelines(self.label_tree_text)
            wr.writelines('\n')
            
            wr.writelines('VERTEX { float[3] VertexCoordinates } @'+str(self.decorators.decorator_vertices)+'\n')
            wr.writelines('VERTEX { int GraphLabels } @'+str(self.decorators.decorator_vertex_labels)+'\n')
            wr.writelines('EDGE { int[2] EdgeConnectivity } @'+str(self.decorators.decorator_edge_connectivity)+'\n')
            wr.writelines('EDGE { int NumEdgePoints } @'+str(self.decorators.decorator_num_pts_per_edge)+'\n')
            wr.writelines('EDGE { int GraphLabels } @'+str(self.decorators.decorator_edge_labels)+'\n')
            wr.writelines('POINT { float[3] EdgePointCoordinates } @'+str(self.decorators.decorator_edge_pt_coords)+'\n')
            wr.writelines('POINT { float Radius } @'+str(self.decorators.decorator_edge_pt_radii)+'\n')
            wr.writelines('\n\n')
            wr.writelines('# Data section follows\n')
            
            wr.writelines('@{}\n'.format(self.decorators.decorator_vertices))
            for vert in self.graph_data.vertices:
                #print(vert)
                wr.writelines('{} {} {}\n'.format(vert[0],vert[1],vert[2]))
            
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_vertex_labels))
            for vert_label in self.graph_data.vertex_labels:
                wr.writelines('{}\n'.format(vert_label))
            
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_connectivity))
            for edge_conn in self.graph_data.edge_connectivity:
                wr.writelines('{} {}\n'.format(edge_conn[0],edge_conn[1]))
            
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_num_pts_per_edge))
            for num_pts in self.graph_data.num_pts_per_edge:
                wr.writelines('{}\n'.format(num_pts))
                
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_labels))
            for edge_label in self.graph_data.edge_labels:
                wr.writelines('{}\n'.format(edge_label))
            
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_coords))
            for coord in self.graph_data.edge_pt_coords:
                wr.writelines('{} {} {}\n'.format(coord[0],coord[1],coord[2]))
            
            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_radii))
            for radius in self.graph_data.edge_pt_radii:
                wr.writelines('{}\n'.format(radius))

    def write_spatial_graph_2(self,filename):
        '''
        File includes following sections
        1) DEFINE section for the vertex,edge and point numbers
        2) Label Tree
        3) various tags
        4) Actual data containing verts,edges and poits
        '''
        with open(filename,'w') as wr:
            wr.writelines('# Avizo 3D ASCII 2.0\n\n\n')
            wr.writelines('define VERTEX {}\n'.format(len(self.graph_data.vertices)))
            wr.writelines('define EDGE {}\n'.format(len(self.graph_data.num_pts_per_edge)))
            wr.writelines('define POINT {}\n'.format(len(self.graph_data.edge_pt_coords)))
            wr.writelines('\n')

            wr.writelines(self.label_tree_text)
            wr.writelines('\n')

            wr.writelines('VERTEX { float[3] VertexCoordinates } @'+str(self.decorators.decorator_vertices)+'\n')
            wr.writelines('VERTEX { int LabelGroup } @'+str(self.decorators.decorator_vertex_labels)+'\n')
            wr.writelines('EDGE { int[2] EdgeConnectivity } @'+str(self.decorators.decorator_edge_connectivity)+'\n')
            wr.writelines('EDGE { int NumEdgePoints } @'+str(self.decorators.decorator_num_pts_per_edge)+'\n')
            wr.writelines('EDGE { int LabelGroup } @'+str(self.decorators.decorator_edge_labels)+'\n')
            wr.writelines('POINT { float[3] EdgePointCoordinates } @'+str(self.decorators.decorator_edge_pt_coords)+'\n')
            wr.writelines('POINT { float thickness } @'+str(self.decorators.decorator_edge_pt_radii)+'\n')
            wr.writelines('\n\n')
            wr.writelines('# Data section follows\n')

            wr.writelines('@{}\n'.format(self.decorators.decorator_vertices))
            for vert in self.graph_data.vertices:
                #print(vert)
                wr.writelines('{} {} {}\n'.format(vert[0],vert[1],vert[2]))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_vertex_labels))
            for vert_label in self.graph_data.vertex_labels:
                wr.writelines('{}\n'.format(vert_label))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_connectivity))
            for edge_conn in self.graph_data.edge_connectivity:
                wr.writelines('{} {}\n'.format(edge_conn[0],edge_conn[1]))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_num_pts_per_edge))
            for num_pts in self.graph_data.num_pts_per_edge:
                wr.writelines('{}\n'.format(num_pts))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_labels))
            for edge_label in self.graph_data.edge_labels:
                wr.writelines('{}\n'.format(edge_label))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_coords))
            for coord in self.graph_data.edge_pt_coords:
                wr.writelines('{} {} {}\n'.format(coord[0],coord[1],coord[2]))

            wr.writelines('\n')
            wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_radii))
            for radius in self.graph_data.edge_pt_radii:
                wr.writelines('{}\n'.format(radius))

    def write_sub_spatial_graph(self,graph_data,filename):
        ''' 
        File includes following sections
        1) DEFINE section for the vertex,edge and point numbers
        2) Label Tree
        3) various tags
        4) Actual data containing verts,edges and poits
        '''
        #graph_data = self.get_labelled_subgraph([label_name])
        if graph_data is not None:
            if len(graph_data.vertices) > 0:
                with open(filename,'w') as wr:
                    wr.writelines('# Avizo 3D ASCII 2.0\n\n\n')
                    wr.writelines('define VERTEX {}\n'.format(len(graph_data.vertices)))
                    wr.writelines('define EDGE {}\n'.format(len(graph_data.num_pts_per_edge)))
                    wr.writelines('define POINT {}\n'.format(len(graph_data.edge_pt_coords)))
                    wr.writelines('\n')

                    wr.writelines(self.label_tree_text)
                    wr.writelines('\n')

                    wr.writelines('VERTEX { float[3] VertexCoordinates } @'+str(self.decorators.decorator_vertices)+'\n')
                    wr.writelines('VERTEX { int GraphLabels } @'+str(self.decorators.decorator_vertex_labels)+'\n')
                    wr.writelines('EDGE { int[2] EdgeConnectivity } @'+str(self.decorators.decorator_edge_connectivity)+'\n')
                    wr.writelines('EDGE { int NumEdgePoints } @'+str(self.decorators.decorator_num_pts_per_edge)+'\n')
                    wr.writelines('EDGE { int GraphLabels } @'+str(self.decorators.decorator_edge_labels)+'\n')
                    wr.writelines('POINT { float[3] EdgePointCoordinates } @'+str(self.decorators.decorator_edge_pt_coords)+'\n')
                    wr.writelines('POINT { float Radius } @'+str(self.decorators.decorator_edge_pt_radii)+'\n')
                    wr.writelines('\n\n')
                    wr.writelines('# Data section follows\n')

                    wr.writelines('@{}\n'.format(self.decorators.decorator_vertices))
                    for vert in graph_data.vertices:
                        wr.writelines('{} {} {}\n'.format(vert[0],vert[1],vert[2]))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_vertex_labels))
                    for vert_label in graph_data.vertex_labels:
                        wr.writelines('{}\n'.format(vert_label))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_edge_connectivity))
                    for edge_conn in graph_data.edge_connectivity:
                        wr.writelines('{} {}\n'.format(edge_conn[0],edge_conn[1]))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_num_pts_per_edge))
                    for num_pts in graph_data.num_pts_per_edge:
                        wr.writelines('{}\n'.format(num_pts))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_edge_labels))
                    for edge_label in graph_data.edge_labels:
                        wr.writelines('{}\n'.format(edge_label))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_coords))
                    for coord in graph_data.edge_pt_coords:
                        wr.writelines('{} {} {}\n'.format(coord[0],coord[1],coord[2]))

                    wr.writelines('\n')
                    wr.writelines('@{}\n'.format(self.decorators.decorator_edge_pt_radii))
                    for radius in graph_data.edge_pt_radii:
                        wr.writelines('{}\n'.format(radius))
            else:
                print('Label not found !!!')
        else:
            print('Label not found !!!')

    '''Convert Amira spatial graph into networkx format'''
    def convert_nx(self):
        g = nx.DiGraph()

        edges = []
        for edge_conn in self.graph_data.edge_connectivity:
            edges.append((edge_conn[0], edge_conn[1]))
        g.add_edges_from(edges)

        edgeIndex = 0
        edgePointIdx = 0
        radii = self.graph_data.edge_pt_radii
        coords = self.graph_data.edge_pt_coords
        edgeLabels = self.graph_data.edge_labels
        for num_pts in self.graph_data.num_pts_per_edge:
            edge = edges[edgeIndex]
            edgeLabelIdx = edgeLabels[edgeIndex]
            edgeLabel = self.getLabelName(edgeLabelIdx)
            edgePoints = []
            for k in range(0,num_pts):
                coord = coords[edgePointIdx]
                radius = radii[edgePointIdx]
                edgePoint = [coord[0], coord[1], coord[2], radius]
                edgePoints.append(edgePoint)
                edgePointIdx += 1
            edgeIndex +=1
            g.edges[edge[0], edge[1]]["label"] = edgeLabel
            g.edges[edge[0], edge[1]]["points"] = np.asarray(edgePoints)

        return g


    def write_all_sub_spatial_graphs(self,filename):
        ##print(filename)
        self.write_sub_spatial_graph(self.pia,filename[0:len(filename)-1]+'_pia.am')
        self.write_sub_spatial_graph(self.wm,filename[0:len(filename)-1]+'_wm.am')
        self.write_sub_spatial_graph(self.bvs,filename[0:len(filename)-1]+'_bvs.am')
        self.write_sub_spatial_graph(self.barrels.all_rows_graphdata,filename[0:len(filename)-1]+'_barrels.am')
        self.write_sub_spatial_graph(self.neuron.all_neurites_subgraphdata,filename[0:len(filename)-1]+'_neuron.am')
        
    def write_neuron_sub_spatial_graph(self,filename):
        self.write_sub_spatial_graph(self.neuron.all_neurites_subgraphdata,filename[0:len(filename)-1]+'_neuron.am')
        
    def write_all_section_spatial_graphs(self,filename,reverse_direction=False):
        # check the direction of the spatial graph
        # make sure that the smallest section has the least section number
        #if len(self.section_subgraphs[0].edge_pt_coords) < \
        #    len(self.section_subgraphs[self.number_z_sections-1].edge_pt_coords)
        if self.compute_edge_length:
            # decide the z direciton based on the edge lengths of first and last pia edges
            if self.section_subgraphs[0].length < self.section_subgraphs[self.number_z_sections-1].length:
                # sec 1 is the smallest write as it is
                for i in range(self.number_z_sections):
                    self.write_sub_spatial_graph(self.section_subgraphs[i],filename+'{}.am'.format(i+1))
            else:
                # should reverse the oder of writing
                ##print(self.number_z_sections)
                for i in range(self.number_z_sections):
                    ##print(filename+'{}.am'.format(self.number_z_sections-i))
                    self.write_sub_spatial_graph(self.section_subgraphs[i],filename+'{}.am'.format(self.number_z_sections-i))
        else:
            if not self.reverse_direction:
                # sec 1 is the smallest write as it is
                for i in range(self.number_z_sections):
                    self.write_sub_spatial_graph(self.section_subgraphs[i],filename+'{}.am'.format(i+1))
            else:
                # should reverse the oder of writing
                ##print(self.number_z_sections)
                for i in range(self.number_z_sections):
                    ##print(filename+'{}.am'.format(self.number_z_sections-i))
                    self.write_sub_spatial_graph(self.section_subgraphs[i],filename+'{}.am'.format(self.number_z_sections-i))

    def write_section_spatial_graphs_as_single_graph(self,filename,limit=None):
        merged_graph = SpatialGraphData()
        if limit is None:
            for i in range(self.number_z_sections):
                merged_graph = self.combine_subgraphs([merged_graph,sp.section_subgraphs[i]])
        else:
            for i in range(self.number_z_sections):
                merged_graph = self.combine_subgraphs([merged_graph,sp.section_subgraphs[i]])
        self.write_sub_spatial_graph(merged_graph,filename)
        
    def set_z_coord(self,z,to_neuron=False):
        #self.graph_data = self.set_z_coord_to_subgraph(self.graph_data,z)
        if to_neuron:
            self.neuron.all_neurites_subgraphdata = self.set_z_coord_to_subgraph(self.neuron.all_neurites_subgraphdata,z,translate=True)
        self.pia = self.set_z_coord_to_subgraph(self.pia,z)
        self.wm = self.set_z_coord_to_subgraph(self.wm,z)
        self.bvs = self.set_z_coord_to_subgraph(self.bvs,z)
        self.barrels.all_rows_graphdata = self.set_z_coord_to_subgraph(self.barrels.all_rows_graphdata,z)
        
        self.update_graph_data()
            
    def set_z_coord_to_subgraph(self,graph_data,z,translate=False):
        for i in range(len(graph_data.vertices)):
            pt = graph_data.vertices[i]
            if translate:
                if len(graph_data.vertices[i]) > 3:
                    graph_data.vertices[i] = [pt[0],pt[1],pt[2]-(self.pia.vertices[0][2]-z),0]
                else:
                    graph_data.vertices[i] = [pt[0],pt[1],pt[2]-(self.pia.vertices[0][2]-z)]
            else:
                if len(graph_data.vertices[i]) > 3:
                    graph_data.vertices[i] = [pt[0],pt[1],z,0]
                else:
                    graph_data.vertices[i] = [pt[0],pt[1],z]
        for i in range(len(graph_data.edge_pt_coords)):
            pt = graph_data.edge_pt_coords[i]
            if translate:
                if len(graph_data.edge_pt_coords[i]) > 3:
                    graph_data.edge_pt_coords[i] = [pt[0],pt[1],pt[2]-(self.pia.vertices[0][2]-z),0]
                else:
                    graph_data.edge_pt_coords[i] = [pt[0],pt[1],pt[2]-(self.pia.vertices[0][2]-z)]
            else:
                if len(graph_data.edge_pt_coords[i]) > 3:
                    graph_data.edge_pt_coords[i] = [pt[0],pt[1],z,0]
                else:
                    graph_data.edge_pt_coords[i] = [pt[0],pt[1],z]
            
        return graph_data
    
    def set_scaling(self,scaling,to_neuron=False):
        #self.graph_data = self.set_z_coord_to_subgraph(self.graph_data,z)
        self.pia = self.set_scaling_to_subgraph(self.pia,scaling)
        self.wm = self.set_scaling_to_subgraph(self.wm,scaling)
        self.bvs = self.set_scaling_to_subgraph(self.bvs,scaling)
        self.barrels.all_rows_graphdata = self.set_scaling_to_subgraph(self.barrels.all_rows_graphdata,scaling)
        if to_neuron:
            self.neuron.all_neurites_subgraphdata = self.set_scaling_to_subgraph(self.neuron.all_neurites_subgraphdata,scaling)
        self.update_graph_data()
            
    def set_scaling_to_subgraph(self,graph_data,scaling):
        ##print(scaling)
        for i in range(len(graph_data.vertices)):
            pt = graph_data.vertices[i]
            graph_data.vertices[i] = [pt[0]*scaling[0],pt[1]*scaling[1],pt[2]*scaling[2]]
        for i in range(len(graph_data.edge_pt_coords)):
            pt = graph_data.edge_pt_coords[i]
            graph_data.edge_pt_coords[i] = [pt[0]*scaling[0],pt[1]*scaling[1],pt[2]*scaling[2]]
        return graph_data
    
    def update_graph_data(self,):
        #list_to_be_sent = [self.pia,self.wm,\
        #                   self.bvs,self.barrels.all_rows_graphdata,\
        #                   self.neuron.all_neurites_subgraphdata,\
        #                   self.rabies_landmarks,self.rabies_landmarks_centroid,\
        #                   self.rabies_S1_landmarks,self.rabies_S1_landmarks_centroid,\
        #                   self.rabies_M1_landmarks,self.rabies_M1_landmarks_centroid,\
        #                   self.rabies_M1_l5_landmarks,self.rabies_M1_L5_landmarks_centroid,\
        #                   self.pia_wm_axis_field,self.wm_pia_axis_field,\
        #                  self.S1_axis_field,self.M1_axis_field]
                          
        
        #for subgraph in self.labelled_subgraphs_list:
        #    list_to_be_sent.append(subgraph)
        #print(len(self.labelled_subgraphs_list))
        if len(self.labelled_subgraphs_list) > 0:
            self.graph_data = self.combine_subgraphs(self.labelled_subgraphs_list)
        else:
            #print('updating the subgraphs')
            self.graph_data = self.combine_subgraphs([self.pia,self.wm,\
                                                       self.bvs,self.barrels.all_rows_graphdata,\
                                                       self.neuron.all_neurites_subgraphdata])

    def set_label_id_in_label_tree_text(self,label_name,label_id):
        for i in range(len(self.label_tree_text)):
            line = self.label_tree_text[i]
            if label_name == line[0:line.rfind('{')].rstrip().lstrip():
                index = self.label_tree_text[i+2].rfind('Id')
                self.label_tree_text[i+2] = self.label_tree_text[i+2][0:index+3] + label_id
                break
                
    def add_axis_field(self,verts_from=[],verts_to=[],edge_list=[],label_name='',label_id=0,label_color=[0,0,0]):
        if len(verts_from) > 0 or len(edge_list)>0:
            sub_graph = SpatialGraphData() 

            sub_graph.label_name = label_name
            sub_graph.label_id = label_id
            sub_graph.label_color = label_color

            if len(verts_from) > 0:
                for i in range(len(verts_from)):
                    sub_graph.add_edge(verts_from[i],verts_to[i],label_id)
            elif len(edge_list)>0:
                for i in range(len(edge_list)):
                    sub_graph.add_edge(edge_list[i][0],edge_list[i][1],label_id)

            self.labelled_subgraphs_list.append(sub_graph)
            self.update_graph_data()

            self.add_label_to_label_tree(label_name,label_id,label_color)
        
    def add_label_to_label_tree(self,label_name,label_id,label_color):
        # Need to insert the new graph label at the next line to GraphLabels field
        insert_index = []
        for i in range(len(self.label_tree_text)):
            if self.label_tree_text[i].find('GraphLabels ') >-1 :
                insert_index.append(i)
        index = max(insert_index) +1

        self.label_tree_text.insert(index,'        ' + str(label_name) + '{' + '\n')
        self.label_tree_text.insert(index+1,'            Color {} {} {},\n'.format(label_color[0],label_color[1],label_color[2]))
        self.label_tree_text.insert(index+2,'            Id {}\n'.format(label_id))
        self.label_tree_text.insert(index+3,'        }\n')
        
        self.label_tree = []
        self.label_tree_depth_level = 0
        self.label_tree.append(Node(['Root',0,[0,0,0]]))
        # parse the label tree text to find label id and color for the leaf nodes
        self.parse_label_tree(0)
        
    def get_label_id_for_label_name(self,label_name,):
        index = 0
        for i in range(len(self.label_tree_text)):
            if self.label_tree_text[i].find(label_name) >-1 :
                index  = i
                break
        id_line = self.label_tree_text[index+2]
        
        label_id = int(id_line.split()[1])
        return label_id

    def loadLabelMap(self):
        with open(os.path.join(os.path.dirname(__file__),"LabelTree.json")) as f:
            labelTree = json.load(f)
            revMap = {}
            for label in labelTree:
                for labelName, labelId in label.items():
                    revMap[labelId] = labelName
            self.labelMap = revMap

    def getLabelName(self, labelIndex):
        if(not self.labelMap):
            self.loadLabelMap()

        if(labelIndex not in self.labelMap.keys()):
            raise RuntimeError("Unknown label index: {}".format(labelIndex))
        else:
            return self.labelMap[labelIndex]