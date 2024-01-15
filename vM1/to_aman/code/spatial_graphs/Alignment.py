import vtk
from spatial_graphs.AmiraSpatialGraph import AmiraSpatialGraph
from spatial_graphs.Surfaces import Surface
from spatial_graphs.Landmarks import Landmarks
from spatial_graphs.Vectors import Vectors

import glob
import os
import shutil
import pathlib
import numpy as np

class Alignment:
    '''
    This class will align input surface to the reference surface
    ICP algorithm used for alignment
    '''
    def __init__(self,fixed=None,moving=None,mode='rigid',nr_iterations=10000,align_by_pca=False,lhs=True):
        self.fixed = fixed
        self.moving = moving
        self.mode = mode
        self.iterations = nr_iterations
        self.txmat = []
        self.lhs = lhs
        if fixed is not None and moving is not None:
            self.icp_transform = self.icp_transform(fixed,moving,mode,nr_iterations,align_by_pca,lhs)

    def get_vtk_transformation_matrix(self,):
        return self.icp_transform.GetMatrix()

    def get_icp_transform(self,):
        return self.icp_transform

    def get_transformation_matrix(self,):
        tr_matrix = []
        for i in range(4):
            for j in range(4):
                tr_matrix.append(self.icp_transform.GetMatrix().GetElement(j,i))
        self.txmat=tr_matrix
        return tr_matrix

    def set_transformation_matrix(self,txmat):
        self.txmat = txmat

    def icp_transform(self,fixed,moving,mode,nr_iterations,align_by_pca=False,lhs=True):
        if align_by_pca:
            ref_pca_pts = self.get_pca_for_surf(fixed,ref_coronal=True,lhs=lhs)
            pca_pts = self.get_pca_for_surf(moving,ref_coronal=False,lhs=lhs)
            icp,txmat = Landmarks(pts=pca_pts).align_landmarks(ref_pca_pts)
            #moving.apply_icp_transform(icp)
        else:
            icp = vtk.vtkIterativeClosestPointTransform()
            icp.SetSource(moving)
            icp.SetTarget(fixed)
            #icp.SetMaximumMeanDistance(0.00001)

            if mode == 'rigid':
                # rotation only
                icp.GetLandmarkTransform().SetModeToRigidBody()
            elif mode == 'similarity':
                # rotation and scaling
                icp.GetLandmarkTransform().SetModeToSimilarity()
            else:
                icp.GetLandmarkTransform().SetModeToAffine()

            icp.SetMaximumNumberOfIterations(nr_iterations)
            #if not align_by_pca_first:
            icp.StartByMatchingCentroidsOn()
            icp.Modified()
            icp.SetMaximumMeanDistance(0.0000001)
            icp.SetMeanDistanceModeToRMS()

            icp.Update()

        return icp

    def get_desired_direction(self,pca_component,positive=True):
        max_ind = np.argmax(abs(pca_component))
        if positive:
            if pca_component[max_ind]<0:
                end_pt = Vectors().create_pt_along_vector_at_given_distance(-1000,[0,0,0],pca_component)
            else:
                end_pt = Vectors().create_pt_along_vector_at_given_distance(1000,[0,0,0],pca_component)
        else:
            if pca_component[max_ind]<0:
                end_pt = Vectors().create_pt_along_vector_at_given_distance(1000,[0,0,0],pca_component)
            else:
                end_pt = Vectors().create_pt_along_vector_at_given_distance(-1000,[0,0,0],pca_component)
        return end_pt

    def get_pca_for_surf(self,surf,output_path=None,ref_coronal=False,lhs=True):
        pca_axes_pts = []
        tr_mat,rot_mat,translation,pca_components,tr_pts = \
                Vectors().get_pca_transformation_matrix(Surface(polydata=surf).get_surface_pts(),translate_to_center=True)
        centroid = Surface(polydata=surf).get_center_of_mass()
        # fix directions based on lhs and rhs
        pca_axes_pts.append(centroid)
        for i in range(len(pca_components)):
            pca_component = pca_components[i]
            max_ind = np.argmax(abs(pca_component))

            if lhs:
                if ref_coronal and i==1:
                    end_pt = self.get_desired_direction(pca_component,positive=False)
                else:
                    end_pt = self.get_desired_direction(pca_component,positive=True)

            else:# rhs
                if ref_coronal and i==2:
                    end_pt = self.get_desired_direction(pca_component,positive=False)
                else:
                    end_pt = self.get_desired_direction(pca_component,positive=True)

            pca_axes_pts.append(end_pt+centroid)
        return pca_axes_pts

    def transform_folder(self,input_path,output_path,exp=None,icp=None,txmat=None,txfilename=None,apply_polydata_transform=False,hem_specific_alignment=False):
        if icp is not None:
            self.icp_transform = icp
        if txmat is not None:
            self.txmat = txmat
        if txfilename is not None:
            self.write_transformation_matrix(txfilename)

        for file in glob.glob(input_path+'/*'):
            if hem_specific_alignment is False or \
                    (hem_specific_alignment==True and self.lhs==True and (os.path.basename(file).find('lhs')>=0)) or \
                            (hem_specific_alignment==True and self.lhs==False and (os.path.basename(file).find('rhs')>=0)):
                if os.path.isfile(file) and len(os.path.basename(file).split('.'))>1:
                    #print(os.path.basename(file).split('.'))
                    if os.path.basename(file).split('.')[1] == 'vtk':
                        surf = Surface(filename=file)
                        if apply_polydata_transform:
                            txed_surf = self.apply_polydata_transformation(surf.surface,txmat)
                            Surface(polydata=txed_surf).write_surface_mesh(output_path+'/'+os.path.basename(file))
                        else:
                            surf.apply_icp_transform(self.icp_transform)
                            surf.write_surface_mesh(output_path+'/'+os.path.basename(file))
                    elif os.path.basename(file).split('.')[1] == 'am':
                        if os.path.basename(file).find('Axis_Field') > 0 or os.path.basename(file).find('axis_field') > 0:
                            sg = AmiraSpatialGraph(file,generic_graph=True)
                        else:
                            sg = AmiraSpatialGraph(file)
                        sg.apply_transformation(self.txmat)
                        sg.write_spatial_graph(output_path+'/'+os.path.basename(file))
                    elif os.path.basename(file).split('.')[1] == 'landmarksAscii' or \
                            os.path.basename(file).split('.')[1] == 'LandmarksAscii' or \
                            os.path.basename(file).split('.')[1] == 'landmarkAscii' or \
                            os.path.basename(file).split('.')[1] == 'LandmarkAscii' :
                        l = Landmarks(file)
                        l.apply_transformation(self.txmat)
                        l.write_landmarks(output_path+'/'+os.path.basename(file))

    def transform_folders(self,input_path,output_path,exp_name=None,icp=None,txmat=None,txfilename=None,apply_polydata_transform=False,hem_specific_alignment=False):
        if icp is not None:
            self.icp_transform = icp
        if txmat is not None:
            self.txmat = txmat
        if txfilename is not None:
            self.write_transformation_matrix(txfilename)
        for folder in glob.glob(input_path+'*'):
            #if os.path.basename(folder).find('Align') < 0:
            if os.path.isdir(folder):
                output_path_final = output_path + os.path.basename(folder)
                print(output_path_final)
                #if os.path.exists(output_path_final):
                #    shutil.rmtree(output_path_final)
                #os.makedirs(output_path_final)
                pathlib.Path(output_path_final).mkdir(exist_ok=True)

                for file in glob.glob(folder+'/*'):
                    if (exp_name is not None and os.path.basename(file).find(exp_name) >= 0) or exp_name is None:
                        if hem_specific_alignment is False or \
                            (hem_specific_alignment==True and self.lhs==True and (os.path.basename(file).find('lhs')>=0)) or \
                                    (hem_specific_alignment==True and self.lhs==False and (os.path.basename(file).find('rhs')>=0)):
                            if os.path.basename(file).split('.')[-1] == 'vtk':
                                #print(file)
                                surf = Surface(filename=file)
                                if apply_polydata_transform:
                                    txed_surf = self.apply_polydata_transformation(surf.surface,txmat)
                                    Surface(polydata=txed_surf).write_surface_mesh(output_path_final+'/'+os.path.basename(file))
                                else:
                                    surf.apply_icp_transform(self.icp_transform)
                                    surf.write_surface_mesh(output_path_final+'/'+os.path.basename(file))
                            elif os.path.basename(file).split('.')[-1] == 'am':
                                #print(file)
                                if os.path.basename(file).find('Axis') >= 0 or os.path.basename(file).find('axis') >= 0\
                                        or os.path.basename(file).find('PCA') >= 0 or os.path.basename(file).find('pca') >= 0:
                                    sg = AmiraSpatialGraph(file,generic_graph=True)
                                else:
                                    sg = AmiraSpatialGraph(file)
                                sg.apply_transformation(self.txmat)
                                sg.write_spatial_graph(output_path_final+'/'+os.path.basename(file))
                            elif os.path.basename(file).split('.')[-1] == 'landmarksAscii' or \
                                    os.path.basename(file).split('.')[-1] == 'LandmarksAscii' or \
                                    os.path.basename(file).split('.')[-1] == 'landmarkAscii' or \
                                    os.path.basename(file).split('.')[-1] == 'LandmarkAscii' :
                                #print(file)
                                l = Landmarks(file)
                                l.apply_transformation(self.txmat)
                                l.write_landmarks(output_path_final+'/'+os.path.basename(file))

    def transform_folders_into_semicoronal(self,input_path,output_path,):
        '''Rotate things by 45 degree around y axis'''

        for folder in glob.glob(input_path+'*'):
            if os.path.isdir(folder):
                output_path_final = output_path + os.path.basename(folder)
                pathlib.Path(output_path_final).mkdir(exist_ok=True)
                for file in glob.glob(folder+'/*'):
                    if os.path.basename(file).startswith('MG48'):
                        txmat = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
                    elif os.path.basename(file).find('lhs') > 0:
                        txmat = [0.7071,0,-0.7071,0, 0,1,0,0, 0.7071,0,0.7071,0, 0,0,0,1]
                    elif os.path.basename(file).find('rhs') > 0:
                        txmat = [0.7071,0,0.7071,0, 0,1,0,0, -0.7071,0,0.7071,0, 0,0,0,1]
                    #if (exp_name is not None and os.path.basename(file).find(exp_name) >= 0) or exp_name is None:
                    if os.path.basename(file).split('.')[-1] == 'vtk':
                        surf = Surface(filename=file)
                        txed_surf = self.apply_polydata_transformation(surf.surface,txmat)
                        Surface(polydata=txed_surf).write_surface_mesh(output_path_final+'/'+os.path.basename(file))

                    elif os.path.basename(file).split('.')[-1] == 'am':
                        if os.path.basename(folder) == 'Axis_Fields':
                            sg = AmiraSpatialGraph(file,generic_graph=True)
                        else:
                            sg = AmiraSpatialGraph(file)
                        sg.apply_transformation(txmat)
                        sg.write_spatial_graph(output_path_final+'/'+os.path.basename(file))
                    elif os.path.basename(file).split('.')[-1] == 'landmarksAscii' or \
                            os.path.basename(file).split('.')[-1] == 'LandmarksAscii' or \
                            os.path.basename(file).split('.')[-1] == 'landmarkAscii' or \
                            os.path.basename(file).split('.')[-1] == 'LandmarkAscii' :
                        l = Landmarks(file)
                        l.apply_transformation(txmat)
                        l.write_landmarks(output_path_final+'/'+os.path.basename(file))

    def transform_rhs_into_lhs(self,input_path,output_path,):
        '''Rotate things by 45 degree around y axis'''

        for folder in glob.glob(input_path+'*'):
            if os.path.isdir(folder):
                output_path_final = output_path + os.path.basename(folder)
                pathlib.Path(output_path_final).mkdir(exist_ok=True)
                for file in glob.glob(folder+'/*'):

                    if os.path.basename(file).find('rhs') > 0:
                        txmat = [-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
                    else:
                        txmat = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
                    #if (exp_name is not None and os.path.basename(file).find(exp_name) >= 0) or exp_name is None:
                    if os.path.basename(file).split('.')[-1] == 'vtk':
                        surf = Surface(filename=file)
                        txed_surf = self.apply_polydata_transformation(surf.surface,txmat)
                        Surface(polydata=txed_surf).write_surface_mesh(output_path_final+'/'+os.path.basename(file))

                    elif os.path.basename(file).split('.')[-1] == 'am':
                        if os.path.basename(folder) == 'Axis_Fields':
                            sg = AmiraSpatialGraph(file,generic_graph=True)
                        else:
                            sg = AmiraSpatialGraph(file)
                        sg.apply_transformation(txmat)
                        sg.write_spatial_graph(output_path_final+'/'+os.path.basename(file))
                    elif os.path.basename(file).split('.')[-1] == 'landmarksAscii' or \
                            os.path.basename(file).split('.')[-1] == 'LandmarksAscii' or \
                            os.path.basename(file).split('.')[-1] == 'landmarkAscii' or \
                            os.path.basename(file).split('.')[-1] == 'LandmarkAscii' :
                        l = Landmarks(file)
                        l.apply_transformation(txmat)
                        l.write_landmarks(output_path_final+'/'+os.path.basename(file))



    def transform_surfaces_from_txmat(self,input_path,output_path,txmat):
        for file in glob.glob(input_path+'*'):
            if not os.path.isdir(file) and  os.path.basename(file).split('.')[-1] == 'vtk':
                print(file)
                txed_surf = self.apply_polydata_transformation(Surface(filename=file).surface,txmat)
                Surface(polydata=txed_surf).write_surface_mesh(output_path+os.path.basename(file))

        return


    def transform_files(self,input_path,output_path,exp_name,icp=None,txmat=None,txfilename=None,write_neuron_seperately=True):
        ret_list = []
        if icp is not None:
            self.icp_transform = icp
        if txmat is not None:
            self.txmat = txmat
        if txfilename is not None:
            self.write_transformation_matrix(txfilename)
        #for file in glob.glob(input_path+exp_name+'*'):
        if os.path.exists(input_path+exp_name+'.am'):
            sg = AmiraSpatialGraph(input_path+exp_name+'.am')
            sg.apply_transformation(txmat)
            sg.write_spatial_graph(output_path+'/'+exp_name+'.am')
            sg = AmiraSpatialGraph(output_path+'/'+exp_name+'.am')
            ret_list.append(sg)
            if write_neuron_seperately:
                sg_neuron = AmiraSpatialGraph()
                sg_neuron.graph_data = sg.combine_subgraphs([sg_neuron.graph_data,sg.get_labelled_subgraph(['ApicalDendrite']),\
                                                       sg.get_labelled_subgraph(['BasalDendrite']),sg.get_labelled_subgraph(['Axon']),\
                                                        sg.get_labelled_subgraph(['Soma'])])
                sg_neuron.write_spatial_graph(output_path+'/'+exp_name+'_reg_neuron.am')
                sg_neuron = AmiraSpatialGraph(output_path+'/'+exp_name+'_reg_neuron.am')
                ret_list.append(sg_neuron)
        else:
            ret_list.append(AmiraSpatialGraph())
            ret_list.append(AmiraSpatialGraph())

        if os.path.exists(input_path+exp_name+'_pia.vtk'):
            pia = Surface(filename=input_path+exp_name+'_pia.vtk')
            pia.apply_icp_transform(icp)
            pia.write_surface_mesh(output_path+'/'+exp_name+'_pia.vtk')
            ret_list.append(pia)
        else:
            ret_list.append(Surface())

        if os.path.exists(input_path+exp_name+'_WM.vtk'):
            wm = Surface(filename=input_path+exp_name+'_WM.vtk')
            wm.apply_icp_transform(icp)
            wm.write_surface_mesh(output_path+'/'+exp_name+'_WM.vtk')
            ret_list.append(wm)
        else:
            ret_list.append(Surface())

        if os.path.exists(input_path+exp_name+'_vM1_pia.vtk'):
            pia = Surface(filename=input_path+exp_name+'_vM1_pia.vtk')
            pia.apply_icp_transform(icp)
            pia.write_surface_mesh(output_path+'/'+exp_name+'_vM1_pia.vtk')
            ret_list.append(pia)
        else:
            ret_list.append(Surface())

        if os.path.exists(input_path+exp_name+'_vM1_WM.vtk'):
            wm = Surface(filename=input_path+exp_name+'_vM1_WM.vtk')
            wm.apply_icp_transform(icp)
            wm.write_surface_mesh(output_path+'/'+exp_name+'_vM1_WM.vtk')
            ret_list.append(wm)
        else:
            ret_list.append(Surface())

        return ret_list

    def write_transformation_matrix(self,filename):
        ''' transformation matrix in the row first from'''
        with open(filename,'w') as fp:
            tr_mat = np.reshape(self.txmat,[4,4])
            for i in range(len(tr_mat)):
                for j in range(len(tr_mat[i])):
                    fp.write('{}'.format(tr_mat[i][j]) + ', ')
            fp.write('\n')

    def read_transformation_matrix(self,filename):
            ''' transformation matrix in the row first from'''
            with open(filename,'r') as fp:
                tr_mat = []
                elements = fp.readline().split('\n')[0].split(', ')
                #print(elements)

                for i in range(len(elements)-1):
                    tr_mat.append(float(elements[i]))
                return tr_mat

    def apply_polydata_transformation(self,polydata,trmat):
        vtk_tr = vtk.vtkTransform()
        mat = vtk.vtkMatrix4x4()
        k=0
        for i in range(4):
            for j in range(4):
                mat.SetElement(j,i,trmat[k])
                k=k+1

        vtk_tr.SetMatrix(mat)

        ptr = vtk.vtkTransformPolyDataFilter()
        ptr.SetInputData(polydata)
        ptr.SetTransform(vtk_tr)
        ptr.Update()

        return ptr.GetOutput()
