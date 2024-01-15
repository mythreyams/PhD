import os
import numpy as np
from sklearn.decomposition import PCA
from scipy.spatial import distance
import vtk
#from spatial_graphs.Surfaces import Surface
#from spatial_graphs.AmiraSpatialGraph import AmiraSpatialGraph

class Vectors():
    ''' This class contains all necessary functions related to vectors.. Like finding angles etc'''
    def __intit__(self,):
        from_verts = []
        to_verts = []

    def get_angle_between_vectors(self,a,b,ignore_opposite_direction=True,keep_negative_angles=False):
        #print(a)
        #print(b)
        center_a, unit_vec_a = self.get_center_and_unit_vec(a[0],a[1])
        center_b, unit_vec_b = self.get_center_and_unit_vec(b[0],b[1])
        if keep_negative_angles:
            cos_theta = (np.arccos(np.dot(unit_vec_a,unit_vec_b))*180/np.pi)
        else:
            cos_theta = np.abs(np.arccos(np.dot(unit_vec_a,unit_vec_b))*180/np.pi)
        if ignore_opposite_direction:
            if cos_theta > 90:
                cos_theta = abs(180 - cos_theta)
        return cos_theta

    def get_angle_between_unit_vectors(self,unit_vec_a,unit_vec_b,ignore_opposite_direction=True):
        #print(a)
        #print(b)

        cos_theta = np.abs(np.arccos(np.dot(unit_vec_a,unit_vec_b))*180/np.pi)
        if ignore_opposite_direction:
            if cos_theta > 90:
                cos_theta = abs(180 - cos_theta)
        return cos_theta


    def get_vec_lengths(self,vecs):
        lengths = []
        for vec in vecs:
            a = vec[0]
            b = vec[1]
            ab = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
            #lengths.append(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))
            lengths.append(np.linalg.norm(ab))

        return lengths,np.array(lengths).mean(),np.array(lengths).std()

    def get_vec_length(self,vec):
        a = vec[0]
        b = vec[1]
        ab = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        #lengths.append(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))
        return (np.linalg.norm(ab))

    def get_negative_vector(self,vec):
        length = self.get_vec_length((vec))
        a = vec[0]
        b = vec[1]
        n = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        n_unit = n / np.linalg.norm(n)
        neg_vec = (length * -n_unit)

        return neg_vec

    def create_pt_along_vector_at_given_distance(self,dist,a,b):
        # find a unit vector along edge. add the distance given
        #ap = [pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]]
        n = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        n_unit = n / np.linalg.norm(n)
        ac = (dist * n_unit)
        c = ac + a

        return c

    def create_pt_along_unit_vector_at_given_distance(self,dist,a,n_unit):
        # find a unit vector along edge. add the distance given
        #ap = [pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]]
        ac = (dist * n_unit)
        c = ac + a

        return c

    def get_center_and_unit_vec(self,a,b):
        # find a unit vector along edge. add the distance given

        #ap = [pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]]
        n = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        n_unit = n / np.linalg.norm(n)
        center = [(b[0] + a[0])/2, (b[1] + a[1])/2,( b[2] + a[2])/2]


        return center,n_unit

    def get_unit_vec(self,a,b):
        # find a unit vector along edge. add the distance given

        #ap = [pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]]
        n = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        n_unit = n / np.linalg.norm(n)
        #center = [(b[0] + a[0])/2, (b[1] + a[1])/2,( b[2] + a[2])/2]

        return n_unit

    def get_avg_vectors(self,axes):
        unit_vecs = []
        for axis in axes:
            unit_vecs.append(Vectors().get_unit_vec(axis[0],axis[1]))
        avg_uv = np.array(unit_vecs).mean(axis=0)
        return avg_uv

    def get_distance_of_pt_from_line(self,pt,edge):
        # distance from a pt to a vector
        # Let 'p' be the pt and 'a' be one pt in the vector
        # 'p' and 'a' can be considered as vectors from origin
        # let n be a unit vector along the edge 'ab'
        # Now p-a defines a vector from a to p
        # therefore p-a.n would be the projection of 'ap' vector on 'n' called 'ac'
        # Now the distance between pt 'p' and vector 'n' is given as :
        # norm( ac - ap )

        a = edge[0]
        b = edge[1]
        ap = [pt[0] - a[0], pt[1] - a[1], pt[2] - a[2]]
        n = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
        n_unit = n / np.linalg.norm(n)
        ac = (np.dot(ap,n_unit) * n_unit)
        c = ac+a
        bc = c - b
        pc = ac - ap
        dist = np.linalg.norm( pc )

        return dist, np.linalg.norm( ac ), np.linalg.norm( bc ), c

    def get_distance_of_pt_from_plane(self,pt,vector):
        pl = vtk.vtkPlane()
        unit_normal = self.get_unit_vec(vector[0],vector[1])
        pl.SetNormal(unit_normal)
        pl.SetOrigin(vector[0])
        proj_pt = [0,0,0]
        pl.ProjectPoint(pt,vector[0],unit_normal,proj_pt)
        return distance.euclidean(pt,proj_pt)

    def get_projection_pts_to_plane(self,pts,vector,center):
        projection_pts = []
        pl = vtk.vtkPlane()
        unit_normal = self.get_unit_vec(vector[0],vector[1])
        pl.SetNormal(unit_normal)
        pl.SetOrigin(center)
        for pt in pts:
            proj_pt = [0,0,0]
            pl.ProjectPoint(pt,center,unit_normal,proj_pt)
            projection_pts.append(proj_pt)
        return projection_pts

    def get_nearest_pt_to_axis(self,axis,pt_list):
        min_dist = 99999
        min_edge = [[0,0,0],[1,1,1]]
        min_dist_from_wm = 99999
        min_dist_from_pia = 99999
        pt_on_verctor = []
        for pt in pt_list:

            dist,dist_from_wm,dist_from_pia, wm_dist_vector = self.get_distance_of_pt_from_line(pt,axis)
            if dist < min_dist:
                min_dist = dist
                min_edge = axis
                min_dist_from_wm = dist_from_wm
                min_dist_from_pia = dist_from_pia
                pt_on_vector = wm_dist_vector

        return min_edge,min_dist, min_dist_from_wm, min_dist_from_pia, pt_on_vector

    def get_farthest_pt_to_axis(self,axis,pt_list,ignore_negative_direction=False):
        max_dist = 0
        min_edge = [[0,0,0],[1,1,1]]
        min_dist_from_wm = 99999
        min_dist_from_pia = 99999
        pt_on_verctor = []
        farthest_pt = []
        for pt in pt_list:
            dist,dist_from_wm,dist_from_pia, wm_dist_vector = self.get_distance_of_pt_from_line(pt,axis)
            if dist > max_dist:
                if ignore_negative_direction:
                    # allow only if the prijection pt lies on the positive side of the axis
                    # first extend the axis so that we dont fall short
                    cyl = Surface().get_cylinder(50,axis)
                    cylinder_srf = Surface(polydata=cyl).create_delunay_surface_3d(return_hull=True,make_cube=False)
                    valid_pts,invalid_pts = Surface().get_landmarks_within_given_surface(cylinder_srf,[wm_dist_vector])
                    if len(invalid_pts) > 0:
                            #Surface(polydata=cylinder_srf).write_surface_mesh(output_root+'neg_cylinder.vtk')
                            #and wm_dist_vector[2] > axis[0][2] and wm_dist_vector[2] < end_pt[2]:
                        max_dist = dist
                        min_edge = axis
                        min_dist_from_wm = dist_from_wm
                        min_dist_from_pia = dist_from_pia
                        pt_on_vector = wm_dist_vector
                        farthest_pt = pt
                else:
                    max_dist = dist
                    min_edge = axis
                    min_dist_from_wm = dist_from_wm
                    min_dist_from_pia = dist_from_pia
                    pt_on_vector = wm_dist_vector
                    farthest_pt = pt

        return min_edge,max_dist, min_dist_from_wm, min_dist_from_pia, pt_on_vector, farthest_pt

    def get_farthest_pt_to_plane(self,axis,pt_list,):
        max_dist = 0
        min_edge = [[0,0,0],[1,1,1]]
        min_dist_from_wm = 99999
        min_dist_from_pia = 99999
        pt_on_verctor = []
        farthest_pt = []
        for pt in pt_list:
            dist = self.get_distance_of_pt_from_plane(pt,axis)
            if dist > max_dist:
                    max_dist = dist
                    farthest_pt = pt

        return farthest_pt,max_dist

    def get_nearest_pts_to_axis_field_within_radius(self,axis_field,pt_list,radius):
        final_pts_list = []
        for axis in axis_field:
            selected_pts = self.get_nearest_pts_to_axis_within_radius(axis,pt_list,radius,return_pts_only=True)
            for pt in selected_pts:
                if not self.check_if_landmark_already_present(final_pts_list,pt):
                    final_pts_list.append(pt)
        return final_pts_list

    def get_nearest_pts_to_axis_within_radius(self,axis,pt_list,radius,return_pts_only=False):
        min_dist = 99999
        min_edge = [[0,0,0],[1,1,1]]
        min_dists_from_wm = []
        min_dists_from_pia = []
        pts_on_vector = []
        selected_pts = []
        for pt in pt_list:
            #print(pt,axis)
            dist,dist_from_wm,dist_from_pia, wm_dist_vector = self.get_distance_of_pt_from_line(pt,axis)

            if dist < radius:
                #print(wm_dist_vector)
                selected_pts.append(pt)
                min_dists_from_wm.append(dist_from_wm)
                min_dists_from_pia.append(dist_from_pia)
                pts_on_vector.append(wm_dist_vector)

        if return_pts_only:
            return selected_pts
        else:
            return pts_on_vector,min_dists_from_wm,min_dists_from_pia

    def check_if_landmark_already_present(self,pt_list,pt):
        #
        for tmp_pt in pt_list:
            if tmp_pt[0] == pt[0] and tmp_pt[1] == pt[1] and tmp_pt[2] == pt[2]:
                return True
        return False

    def check_if_axis_already_present(self,axis_list,axis):
        #
        for tmp_axis in axis_list:
            if tmp_axis[0][0] == axis[0][0] and tmp_axis[0][1] == axis[0][1] and tmp_axis[0][2] == axis[0][2]\
                and tmp_axis[1][0] == axis[1][0] and tmp_axis[1][1] == axis[1][1] and tmp_axis[1][2] == axis[1][2]:
                return True
        return False

    def append(self,vectors1,vectors2,allow_duplicates=False):
        vectors = []

        for vec in vectors1:
            vectors.append(vec)

        for vec in vectors2:
            if not allow_duplicates and not self.check_if_axis_already_present(vectors,vec):
               vectors.append(vec)
        return vectors

    def get_nearest_axes_to_ref_axis_within_radius(self,ref_axis,axis_list,radius,positive_direction=False):

        nearby_axes = []

        for axis in axis_list:
            axis_mid_pt = (np.array(axis[0]) + np.array(axis[1]))/2
            dist,dist_from_wm,dist_from_pia, pt_on_vector = self.get_distance_of_pt_from_line(axis_mid_pt,ref_axis)
            if positive_direction:
                if dist < radius and pt_on_vector[0] > ref_axis[0][0]:
                    nearby_axes.append(axis)
            else:
                if dist < radius:
                    nearby_axes.append(axis)
        return nearby_axes

    def get_nearest_pts_to_ref_axis_within_radius(self,ref_axis,pt_list,radius,axis_list=[],positive_direction=True):

        dist_from_origin = 0
        nearby_pt_inds = []
        pt_farthest_from_origin_ind = 0

        for i in range(len(pt_list)):
            pt = pt_list[i]
            dist,dist_from_wm,dist_from_pia, pt_on_vector = self.get_distance_of_pt_from_line([pt[0],pt[1],pt[2]],ref_axis)
            if positive_direction:
                if dist < radius and pt_on_vector[0] > ref_axis[0][0]:
                    nearby_pt_inds.append(i)
                    if dist_from_wm > dist_from_origin:
                        dist_from_origin = dist_from_wm
                        pt_farthest_from_origin_ind = i
            else:
                if dist < radius and pt_on_vector[0] < ref_axis[0][0]:
                    nearby_pt_inds.append(i)
                    if dist_from_wm > dist_from_origin:
                        dist_from_origin = dist_from_wm
                        pt_farthest_from_origin_ind = i
                # optionally also send the nearest axis to this pt

        return nearby_pt_inds,pt_farthest_from_origin_ind

    def get_nearest_axes_to_pts_within_radius(self,axis_list,pt_list,radius,axis_validation_distance=None,axis_validation_xy_distance=None,final_validation_distance=5000):

        min_dists_from_wm = []
        min_dists_from_pia = []
        nearby_axes = []

        valid_axis_list = []
        if axis_validation_distance is not None:
            # need to first select the set of axes to pick from
            # since the get_distance_of_pt_from_line() looks for pt projection to line direction, it can find
            # line that is closest in orientation but far away physically
            if len(pt_list)>1:
                pts_center = np.array(pt_list).mean(axis=0)
            else:
                pts_center = pt_list[0]
            for axis in axis_list:
                if distance.euclidean(pts_center,axis[0]) < axis_validation_distance and \
                        distance.euclidean(pts_center,axis[1]) < axis_validation_distance:
                    valid_axis_list.append(axis)
        elif axis_validation_xy_distance is not None:
            # need to first select the set of axes to pick from
            # since the get_distance_of_pt_from_line() looks for pt projection to line direction, it can find
            # line that is closest in orientation but far away physically
            if len(pt_list)>1:
                pts_center = np.array(pt_list).mean(axis=0)
            else:
                pts_center = pt_list[0]
            for axis in axis_list:
                if distance.euclidean(pts_center[:2],axis[0][:2]) < axis_validation_xy_distance and \
                        distance.euclidean(pts_center[:2],axis[1][:2]) < axis_validation_xy_distance:
                    valid_axis_list.append(axis)
        else:
            valid_axis_list =  axis_list

        for pt in pt_list:
            min_dist = 99999
            min_edge = [[0,0,0],[1,1,1]]
            for axis in valid_axis_list:
                dist,dist_from_wm,dist_from_pia, wm_dist_vector = self.get_distance_of_pt_from_line(pt,axis)
                if dist < radius and not self.check_if_axis_already_present(nearby_axes,axis):
                    if distance.euclidean(pt,axis[0]) < final_validation_distance:
                        #print(wm_dist_vector)
                        nearby_axes.append(axis)

        return nearby_axes

    def get_nearest_axis_to_pt(self,pt,axis_list,axis_validation_distance=None):
        min_dist = 99999
        min_edge = [[0,0,0],[1,1,1]]
        min_dist_from_wm = 99999
        min_dist_from_pia = 99999
        pt_on_vector = []
        valid_axis_list = []
        if axis_validation_distance is not None:
            # need to first select the set of axes to pick from
            # since the get_distance_of_pt_from_line() looks for pt projection to line direction, it can find
            # line that is closest in orientation but far away physically
            for axis in axis_list:
                if distance.euclidean(pt,axis[0]) < axis_validation_distance and \
                        distance.euclidean(pt,axis[1]) < axis_validation_distance:
                    valid_axis_list.append(axis)
        else:
            valid_axis_list =  axis_list
        for axis in valid_axis_list:
            dist,dist_from_wm,dist_from_pia, wm_dist_vector = self.get_distance_of_pt_from_line(pt,axis)
            if dist < min_dist:
                min_dist = dist
                min_edge = axis
                min_dist_from_wm = dist_from_wm
                min_dist_from_pia = dist_from_pia
                pt_on_vector = wm_dist_vector

        return min_edge,min_dist, min_dist_from_wm, min_dist_from_pia, pt_on_vector

    def get_layer_pts_for_given_depth(self,axis_field,depth,realtive_depths=True,invert_axes=False,rel_depth_given=False):
        #
        equidist_pts = []

        if realtive_depths:
            # compute the raltive depth.. by normalizing to the avg vector in the list
            vec_lens,vec_lens_mean,vec_lens_std = self.get_vec_lengths(axis_field)

            for edge in axis_field:
                edge_len, a,b = self.get_vec_lengths([edge])
                if rel_depth_given:
                    # depth given is already a percent of vector length
                    rel_depth = edge_len[0] * depth/100
                else:
                    rel_depth = edge_len[0] * depth/vec_lens_mean
                if invert_axes:
                    equidist_pts.append(self.create_pt_along_vector_at_given_distance(rel_depth,edge[1],edge[0]))
                else:
                    equidist_pts.append(self.create_pt_along_vector_at_given_distance(rel_depth,edge[0],edge[1]))

        else:
            for edge in axis_field:
                if invert_axes:
                    equidist_pts.append(self.create_pt_along_vector_at_given_distance(depth,edge[1],edge[0]))
                else:
                    equidist_pts.append(self.create_pt_along_vector_at_given_distance(depth,edge[0],edge[1]))
        return equidist_pts

    def get_layer_pts_for_given_per_depth(self,axis_field,per_depth,invert_axes=False):
        #
        equidist_pts = []

        for edge in axis_field:
            edge_len, a,b = self.get_vec_lengths([edge])
            rel_depth = edge_len[0] * per_depth / 100
            if invert_axes:
                equidist_pts.append(self.create_pt_along_vector_at_given_distance(rel_depth,edge[1],edge[0]))
            else:
                equidist_pts.append(self.create_pt_along_vector_at_given_distance(rel_depth,edge[0],edge[1]))

        return equidist_pts

    def get_rotation_matrix(self,vec,ref_vec):
        ''' Find rotation matrix to align given vec to ref vec
            R = [cos(theta) -sin(theta) 0
                 sin(theta) cos(theta)  0
                 0          0           1 ]
        '''
        # find unit vectors
        au = vec / np.linalg.norm(vec)
        bu = ref_vec / np.linalg.norm(ref_vec)
        #print(au,bu)
        cos_theta = np.dot(au,bu)
        sin_theta = np.linalg.norm(np.cross(au,bu))

        rot_mat = [[cos_theta, -sin_theta, 0],[sin_theta, cos_theta, 0],[ 0, 0, 1]]

        return np.array(rot_mat)

    def get_rotation_matrix_around_given_axis(self,theta=180,axis_of_rotation='z'):
        ''' Find rotation matrix to align given vec to ref vec
            Rz = [cos(theta) -sin(theta) 0
                 sin(theta) cos(theta)  0
                 0          0           1 ]
        '''
        # get angle in radians
        theta = theta*np.pi / 180
        if axis_of_rotation == 'z':
            rot_mat = [[ np.cos(theta), - np.sin(theta), 0],[ np.sin(theta),  np.cos(theta), 0],[ 0, 0, 1]]

        elif axis_of_rotation == 'y':
            rot_mat = [[ np.cos(theta), 0, np.sin(theta)], [0, 1, 0 ],  [ - np.sin(theta), 0, np.cos(theta)], ]

        else:
            rot_mat = [[1,0,0], [0, np.cos(theta), - np.sin(theta)],[ 0, np.sin(theta),  np.cos(theta)]]

        return np.array(rot_mat)

    def get_pca_transformation_matrix(self,pts,translate_to_center=True,output_file = None):
        pca = PCA(n_components=3,)

        translation = np.array(pts).mean(axis=0)
        pts_entered = pts - translation
        transformed_pca = pca.fit_transform(X=(pts))


        ### Directions:
        # PCA 1 along x axis, 2 along y and 3 along z
        # if the x component of PCA 1 is negative then negate it etc
        # pca.components_ is row wise matrix.. i.e first row is first component etc

        rot_mat = np.zeros([3,3])
        for i in range(3):
            #if pca.components_[i,i] < 0:
            #    # ref is X axis
            #    rot_mat[i,:] =  - pca.components_[i,:]
            #else:
            rot_mat[i,:] =  pca.components_[i,:]

        # get it in column vector form
        rot_mat = np.transpose(rot_mat)
        # add translation
        if translate_to_center is True:
            translation_mat = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[-translation[0],-translation[1],-translation[2],1]]).reshape([4,4])

            rot_mat = np.array([[rot_mat[0][0],rot_mat[0][1],rot_mat[0][2],0],[rot_mat[1][0],rot_mat[1][1],rot_mat[1][2],0],\
                      [rot_mat[2][0],rot_mat[2][1],rot_mat[2][2],0],[0,0,0,1]]).reshape([4,4])
            tr_mat = np.matmul(translation_mat,rot_mat)
        else:
            tr_mat = [rot_mat[0][0],rot_mat[0][1],rot_mat[0][2],0],[rot_mat[1][0],rot_mat[1][1],rot_mat[1][2],0],\
                    [rot_mat[2][0],rot_mat[2][1],rot_mat[2][2],0],[0,0,0,1]

        if output_file is not None:
            pca0 = AmiraSpatialGraph()
            pca0.graph_data.add_edge(translation,np.array(pca.components_[0])*1000)
            pca0.write_spatial_graph(output_file+'_pca0.am')

            pca1 = AmiraSpatialGraph()
            pca1.graph_data.add_edge(translation,np.array(pca.components_[1])*1000)
            pca1.write_spatial_graph(output_file+'_pca1.am')

            pca2 = AmiraSpatialGraph()
            pca2.graph_data.add_edge(translation,np.array(pca.components_[2])*1000)
            pca2.write_spatial_graph(output_file+'_pca2.am')

        return tr_mat,rot_mat,translation,pca.components_,transformed_pca

    def add_pca_axes_to_spatial_graph(self,pca_components,dist,origin, spatial_graph):
        for i in range(3):
            #print(origin,pca_components[i,:])
            #edge = [origin,pca_components[i,:]]
            end_pt = self.create_pt_along_vector_at_given_distance(dist,[0,0,0],pca_components[i,:])

            spatial_graph.add_axis_field(verts_from=[origin],verts_to = [end_pt+origin],label_name='M1_PCA_{}'.format(i),label_id=93+i,label_color=[0,0,0])
            #pca_sg_pca.write_spatial_graph(filename+'_{}.am'.format(i))

    def get_vm1_axis_filed(self,vm1_rabies_hull,axis_field_sg):
        obb_rabies = vtk.vtkOBBTree()
        obb_rabies.SetDataSet(vm1_rabies_hull)
        obb_rabies.BuildLocator()

        # find the vectors within this
        points = vtk.vtkPoints()

        wm_pts = []
        pia_pts = []
        edge_list = []
        for edge in axis_field_sg.edge_list:
            code = obb_rabies.IntersectWithLine(edge[0], edge[1], points, None)
            if points.GetNumberOfPoints() > 0:
                wm_pts.append(edge[0])
                pia_pts.append(edge[1])
                edge_list.append(edge)
        return wm_pts,pia_pts,edge_list

    def get_landmarks_depth_along_axis_field(self,axis_list,pt_list,radius=200):
        ''' 1) Find the axis that is closest to the rabies pt
            2) Find the intersection pt on the axis along with the distances
            from the wm and  pia'''
        min_pts = []
        max_pts = []
        mean_pts = []
        min_dists = []
        max_dists = []
        mean_dists = []
        axes_lens = []
        landmarks_on_axes = []
        l5_upper_depth = []
        l5_lower_depth = []
        avg_l5_depth = []

        for axis in axis_list:
            pts_on_vector,min_dists_wm,min_dists_pia = self.get_nearest_pts_to_axis_within_radius(axis,pt_list,radius)
            #print(self.get_vec_lengths([axis]),min_dists_wm,min_dists_pia)
            if len(pts_on_vector) > 0:
                axis_len,a,b = self.get_vec_lengths([axis])
                axes_lens.append(axis_len[0])

                wm_ind = np.argmin(np.array(min_dists_wm))
                depth_lower = axis_len[0] - min_dists_wm[wm_ind]
                l5_lower_depth.append(depth_lower)

                pia_ind = np.argmin(np.array(min_dists_pia))
                depth_upper = min_dists_pia[pia_ind]
                l5_upper_depth.append(depth_upper)

                avg_l5_depth.append((depth_lower+depth_upper)/2)

                min_pts.append(pts_on_vector[wm_ind])
                max_pts.append(pts_on_vector[pia_ind])
                mean_pts.append(np.array(pts_on_vector).mean(axis=0))

            for pt in pts_on_vector:
                landmarks_on_axes.append(pt)

        return axes_lens, l5_upper_depth,l5_lower_depth,avg_l5_depth,landmarks_on_axes#min_pts, max_pts, mean_pts

