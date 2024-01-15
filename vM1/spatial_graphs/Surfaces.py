import vtk
import os
import numpy as np
from spatial_graphs.AmiraSpatialGraph import AmiraSpatialGraph
from spatial_graphs.Vectors import Vectors
from scipy.spatial import distance
import pandas as pd

class Surface(AmiraSpatialGraph,Vectors):
    ''' This class contains all the functions to handle surfaces'''

    def __init__(self,filename='',pts = [],polydata=None):
        if len(filename) > 0:
            self.surface = self.read_polydata(filename)
            self.transform_applied = False
        elif len(pts) > 0:
            self.surface = self.convert_points_to_polydata(pts)
            self.transform_applied = False
        elif polydata is not None:
            self.surface = polydata
            self.transform_applied = False
        else:
            self.surface = vtk.vtkPolyData()
            self.transform_applied = False

    def read_polydata(self,filename):
        ##print(os.path.basename(filename)[-3:])
        if os.path.basename(filename)[-3:] == 'vtk':
            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(filename)
            reader.Update()
            return reader.GetOutput()
        elif os.path.basename(filename)[-3:] == 'stl':
            reader = vtk.vtkSTLReader()
            reader.SetFileName(filename)
            reader.Update()
            return reader.GetOutput()

    def divide_surface(self,nr_divisions,surface,output_filename=None):
        surf_devider = vtk.vtkLoopSubdivisionFilter()
        surf_devider.SetNumberOfSubdivisions(nr_divisions)
        surf_devider.SetInputData(surface)
        surf_devider.Update()

        #print(output_filename)
        #print(surf_devider.GetOutput().GetNumberOfCells())
        self.surface = surf_devider.GetOutput()

        if output_filename is not None:
             Surface(polydata=surf_devider.GetOutput()).write_surface_mesh(output_filename)

        return surf_devider.GetOutput()

    def convert_points_to_polydata(self,pts):
        points = vtk.vtkPoints()
        points.SetDataTypeToFloat()
        vertices = vtk.vtkCellArray()

        for i in range(len(pts)):
            ##print(len(pts[i]))
            ids = points.InsertNextPoint([pts[i][0],pts[i][1],pts[i][2]])

            vertices.InsertNextCell(1,[i])
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts(vertices)
        return polydata

    def convert_points_to_vtkpoints(self,pts):
        points = vtk.vtkPoints()
        points.SetDataTypeToFloat()
        for i in range(len(pts)):
            ##print(len(pts[i]))
            ids = points.InsertNextPoint([pts[i][0],pts[i][1],pts[i][2]])
        return points

    def write_surface_mesh(self,filename):
        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(self.surface)
        writer.SetFileName(filename)
        writer.Write()

    def apply_icp_transform(self,icp):
        if self.transform_applied  is False:
            icpTransformFilter = vtk.vtkTransformPolyDataFilter()
            icpTransformFilter.SetInputData(self.surface)
            icpTransformFilter.SetTransform(icp)
            icpTransformFilter.Update()
            self.surface = icpTransformFilter.GetOutput()

    def apply_transformation(self,tr_mat,inverse=False):
        tr_mat_4x4 = np.reshape(tr_mat,[4,4])
        # Set the matrix
        mat44 = vtk.vtkMatrix4x4()
        for i in range(4):
            for j in range(4):
                mat44.SetElement(i,j,tr_mat_4x4[j,i])
        tr = vtk.vtkTransform()
        tr.SetMatrix(mat44)
        if inverse:
            tr.Inverse()
        tr.Update()
        # apply this transformation to the surface polydata
        icpTransformFilter = vtk.vtkTransformPolyDataFilter()
        icpTransformFilter.SetTransform(tr)
        icpTransformFilter.SetInputData(self.surface)
        icpTransformFilter.Update()
        self.surface = icpTransformFilter.GetOutput()

    def apply_rotation(self,rot_mat=[],translation=None):


        # now set the 4x4 transformation matrix
        tr_mat = [rot_mat[0][0],rot_mat[0][1],rot_mat[0][2],0],[rot_mat[1][0],rot_mat[1][1],rot_mat[1][2],0],\
                 [rot_mat[2][0],rot_mat[2][1],rot_mat[2][2],0],[0,0,0,1]
        tr_mat = np.reshape(tr_mat,[4,4])
        # Set the matrix
        mat44 = vtk.vtkMatrix4x4()
        for i in range(4):
            for j in range(4):
                mat44.SetElement(i,j,tr_mat[i,j])
        tr = vtk.vtkTransform()
        tr.SetMatrix(mat44)
        tr.Update()
        # apply this transformation to the surface polydata
        icpTransformFilter = vtk.vtkTransformPolyDataFilter()
        icpTransformFilter.SetTransform(tr)
        icpTransformFilter.SetInputData(self.surface)
        icpTransformFilter.Update()
        self.surface = icpTransformFilter.GetOutput()

        if translation is not None:
            # Set translation as transformation and apply to polydata
            translator = vtk.vtkTransform()
            translator.Translate(-translation)
            TransformFilter = vtk.vtkTransformPolyDataFilter()
            TransformFilter.SetInputData(self.surface)
            TransformFilter.SetTransform(translator)
            TransformFilter.Update()
            self.surface = TransformFilter.GetOutput()

    def get_surface(self,):
        return self.surface

    def append(self,surf2,surf1=None,return_surf=False):
        appender = vtk.vtkAppendPolyData()
        if surf1 is None:
            appender.AddInputData(self.surface)
        else:
            appender.AddInputData(surf1)
        appender.AddInputData(surf2)
        appender.Update()
        self.surface = appender.GetOutput()
        if return_surf:
            return appender.GetOutput()

    def get_surface_pts(self,):
        pts = []
        for i in range(self.surface.GetNumberOfPoints()):
            pts.append(self.surface.GetPoints().GetPoint(i))

        return pts

    def decimate_surface(self,reduction_ratio):
        deci = vtk.vtkDecimatePro()
        deci.SetInputData(self.surface)
        deci.SetTargetReduction(reduction_ratio)
        deci.Update()

        self.surface = deci.GetOutput()

    def get_ray_cast_surface_pts(self,center=[0,0,0],surf_list = [],theta_res=6,phi_res=6,start_phi=0,end_phi=90, \
                                 start_theta=0,end_theta=360,get_rays_sg = False,get_ray_theta=False,get_ray_phi=False,\
                                 get_rays_for_given_theta=None,get_rays_for_given_phi=None,get_boundary_rays=False,\
                                 validation_nr_pts=0):
        '''Create spherical core from which to cast ray'''
        #center = [0,0,-4000]
        pia_sphere = vtk.vtkSphereSource()
        pia_sphere.SetCenter(center)
        pia_sphere.SetRadius(1)
        pia_sphere.SetThetaResolution(theta_res)
        pia_sphere.SetPhiResolution(phi_res)
        pia_sphere.SetStartPhi(start_phi)
        pia_sphere.SetEndPhi(end_phi)
        pia_sphere.SetStartTheta(start_theta)
        pia_sphere.SetEndTheta(end_theta)
        pia_sphere.SetLatLongTessellation(True)
        pia_sphere.Update()

        pia_cell_centers = vtk.vtkCellCenters()
        pia_cell_centers.SetInputData(pia_sphere.GetOutput())
        pia_cell_centers.Update()
        pia_centroid_cell_centers = pia_cell_centers.GetOutput()
        ##print('pia_centroid_cell_centers {}'.format(pia_centroid_cell_centers))
        pia_sphere_normals = vtk.vtkPolyDataNormals()
        pia_sphere_normals.SetInputConnection(pia_sphere.GetOutputPort())
        pia_sphere_normals.ComputePointNormalsOff()
        pia_sphere_normals.ComputeCellNormalsOn()
        pia_sphere_normals.SplittingOff()
        pia_sphere_normals.FlipNormalsOff()
        pia_sphere_normals.AutoOrientNormalsOn()
        pia_sphere_normals.Update()
        pia_source_normals = pia_sphere_normals.GetOutput().GetCellData().GetNormals()


        obbPia = vtk.vtkOBBTree()
        obbPia.SetDataSet(self.surface)
        obbPia.BuildLocator()

        surf_instersecting_pts = []
        rays_sg = AmiraSpatialGraph()
        #rays_for_given_theta_sg = AmiraSpatialGraph()
        #rays_for_given_phi_sg = AmiraSpatialGraph()
        #rays_theta = []
        #rays_phi = []
        idx =-1

        #min_theta = 10000
        #max_theta = 0
        #min_phi = 10000
        #max_phi = 0

        #min_rays_theta_sg = AmiraSpatialGraph()
        #min_rays_phi_sg = AmiraSpatialGraph()
        #max_rays_theta_sg = AmiraSpatialGraph()
        #max_rays_phi_sg = AmiraSpatialGraph()

        rays_per_phi_angle_list = []
        # loop through rays and find the intersecting for on the surface
        ##print(pia_centroid_cell_centers.GetNumberOfPoints())
        phi_step = int((end_phi-start_phi)/phi_res)
        theta_step = int((end_theta-start_theta)/theta_res)
        for phi in range(start_phi,end_phi-phi_step,phi_step):
            rays_for_given_phi_sg = AmiraSpatialGraph()
            #angles_theta = []
            for theta in range(start_theta,end_theta,theta_step):

        #for idx in range(pia_centroid_cell_centers.GetNumberOfPoints()):
                idx = idx + 1
                point = pia_centroid_cell_centers.GetPoint(idx)
                ##print(point)
                normal = pia_source_normals.GetTuple(idx)

                # Calculate the 'target' of the ray based on 'RayCastLength'
                pointRayTarget = list(np.array(point) + 10000*np.array(normal))

                if len(surf_list) == 0:
                    points = vtk.vtkPoints()
                    code = obbPia.IntersectWithLine(pointRayTarget, point, points, None)
                    if points.GetNumberOfPoints() > validation_nr_pts:
                        #pts = []
                        #for i in range(points.GetNumberOfPoints()):
                        #    pts.append(points.GetPoint(i))
                        surf_instersecting_pts.append(points.GetPoint(0))
                else:
                    # need to iterate through the list of surfaces in order to find average
                    intersection_pt_list = []
                    for surf in surf_list:

                        obb = vtk.vtkOBBTree()
                        obb.SetDataSet(surf)
                        obb.BuildLocator()
                        points = vtk.vtkPoints()
                        code = obb.IntersectWithLine(pointRayTarget, point, points, None)
                        if points.GetNumberOfPoints() > 0:
                            intersection_pt_list.append(points.GetPoint(0))
                    if len(intersection_pt_list)>validation_nr_pts:
                        if len(intersection_pt_list) > 1:
                            surf_instersecting_pts.append(np.array(intersection_pt_list).mean(axis=0))
                        elif len(intersection_pt_list)==1:
                            surf_instersecting_pts.append(intersection_pt_list[0])


                    #if get_rays_sg:
                        #rays_sg.add_edge(center,list(np.array(points.GetPoint(0)) + 10000*np.array(normal)))

                    #if get_ray_theta:
                    #    rays_theta.append(theta)

                    #if get_ray_phi:
                    #    rays_phi.append(phi)

                    #if get_rays_for_given_theta is not None:
                    #    if get_rays_for_given_theta == theta:
                    #        rays_for_given_theta_sg.add_edge(center,list(np.array(points.GetPoint(0)) + 10000*np.array(normal)))

                    #if get_rays_for_given_phi is not None:
                    #    if get_rays_for_given_phi == phi:

                    #rays_for_given_phi_sg.add_edge(center,list(np.array(points.GetPoint(0)) + 10000*np.array(normal)))

            #if len(rays_for_given_phi_sg.graph_data.edge_list) > 0:
            #    rays_per_phi_angle_list.append(rays_for_given_phi_sg)

        return_list = []

        return_list.append(surf_instersecting_pts)

        if get_rays_sg:
            return_list.append(rays_sg)

        if get_ray_theta:
            return_list.append(rays_theta)

        if get_ray_phi:
            return_list.append(rays_phi)

        if get_rays_for_given_theta is not None:
            return_list.append(rays_for_given_theta_sg)

        if get_rays_for_given_phi is not None:
            return_list.append(rays_for_given_phi_sg)


        if get_boundary_rays:
            #return_list.append(min_rays_theta_sg)
            #return_list.append(max_rays_theta_sg)
            return_list.append(rays_per_phi_angle_list[0])
            return_list.append(rays_per_phi_angle_list[len(rays_per_phi_angle_list)-1])

        return return_list

    def get_ray_cast_surface_pts_avg(self,center=[0,0,0],theta_res=6,phi_res=6,start_phi=0,end_phi=90, \
                                 start_theta=0,end_theta=360,validation_nr_pts=0):
        '''Create spherical core from which to cast ray'''
        #center = [0,0,-4000]
        pia_sphere = vtk.vtkSphereSource()
        pia_sphere.SetCenter(center)
        pia_sphere.SetRadius(1)
        pia_sphere.SetThetaResolution(theta_res)
        pia_sphere.SetPhiResolution(phi_res)
        pia_sphere.SetStartPhi(start_phi)
        pia_sphere.SetEndPhi(end_phi)
        pia_sphere.SetStartTheta(start_theta)
        pia_sphere.SetEndTheta(end_theta)
        pia_sphere.SetLatLongTessellation(True)
        pia_sphere.Update()

        pia_cell_centers = vtk.vtkCellCenters()
        pia_cell_centers.SetInputData(pia_sphere.GetOutput())
        pia_cell_centers.Update()
        pia_centroid_cell_centers = pia_cell_centers.GetOutput()
        ##print('pia_centroid_cell_centers {}'.format(pia_centroid_cell_centers))
        pia_sphere_normals = vtk.vtkPolyDataNormals()
        pia_sphere_normals.SetInputConnection(pia_sphere.GetOutputPort())
        pia_sphere_normals.ComputePointNormalsOff()
        pia_sphere_normals.ComputeCellNormalsOn()
        pia_sphere_normals.SplittingOff()
        pia_sphere_normals.FlipNormalsOff()
        pia_sphere_normals.AutoOrientNormalsOn()
        pia_sphere_normals.Update()
        pia_source_normals = pia_sphere_normals.GetOutput().GetCellData().GetNormals()

        obb_pia = vtk.vtkOBBTree()
        obb_pia.SetDataSet(self.surface)
        obb_pia.BuildLocator()

        avg_pts = []
        idx =-1

        phi_step = int((end_phi-start_phi)/phi_res)
        theta_step = int((end_theta-start_theta)/theta_res)
        for phi in range(start_phi,end_phi-phi_step,phi_step):
            for theta in range(start_theta,end_theta,theta_step):
                idx = idx + 1
                point = pia_centroid_cell_centers.GetPoint(idx)

                normal = pia_source_normals.GetTuple(idx)

                # Calculate the 'target' of the ray based on 'RayCastLength'
                pointRayTarget = list(np.array(point) + 10000*np.array(normal))

                # need to iterate through the list of surfaces in order to find average
                intersection_pt_list = []

                points = vtk.vtkPoints()
                code = obb_pia.IntersectWithLine(pointRayTarget, point, points, None)
                intersection_pt_list = []
                for pt_ind in range(points.GetNumberOfPoints()):
                    intersection_pt_list.append(points.GetPoint(pt_ind))
                if len(intersection_pt_list) > validation_nr_pts:
                    avg_pts.append(np.array(intersection_pt_list).mean(axis=0))

        return avg_pts

    def clip_surface(self, path,pattern, op_filename = None):

        # read the each ray from ray burst to get last z contour
        zmin_pts = []
        for file in glob.glob(path+pattern):
            l = Landmarks(file)
            zmin_pts.append(sorted(l.pts,key=lambda x:x[2])[0])

        # find extreme pts of the last z contour to get make a square
        center = np.mean(np.array(zmin_pts),axis=0)
        xmin_ind = np.argmin(np.array(zmin_pts),axis=0)[0]
        xmax_ind = np.argmax(np.array(zmin_pts),axis=0)[0]
        xmin = zmin_pts[xmin_ind]
        xmax = zmin_pts[xmax_ind]

        ymin_ind = np.argmin(np.array(zmin_pts),axis=0)[1]
        ymax_ind = np.argmax(np.array(zmin_pts),axis=0)[1]
        ymin = zmin_pts[ymin_ind]
        ymax = zmin_pts[ymax_ind]

        # make a sqaure using the pts above
        square = vtk.vtkPlaneSource()
        square.SetCenter(center)
        square.SetPoint1(ymin)
        square.SetPoint2(ymax)
        square.Update()
        #Surface(polydata=square.GetOutput()).write_surface_mesh(path+'square.vtk')

        # use this sqaure as the cutting plane
        pl = vtk.vtkPlane()
        pl.SetOrigin(center)
        pl.SetNormal(square.GetNormal())

        # Cut the surface using the cutting plane
        ex = vtk.vtkExtractPolyDataGeometry()
        ex.SetImplicitFunction(pl)
        ex.SetInputData(self.surface)
        ex.SetExtractInside(False)
        ex.Update()

        if op_filename is not None:
            Surface(polydata=ex.GetOutput()).write_surface_mesh(op_filename)

    def clip_surface_at_given_z(self, z_offset, z_limit=None, output_filename = None,check_correct_clipped_part=False):

        # use this sqaure as the cutting plane
        lower_bound = [self.surface.GetBounds()[0],self.surface.GetBounds()[2],self.surface.GetBounds()[4]]
        upper_bound = [self.surface.GetBounds()[1],self.surface.GetBounds()[3],self.surface.GetBounds()[5]]
        #center = [(lower_bound[0]+upper_bound[0])/2,(lower_bound[1]+upper_bound[1])/2,(lower_bound[2]+upper_bound[2])/2]
        if z_limit is not None:
            z = z_limit
        else:
            z = upper_bound[2]+z_offset

        pt1 = [lower_bound[0],lower_bound[1],z]
        pt2 = [upper_bound[0],upper_bound[1],z]
        pt3 = [lower_bound[0],upper_bound[1],z]

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,0)
        triangle.GetPointIds().SetId(1,1)
        triangle.GetPointIds().SetId(2,2)

        points = vtk.vtkPoints()
        points.InsertNextPoint(pt1)
        points.InsertNextPoint(pt2)
        points.InsertNextPoint(pt3)

        triangles = vtk.vtkCellArray()
        triangles.InsertNextCell(triangle )

        tri_min = vtk.vtkPolyData()
        tri_min.SetPoints(points)
        tri_min.SetPolys(triangles)

        #Surface(polydata=tri_min).write_surface_mesh(output_path+exp_name+'tri_min.vtk')

        norm_gen = vtk.vtkPolyDataNormals()
        norm_gen.SetInputData(tri_min)
        norm_gen.ComputeCellNormalsOn()
        norm_gen.Update()
        pdata_with_norms = norm_gen.GetOutput()
        arr = pdata_with_norms.GetCellData().GetNormals()
        norms_min = [arr.GetComponent(0,0),arr.GetComponent(0,1),arr.GetComponent(0,2)]

        plane = vtk.vtkPlane()
        plane.SetNormal(norms_min)
        plane.SetOrigin(pt1)
        #Surface(polydata=tri).write_surface_mesh(output_path+'triangle.vtk')

        clipper = vtk.vtkClipPolyData()
        clipper.SetClipFunction(plane)
        clipper.SetInputData(self.surface)
        clipper.SetGenerateClippedOutput(True)
        clipper.Update()

        clipped_surf = None
        #print(clipper.GetOutput().GetNumberOfCells() , clipper.GetClippedOutput().GetNumberOfCells())
        # select the larger of the two pieces as we are only trimming the bottom or top
        if check_correct_clipped_part:
            # need to check which is the correct clipped output
            if clipper.GetOutput().GetNumberOfCells() > clipper.GetClippedOutput().GetNumberOfCells():
                clipped_surf = clipper.GetOutput()
            else:
                clipped_surf = clipper.GetClippedOutput()
        else:
            clipped_surf = clipper.GetClippedOutput()
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(clipped_surf)
        cleaner.Update()

        connectivityfilt = vtk.vtkPolyDataConnectivityFilter()
        connectivityfilt.SetInputData(cleaner.GetOutput())
        connectivityfilt.SetExtractionModeToLargestRegion()
        connectivityfilt.Update()

        self.surface = connectivityfilt.GetOutput()
        if output_filename is not None:
            Surface(polydata=connectivityfilt.GetOutput()).write_surface_mesh(output_filename)

        self.surface = connectivityfilt.GetOutput()

        return connectivityfilt.GetOutput()

    def clip_surface_at_given_x(self, x_offset, x_limit=None, output_filename = None):

        # use this sqaure as the cutting plane
        lower_bound = [self.surface.GetBounds()[0],self.surface.GetBounds()[2],self.surface.GetBounds()[4]]
        upper_bound = [self.surface.GetBounds()[1],self.surface.GetBounds()[3],self.surface.GetBounds()[5]]
        #center = [(lower_bound[0]+upper_bound[0])/2,(lower_bound[1]+upper_bound[1])/2,(lower_bound[2]+upper_bound[2])/2]
        if x_limit is not None:
            x = x_limit
        else:
            x = upper_bound[0]+x_offset

        pt1 = [x,upper_bound[1],upper_bound[2]]
        pt2 = [x,lower_bound[1],lower_bound[2]]
        pt3 = [x,lower_bound[1],upper_bound[2]]

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,0)
        triangle.GetPointIds().SetId(1,1)
        triangle.GetPointIds().SetId(2,2)

        points = vtk.vtkPoints()
        points.InsertNextPoint(pt1)
        points.InsertNextPoint(pt2)
        points.InsertNextPoint(pt3)

        triangles = vtk.vtkCellArray()
        triangles.InsertNextCell(triangle )

        tri_min = vtk.vtkPolyData()
        tri_min.SetPoints(points)
        tri_min.SetPolys(triangles)

        #Surface(polydata=tri_min).write_surface_mesh(output_path+exp_name+'tri_min.vtk')

        norm_gen = vtk.vtkPolyDataNormals()
        norm_gen.SetInputData(tri_min)
        norm_gen.ComputeCellNormalsOn()
        norm_gen.Update()
        pdata_with_norms = norm_gen.GetOutput()
        arr = pdata_with_norms.GetCellData().GetNormals()
        norms_min = [arr.GetComponent(0,0),arr.GetComponent(0,1),arr.GetComponent(0,2)]

        plane = vtk.vtkPlane()
        plane.SetNormal(norms_min)
        plane.SetOrigin(pt1)
        #Surface(polydata=tri).write_surface_mesh(output_path+'triangle.vtk')

        clipper = vtk.vtkClipPolyData()
        clipper.SetClipFunction(plane)
        clipper.SetInputData(self.surface)
        clipper.SetGenerateClippedOutput(True)
        clipper.Update()

        clipped_surf = None
        #print(clipper.GetOutput().GetNumberOfCells() , clipper.GetClippedOutput().GetNumberOfCells())
        # select the larger of the two pieces as we are only trimming the bottom or top
        if clipper.GetOutput().GetNumberOfCells() > clipper.GetClippedOutput().GetNumberOfCells():
            clipped_surf = clipper.GetOutput()
        else:
            clipped_surf = clipper.GetClippedOutput()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(clipped_surf)
        cleaner.Update()

        connectivityfilt = vtk.vtkPolyDataConnectivityFilter()
        connectivityfilt.SetInputData(cleaner.GetOutput())
        connectivityfilt.SetExtractionModeToLargestRegion()
        connectivityfilt.Update()

        if output_filename is not None:
            Surface(polydata=connectivityfilt.GetOutput()).write_surface_mesh(output_filename)

        self.surface = connectivityfilt.GetOutput()

        return connectivityfilt.GetOutput()

    def clip_surface_at_given_y(self, y_offset, y_limit=None, output_filename = None,get_top=False):

        # use this sqaure as the cutting plane
        lower_bound = [self.surface.GetBounds()[0],self.surface.GetBounds()[2],self.surface.GetBounds()[4]]
        upper_bound = [self.surface.GetBounds()[1],self.surface.GetBounds()[3],self.surface.GetBounds()[5]]
        #center = [(lower_bound[0]+upper_bound[0])/2,(lower_bound[1]+upper_bound[1])/2,(lower_bound[2]+upper_bound[2])/2]
        if y_limit is not None:
            y = y_limit
        else:
            y = upper_bound[0]+y_offset

        pt1 = [upper_bound[0],y,upper_bound[2]]
        pt2 = [lower_bound[0],y,lower_bound[2]]
        pt3 = [lower_bound[0],y,upper_bound[2]]

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,0)
        triangle.GetPointIds().SetId(1,1)
        triangle.GetPointIds().SetId(2,2)

        points = vtk.vtkPoints()
        points.InsertNextPoint(pt1)
        points.InsertNextPoint(pt2)
        points.InsertNextPoint(pt3)

        triangles = vtk.vtkCellArray()
        triangles.InsertNextCell(triangle )

        tri_min = vtk.vtkPolyData()
        tri_min.SetPoints(points)
        tri_min.SetPolys(triangles)

        #Surface(polydata=tri_min).write_surface_mesh(output_path+exp_name+'tri_min.vtk')

        norm_gen = vtk.vtkPolyDataNormals()
        norm_gen.SetInputData(tri_min)
        norm_gen.ComputeCellNormalsOn()
        norm_gen.Update()
        pdata_with_norms = norm_gen.GetOutput()
        arr = pdata_with_norms.GetCellData().GetNormals()
        norms_min = [arr.GetComponent(0,0),arr.GetComponent(0,1),arr.GetComponent(0,2)]

        plane = vtk.vtkPlane()
        plane.SetNormal(norms_min)
        plane.SetOrigin(pt1)
        #Surface(polydata=tri).write_surface_mesh(output_path+'triangle.vtk')

        clipper = vtk.vtkClipPolyData()
        clipper.SetClipFunction(plane)
        clipper.SetInputData(self.surface)
        clipper.SetGenerateClippedOutput(True)
        clipper.Update()

        clipped_surf = None
        #print(clipper.GetOutput().GetNumberOfCells() , clipper.GetClippedOutput().GetNumberOfCells())
        # select the larger of the two pieces as we are only trimming the bottom or top
        if get_top:
            clipped_surf = clipper.GetClippedOutput()
        else:
            if clipper.GetOutput().GetNumberOfCells() > clipper.GetClippedOutput().GetNumberOfCells():
                clipped_surf = clipper.GetOutput()
            else:
                clipped_surf = clipper.GetClippedOutput()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(clipped_surf)
        cleaner.Update()

        connectivityfilt = vtk.vtkPolyDataConnectivityFilter()
        connectivityfilt.SetInputData(cleaner.GetOutput())
        connectivityfilt.SetExtractionModeToLargestRegion()
        connectivityfilt.Update()

        if output_filename is not None:
            Surface(polydata=connectivityfilt.GetOutput()).write_surface_mesh(output_filename)

        self.surface = connectivityfilt.GetOutput()

        return connectivityfilt.GetOutput()

    def get_distance(self,a,b):
        return np.linalg.norm( np.array(b)-np.array(a) )

    def get_euclidean_distance(self,a,b):
        return distance.euclidean(a,b)

    def get_vector_intersection_pt(self,start_vert,end_vert,extrapolation_len=10000):
        '''for the given vector find its intersection with the surface as well as depth from it'''
        intersec_pt = []
        dist = 0

        obbsurf = vtk.vtkOBBTree()
        obbsurf.SetDataSet(self.surface)
        obbsurf.BuildLocator()

        if extrapolation_len != 0:
            new_end_vert = self.create_pt_along_vector_at_given_distance(extrapolation_len,start_vert,end_vert)
            #new_start_vert = self.create_pt_along_vector_at_given_distance(-extrapolation_len,start_vert,end_vert)
        else:
            #new_start_vert = start_vert
            new_end_vert = end_vert

        intersection_points = vtk.vtkPoints()
        intersection_cellidlist = vtk.vtkIdList()

        code = obbsurf.IntersectWithLine(new_end_vert, start_vert, intersection_points, intersection_cellidlist)
        #print(new_end_vert,start_vert,code)
        if intersection_points.GetNumberOfPoints() > 0:
            intersec_pt = intersection_points.GetPoint(0)
            if len(intersec_pt) > 0:
                dist = self.get_distance(start_vert,intersec_pt)

        return intersec_pt,dist

    def get_vector_intersection_pt_from_obb_tree(self,obbsurf,start_vert,end_vert,extrapolation_len=10000):
        '''for the given vector find its intersection with the surface as well as depth from it'''
        intersec_pt = []
        dist = 0

        #obbsurf = vtk.vtkOBBTree()
        #obbsurf.SetDataSet(self.surface)
        #obbsurf.BuildLocator()

        if extrapolation_len != 0:
            new_end_vert = self.create_pt_along_vector_at_given_distance(extrapolation_len,start_vert,end_vert)
            #new_start_vert = self.create_pt_along_vector_at_given_distance(-extrapolation_len,start_vert,end_vert)
        else:
            #new_start_vert = start_vert
            new_end_vert = end_vert

        intersection_points = vtk.vtkPoints()
        intersection_cellidlist = vtk.vtkIdList()

        code = obbsurf.IntersectWithLine(new_end_vert, start_vert, intersection_points, intersection_cellidlist)
        #print(new_end_vert,start_vert,code)
        if intersection_points.GetNumberOfPoints() > 0:
            intersec_pt = intersection_points.GetPoint(0)
            if len(intersec_pt) > 0:
                dist = self.get_distance(start_vert,intersec_pt)

        return intersec_pt,dist

    def is_neuron_sticking_out(self,sg,):
        linear_ind = 0
        sticking_out_pts = []
        need_to_check_further = False
        # first check if all verts are inside
        soma_center = soma_center = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)
        for vertex in sg.neuron.all_neurites_subgraphdata.vertices:
            pia_proj_pt_neuron,soma_depth = self.get_vector_intersection_pt(vertex,soma_center,extrapolation_len=0)
            if len(pia_proj_pt_neuron) > 0:
                return True
        return False
        for i in range(len(sg.neuron.all_neurites_subgraphdata.edge_connectivity)):
            #vert_from = sg.neuron.all_neurites_subgraphdata.vertices[(sg.neuron.all_neurites_subgraphdata.edge_connectivity[i][0])]
            #vert_to = sg.neuron.all_neurites_subgraphdata.vertices[(sg.neuron.all_neurites_subgraphdata.edge_connectivity[i][1])]
            for j in range(sg.neuron.all_neurites_subgraphdata.num_pts_per_edge[i]-1):
                from_edge_pt = sg.neuron.all_neurites_subgraphdata.edge_pt_coords[linear_ind+j]
                to_edge_pt = sg.neuron.all_neurites_subgraphdata.edge_pt_coords[linear_ind+j+1]
                intersection_pt,dist = self.get_vector_intersection_pt(from_edge_pt[0:3],to_edge_pt[0:3],extrapolation_len=0)
                if len(intersection_pt) > 0:
                        return True
            linear_ind = linear_ind +j
#         for edge in sg.neuron.all_neurites_subgraphdata.edge_list:
#             for i in range(len(edge)-1):
#                 intersection_pt,dist = self.get_vector_intersection_pt(edge[i+1][0:3],edge[i][0:3],extrapolation_len=0)
            #print(intersection_pt)

        return False

    def create_axis_field(self,target_surface,op_sg_name = None,flip_normals=True, check_source_double_touch = False,return_axes=False, max_allowed_axis_length=5000):
        '''
        from the self.surface create rays till they instersect the given surface
        optionally write out the axis field spatial graph
        '''
        obbPia = vtk.vtkOBBTree()
        obbPia.SetDataSet(target_surface.surface)
        obbPia.BuildLocator()

        obbwm = vtk.vtkOBBTree()
        obbwm.SetDataSet(self.surface)
        obbwm.BuildLocator()

        source_cell_centers = vtk.vtkCellCenters()
        source_cell_centers.SetInputData(self.surface)
        source_cell_centers.Update()
        source_centroid_cell_centers = source_cell_centers.GetOutput()

        source_sphere_normals = vtk.vtkPolyDataNormals()
        source_sphere_normals.SetInputData(self.surface)
        source_sphere_normals.ComputePointNormalsOff()
        source_sphere_normals.ComputeCellNormalsOn()
        source_sphere_normals.SplittingOff()
        if flip_normals is True:
            source_sphere_normals.FlipNormalsOn()
        else:
            source_sphere_normals.FlipNormalsOff()
        source_sphere_normals.AutoOrientNormalsOn()
        source_sphere_normals.Update()
        source_source_normals = source_sphere_normals.GetOutput().GetCellData().GetNormals()

        dest_cell_centers = vtk.vtkCellCenters()
        dest_cell_centers.SetInputData(target_surface.surface)
        dest_cell_centers.Update()
        dest_centroid_cell_centers = dest_cell_centers.GetOutput()

        dest_sphere_normals = vtk.vtkPolyDataNormals()
        dest_sphere_normals.SetInputData(target_surface.surface)
        dest_sphere_normals.ComputePointNormalsOff()
        dest_sphere_normals.ComputeCellNormalsOn()
        dest_sphere_normals.SplittingOff()
        if flip_normals is True:
            dest_sphere_normals.FlipNormalsOff()
        else:
            dest_sphere_normals.FlipNormalsOn()
        #dest_sphere_normals.AutoOrientNormalsOn()
        dest_sphere_normals.Update()
        dest_source_normals = dest_sphere_normals.GetOutput().GetCellData().GetNormals()

        axis_field_sg_source = AmiraSpatialGraph()
        axis_field_sg_dest = AmiraSpatialGraph()
        from_vert = []
        to_vert = []
        prev_from_vert = []
        prev_to_vert = []
        for idx in range(source_centroid_cell_centers.GetNumberOfPoints()):

            source_point = source_centroid_cell_centers.GetPoint(idx)
            ##print(point)
            source_normal = source_source_normals.GetTuple(idx)

            # Calculate the 'target' of the ray based on 'RayCastLength'
            point_target = list(np.array(source_point) + 5000*np.array(source_normal))

            ##print('are we here')
            intersection_points = vtk.vtkPoints()
            intersection_cellidlist = vtk.vtkIdList()
            code = obbPia.IntersectWithLine(point_target, source_point, intersection_points, intersection_cellidlist)

            if intersection_points.GetNumberOfPoints() >= 1:
                from_vert = source_point
                if intersection_points.GetNumberOfPoints() > 1:
                    # get the last point
                    to_vert = intersection_points.GetPoint(intersection_points.GetNumberOfPoints()-1)
                else:
                    to_vert = intersection_points.GetPoint(0)

                edge_len,edge_len_mean,edge_len_std = self.get_vec_lengths([[from_vert,to_vert]])
                if edge_len[0] < max_allowed_axis_length:
                    if check_source_double_touch:
                        # check if it intersects the source surface
                        intersection_points_wm = vtk.vtkPoints()
                        intersection_cellidlist_wm = vtk.vtkIdList()
                        wm_pt = Vectors().create_pt_along_vector_at_given_distance(100,from_vert,to_vert)
                        #new_edge = [wm_pt,edge[1]]
                        code = obbwm.IntersectWithLine(point_target, wm_pt, intersection_points_wm, intersection_cellidlist_wm)
                        if intersection_points_wm.GetNumberOfPoints() <= 0:
                            axis_field_sg_source.graph_data.add_edge(from_vert,to_vert)
                    else:
                        axis_field_sg_source.graph_data.add_edge(from_vert,to_vert)

        if op_sg_name is not None:
            axis_field_sg_source.write_spatial_graph(op_sg_name)
        if return_axes:
            return axis_field_sg_source.graph_data.edge_list

    def fill_holes(self,surf,hole_size=10000):
        filler = vtk.vtkFillHolesFilter()
        filler.SetHoleSize(hole_size)
        filler.SetInputData(surf)
        filler.Update()
        return filler.GetOutput()

    def create_delunay_surface_2d(self,return_hull=False,output_filename = None,alpha=0,offset=1):
        delu2d = vtk.vtkDelaunay2D()
        delu2d.SetInputData(self.surface)
        delu2d.SetAlpha(alpha)
        delu2d.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE )
        delu2d.SetOffset(offset)
        #delu2d.SetProjectionPlaneMode(vtk.VTK_DELAUNAY_XY_PLANE )
        delu2d.Update()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(delu2d.GetOutput())
        cleaner.Update()

        surf = vtk.vtkDataSetSurfaceFilter()
        surf.SetInputData(cleaner.GetOutput())
        surf.Update()

        filled_surf = surf.GetOutput()#self.fill_holes(surf.GetOutput())
        if output_filename is not None:
            Surface(polydata=filled_surf).write_surface_mesh(output_filename)

        if return_hull:
            return filled_surf

    def create_delunay_surface_3d(self,surf=None,return_hull=False,output_filename = None,make_cube=False,bounding_offset=0,bounds=None):
        if make_cube:
            hull = vtk.vtkDelaunay3D()
            if surf is None:
                hull.SetInputData(self.surface)
            else:
                hull.SetInputData(surf)
            hull.Update()

            if bounds is None:
                bounds = hull.GetOutput().GetBounds()

            x_low = bounds[0] - bounding_offset
            x_up = bounds[1] + bounding_offset
            y_low = bounds[2] - bounding_offset
            y_up = bounds[3] + bounding_offset

            cube = vtk.vtkCubeSource()
            cube.SetBounds([x_low,x_up,y_low,y_up,bounds[4],bounds[5]])
            cube.SetCenter(hull.GetOutput().GetCenter())
            cube.Update()

            surf = vtk.vtkDataSetSurfaceFilter()
            surf.SetInputData(cube.GetOutput())
            surf.Update()
        else:

            hull = vtk.vtkDelaunay3D()
            hull.SetInputData(self.surface)
            hull.Update()

            surf = vtk.vtkDataSetSurfaceFilter()
            surf.SetInputData(hull.GetOutput())
            surf.Update()

        self.surface = surf.GetOutput()

        if output_filename is not None:
            Surface(polydata=surf.GetOutput()).write_surface_mesh(output_filename)
        #Surface(polydata=surf.GetOutput()).write_surface_mesh(output_path+'new_hull1.vtk')

        if return_hull:
            return surf.GetOutput()

    def get_center_of_mass(self,):
        com = vtk.vtkCenterOfMass()
        com.SetInputData(self.surface)
        com.Update()

        return com.GetCenter()

    def get_surface_properties(self,surf=None,prop='area'):
        props = vtk.vtkMassProperties()
        if surf is not None:
            props.SetInputData(surf)
        else:
            props.SetInputData(self.surface)
        props.Update()

        if prop == 'area':
            return props.GetSurfaceArea()
        elif prop == 'volume':
            return props.GetVolume()

    def get_axis_field_within_instersecting_surface_with_trimming(self,pia_surface,wm_surface,edge_list,vec_extension_offset=500,):
        obb_tree_pia = vtk.vtkOBBTree()
        obb_tree_pia.SetDataSet(pia_surface)
        obb_tree_pia.BuildLocator()

        obb_tree_wm = vtk.vtkOBBTree()
        obb_tree_wm.SetDataSet(wm_surface)
        obb_tree_wm.BuildLocator()

        # find the vectors within this
        vm1_axis_field_sg = AmiraSpatialGraph()
        missing_axis_field_sg = AmiraSpatialGraph()
        wm_pts = []
        pia_pts = []
        #edge_list = []
        for edge in edge_list:
            #print(edge)
            #tmp_start = self.create_pt_along_vector_at_given_distance(1000, edge[0], edge[1])
            #tmp_end = self.create_pt_along_vector_at_given_distance(1000, edge[1], edge[0])
            edge_len = Vectors().get_vec_length(edge)
            start_vert = Vectors().create_pt_along_vector_at_given_distance(edge_len+vec_extension_offset*edge_len, edge[1], edge[0])
            end_vert = Vectors().create_pt_along_vector_at_given_distance(10000, edge[0], edge[1])

            #vm1_axis_field_sg.graph_data.add_edge(start_vert, end_vert)
            #return vm1_axis_field_sg,wm_pts,pia_pts
            #print(start_vert,end_vert)

            pia_intersection_pt = vtk.vtkPoints()
            wm_intersection_pt = vtk.vtkPoints()

            code1 = obb_tree_pia.IntersectWithLine(end_vert, start_vert, pia_intersection_pt, None)
            code2 = obb_tree_wm.IntersectWithLine(start_vert, end_vert, wm_intersection_pt, None)
            if pia_intersection_pt.GetNumberOfPoints() > 0 and wm_intersection_pt.GetNumberOfPoints() > 0:
                #print(edge)
                pia_pt =  pia_intersection_pt.GetPoint(0)
                wm_pt =  wm_intersection_pt.GetPoint(0)
                #if pia_intersection_pt.GetNumberOfPoints() > 1:
                #    if Vectors().get_vec_length([wm_intersection_pt.GetPoint(0),pia_intersection_pt.GetPoint(0)]) < \
                #    Vectors().get_vec_length([wm_intersection_pt.GetPoint(0),pia_intersection_pt.GetPoint(1)]):
                #        pia_pt = pia_intersection_pt.GetPoint(0)
                #    else:
                #        pia_pt = pia_intersection_pt.GetPoint(1)
                #if wm_intersection_pt.GetNumberOfPoints() > 1:
                #    if Vectors().get_vec_length([pia_intersection_pt.GetPoint(0),wm_intersection_pt.GetPoint(0)]) < \
                #    Vectors().get_vec_length([pia_intersection_pt.GetPoint(0),wm_intersection_pt.GetPoint(1)]):
                #        wm_pt = wm_intersection_pt.GetPoint(0)
                #    else:
                #        wm_pt = wm_intersection_pt.GetPoint(1)
                pia_pts.append(pia_pt)
                wm_pts.append(wm_pt)
                vm1_axis_field_sg.graph_data.add_edge(wm_pt,pia_pt)
            else:
                # for some reason could not clip the axis by the surfaces
                missing_axis_field_sg.graph_data.add_edge(start_vert,end_vert)

        return vm1_axis_field_sg,wm_pts,pia_pts,missing_axis_field_sg

    def get_axis_field_within_instersecting_surface(self,surface,edge_list,trim_edge=False,vec_extension_len=None,):
        obb_tree = vtk.vtkOBBTree()
        obb_tree.SetDataSet(surface)
        obb_tree.BuildLocator()

        # find the vectors within this
        points = vtk.vtkPoints()

        vm1_axis_field_sg = AmiraSpatialGraph()

        wm_pts = []
        pia_pts = []
        #edge_list = []
        for edge in edge_list:
            #print(edge)
            if vec_extension_len is not None:
                start_vert = edge[0]
                end_vert = self.create_pt_along_vector_at_given_distance(vec_extension_len, edge[0], edge[1])
            else:
                start_vert = edge[0]
                end_vert = edge[1]

            code = obb_tree.IntersectWithLine(start_vert, end_vert, points, None)
            if points.GetNumberOfPoints() > 0:
                #print(edge)
                wm_pts.append(edge[0])
                pia_pts.append(edge[1])
                if trim_edge:
                    vm1_axis_field_sg.graph_data.add_edge(edge[0],points.GetPoint(0))
                else:
                    vm1_axis_field_sg.graph_data.add_edge(edge[0],edge[1])
        return vm1_axis_field_sg,wm_pts,pia_pts

    def get_cutting_plane(self,pt1,pt2,pt3):
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,0)
        triangle.GetPointIds().SetId(1,1)
        triangle.GetPointIds().SetId(2,2)

        points = vtk.vtkPoints()
        points.InsertNextPoint(pt1)
        points.InsertNextPoint(pt2)
        points.InsertNextPoint(pt3)

        triangles = vtk.vtkCellArray()
        triangles.InsertNextCell(triangle )

        tri_min = vtk.vtkPolyData()
        tri_min.SetPoints(points)
        tri_min.SetPolys(triangles)

        #Surface(polydata=tri_min).write_surface_mesh(output_path+exp_name+'tri_min.vtk')

        norm_gen = vtk.vtkPolyDataNormals()
        norm_gen.SetInputData(tri_min)
        norm_gen.ComputeCellNormalsOn()
        norm_gen.Update()
        pdata_with_norms = norm_gen.GetOutput()
        arr = pdata_with_norms.GetCellData().GetNormals()
        norms_min = [arr.GetComponent(0,0),arr.GetComponent(0,1),arr.GetComponent(0,2)]

        plane = vtk.vtkPlane()
        plane.SetNormal(norms_min)
        plane.SetOrigin(pt1)

        return plane,tri_min

    def clip_surfaces_using_given_planes(self,pia,wm,plane_right,plane_left,plane_top,plane_bottom,central_axis=[],landmarks=[],order_of_cutting=[True,False,True,False],invert_cell_direction=False):

        #pia,wm = self.append(wm,surf1=pia,return_surf=True)
        pia1,wm1 = self.clip_surfaces(pia,wm,plane_right,central_axis,get_left=order_of_cutting[0])
        pia2,wm2 = self.clip_surfaces(pia1,wm1,plane_left,central_axis,get_left=order_of_cutting[1])
        pia3,wm3 = self.clip_surfaces(pia2,wm2,plane_top,central_axis,get_left=order_of_cutting[2])
        pia4,wm4 = self.clip_surfaces(pia3,wm3,plane_bottom,central_axis,get_left=order_of_cutting[3])

        appender = vtk.vtkAppendPolyData()
        appender.AddInputData(pia4)
        appender.AddInputData(wm4)
        appender.Update()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(appender.GetOutput())
        cleaner.Update()

        delu3d = Surface(polydata=cleaner.GetOutput())
        hull = delu3d.create_delunay_surface_3d(return_hull=True)

        if wm4.GetNumberOfCells() > 0:
            final = self.crop_hull(hull,wm4,invert_cell_direction)
            #final = self.crop_hull(final1,pia4)
            cleaner = vtk.vtkCleanPolyData()
            cleaner.SetInputData(final)
            cleaner.Update()
            final = cleaner.GetOutput()
        else:
            final = hull
        return final,hull,pia4,wm4

    def clip_surfaces(self,pia,wm,plane,central_axis=[],landmarks=[],get_left=True):
        #print('clipping_surface')
        pia1,pia2 = self.clip_surface(pia,plane)
        wm1,wm2 = self.clip_surface(wm,plane)

        if len(central_axis) > 0:
            print(central_axis)
            # find the right half by checking which one contains the central axis
            print(pia1.GetNumberOfCells(),wm1.GetNumberOfCells())
            print(pia2.GetNumberOfCells(),wm2.GetNumberOfCells())
            if pia1.GetNumberOfCells()>0 :
                pia_pt1,dist= Surface(polydata=pia1).get_vector_intersection_pt(central_axis[0],central_axis[1],)
                #wm_pt1,dist = Surface(polydata=wm1).get_vector_intersection_pt(central_axis[1],central_axis[0],)
            else:
                pia_pt1,dist = Surface(polydata=pia2).get_vector_intersection_pt(central_axis[0],central_axis[1],)
                #wm_pt1,dist = Surface(polydata=wm2).get_vector_intersection_pt(central_axis[1],central_axis[0],)

            if len(pia_pt1) > 0:
                print(pia_pt1)
                return pia1,wm1
            else:
                print(pia_pt1)
                return pia2,wm2
        else:
            if get_left:
                return pia1,wm1
            else:
                return pia2,wm2
        # figure out which part is the correct one
        # Do so by checking the ray intersection for the pts enclising cube
        #if len(landmarks) > 1:
        #    # this is the case when we are trying to create vm1 from given pts
        #    delun3d = Surface(pts=landmarks)
        #    cube = delun3d.create_delunay_surface_3d(return_hull=True,make_cube=True)
        #
        #    edges1 = Surface(polydata=cube).create_axis_field(Surface(polydata=pia1),return_axes=True)
        #    edges2 = Surface(polydata=cube).create_axis_field(Surface(polydata=pia2),return_axes=True)


            #print(len(edges1),len(edges2))

        #    if len(edges1) > len(edges2):
        #        return pia1,wm1
        #    else:
        #        return pia2,wm2
        #else:
            # determine the right clippd part by checking that the center pt remains in the part
            # this is used for getting the avg vM1, given the pcas
        #if get_left:
        #    return pia1,wm1
        #else:
        #    return pia2,wm2
            #print('trying to validate using method 2')
            #pia_wm1 = Surface(polydata=pia1)
            #pia_wm1.append(wm1)
            #print(pia_wm1.surface.GetNumberOfCells())
            #pia_wm1_hull = pia_wm1.create_delunay_surface_3d(return_hull=True)
            #print(pia_wm1_hull.GetNumberOfCells())
            #pts_within_hull = self.get_landmarks_within_given_surface(pia_wm1.surface,pts_to_be_used=landmarks)
            #print(len(pts_within_hull))
            #if len(pts_within_hull) > 0:
            #    return pia1,wm1
            #else:
            #    return pia2,wm2

    def get_landmarks_within_given_surface(self,surface,pts_to_be_used = None,output_file = None):

        if pts_to_be_used is not None:
            pts_being_used = pts_to_be_used
        else:
            pts_being_used = self.pts

        sel = vtk.vtkSelectEnclosedPoints()
        sel.SetInputData(Surface(pts=pts_being_used).convert_points_to_polydata(pts_being_used))
        sel.SetSurfaceData(surface)
        sel.SetInsideOut(True)
        sel.Update()

        valid_pts = []
        invalid_pts = []
        arr =  sel.GetOutput().GetPointData().GetArray("SelectedPoints")
        for i in range((arr.GetNumberOfTuples())):
            if (arr.GetComponent(i,0)) == 1:
                valid_pts.append(pts_being_used[i])
            else:
                invalid_pts.append(pts_being_used[i])

        if output_file is not None:
            self.write_landmarks(output_file,pts_to_be_written=valid_pts)

        return valid_pts,invalid_pts

    def clean_polydata(self,surf):
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(surf)
        cleaner.Update()

        return cleaner.GetOutput()

    def get_largest_connected_polydata(self,polydata):
        connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
        connectivityFilter.SetInputData(polydata)
        connectivityFilter.SetExtractionModeToLargestRegion()
        connectivityFilter.Update()
        return connectivityFilter.GetOutput()

    def clip_surface(self,surface,plane):
        clipper = vtk.vtkClipPolyData()
        clipper.SetClipFunction(plane)
        clipper.SetInputData(surface)
        clipper.SetGenerateClippedOutput(True)

        clipper.Update()
        #Surface(polydata=clipper.GetOutput()).write_surface_mesh('clipped_pia1.vtk')
        #Surface(polydata=clipper.GetClippedOutput()).write_surface_mesh('clipped_other_pia1.vtk')
        return self.get_largest_connected_polydata(clipper.GetOutput()),self.get_largest_connected_polydata(clipper.GetClippedOutput())

    def crop_hull(self,surface1,surface2,invert_cell_direction=False):
        '''This filter trims the hull (surf 1) with the surf 2'''

        sd = vtk.vtkDistancePolyDataFilter()
        sd.SetInputData(0,surface1)
        sd.SetInputData(1,surface2)
        sd.Update()

        ids = vtk.vtkIdTypeArray()
        ids.SetNumberOfComponents(1)

        if invert_cell_direction:
            for i in range(surface1.GetNumberOfPoints()):
                if sd.GetOutput().GetPointData().GetScalars().GetComponent(i,0) >0:
                    #print(i)
                    ids.InsertNextValue(i)
        else:
            for i in range(surface1.GetNumberOfPoints()):
                if sd.GetOutput().GetPointData().GetScalars().GetComponent(i,0) <0:
                    #print(i)
                    ids.InsertNextValue(i)

        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT);
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES);
        selectionNode.SetSelectionList(ids);
        selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1);

        selection = vtk.vtkSelection()
        selection.AddNode(selectionNode)

        extractSelection  = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0,surface1);
        extractSelection.SetInputData(1, selection);
        extractSelection.Update();

        surf = vtk.vtkDataSetSurfaceFilter()
        surf.SetInputData(extractSelection.GetOutput())
        surf.Update()

        appender = vtk.vtkAppendPolyData()
        appender.AddInputData(surf.GetOutput())
        appender.AddInputData(surface2)
        appender.Update()

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(appender.GetOutput())
        cleaner.Update()

        return cleaner.GetOutput()
        #return surf.GetOutput()

    def clip_surface_between_left_right_planes(self,plane1,plane2,):
        # the input planes must be left and then right,
        # since the cutter first output is always the left and then right
        clipped_1_left,clipped_1_right = self.clip_surface(self.surface,plane1)
        clipped_2_left,clipped_2_right = self.clip_surface(clipped_1_right,plane2)

    #     Surface(polydata=clipped_2_left).write_surface_mesh(output_path+exp_name+'clipped_2_left.vtk')
    #     Surface(polydata=clipped_2_right).write_surface_mesh(output_path+exp_name+'clipped_2_right.vtk')


        return clipped_2_left

    def clip_surface_between_top_bottom_planes(self,plane1,plane2):
        # the input planes must be left and then right,
        # since the cutter first output is always the left and then right
        clipped_1_bottom,clipped_1_top = self.clip_surface(self.surface,plane1)
        clipped_2_bottom,clipped_2_top = self.clip_surface(clipped_1_top,plane2)

    #     Surface(polydata=clipped_1_top).write_surface_mesh(output_path+exp_name+'clipped_1_top.vtk')
    #     Surface(polydata=clipped_1_bottom).write_surface_mesh(output_path+exp_name+'clipped_1_bottom.vtk')
    #     Surface(polydata=clipped_2_top).write_surface_mesh(output_path+exp_name+'clipped_2_top.vtk')
    #     Surface(polydata=clipped_2_bottom).write_surface_mesh(output_path+exp_name+'clipped_2_bottom.vtk')

        return clipped_2_bottom

    def create_surface_from_unorganized_points(self,pts,sample_resolution=200,clip_offset=0,limit=None,output_filename=None,clipping_plane='z',hole_size=1000):

        surf = Surface(pts=pts)

        recon = vtk.vtkSurfaceReconstructionFilter()
        recon.SetSampleSpacing(sample_resolution)
        #recon.SetNeighborhoodSize(sample_resolution)
        recon.SetInputData(surf.surface)
        recon.Update()

    #     cont = vtk.vtkContourFilter()
    #     cont.SetInputData(recon.GetOutput())
    #     cont.SetValue(0,0.0)
    #     cont.Update()

        mc = vtk.vtkMarchingContourFilter()
        mc.SetInputData(recon.GetOutput())
        #mc.SetNumberOfContours(0)
        mc.SetValue(0,0)
        mc.Update()

        rev = vtk.vtkReverseSense()
        rev.SetInputData(mc.GetOutput())
        rev.ReverseCellsOn()
        rev.ReverseNormalsOn()
        rev.Update()

        if clipping_plane == 'z':
            clipped_surf = Surface(polydata=rev.GetOutput()).clip_surface_at_given_z(z_offset=clip_offset,z_limit=limit)
        else:
            clipped_surf = Surface(polydata=rev.GetOutput()).clip_surface_at_given_x(x_offset=clip_offset,x_limit=limit)
        #divided_surf = self.divide_surface(1,clipped_surf)

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputData(clipped_surf)
        cleaner.Update()

        hole_filled = self.fill_holes(cleaner.GetOutput(),hole_size=hole_size)

        smoothed_surf = self.smooth_surface(hole_filled)

        if output_filename:
            Surface(polydata=smoothed_surf).write_surface_mesh(output_filename)

        return smoothed_surf

    def smooth_surface(self,surface,num_iter=20,pass_band=0.01):
        smoother = vtk.vtkWindowedSincPolyDataFilter()
        smoother.BoundarySmoothingOff()
        smoother.FeatureEdgeSmoothingOff()
        smoother.NormalizeCoordinatesOn()
        smoother.SetNumberOfIterations(20)
        smoother.SetPassBand(0.01)
        smoother.SetInputData(surface)
        smoother.Update()

        return smoother.GetOutput()

    def find_max_pca(self,pca,axis=0):
        correct_pca = []
        max_component = 0
        correct_ind = 0
        for i in range(len(pca)):
            if abs(pca[i][axis]) > max_component:
                max_component = abs(pca[i][axis])
                correct_pca = pca[i]
                correct_ind = i
        #print(correct_pca)
        if correct_pca[axis] < 0:
            correct_pca = -correct_pca
        other_pcas = []
        for i in range(len(pca)):
            if i != correct_ind:
                other_pcas.append(pca[i])
        #print(other_pcas)
        return correct_pca,other_pcas

    def sort_pca_as_desired(self,pca,pts,output_path_vm1_axis,exp_name):
        # Desired direction
        # pca0 is one closest to x axis
        delunay_3d = Surface(pts=pts)
        rabies_hull_extended = delunay_3d.create_delunay_surface_3d(return_hull=True)
        center = np.array(pts).mean(axis=0)
        pcas = [pca[0,:],pca[1,:],pca[2,:]]
        pca0,other_pcas = self.find_max_pca(pcas,axis=0)
        pca1,other_pcas = self.find_max_pca(other_pcas,axis=2)
        pca2,other_pcas = self.find_max_pca(other_pcas,axis=1)

        corrected_pcas_list = [pca0,pca1,pca2]

        if exp_name == 'MG50_lhs':
            corrected_pcas_list = [corrected_pcas_list[2],corrected_pcas_list[0],corrected_pcas_list[1]]

        #print(corrected_pcas_list)
        pcas=[]
        for i in range(len(corrected_pcas_list)):
            pca_component = corrected_pcas_list[i]
            end_pt = Vectors().create_pt_along_vector_at_given_distance(10000,[0,0,0],pca_component)

            vectors_sg_pca,wm_pts,pia_pts = delunay_3d.get_axis_field_within_instersecting_surface\
                                                (rabies_hull_extended,[[list(center),list(end_pt+center)]],\
                                                 trim_edge=True)
            vectors_sg_pca.write_spatial_graph(output_path_vm1_axis+exp_name+'_pca_{}.am'.format(i))
            pcas.append(vectors_sg_pca.graph_data.edge_list[0])

        return pcas

    def get_pca_from_pts(self,pts,output_path_vm1_axis,exp_name):
        tr_mat,rot_mat,translation,pca_components,tr_pts = Vectors().get_pca_transformation_matrix(pts,translate_to_center=True)
        pcas = self.sort_pca_as_desired(pca_components,pts,output_path_vm1_axis,exp_name)
        #print(pcas)
        return pcas

    def get_vm1(self,pia_surf,wm_surf,pts=[],given_pcas=[],output_path_vm1_axis=None,exp_name=None):
        pcas = []
        if len(pts) > 0:
            pcas = self.get_pca_from_pts(pts,output_path_vm1_axis,exp_name)
            #print(pcas)
        else:
            pcas = given_pcas
        x_max = Vectors().get_vec_length(pcas[0])
        y_max = Vectors().get_vec_length(pcas[1])

        cone_center = Vectors().create_pt_along_vector_at_given_distance(x_max,pcas[2][0],pcas[2][1])

        x = Vectors().create_pt_along_vector_at_given_distance(x_max,pcas[0][0],pcas[0][1])
        y = Vectors().create_pt_along_vector_at_given_distance(y_max,pcas[1][0],pcas[1][1])
        a,x_delta = Vectors().get_center_and_unit_vec(pcas[0][0],pcas[0][1])
        b,y_delta = Vectors().get_center_and_unit_vec(pcas[1][0],pcas[1][1])
        neg_x = Vectors().create_pt_along_vector_at_given_distance(-x_max,pcas[0][0],pcas[0][1])
        neg_y = Vectors().create_pt_along_vector_at_given_distance(-y_max,pcas[1][0],pcas[1][1])
        neg_x_delta = Vectors().create_pt_along_vector_at_given_distance(-100,pcas[0][0],pcas[0][1])
        neg_y_delta = Vectors().create_pt_along_vector_at_given_distance(-100,pcas[1][0],pcas[1][1])

        plane_right,tri_right = Surface().get_cutting_plane(cone_center,x,x+100*y_delta)
        #Landmarks(pts=[x+100*y_delta]).write_landmarks(output_path_vM1+'b')
        #Landmarks(pts=[x]).write_landmarks(output_path_vM1+'x')
        plane_left,tri_left = Surface().get_cutting_plane(cone_center,neg_x,neg_x+100*y_delta)
        plane_top,tri_top = Surface().get_cutting_plane(cone_center,y,y+100*x_delta)
        plane_bottom,tri_bottom = Surface().get_cutting_plane(cone_center,neg_y,neg_y+100*x_delta)

        #Surface(polydata=tri_right).write_surface_mesh(output_path_avg+'tri_right.vtk')
        #Surface(polydata=tri_left).write_surface_mesh(output_path_avg+'tri_left.vtk')
        #Surface(polydata=tri_top).write_surface_mesh(output_path_avg+'tri_top.vtk')
        #Surface(polydata=tri_bottom).write_surface_mesh(output_path_avg+'tri_bottom.vtk')

        if len(pts) == 0:
            # push the pts to be pca center if pcas are already given
            # this help figure out which clipped parts are desired ones
            pts.append(pcas[0][0])
        vm1_final,vM1_hull,vm1_pia,vm1_wm = Surface().clip_surfaces_using_given_planes(pia_surf,wm_surf,plane_right,plane_left,\
                                                                               plane_top,plane_bottom,pts)

        return pcas,vm1_final,vM1_hull,vm1_pia,vm1_wm

    def get_intersection_plane(self,z_plane_coord):
        cutPlane = vtk.vtkPlane()
        center = self.surface.GetCenter()
        cutPlane.SetOrigin([center[0],center[1],z_plane_coord])
        cutPlane.SetNormal(0, 0, 1)

        cutter = vtk.vtkCutter()
        cutter.SetInputData(self.surface)
        cutter.SetCutFunction(cutPlane)
        cutter.Update()

        stripper = vtk.vtkStripper()
        stripper.SetInputData(cutter.GetOutput())
        stripper.Update()


        return stripper.GetOutput()

    def get_rays(self,center,theta_res=33,phi_res=32,ray_length=10000,ray_core_radius=1,tesselation=False,return_sphere=False):
        pia_sphere = vtk.vtkSphereSource()
        pia_sphere.SetCenter(center)
        pia_sphere.SetRadius(ray_core_radius)
        pia_sphere.SetThetaResolution(theta_res)
        pia_sphere.SetPhiResolution(phi_res)
        pia_sphere.SetStartPhi(0)
        pia_sphere.SetEndPhi(360)
        pia_sphere.SetStartTheta(0)
        pia_sphere.SetEndTheta(360)
        pia_sphere.SetLatLongTessellation(tesselation)
        pia_sphere.Update()

        pia_cell_centers = vtk.vtkCellCenters()
        pia_cell_centers.SetInputData(pia_sphere.GetOutput())
        pia_cell_centers.Update()
        pia_centroid_cell_centers = pia_cell_centers.GetOutput()
        ##print('pia_centroid_cell_centers {}'.format(pia_centroid_cell_centers))
        pia_sphere_normals = vtk.vtkPolyDataNormals()
        pia_sphere_normals.SetInputConnection(pia_sphere.GetOutputPort())
        pia_sphere_normals.ComputePointNormalsOff()
        pia_sphere_normals.ComputeCellNormalsOn()
        pia_sphere_normals.SplittingOff()
        pia_sphere_normals.FlipNormalsOff()
        pia_sphere_normals.AutoOrientNormalsOn()
        pia_sphere_normals.Update()
        pia_source_normals = pia_sphere_normals.GetOutput().GetCellData().GetNormals()

        rays_sg = AmiraSpatialGraph(generic_graph=True)
        phi_step = int((360-0)/theta_res)
        theta_step = int((360-0)/phi_res)
        idx = -1
        # for phi in range(0,360-phi_step,phi_step):
        # for theta in range(0,360-theta_step,theta_step):
        for idx in range(pia_centroid_cell_centers.GetNumberOfPoints()-1):
            idx = idx + 1
            point = pia_centroid_cell_centers.GetPoint(idx)
            ##print(point)
            normal = pia_source_normals.GetTuple(idx)

            # Calculate the 'target' of the ray based on 'RayCastLength'
            pointRayTarget = list(np.array(point) + ray_length*np.array(normal))

            # check if needs to be reveresed by looking for intersection wih the core sphere
            if self.check_ray_intersection([point,pointRayTarget],pia_sphere.GetOutput()):
                pointRayTarget = list(np.array(point) - ray_length*np.array(normal))

            rays_sg.graph_data.add_edge(point,pointRayTarget)

        if return_sphere:
            return rays_sg,pia_sphere.GetOutput()
        else:
            return rays_sg

    def get_aligned_surfaces_list(self,path,surf_res,surf_type='WB',closed_surf=False):
        aligned_pia_wm_surfs = Surface()
        aligned_pia_wm_surfs_lhs = Surface()
        aligned_pia_wm_surfs_rhs = Surface()
        aligned_pias_list_lhs = []
        aligned_wms_list_lhs = []
        aligned_pias_list_rhs = []
        aligned_wms_list_rhs = []

        for exp in ['MG48','MG49','MG50',]:
            for hem in ['lhs','rhs']:
                # read all aligned surfaces
                if surf_type == 'WB':
                    if closed_surf and exp == 'MG48':
                        pia = Surface(path + '/Surfaces/{}_{}_pia_voxel_size_{}.vtk'.format(exp,hem,surf_res))
                        wm = Surface(path + '/Surfaces/{}_{}_WM_voxel_size_{}.vtk'.format(exp,hem,surf_res))
                    else:
                        pia = Surface(path + '/Surfaces/{}_{}_pia_open_bottom_voxel_size_{}.vtk'.format(exp,hem,surf_res))
                        wm = Surface(path + '/Surfaces/{}_{}_WM_open_bottom_voxel_size_{}.vtk'.format(exp,hem,surf_res))
                elif surf_type == 'vS1':
                    pia = Surface(path + '/vS1_Ref_Surfaces/{}_{}_vs1_pia_surf_res_{}_voxel_size_500_cutoff_0.05.vtk'.format(exp,hem,surf_res))
                    wm = Surface(path + '/vS1_Ref_Surfaces/{}_{}_vs1_wm_surf_res_{}_voxel_size_500_cutoff_0.05.vtk'.format(exp,hem,surf_res))
                else:
                    pia = Surface(path + '/vM1_Ref_Surfaces/{}_{}_vm1_pia_surf_res_{}_voxel_size_500_cutoff_0.05.vtk'.format(exp,hem,surf_res))
                    wm = Surface(path + '/vM1_Ref_Surfaces/{}_{}_vm1_wm_surf_res_{}_voxel_size_500_cutoff_0.05.vtk'.format(exp,hem,surf_res))

                aligned_pia_wm_surfs.append(pia.surface)
                aligned_pia_wm_surfs.append(wm.surface)

                if hem== 'lhs':
                    aligned_pias_list_lhs.append(pia.surface)
                    aligned_wms_list_lhs.append(wm.surface)
                    aligned_pia_wm_surfs_lhs.append(pia.surface)
                    aligned_pia_wm_surfs_lhs.append(wm.surface)
                else:
                    aligned_pias_list_rhs.append(pia.surface)
                    aligned_wms_list_rhs.append(wm.surface)
                    aligned_pia_wm_surfs_rhs.append(pia.surface)
                    aligned_pia_wm_surfs_rhs.append(wm.surface)

        centroid = aligned_pia_wm_surfs.get_center_of_mass()
        centroid_lhs = aligned_pia_wm_surfs_lhs.get_center_of_mass()
        centroid_rhs = aligned_pia_wm_surfs_rhs.get_center_of_mass()

        return  centroid,centroid_lhs,centroid_rhs,aligned_pias_list_lhs,aligned_wms_list_lhs,aligned_pias_list_rhs,aligned_wms_list_rhs

    def get_obb_tree_list(self,surf_list):
        obb_tree_list = []
        for surf in surf_list:
            obb = vtk.vtkOBBTree()
            obb.SetDataSet(surf)
            obb.BuildLocator()
            obb_tree_list.append(obb)
        return obb_tree_list

    def get_avg_pts(self,rays_sg,obb_tree_list_pia,obb_tree_list_wm,create_error_ellisoids=False,\
                    alignment_type='',hem='',aligned_surf_res='',surf_res='',surf_type='WB'):
        avg_pts_pia = []
        avg_pts_wm = []
        wm_pia_rays_sg = AmiraSpatialGraph(generic_graph=True)
        error_ellipsoids = Surface()
        pia_stds = []
        wm_stds = []
        dists_pia_3d = []
        dists_wm_3d = []
        df_cols = ['Alignment_Type','Surface_Type','Hemisphere','Aligned_Using_Surf_Resolution','Present_Surf_Resolution',\
                                   'Pia_Mean_Std_X','Pia_Mean_Std_Y','Pia_Mean_Std_Z',\
                                   'WM_Mean_Std_X','WM_Mean_Std_Y','WM_Mean_Std_Z',\
                                   'Pia_Mean_Dist_3D','Pia_STD_Dist_3D','Pia_CoV_Dist_3D','Pia_SEM_Dist_3D',\
                                    'WM_Mean_Dist_3D','WM_STD_Dist_3D','WM_CoV_Dist_3D','WM_SEM_Dist_3D']
        df = pd.DataFrame(columns=df_cols)
        for edge in rays_sg.graph_data.edge_list:
            intersected_pts_pia = []
            intersected_pts_wm = []
            for i in range(len(obb_tree_list_pia)):
                obb_pia = obb_tree_list_pia[i]
                obb_wm = obb_tree_list_wm[i]
                points_pia = vtk.vtkPoints()
                points_wm = vtk.vtkPoints()
                code_pia = obb_pia.IntersectWithLine(edge[1], edge[0], points_pia, None)
                code_wm = obb_wm.IntersectWithLine(edge[1], edge[0], points_wm, None)
                if points_pia.GetNumberOfPoints() > 0:
                    intersected_pts_pia.append(points_pia.GetPoint(0))
                if points_wm.GetNumberOfPoints() > 0:
                    intersected_pts_wm.append(points_wm.GetPoint(0))

            avg_pt_pia = []
            avg_pt_wm = []
            std_pia = []
            std_wm = []

            if len(intersected_pts_pia) > 1:
                avg_pt_pia = np.array(intersected_pts_pia).mean(axis=0)
                std_pia = np.array(intersected_pts_pia).std(axis=0)
                pia_stds.append(std_pia)
                avg_pts_pia.append(avg_pt_pia)
                dists_3d =  []
                for pt in intersected_pts_pia:
                    dists_3d.append(self.get_euclidean_distance(pt,avg_pt_pia))
                dists_pia_3d.append(np.array(dists_3d).mean(axis=0))

            elif len(intersected_pts_pia) > 0:
                avg_pt_pia = intersected_pts_pia[0]
                avg_pts_pia.append(avg_pt_pia)

            if len(intersected_pts_wm) > 1:
                avg_pt_wm = np.array(intersected_pts_wm).mean(axis=0)
                std_wm = np.array(intersected_pts_wm).std(axis=0)
                wm_stds.append(std_wm)
                avg_pts_wm.append(avg_pt_wm)
                dists_3d =  []
                for pt in intersected_pts_wm:
                    dists_3d.append(self.get_euclidean_distance(pt,avg_pt_wm))
                dists_wm_3d.append(np.array(dists_3d).mean(axis=0))

            elif len(intersected_pts_wm) > 0:
                avg_pt_wm = intersected_pts_wm[0]
                avg_pts_wm.append(avg_pt_wm)

            if len(avg_pt_pia) > 0 and len(avg_pt_wm) > 0:
                wm_pia_rays_sg.graph_data.add_edge(avg_pt_wm,avg_pt_pia)
                if len(intersected_pts_pia) > 1 and len(intersected_pts_wm) > 1 and create_error_ellisoids==True:
                    # create arror ellipsoids
                    pia_ellipsoid = Surface(polydata=self.get_ellipsoid(std_pia[0],std_pia[1],std_pia[2]))
                    wm_ellipsoid =  Surface(polydata=self.get_ellipsoid(std_wm[0],std_wm[1],std_wm[2]))
                    pia_ellipsoid_txed,wm_ellipsoid_txed = self.transform_error_ellipsoid(pia_ellipsoid,wm_ellipsoid,[avg_pt_wm,avg_pt_pia])
                    #print(pia_ellipsoid_txed.surface.GetNumberOfCells(),wm_ellipsoid_txed.surface.GetNumberOfCells())
                    error_ellipsoids.append(pia_ellipsoid_txed.surface)
                    error_ellipsoids.append(wm_ellipsoid_txed.surface)

        # fill df
        df['Alignment_Type'] = [alignment_type]
        df['Surface_Type'] = [surf_type]
        df['Hemisphere'] = [hem]
        df['Aligned_Using_Surf_Resolution'] = [aligned_surf_res]
        df['Present_Surf_Resolution'] = [surf_res]
        df['Pia_Mean_Std_X'] = [np.array(pia_stds).mean(axis=0)[0]]
        df['Pia_Mean_Std_Y'] = [np.array(pia_stds).mean(axis=0)[1]]
        df['Pia_Mean_Std_Z'] =[ np.array(pia_stds).mean(axis=0)[2]]
        df['WM_Mean_Std_X'] = [np.array(wm_stds).mean(axis=0)[0]]
        df['WM_Mean_Std_Y'] = [np.array(wm_stds).mean(axis=0)[1]]
        df['WM_Mean_Std_Z'] = [np.array(wm_stds).mean(axis=0)[2]]
        df['Pia_Mean_Dist_3D'] = [np.array(dists_pia_3d).mean(axis=0)]
        df['Pia_STD_Dist_3D'] = [np.array(dists_pia_3d).std(axis=0)]
        df['Pia_CoV_Dist_3D'] = [np.array(dists_pia_3d).std(axis=0) / np.array(dists_pia_3d).mean(axis=0)]
        df['Pia_SEM_Dist_3D'] = [np.array(dists_pia_3d).std(axis=0) / np.sqrt(len(dists_pia_3d))]
        df['WM_Mean_Dist_3D'] = [np.array(dists_wm_3d).mean(axis=0)]
        df['WM_STD_Dist_3D'] = [np.array(dists_wm_3d).std(axis=0)]
        df['WM_CoV_Dist_3D'] = [np.array(dists_wm_3d).std(axis=0) / np.array(dists_wm_3d).mean(axis=0)]
        df['WM_SEM_Dist_3D'] = [np.array(dists_wm_3d).std(axis=0) / np.sqrt(len(dists_wm_3d))]

        return avg_pts_pia,avg_pts_wm,wm_pia_rays_sg,error_ellipsoids,df

    def get_surf_center_errors(self,aligned_pias_list_lhs,aligned_wms_list_lhs):
        centers_lhs = []
        surfs_lhs = Surface()
        dists = []
        for i in range(len(aligned_pias_list_lhs)):
            pia_wm = Surface()
            pia_wm.append(aligned_pias_list_lhs[i])
            pia_wm.append(aligned_wms_list_lhs[i])
            center = pia_wm.get_center_of_mass()
            centers_lhs.append(center)
            surfs_lhs.append(pia_wm.surface)
        centroid_lhs = np.array(centers_lhs).mean(axis=0)
        for i in range(len(centers_lhs)):
            dists.append(distance.euclidean(centers_lhs[i],centroid_lhs))
        #print(dists)
        return np.array(dists).mean(axis=0),np.array(dists).std(axis=0)

    def get_ellipsoid(self,rad_x,rad_y,rad_z):
        #print(rad_x,rad_y,rad_z)
        ell = vtk.vtkParametricEllipsoid()
        ell.SetXRadius(rad_x)
        ell.SetYRadius(rad_y)
        ell.SetZRadius(rad_z)

        pm = vtk.vtkParametricFunctionSource()
        pm.SetParametricFunction(ell)
        pm.Update()
        #print(pm.GetOutput().GetNumberOfCells())
        return pm.GetOutput()

    def transform_error_ellipsoid(self,pia_ellipsoid,wm_ellipsoid,edge):
        ''' Align error ellipsoids to its axis.. i.e translate and rotate from origin and y axis to wm-pia avg axis '''
        # do it for wm
        pts = []
        ref_pts = []
        ref_pts.append(edge[0])
        #ref_pts.append(edge[1])
        pts.append([0,0,0])
        #pts.append([0,-Vectors().get_vec_length(edge),0])
        icp,txmat = self.align_landmarks(ref_pts,pts)
        wm_ellipsoid.apply_icp_transform(icp)

        # for pia
        pts = []
        ref_pts = []
        ref_pts.append(edge[1])
        #pts.append([0,Vectors().get_vec_length(edge),0])
        pts.append([0,0,0])
        icp,txmat = self.align_landmarks(ref_pts,pts)
        pia_ellipsoid.apply_icp_transform(icp)

        return pia_ellipsoid,wm_ellipsoid

    def align_landmarks(self,ref_pts,pts,mode=0,):
        ''' This function performs rigid(+scaling) transformation of the input landmarks to register them to
        the given target landmarks and apply the transformation to other landmarks / surfaces provided'''

        # First convert the given landmarks into polydata in order to perform icp
        #ref_landmarks = Landmarks(pts=ref_matching_landmarks)
        #landmarks = Landmarks(pts=matching_landmarks)

        ##print(pt_list)
        icp = vtk.vtkLandmarkTransform()
        icp.SetSourceLandmarks(self.convertPointsToVtkPoints(pts))
        icp.SetTargetLandmarks(self.convertPointsToVtkPoints(ref_pts))
        if mode == 0:
            icp.SetModeToRigidBody()
        elif mode == 1:
            icp.SetModeToSimilarity()
        else:
            icp.SetModeToAffine()
        icp.Modified()
        icp.Update()

        tr_matrix = []
        for i in range(4):
            for j in range(4):
                tr_matrix.append(icp.GetMatrix().GetElement(j,i))


        return icp,tr_matrix

    def convertPointsToVtkPoints(self,pts = []):
        plain_pts = []
        if len(pts) > 0:
            plain_pts = pts
        else:
            plain_pts = self.pts
        vtkpoints = vtk.vtkPoints()
        #points = axon_front_contour_2D + axon_back_contour_2D + dend_front_contour_2D + dend_back_contour_2D
        for i in range(len(plain_pts)):
            vtkpoints.InsertPoint(i,plain_pts[i][0],plain_pts[i][1],plain_pts[i][2])

        return vtkpoints

    def check_ray_intersection(self,ray,surf):
        obb_tree = vtk.vtkOBBTree()
        obb_tree.SetDataSet(surf)
        obb_tree.BuildLocator()
        points = vtk.vtkPoints()
        code = obb_tree.IntersectWithLine(ray[1], ray[0], points, None)
        if points.GetNumberOfPoints()>0:
            return True
        else:
            return False

    def get_intersection_cube(self,surf_list):
        bounds = []
        for surf in surf_list:
            bounds.append(surf.surface.GetBounds())
            bounds.append(surf.surface.GetBounds())
            bounds.append(surf.surface.GetBounds())

        min_bounds = np.array(bounds).min(axis=0)
        max_bounds = np.array(bounds).max(axis=0)

        inner_cube = vtk.vtkCubeSource()
        inner_cube.SetBounds([max_bounds[0],min_bounds[1],max_bounds[2],min_bounds[3],max_bounds[4],min_bounds[5]])
        inner_cube.Update()
        inner_cube_diagonal = distance.euclidean([max_bounds[0],max_bounds[2],max_bounds[4]],[min_bounds[1],min_bounds[3],min_bounds[5]])
        inner_vol = inner_cube_diagonal*inner_cube_diagonal*inner_cube_diagonal

        outer_cube = vtk.vtkCubeSource()
        outer_cube.SetBounds([min_bounds[0],max_bounds[1],min_bounds[2],max_bounds[3],min_bounds[4],max_bounds[5]])
        outer_cube.Update()
        outer_cube_diagonal = distance.euclidean([min_bounds[0],min_bounds[2],min_bounds[4]],[max_bounds[1],max_bounds[3],max_bounds[5]])
        outer_vol = outer_cube_diagonal*outer_cube_diagonal*outer_cube_diagonal
        return outer_cube.GetOutput(),inner_cube.GetOutput(),inner_vol/outer_vol

    def get_distance_between_centroids(self,data_list):
        centers = []
        for surf in data_list:
            centers.append(surf.get_center_of_mass())

        centroid = np.array(centers).mean(axis=0)

        dists = []
        for center in centers:
            #print(center,centroid)
            dists.append(distance.euclidean(center,centroid))

        return dists

    def get_surface_volumes(self,surf_list):
        volumes = []
        for surf in surf_list:
            volumes.append(surf.get_surface_properties(prop='volume'))
        return volumes

    def get_delunay_volume(self,surf_list):
        all_surfs = Surface()
        for surf in surf_list:
            all_surfs.append(surf.surface)
        all_surfs.create_delunay_surface_3d()
        vol = all_surfs.get_surface_properties(prop='volume')
        return vol

    def get_contour_containing_vertex(self,vertex,edge_list):
        for edge in edge_list:
            for pt in edge:
                if pt == vertex:
                    return edge

    def get_surface_from_contours(self,sg,axis=1,output_filename=None):
        sorted_verts =  sorted( sg.graph_data.vertices, key=lambda coord: coord[axis] )

        sorted_edge_list = []
        for vert in sorted_verts:
            sorted_edge_list.append(self.get_contour_containing_vertex(vert,sg.graph_data.edge_list))

        surf = Surface()
        for i in range(len(sorted_edge_list)-1):
            pts = []
            for pt in sorted_edge_list[i]:
                pts.append(pt)
            for pt in sorted_edge_list[i+1]:
                pts.append(pt)

            hull = Surface(pts=pts).create_delunay_surface_3d(return_hull=True)
            surf.append(hull)

        if output_filename is not None:
            surf.write_surface_mesh(output_filename)

        return surf
