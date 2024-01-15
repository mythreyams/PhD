import vtk
import numpy as np
from spatial_graphs.Surfaces import Surface
from spatial_graphs.Vectors import Vectors
from spatial_graphs.AmiraSpatialGraph import AmiraSpatialGraph,Contours
from scipy.spatial import distance
#import SimpleITK as sitk
from scipy.spatial import distance_matrix

class Landmarks:
    ''' This class contains functions to handle a list of 3d coordinates in the amira format'''

    def __init__(self,filename='',pts=[],axis_directions=[1,1,1]):
        self.transformation_matrix = []
        self.transformation_applied = False
        self.axis_directions = axis_directions
        if len(filename)>0:
            self.pts = self.read_landmarks(filename)
        elif len(pts)>0:
            self.pts = pts
        else:
            self.pts = []

    def read_landmarks(self,filename):
        manual_landmark_list = []
        lines = []
        with open(filename, 'r') as csb:
            lines = csb.readlines()
        number = 0

        # for each manual landmark
        for line in lines:
            #print line
            if line.startswith("@1"):
                number = 1
                continue
            if number == 1 and line.isspace() == False:
                #print line.split()
                pt = list(map(float,line.split()))
                pt = [pt[0]*self.axis_directions[0],pt[1]*self.axis_directions[1],pt[2]*self.axis_directions[2]]
                manual_landmark_list.append(pt)

        return manual_landmark_list

    def apply_transformation(self,transformation_matrix=[],z_desired=None,trans_2d_only=False,inverse=False):

        if len(transformation_matrix) > 0:
            self.transformation_matrix = transformation_matrix
        # transformation matrix 4x4 as a row first list
        if len(self.pts)>0:
            self.transformation_applied = True
            tr_mat = np.reshape(self.transformation_matrix,[4,4])
            if z_desired is not None:
                # to translate in z to the desired position for section graphs
                tr_mat[2,3] = z_desired - self.pts[0][2]
            if trans_2d_only:
                tr_mat[2,2] = 1
                tr_mat[0,2] = 0
                tr_mat[1,2] = 0
                tr_mat[2,0] = 0
                tr_mat[2,1] = 0

            if inverse:
                self.pts = np.matmul((self.add_ones_column(self.pts)),np.transpose(tr_mat))
            else:
                self.pts = np.matmul((self.add_ones_column(self.pts)),tr_mat)

    def apply_translation(self,translation):
        self.pts = self.pts - translation

    def get_cylinder(self,radius,edge):
        center = (np.array(edge[0])+np.array(edge[1]))/2
        height = Vectors().get_vec_length(edge)
        cylinderSource = vtk.vtkCylinderSource()
        cylinderSource.SetCenter(center)
        cylinderSource.SetRadius(radius)
        cylinderSource.SetHeight(height)
        cylinderSource.SetResolution(10)
        cylinderSource.Update()

        pts = []
        ref_pts = []
        pts.append(center)
        end_pt = [center[0],cylinderSource.GetOutput().GetBounds()[2],center[2] ]
        pts.append(end_pt)

        ref_pts.append(center)
        end_pt = Vectors().get_unit_vec(center,edge[1])*height/2 + center
        ref_pts.append(end_pt)
        #Landmarks(pts=pts).write_landmarks(output_path_neun_profiles+'pts')
        #Landmarks(pts=ref_pts).write_landmarks(output_path_neun_profiles+'ref_pts')

        landmarks = Landmarks(pts=pts)
        icp,txmat = landmarks.align_landmarks(ref_pts)

        cyl_surf = Surface(polydata=cylinderSource.GetOutput())
        cyl_surf.apply_icp_transform(icp)
        #print(txmat)

        #Surface(polydata=cylinderSource.GetOutput()).write_surface_mesh(output_path_neun_profiles+'_cylinder.vtk')

        return cyl_surf.surface

    def apply_rotation(self,rot_mat=[],translation=None,trans_2d_only=False):

        rot_mat = np.reshape(rot_mat,[3,3])

        tmp = np.matmul(rot_mat,np.transpose(self.pts))
        self.pts = np.transpose(tmp)

        if translation is not None:
            self.pts = self.pts - translation

    def add_ones_column(self,data):
        a = np.array(data).transpose()
        b = [a[0],a[1],a[2],np.ones(len(data))]
        c = np.array(b).transpose()
        return c

    def write_landmarks(self,filename,pts_to_be_written=None,pts_2d=False):
        if pts_to_be_written is not None:
            pts_being_written = pts_to_be_written
        else:
            pts_being_written = self.pts
        with open(filename, 'w') as file:
            file.write("# Avizo 3D ASCII 2.0 \n")
            file.write("define Markers ")
            file.write("{}".format(len(pts_being_written)))
            file.write("\n")
            file.write("Parameters {\n")
            file.write("    NumSets 1,\n")
            file.write("    ContentType \"LandmarkSet\"\n")
            file.write("}\n")

            file.write("Markers { float[3] Coordinates } @1\n")

            file.write("# Data section follows\n")
            file.write("@1\n")

            if len(pts_being_written)!=0:
                for landmark in pts_being_written:
                    ##print(landmark[0][0])
                    file.write("{}".format(landmark[0]))
                    file.write(" ")
                    file.write("{}".format(landmark[1]))
                    file.write(" ")
                    if pts_2d:
                        file.write("{}".format(0))
                    else:
                        file.write("{}".format(landmark[2]))
                    file.write("\n")

    def append_landmarks(self,pts):
        for pt in pts:
            self.pts.append(pt)

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

    def align_landmarks(self,ref_pts,mode=0,):
        ''' This function performs rigid(+scaling) transformation of the input landmarks to register them to
        the given target landmarks and apply the transformation to other landmarks / surfaces provided'''

        # First convert the given landmarks into polydata in order to perform icp
        #ref_landmarks = Landmarks(pts=ref_matching_landmarks)
        #landmarks = Landmarks(pts=matching_landmarks)

        ##print(pt_list)
        icp = vtk.vtkLandmarkTransform()
        icp.SetSourceLandmarks(self.convertPointsToVtkPoints(self.pts))
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

    def set_z_coord(self,z):
        for i in range(len(self.pts)):
            if len(self.pts[i]) > 3:
                self.pts[i] = [self.pts[i][0],self.pts[i][1],z,0]
            else:
                self.pts[i] = [self.pts[i][0],self.pts[i][1],z]

    def add_z_offset(self,z):
            for i in range(len(self.pts)):
                if len(self.pts[i]) > 3:
                    self.pts[i] = [self.pts[i][0],self.pts[i][1],z+self.pts[i][2],0]
                else:
                    self.pts[i] = [self.pts[i][0],self.pts[i][1],z+self.pts[i][2]]

    def set_axis_directions(self,axis=[1,1,1]):
        for i in range(len(self.pts)):
            self.pts[i] = [self.pts[i][0]*axis[0],self.pts[i][1]*axis[1],self.pts[i][2]*axis[2]]

    def get_landmarks_in_bounding_box(self,bounds,pts_to_be_used=None):
        cropped_pts = []

        if pts_to_be_used is not None:
            pts_to_being_used = pts_to_be_used
        else:
            pts_to_being_used = self.pts
        for pt in pts_to_being_used:
            if pt[0]>= bounds[0] and pt[0]<= bounds[1] \
                    and pt[1]>= bounds[2] and pt[1]<= bounds[3] \
                    and pt[2]>= bounds[4] and pt[2]<= bounds[5]:

                cropped_pts.append(pt)

        #self.pts = cropped_pts
        return cropped_pts

    def get_contour_spatial_graph(self,binim,dim=2,translation=[0,0,0],scaling=1,voxel_sze=500):
        if dim==0:
            upper_bound = [0,1,1]
            lower_bound = [0,1,1]
        elif dim==1:
            upper_bound = [1,0,1]
            lower_bound = [1,0,1]
        else:
            upper_bound = [1,1,0]
            lower_bound = [1,1,0]

        binim = sitk.ConstantPad(sitk.Cast(binim,sitk.sitkUInt8),lower_bound,upper_bound)
        #sitk.WriteImage(contim,output_root+'contim.tif')
        cont = sitk.BinaryContourImageFilter()
        cont.SetBackgroundValue(255)
        cont.SetForegroundValue(0)
        cont.SetFullyConnected(False)
        contim = sitk.InvertIntensity(cont.Execute(sitk.RescaleIntensity(binim,0,255)))

        npiminds = np.where(sitk.GetArrayFromImage(contim)>0)
        if dim==0:
            #yz projection im
            contour_coords = np.transpose([npiminds[2],npiminds[1],npiminds[0]])
        elif dim==1:
            #xz
            contour_coords = np.transpose([npiminds[2],npiminds[1],npiminds[0]])
        else:
            #xy
            contour_coords = np.transpose([npiminds[1],npiminds[0],np.zeros(len(npiminds[0]))])
        sg = AmiraSpatialGraph(generic_graph=True)
        sg.graph_data.add_contour(contour_coords)
        sg.graph_data.apply_scaling([scaling,scaling,scaling])
        sg.graph_data.apply_transformation([[1,0,0,0], [0,1,0,0], [0,0,1,0], \
                                            [translation[0]-lower_bound[0]*voxel_sze,translation[1]-lower_bound[1]*voxel_sze,translation[2]-lower_bound[2]*voxel_sze,1]])

        return sg

    def get_size_of_region(self,connim,region_id):
        size = 0
        x, y, z = connim.GetDimensions()
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    val = connim.GetScalarComponentAsFloat(i, j, k, 0)
                    if val == region_id:
                        size = size + 1
        return size

    def get_binary_image(self,densim,threshold):

        thre = vtk.vtkImageThreshold()
        thre.SetInputData(densim)
        thre.SetInValue(255)
        thre.SetOutValue(0)
        thre.ThresholdByUpper(threshold)
        thre.Update()

        imconn = vtk.vtkImageConnectivityFilter()
        imconn.SetInputData(thre.GetOutput())
        #imconn.SetExtractionModeToLargestRegion()
        imconn.SetLabelScalarTypeToUnsignedChar()
        imconn.SetLabelModeToSizeRank()
        imconn.Update()

        label_im = imconn.GetOutput()
        bin_im = imconn.GetOutput()
        # iterate through densim and find the label of the peak density
        max_dens =  densim.GetScalarRange()[1]
        num_regions = imconn.GetExtractedRegionSizes()
        x,y,z = densim.GetDimensions()
        max_dens_region_id = -1
        max_dens_region_ids = []
        peak_dens_inds = []
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    val = densim.GetScalarComponentAsFloat(i,j,k,0)
                    if val == max_dens:
                        max_dens_region_ids.append(label_im.GetScalarComponentAsFloat(i,j,k,0))
                        print('region id ',label_im.GetScalarComponentAsFloat(i,j,k,0), 'max dens ',val)
                        #peak_dens_inds.append([i,j,k])

        peak_region_size = 0
        # in case there are multiple regions with peak densities, pick the one with max size
        for region_id in max_dens_region_ids:
            size = self.get_size_of_region(label_im,region_id)
            if size > peak_region_size:
                peak_region_size = size
                max_dens_region_id = region_id

        #print(max_dens_region_id)
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    val = label_im.GetScalarComponentAsFloat(i,j,k,0)
                    if val == max_dens_region_id:
                        peak_dens_inds.append([i, j, k])
                        bin_im.SetScalarComponentFromFloat(i,j,k,0,255)
                    else:
                        bin_im.SetScalarComponentFromFloat(i,j,k,0,0)
        peak_loc = np.array(peak_dens_inds).mean(axis=0)

        return thre.GetOutput(),bin_im,peak_loc

    def read_vtk_image(self,filename,orientation_type):
        reader = vtk.vtkTIFFReader()
        reader.SetFileName(filename)
        reader.SetOrientationType(orientation_type)
        reader.Update()
        return reader.GetOutput()

    def write_vtk_image(self,filename,im):
        writer = vtk.vtkTIFFWriter()
        writer.SetFileName(filename)
        writer.SetInputData(im)
        writer.Update()

    def itk_fill_hole(self,binim,radius):
        holefilled = sitk.BinaryDilateImageFilter()
        holefilled.SetKernelRadius(radius)
        holefilled.SetBackgroundValue(0)
        holefilled.SetForegroundValue(255)
        holefilled.SetKernelType(sitk.sitkBox)
        holefilledim = holefilled.Execute(binim)

        holefilled = sitk.BinaryErodeImageFilter()
        holefilled.SetKernelRadius(radius)
        holefilled.SetBackgroundValue(0)
        holefilled.SetForegroundValue(255)
        holefilled.SetKernelType(sitk.sitkBox)
        holefilledim = holefilled.Execute(holefilledim)

        #sitk.WriteImage(holefilledim,output_root+'holefilled.tif')
        return holefilledim

    def found_pt_at_voxel(self,x,y,z,im,pt_list):
        for pt in pt_list:
            ptx = (pt[0]-im.GetOrigin()[0])/im.GetSpacing()[0]
            pty = (pt[1]-im.GetOrigin()[1])/im.GetSpacing()[1]
            ptz = (pt[2]-im.GetOrigin()[2])/im.GetSpacing()[2]
            #print(ptx,pty,ptz)
            #print(im.GetOrigin())
            #ind = list(map(int,[ptx,pty,ptz]))
            if ptx >= x and ptx < x+1 and pty >= y and pty < y+1 and ptz >= z and ptx < z+1:
                return True
        return False

    def fill_surface_holes(self,surf,hole_size=10000):
        holefiller = vtk.vtkFillHolesFilter()
        holefiller.SetInputData(surf)
        holefiller.SetHoleSize(10000)
        holefiller.Update()
        return holefiller.GetOutput()

    def get_surface_from_image_volume(self,im,fill_hole=True):
        surf = vtk.vtkMarchingCubes()
        surf.SetInputData(im)
        surf.SetValue(0,1)
        surf.Update()
        if fill_hole:
            surface = self.fill_surface_holes(surf.GetOutput())
        else:
            surface = surf.GetOutput()
        #surface = Surface(polydata=surf.GetOutput())
        return Surface(polydata=surface)

    def increase_im_res(self,im,voxel_size1,voxel_size2):
        scale = int(voxel_size1/voxel_size2)
        arr = sitk.GetArrayFromImage(im)
        arr_50= np.zeros([arr.shape[0]*scale,arr.shape[1]*scale,arr.shape[2]*scale],dtype=np.uint8)
        print(arr_50.shape)
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                for k in range(arr.shape[2]):
                    arr_50[i*scale:((i+1)*scale),j*scale:((j+1)*scale),k*scale:((k+1)*scale)] = arr[i,j,k]
                    #print(arr[i,j,k])
        newim = sitk.GetImageFromArray(arr_50)
        return newim

    def create_contour_images(self,im,output_filename,lower_corner,voxel_size,spacing=[0,0,0],dens_im=False,holefillrad=None,padding=1):
        ''' Create XY, XZ and YZ binary projection images.
            Create contour images. Make them thick so that contour images can be visualized as volume using
            volren in Amira
        '''
        if dens_im:
            bin_xy = sitk.MaximumProjection(im,2)
            bin_xz = sitk.MaximumProjection(im,1)
            bin_yz = sitk.MaximumProjection(im,0)
            npim_xy = sitk.GetArrayFromImage(bin_xy)
            #print(npim_xy.shape)
            thickim_xy = np.zeros([npim_xy.shape[0]*2,npim_xy.shape[1],npim_xy.shape[2]],dtype=np.uint8)
            thickim_xy[0:1,:,:] = npim_xy
            thickim_xy[1:2,:,:] = npim_xy

            npim_xz = sitk.GetArrayFromImage(bin_xz)
            thickim_xz = np.zeros([npim_xz.shape[0],npim_xz.shape[1]*2,npim_xz.shape[2]],dtype=np.uint8)
            thickim_xz[:,0:1,:] = npim_xz
            thickim_xz[:,1:2,:] = npim_xz

            npim_yz = sitk.GetArrayFromImage(bin_yz)
            thickim_yz = np.zeros([npim_yz.shape[0],npim_yz.shape[1],npim_yz.shape[2]*2],dtype=np.uint8)
            thickim_yz[:,:,0:1] = npim_yz
            thickim_yz[:,:,1:2] = npim_yz

            sitk.WriteImage(sitk.GetImageFromArray(thickim_xy),output_filename+'_pos_x_{}_y_{}_z_{}_dens_XY.tif'.\
                            format(int(lower_corner[0]),int(lower_corner[1]),int(lower_corner[2])))
            sitk.WriteImage(sitk.GetImageFromArray(thickim_xz),output_filename+'_pos_x_{}_y_{}_z_{}_dens_XZ.tif'.\
                            format(int(lower_corner[0]),int(lower_corner[1]),int(lower_corner[2])))
            sitk.WriteImage(sitk.GetImageFromArray(thickim_yz),output_filename+'_pos_x_{}_y_{}_z_{}_dens_YZ.tif'.\
                            format(int(lower_corner[0]),int(lower_corner[1]),int(lower_corner[2])))

        else:
            # pad to close the boundary
            # Recast binary image with smaller voxel size for a tighter fit
            scale = int(voxel_size/50)
            if holefillrad is not None:
                im = self.itk_fill_hole(sitk.Cast(im,sitk.sitkUInt8),radius=holefillrad)

            highresim = self.increase_im_res(im,voxel_size,50)
            bin_xy = sitk.MaximumProjection(highresim,2)
            bin_xz = sitk.MaximumProjection(highresim,1)
            bin_yz = sitk.MaximumProjection(highresim,0)

            paddedim_xy = sitk.ConstantPad(sitk.Cast(bin_xy,sitk.sitkUInt8),[padding,padding,0],[padding,padding,0])
            paddedim_xz = sitk.ConstantPad(sitk.Cast(bin_xz,sitk.sitkUInt8),[padding,0,padding],[padding,0,padding])
            paddedim_yz = sitk.ConstantPad(sitk.Cast(bin_yz,sitk.sitkUInt8),[0,padding,padding],[0,padding,padding])

            #holefilledim_xy = self.itk_fill_hole(paddedim_xy,radius=int(scale/2))
            #holefilledim_xz = self.itk_fill_hole(paddedim_xz,radius=int(scale/2))
            #holefilledim_yz = self.itk_fill_hole(paddedim_yz,radius=int(scale/2))

            # Find contour
            contourim_xy = sitk.BinaryContour(paddedim_xy,fullyConnected=True,backgroundValue=255,foregroundValue=0)
            contourim_inv_xy = sitk.InvertIntensity(contourim_xy)

            contourim_xz = sitk.BinaryContour(paddedim_xz,fullyConnected=True,backgroundValue=255,foregroundValue=0)
            contourim_inv_xz = sitk.InvertIntensity(contourim_xz)

            contourim_yz = sitk.BinaryContour(paddedim_yz,fullyConnected=True,backgroundValue=255,foregroundValue=0)
            contourim_inv_yz = sitk.InvertIntensity(contourim_yz)

            # convert 2D image to 3D
            npim_xy = sitk.GetArrayFromImage(contourim_inv_xy)
            #print(npim_xy.shape)
            thickim_xy = np.zeros([npim_xy.shape[0]*2,npim_xy.shape[1],npim_xy.shape[2]],dtype=np.uint8)
            thickim_xy[0:1,:,:] = npim_xy
            thickim_xy[1:2,:,:] = npim_xy

            npim_xz = sitk.GetArrayFromImage(contourim_inv_xz)
            thickim_xz = np.zeros([npim_xz.shape[0],npim_xz.shape[1]*2,npim_xz.shape[2]],dtype=np.uint8)
            thickim_xz[:,0:1,:] = npim_xz
            thickim_xz[:,1:2,:] = npim_xz

            npim_yz = sitk.GetArrayFromImage(contourim_inv_yz)
            thickim_yz = np.zeros([npim_yz.shape[0],npim_yz.shape[1],npim_yz.shape[2]*2],dtype=np.uint8)
            thickim_yz[:,:,0:1] = npim_yz
            thickim_yz[:,:,1:2] = npim_yz

            sitk.WriteImage(sitk.GetImageFromArray(thickim_xy),output_filename+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_Contour_XY.tif'.\
                            format(int(lower_corner[0]-(voxel_size/2)),int(lower_corner[1]-(voxel_size/2)),int(lower_corner[2]-(voxel_size/2)),int(spacing[0]/scale),int(spacing[1]/scale),int(spacing[2]/scale)))
            sitk.WriteImage(sitk.GetImageFromArray(thickim_xz),output_filename+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}Contour_XZ.tif'.\
                            format(int(lower_corner[0]-(voxel_size/2)),int(lower_corner[1]-(voxel_size/2)),int(lower_corner[2]-(voxel_size/2)),int(spacing[0]/scale),int(spacing[1]/scale),int(spacing[2]/scale)))
            sitk.WriteImage(sitk.GetImageFromArray(thickim_yz),output_filename+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}Contour_YZ.tif'.\
                            format(int(lower_corner[0]-(voxel_size/2)),int(lower_corner[1]-(voxel_size/2)),int(lower_corner[2]-(voxel_size/2)),int(spacing[0]/scale),int(spacing[1]/scale),int(spacing[2]/scale)))

        return

    def get_high_voxel_centers(self,binim):
        voxel_centers = []
        for x in range(binim.GetDimensions()[0]):
            for y in range(binim.GetDimensions()[1]):
                for z in range(binim.GetDimensions()[2]):

                    if binim.GetScalarComponentAsFloat(x,y,z,0)>0:
                        np_ind = np.array([x,y,z])
                        voxel_center = (np_ind)*binim.GetSpacing() + binim.GetOrigin()
                        voxel_centers.append(voxel_center)
        return voxel_centers

    def get_pts_within_voxel_radius(self,pt_list,radius,voxel_centers):
        selected_pts  = []
        for pt in pt_list:
            for voxel_center in voxel_centers:
                if distance.euclidean(voxel_center,pt)<=radius:
                    selected_pts.append(pt)
                    break
        return selected_pts

    def found_in_list(self,pt_list,pt):
        for pt1 in pt_list:
            if pt1[0]==pt[0] and pt1[1]==pt[1] and pt1[2]==pt[2]:
                return True
        return False

    def combine_landmarks_without_duplication(self,landmarks1,landmarks2):
        combined_landmarks = []
        for pt in landmarks1:
            combined_landmarks.append(pt)

        for pt in landmarks2:
            if not self.found_in_list(combined_landmarks,pt):
                combined_landmarks.append(pt)

        return combined_landmarks

    def get_density_image(self,output_file,voxel_size=400):
        pts_pdata = Surface(pts=self.pts)

        # 1. Compute density volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.5)
        radius = voxel_size*np.sqrt(3)/2
        densfilt.SetRadius(radius)
        densfilt.SetDensityFormToNumberOfPoints()
        bounds = pts_pdata.surface.GetBounds()
        x = abs(bounds[1] - bounds[0])
        y = abs(bounds[3] - bounds[2])
        z = abs(bounds[5] - bounds[4])
        lower_corner = [round(bounds[0]-voxel_size/2),round(bounds[2]-voxel_size/2),round(bounds[4]-voxel_size/2)]
        densfilt.SetSampleDimensions(int(np.ceil(x/voxel_size))+1,int(np.ceil(y/voxel_size))+1,int(np.ceil(z/voxel_size))+1)
        densfilt.SetInputData(pts_pdata.surface)
        #densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.SetModelBounds(bounds[0]-voxel_size/2,bounds[1]+voxel_size/2,\
                                bounds[2]-voxel_size/2,bounds[3]+voxel_size/2,\
                                bounds[4]-voxel_size/2,bounds[5]+voxel_size/2)
        densfilt.Update()
        densimage = densfilt.GetOutput()
        #print(densimage.GetScalarRange())
        imspacing = [int(densimage.GetSpacing()[0]),int(densimage.GetSpacing()[1]),int(densimage.GetSpacing()[2])]
        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]),densimage)
        densimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]))
        self.create_contour_images(densimage_itk,output_file+'_density_image_',lower_corner,voxel_size,dens_im=True)

        return densimage

    def get_binary_and_contour_images_from_high_voxels(self,high_voxels,output_file,voxel_size=400,cutoff=0.1,holefillrad=None):
        pts_pdata = Surface(pts=self.pts)
        #
        surf = Surface(output_file+'_surface.vtk')

        # 1. Compute density volume image to create the volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.5)
        radius = voxel_size*np.sqrt(3)/2
        densfilt.SetRadius(radius)
        densfilt.SetDensityFormToNumberOfPoints()
        bounds = pts_pdata.surface.GetBounds()
        x_range = abs(bounds[1] - bounds[0])
        y_range = abs(bounds[3] - bounds[2])
        z_range = abs(bounds[5] - bounds[4])
        lower_corner = [round(bounds[0]-voxel_size/2),round(bounds[2]-voxel_size/2),round(bounds[4]-voxel_size/2)]
        densfilt.SetSampleDimensions(int(np.ceil(x_range/voxel_size))+1,int(np.ceil(y_range/voxel_size))+1,int(np.ceil(z_range/voxel_size))+1)
        densfilt.SetInputData(pts_pdata.surface)
        #densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.SetModelBounds(bounds[0]-voxel_size/2,bounds[1]+voxel_size/2,\
                                bounds[2]-voxel_size/2,bounds[3]+voxel_size/2,\
                                bounds[4]-voxel_size/2,bounds[5]+voxel_size/2)
        densfilt.Update()
        densimage = densfilt.GetOutput()
        imspacing = [int(densimage.GetSpacing()[0]),int(densimage.GetSpacing()[1]),int(densimage.GetSpacing()[2])]
        binim = densimage

        # clear out the binim before assigning high voxels
        for x in range(densimage.GetDimensions()[0]):
            for y in range(densimage.GetDimensions()[1]):
                for z in range(densimage.GetDimensions()[2]):
                    binim.SetScalarComponentFromFloat(x,y,z,0,0)

        # surf = Surface(output_file+'_surface.vtk')
        # voxels = []
        # for x in range(densimage.GetDimensions()[0]):
        #     for y in range(densimage.GetDimensions()[1]):
        #         for z in range(densimage.GetDimensions()[2]):
        #             ptx = x*densimage.GetSpacing()[0] + densimage.GetOrigin()[0]
        #             pty = y*densimage.GetSpacing()[1] + densimage.GetOrigin()[1]
        #             ptz = z*densimage.GetSpacing()[2] + densimage.GetOrigin()[2]
        #             voxels.append([ptx,pty,ptz])
        # Landmarks(pts=voxels).write_landmarks((output_file + 'voxel_centers.landmarksascii'))
        # selected_pts = Landmarks(pts=voxels).get_landmarks_within_given_surface_using_vtk(surf.surface)
        # Landmarks(pts=selected_pts).write_landmarks((output_file + 'voxel_centers_within_surface.landmarksascii'))

        selected_pts = high_voxels
        if len(selected_pts)==0 or cutoff==1:
            com = surf.get_center_of_mass()
            print('single voxel case')
            # Set the center of surface as the high voxel
            ptx = round((com[0] - densimage.GetOrigin()[0]) / densimage.GetSpacing()[0])
            pty = round((com[1] - densimage.GetOrigin()[1]) / densimage.GetSpacing()[1])
            ptz = round((com[2] - densimage.GetOrigin()[2]) / densimage.GetSpacing()[2])
            binim.SetScalarComponentFromFloat(ptx, pty, ptz, 0, 255)
        else:
            for pt in selected_pts:
                ptx = round((pt[0]-densimage.GetOrigin()[0])/densimage.GetSpacing()[0])
                pty = round((pt[1]-densimage.GetOrigin()[1])/densimage.GetSpacing()[1])
                ptz = round((pt[2]-densimage.GetOrigin()[2])/densimage.GetSpacing()[2])
                binim.SetScalarComponentFromFloat(ptx,pty,ptz,0,255)

        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
                             format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])),binim)

        binimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
                             format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])))

        if cutoff==0:
            self.create_contour_images(binimage_itk,output_file+'_binary_image_',[lower_corner[0]-20,lower_corner[1]-20,lower_corner[2]-20],voxel_size,\
                                   spacing=[int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])],holefillrad=holefillrad,padding=20)
        else:
            self.create_contour_images(binimage_itk,output_file+'_binary_image_',lower_corner,voxel_size,\
                                   spacing=[int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])],holefillrad=holefillrad,padding=1)

    def get_binary_and_contour_images_from_pts(self,output_file,voxel_size=400,cutoff=0.1,holefillrad=None):
        pts_pdata = Surface(pts=self.pts)
        #
        surf = Surface(output_file+'_surface.vtk')

        # 1. Compute density volume image to create the volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.5)
        radius = voxel_size*np.sqrt(3)/2
        densfilt.SetRadius(radius)
        densfilt.SetDensityFormToNumberOfPoints()
        bounds = pts_pdata.surface.GetBounds()
        x_range = abs(bounds[1] - bounds[0])
        y_range = abs(bounds[3] - bounds[2])
        z_range = abs(bounds[5] - bounds[4])
        lower_corner = [round(bounds[0]-voxel_size/2),round(bounds[2]-voxel_size/2),round(bounds[4]-voxel_size/2)]
        densfilt.SetSampleDimensions(int(np.ceil(x_range/voxel_size))+1,int(np.ceil(y_range/voxel_size))+1,int(np.ceil(z_range/voxel_size))+1)
        densfilt.SetInputData(pts_pdata.surface)
        #densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.SetModelBounds(bounds[0]-voxel_size/2,bounds[1]+voxel_size/2,\
                                bounds[2]-voxel_size/2,bounds[3]+voxel_size/2,\
                                bounds[4]-voxel_size/2,bounds[5]+voxel_size/2)
        densfilt.Update()
        densimage = densfilt.GetOutput()
        imspacing = [int(densimage.GetSpacing()[0]),int(densimage.GetSpacing()[1]),int(densimage.GetSpacing()[2])]
        binim = densimage

        # clear out the binim before assigning high voxels
        for x in range(densimage.GetDimensions()[0]):
            for y in range(densimage.GetDimensions()[1]):
                for z in range(densimage.GetDimensions()[2]):
                    binim.SetScalarComponentFromFloat(x,y,z,0,0)

        surf = Surface(output_file+'_surface.vtk')
        voxels = []
        for x in range(densimage.GetDimensions()[0]):
            for y in range(densimage.GetDimensions()[1]):
                for z in range(densimage.GetDimensions()[2]):
                    ptx = x*densimage.GetSpacing()[0] + densimage.GetOrigin()[0]
                    pty = y*densimage.GetSpacing()[1] + densimage.GetOrigin()[1]
                    ptz = z*densimage.GetSpacing()[2] + densimage.GetOrigin()[2]
                    voxels.append([ptx,pty,ptz])
        Landmarks(pts=voxels).write_landmarks((output_file + 'voxel_centers.landmarksascii'))
        selected_pts = Landmarks(pts=voxels).get_landmarks_within_given_surface_using_vtk(surf.surface)
        Landmarks(pts=selected_pts).write_landmarks((output_file + 'voxel_centers_within_surface.landmarksascii'))

        if len(selected_pts)==0 or cutoff==1:
            com = surf.get_center_of_mass()
            print('single voxel case')
            # Set the center of surface as the high voxel
            ptx = round((com[0] - densimage.GetOrigin()[0]) / densimage.GetSpacing()[0])
            pty = round((com[1] - densimage.GetOrigin()[1]) / densimage.GetSpacing()[1])
            ptz = round((com[2] - densimage.GetOrigin()[2]) / densimage.GetSpacing()[2])
            binim.SetScalarComponentFromFloat(ptx, pty, ptz, 0, 255)
        else:
            for pt in selected_pts:
                ptx = round((pt[0]-densimage.GetOrigin()[0])/densimage.GetSpacing()[0])
                pty = round((pt[1]-densimage.GetOrigin()[1])/densimage.GetSpacing()[1])
                ptz = round((pt[2]-densimage.GetOrigin()[2])/densimage.GetSpacing()[2])
                binim.SetScalarComponentFromFloat(ptx,pty,ptz,0,255)

        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
                             format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])),binim)

        binimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
                             format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])))

        self.create_contour_images(binimage_itk,output_file+'_binary_image_',lower_corner,voxel_size,\
                                   spacing=[int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])],holefillrad=holefillrad)

        # #print('Number of cells {}'.format(surf.surface.GetNumberOfCells()))
        # # if the area covered by landmarks is less than a voxel of 400 then use the whole voxel. otherwise make an iso contour
        # if (x_range/voxel_size<1) and (y_range/voxel_size<1) and (z_range/voxel_size<1):
        #     print('only one voxel is present')
        #     for x in range(densimage.GetDimensions()[0]):
        #         for y in range(densimage.GetDimensions()[1]):
        #             for z in range(densimage.GetDimensions()[2]):
        #                 # if x>0 and x<densimage.GetDimensions()[0]-1 and\
        #                 #     y>0 and y<densimage.GetDimensions()[1]-1 and\
        #                 #     z>0 and z<densimage.GetDimensions()[2]-1 :
        #                 #     densimage.SetScalarComponentFromFloat(x,y,z,0,255)
        #                 # else:
        #                 binim.SetScalarComponentFromFloat(x,y,z,0,255)
        # else:
        #     for pt in self.pts:
        #         ptx = round((pt[0]-densimage.GetOrigin()[0])/densimage.GetSpacing()[0])
        #         pty = round((pt[1]-densimage.GetOrigin()[1])/densimage.GetSpacing()[1])
        #         ptz = round((pt[2]-densimage.GetOrigin()[2])/densimage.GetSpacing()[2])
        #         binim.SetScalarComponentFromFloat(ptx,pty,ptz,0,255)
        #     # for x in range(densimage.GetDimensions()[0]):
        #     #     for y in range(densimage.GetDimensions()[1]):
        #     #         for z in range(densimage.GetDimensions()[2]):
        #     #             #ind = list(map(int,[x,y,z]))
        #     #             if self.found_pt_at_voxel(x,y,z,densimage,self.pts):
        #     #                 print(x,y,z)
        #     #                 densimage.SetScalarComponentFromFloat(x,y,z,0,255)
        #     #             else: #== M1_region_id+1:
        #     #                 #print('removing {},{},{}'.format(x,y,z))
        #     #                 densimage.SetScalarComponentFromFloat(x,y,z,0,0)
        # surf = self.get_surface_from_image_volume(binim)
        # surf.write_surface_mesh(output_file+'tmpsurf.vtk')
        #
        # voxels = []
        # for x in range(densimage.GetDimensions()[0]):
        #     for y in range(densimage.GetDimensions()[1]):
        #         for z in range(densimage.GetDimensions()[2]):
        #             ptx = x*densimage.GetSpacing()[0] + densimage.GetOrigin()[0]
        #             pty = y*densimage.GetSpacing()[1] + densimage.GetOrigin()[1]
        #             ptz = z*densimage.GetSpacing()[2] + densimage.GetOrigin()[2]
        #             voxels.append([ptx,pty,ptz])
        # selected_pts = self.get_landmarks_within_given_surface_using_vtk(surf.surface)
        # Landmarks(pts=selected_pts).write_landmarks(output_file+'pts')
        # for pt in selected_pts:
        #     ptx = round((pt[0]-densimage.GetOrigin()[0])/densimage.GetSpacing()[0])
        #     pty = round((pt[1]-densimage.GetOrigin()[1])/densimage.GetSpacing()[1])
        #     ptz = round((pt[2]-densimage.GetOrigin()[2])/densimage.GetSpacing()[2])
        #     binim.SetScalarComponentFromFloat(ptx,pty,ptz,0,255)

        #self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
        #                         format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]),densimage)
        #self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
        #                     format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])),binim)
        #self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),binim)

        #densimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
        #                         format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]))
        #binimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_binary_image.tif'.\
        #                     format(lower_corner[0],lower_corner[1],lower_corner[2],int(binim.GetSpacing()[0]),int(binim.GetSpacing()[1]),int(binim.GetSpacing()[2])))
        #connimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image.tif'.\
        #                     format(lower_corner[0],lower_corner[1],lower_corner[2]))

        #connim_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))

        #self.create_contour_images(densimage_itk,output_file+'_density_image_',lower_corner,voxel_size,dens_im=True)
        #self.create_contour_images(connimage_itk,output_file+'_binary_image_',lower_corner,voxel_size)

    def get_density_cluster(self,output_file,voxel_size=500,dens_cutoff_threshold=0.1,write_density_images=False,min_threshold=0):

        pts_pdata = Surface(pts=self.pts)

        # 1. Compute density volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.5)
        radius = voxel_size*np.sqrt(3)/2
        densfilt.SetRadius(radius)
        densfilt.SetDensityFormToNumberOfPoints()
        bounds = pts_pdata.surface.GetBounds()
        x = abs(bounds[1] - bounds[0])
        y = abs(bounds[3] - bounds[2])
        z = abs(bounds[5] - bounds[4])
        lower_corner = [round(bounds[0]-voxel_size/2),round(bounds[2]-voxel_size/2),round(bounds[4]-voxel_size/2)]
        densfilt.SetSampleDimensions(int(np.ceil(x/voxel_size))+1,int(np.ceil(y/voxel_size))+1,int(np.ceil(z/voxel_size))+1)
        densfilt.SetInputData(pts_pdata.surface)
        #densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.SetModelBounds(bounds[0]-voxel_size/2,bounds[1]+voxel_size/2,\
                                bounds[2]-voxel_size/2,bounds[3]+voxel_size/2,\
                                bounds[4]-voxel_size/2,bounds[5]+voxel_size/2)
        densfilt.Update()
        densimage = densfilt.GetOutput()
        imspacing = [int(densimage.GetSpacing()[0]),int(densimage.GetSpacing()[1]),int(densimage.GetSpacing()[2])]
        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]),densimage)

        # 2. Binarize the image based on density cut off (10% of max density)
        max_dens =  densimage.GetScalarRange()[1]
        min_dens =  densimage.GetScalarRange()[0]
        if dens_cutoff_threshold is 0:
            # the cutoff given is 0. so need to use the threshold to be the least non zero density value
            densimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]))
            densimage_np = sitk.GetArrayFromImage(densimage_itk)
            values = densimage_np[np.where(densimage_np>0)]
            min_value = np.min(values)
            threshold = min_value
        else:
            threshold = max_dens*dens_cutoff_threshold

        # binarize the density image for given thresold and find find largest connected image
        binim,connim,peak_loc = self.get_binary_image(densimage,threshold)

        # write location of peak density
        np_peak_loc = np.array(peak_loc)
        peak_loc_voxel_center = (np_peak_loc)*binim.GetSpacing() + binim.GetOrigin()
        l = Landmarks(pts=[peak_loc_voxel_center])
        #l.apply_transformation([[1,0,0,0], [0,1,0,0], [0,0,1,0], [lower_corner[0],lower_corner[1],lower_corner[2],1]])
        l.write_landmarks(output_file+'_peak_dens_centroid.landmarksAscii')

        # write 3D density and binary images
        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]),densimage)
        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),connim)
        self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),binim)

        # 3. Create surface from the binary image and select the landmarks
        # Create surfaces from the binary volume
        surface=self.get_surface_from_image_volume(connim)
        surface.write_surface_mesh(output_file+'_surface.vtk')
        surface_raw=self.get_surface_from_image_volume(binim)
        surface_raw.write_surface_mesh(output_file+'_surface_raw.vtk')

        # capture landmarks within the surface
        high_voxel_centers_raw = self.get_high_voxel_centers(binim)
        high_voxel_centers = self.get_high_voxel_centers(connim)
        Landmarks(pts=high_voxel_centers_raw).write_landmarks(output_file+'_high_voxel_centers_raw.LandmarkAscii')
        Landmarks(pts=high_voxel_centers).write_landmarks(output_file+'_high_voxel_centers.LandmarkAscii')

        # get the border landmarks that are convered by the surface but might be missed by the other method
        selected_pts1 = self.get_landmarks_within_given_surface_using_vtk(surface.surface)
        selected_pts_raw1 = self.get_landmarks_within_given_surface_using_vtk(surface_raw.surface)

        selected_pts2 = self.get_pts_within_voxel_radius(self.pts,radius,high_voxel_centers)
        selected_pts_raw2 = self.get_pts_within_voxel_radius(self.pts,radius,high_voxel_centers_raw)

        selected_pts = self.combine_landmarks_without_duplication(selected_pts1,selected_pts2)
        selected_pts_raw = self.combine_landmarks_without_duplication(selected_pts_raw1,selected_pts_raw2)

        self.write_landmarks(output_file,selected_pts)
        self.write_landmarks(output_file+'_raw.LandmarkAscii',selected_pts_raw)

        centroid_raw = np.array(selected_pts_raw).mean(axis=0)
        centroid = np.array(selected_pts).mean(axis=0)
        self.write_landmarks(output_file+'_raw_Centroid.LandmarkAscii',[centroid_raw])
        self.write_landmarks(output_file+'_Centroid.LandmarkAscii',[centroid])

        # 4. write projection and contour images for iso contours visualization if needed
        if write_density_images:
            densimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_spacing_x_{}_y_{}_z_{}_density_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2],imspacing[0],imspacing[1],imspacing[2]))
            binimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_raw.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2]))
            connimage_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image.tif'.\
                                 format(lower_corner[0],lower_corner[1],lower_corner[2]))

            #connim_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))

            self.create_contour_images(densimage_itk,output_file+'_density_image_',lower_corner,voxel_size,dens_im=True)
            self.create_contour_images(binimage_itk,output_file+'_binary_image_raw_',lower_corner,voxel_size)
            self.create_contour_images(connimage_itk,output_file+'_binary_image_',lower_corner,voxel_size)

        return selected_pts,max_dens,dens_cutoff_threshold,max_dens*dens_cutoff_threshold
        # Create contour spatial graph
            #xy_densim = self.read_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xy.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),2)
            #yz_densim = self.read_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_yz.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),1)
            #xz_densim = self.read_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xz.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),2)

            #binim_xy,connim_xy,peak_loc_xy = self.get_binary_image(xy_densim,threshold)
            #binim_yz,connim_yz,peak_loc_yz = self.get_binary_image(yz_densim,threshold)
            #binim_xz,connim_xz,peak_loc_xz = self.get_binary_image(xz_densim,threshold)

            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xy_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),connim_xy)
            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_yz_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),connim_yz)
            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xz_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),connim_xz)
            #
            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xy_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),binim_xy)
            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_yz_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),binim_yz)
            # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_max_xz_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),binim_xz)
            #
            # connim_xy=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_xy_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            # connim_yz=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_yz_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            # connim_xz=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_xz_bin.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            #
            # binim_xy_raw=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_xy_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            # binim_yz_raw=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_yz_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            # binim_xz_raw=sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_max_xz_bin_raw.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
            #
            # holefilledim_xy = self.itk_fill_hole(connim_xy)
            # holefilledim_yz = self.itk_fill_hole(connim_yz)
            # holefilledim_xz = self.itk_fill_hole(connim_xz)
            #
            # contour_xy_sg = self.get_contour_spatial_graph(holefilledim_xy,dim=2,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            # contour_yz_sg = self.get_contour_spatial_graph(holefilledim_yz,dim=0,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            # contour_xz_sg = self.get_contour_spatial_graph(holefilledim_xz,dim=1,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            #
            # contour_xy_sg_raw = self.get_contour_spatial_graph(binim_xy_raw,dim=2,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            # contour_yz_sg_raw = self.get_contour_spatial_graph(binim_yz_raw,dim=0,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            # contour_xz_sg_raw = self.get_contour_spatial_graph(binim_xz_raw,dim=1,translation=[lower_corner[0],lower_corner[1],lower_corner[2]],scaling=voxel_size)
            #
            # contour_xy_sg.write_spatial_graph(output_file+'_contour_xy.am')
            # contour_yz_sg.write_spatial_graph(output_file+'_contour_yz.am')
            # contour_xz_sg.write_spatial_graph(output_file+'_contour_xz.am')
            #
            # contour_xy_sg_raw.write_spatial_graph(output_file+'_contour_xy_raw.am')
            # contour_yz_sg_raw.write_spatial_graph(output_file+'_contour_yz_raw.am')
            # contour_xz_sg_raw.write_spatial_graph(output_file+'_contour_xz_raw.am')
            #
            # print('max density is {}'.format(max_dens))
            # print('cutoff is {}'.format(dens_cutoff_threshold))
            # print('threshold {}'.format(max_dens*dens_cutoff_threshold))
        # valid_cluster_raw = []
        # for pt in self.pts:
        #     x = (pt[0]-binim.GetOrigin()[0])/binim.GetSpacing()[0]
        #     y = (pt[1]-binim.GetOrigin()[1])/binim.GetSpacing()[1]
        #     z = (pt[2]-binim.GetOrigin()[2])/binim.GetSpacing()[2]
        #     ind = list(map(int,[x,y,z]))
        #
        #     if binim.GetScalarComponentAsFloat(ind[0],ind[1],ind[2],0)  >0 : #== M1_region_id+1:
        #         valid_cluster_raw.append(pt)
        #
        # valid_cluster = []
        # for pt in self.pts:
        #     x = (pt[0]-connim.GetOrigin()[0])/connim.GetSpacing()[0]
        #     y = (pt[1]-connim.GetOrigin()[1])/connim.GetSpacing()[1]
        #     z = (pt[2]-connim.GetOrigin()[2])/connim.GetSpacing()[2]
        #     ind = list(map(int,[x,y,z]))
        #
        #     if connim.GetScalarComponentAsFloat(ind[0],ind[1],ind[2],0)  >0 : #== M1_region_id+1:
        #         valid_cluster.append(pt)
        # # clean up by removing outlier based on distance from the nearest pt
        # if len(cluster_in_bounds)>1:
        #     npl = np.array(cluster_in_bounds)
        #     dist = distance_matrix(npl,npl)
        #     for i in range(len(dist)):
        #         dist[i,i] =100000
        #     outlier_pts = npl[dist.min(axis=0)>voxel_size]
        #     selected_pts = npl[dist.min(axis=0)<=voxel_size]
        # else:
        #     selected_pts = cluster_in_bounds
        #
        # # make surface for the selected_pts
        # #selected_im_np = np.zeros([binim_itk.GetDepth(),binim_itk.GetHeight(),binim_itk.GetWidth()],dtype=np.uint8)
        # #print(outlier_pts)
        # cleaned_bin_im = vtk.vtkImageData()
        # cleaned_bin_im = connim
        # #for pt in outlier_pts:
        # for x in range(connim.GetDimensions()[1]):
        #     for y in range(connim.GetDimensions()[2]):
        #         for z in range(connim.GetDimensions()[0]):
        #             #x = (pt[0]-connim.GetOrigin()[0])/connim.GetSpacing()[0]
        #             #y = (pt[1]-connim.GetOrigin()[1])/connim.GetSpacing()[1]
        #             #z = (pt[2]-connim.GetOrigin()[2])/connim.GetSpacing()[2]
        #             #ind = list(map(int,[x,y,z]))
        #             if self.found_pt_at_voxel(x,y,z,connim,selected_pts):
        #                 cleaned_bin_im.SetScalarComponentFromFloat(x,y,z,0,255)
        #             else: #== M1_region_id+1:
        #                 #print('removing {},{},{}'.format(x,y,z))
        #                 cleaned_bin_im.SetScalarComponentFromFloat(x,y,z,0,0)
        #
        # self.write_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_cleaned.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),cleaned_bin_im)
        # connim_cleaned_itk = sitk.ReadImage(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_cleaned.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
        # sitk.WriteImage(sitk.ConstantPad(connim_cleaned_itk,[1,1,1],[1,1,1]),output_file+'_pos_x_{}_y_{}_z_{}_binary_image_cleaned_padded.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
        # paddedimcleaned = self.read_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_cleaned_padded.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),0)
        # #selected_im_itk = sitk.GetImageFromArray(selected_im_np)
        # #sitk.WriteImage(selected_im_itk,output_file+'_pos_x_{}_y_{}_z_{}_binary_image_selected.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]))
        # #selected_im_vtk = self.read_vtk_image(output_file+'_pos_x_{}_y_{}_z_{}_binary_image_selected.tif'.format(lower_corner[0],lower_corner[1],lower_corner[2]),0)
        # surf = vtk.vtkMarchingCubes()
        # surf.SetInputData(paddedimcleaned)
        # surf.SetValue(0,1)
        # surf.Update()
        # surface_selected = Surface(polydata=surf.GetOutput())
        # surface_selected.apply_transformation([[1,0,0,0], [0,1,0,0], [0,0,1,0], \
        #                                     [lower_corner[0]-voxel_size,lower_corner[1]-voxel_size,lower_corner[2]-voxel_size,1]])
        #
        # Surface(polydata=surface_selected.surface).write_surface_mesh(output_file+'_surface_selected.vtk')
        #
        # high_voxel_centers_raw = self.get_high_voxel_centers(binim)
        # high_voxel_centers = self.get_high_voxel_centers(connim)
        # Landmarks(pts=high_voxel_centers_raw).write_landmarks(output_file+'_high_voxel_centers_raw.LandmarkAscii')
        # Landmarks(pts=high_voxel_centers).write_landmarks(output_file+'_high_voxel_centers.LandmarkAscii')
        #
        # selected_pts_raw = self.get_pts_within_voxel_radius(self.pts,voxel_size/np.sqrt(2),high_voxel_centers_raw)
        # selected_pts = self.get_pts_within_voxel_radius(self.pts,voxel_size/np.sqrt(2),high_voxel_centers)
        #
        # self.write_landmarks(output_file+'_raw.LandmarkAscii',selected_pts_raw)
        # self.write_landmarks(output_file,selected_pts)

        # high_voxel_centers_raw = self.get_high_voxel_centers(binim)
        # high_voxel_centers = self.get_high_voxel_centers(connim)
        # Landmarks(pts=high_voxel_centers_raw).write_landmarks(output_file+'_high_voxel_centers_raw.LandmarkAscii')
        # Landmarks(pts=high_voxel_centers).write_landmarks(output_file+'_high_voxel_centers.LandmarkAscii')
        #
        # selected_pts_raw = self.get_pts_within_voxel_radius(self.pts,radius,high_voxel_centers_raw)
        # selected_pts = self.get_pts_within_voxel_radius(self.pts,radius,high_voxel_centers)
        # return cluster_in_bounds,max_dens,dens_cutoff_threshold,max_dens*dens_cutoff_threshold


    def get_rostral_most_cluster(self,output_file = None,voxel_size=500,dens_cutoff_threshold=0.15,write_density_images=False):

        #for exp_name in exp_names:
        #pts = Landmarks(output_path_landmarks+'{}_rabies_landmarks.LandmarkAscii'.format(exp_name)).pts
        pts_pdata = Surface(pts=self.pts)
        y_range = pts_pdata.surface.GetBounds()[3] - pts_pdata.surface.GetBounds()[2]

        cleaned_pts = []
        for pt in self.pts:
            if pt[1] < pts_pdata.surface.GetBounds()[2] + y_range/2:
                cleaned_pts.append(pt)

        cleaned_pts = self.pts
        pts_pdata = Surface(pts=cleaned_pts)

        # Compute density volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.25)
        densfilt.SetRadius(voxel_size*1.73/2)
        densfilt.SetDensityFormToNumberOfPoints()

        x = pts_pdata.surface.GetBounds()[1] - pts_pdata.surface.GetBounds()[0]
        y = pts_pdata.surface.GetBounds()[3] - pts_pdata.surface.GetBounds()[2]
        z = pts_pdata.surface.GetBounds()[5] - pts_pdata.surface.GetBounds()[4]

        densfilt.SetSampleDimensions(int(x/voxel_size)+1,int(y/voxel_size)+1,int(z/voxel_size)+1)
        densfilt.SetInputData(pts_pdata.surface)
        densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.Update()
        densimage = densfilt.GetOutput()

        if write_density_images:
            writer = vtk.vtkTIFFWriter()
            writer.SetFileName(output_file+'_density_image.tif')
            writer.SetInputData(densimage)
            writer.Update()

        # binarize the image based on density cut off (10% of max density)
        max_dens =  densimage.GetScalarRange()[1]
        thre = vtk.vtkImageThreshold()
        thre.SetInputData(densimage)
        thre.SetInValue(255)
        thre.SetOutValue(0)
        thre.ThresholdByUpper(max_dens*dens_cutoff_threshold)
        thre.Update()
        binim = thre.GetOutput()

        if write_density_images:
            writer = vtk.vtkTIFFWriter()
            writer.SetFileName(output_file+'_binary_image.tif')
            writer.SetInputData(binim)
            writer.Update()
            print('max density is {}'.format(max_dens))
            print('cutoff is {}'.format(dens_cutoff_threshold))
            print('threshold {}'.format(max_dens*dens_cutoff_threshold))

        valid_cluster = []
        for pt in cleaned_pts:
            x = (pt[0]-binim.GetOrigin()[0])/binim.GetSpacing()[0]
            y = (pt[1]-binim.GetOrigin()[1])/binim.GetSpacing()[1]
            z = (pt[2]-binim.GetOrigin()[2])/binim.GetSpacing()[2]
            ind = list(map(int,[x,y,z]))

            if binim.GetScalarComponentAsFloat(ind[0],ind[1],ind[2],0)  >0 : #== M1_region_id+1:
                valid_cluster.append(pt)

        #Landmarks(pts=valid_cluster).write_landmarks\
        #(output_path_landmarks+'{}_{}_{}_{}.landmarksAscii'.format(exp_name,voxel_size,dens_cutoff_threshold,max_dens))

        ymin = 9999
        ymin_pt = []
        for pt in valid_cluster:
            if pt[1] < ymin:
                ymin_pt = pt
                ymin = pt[1]

        radius = 4000
        final_cluster = []
        for pt in valid_cluster:
            dist = distance.euclidean(pt,ymin_pt)
            if dist < radius:
                final_cluster.append(pt)

        cluster_in_bounds = self.get_landmarks_in_bounding_box(Surface(pts=final_cluster).surface.GetBounds(),cleaned_pts)
        if output_file is not None:
            self.write_landmarks(output_file,cluster_in_bounds)

        return cluster_in_bounds,max_dens,dens_cutoff_threshold,max_dens*dens_cutoff_threshold

    def get_cluster_closest_to_given_pt(self,output_file = None,voxel_size=500,dens_cutoff_threshold=0.15,write_density_images=False,cluster_center=[0,0,0],distance_from_center=1000):

        #for exp_name in exp_names:
        #pts = Landmarks(output_path_landmarks+'{}_rabies_landmarks.LandmarkAscii'.format(exp_name)).pts
        pts_pdata = Surface(pts=self.pts)
        y_range = pts_pdata.surface.GetBounds()[3] - pts_pdata.surface.GetBounds()[2]

        cleaned_pts = self.pts

        pts_pdata = Surface(pts=cleaned_pts)

        # Compute density volume image
        densfilt = vtk.vtkPointDensityFilter()
        #densfilt.SetDensityEstimateToRelativeRadius()
        densfilt.SetDensityEstimateToFixedRadius()
        #densfilt.SetRelativeRadius(0.25)
        densfilt.SetRadius(voxel_size*1.73/2)
        densfilt.SetDensityFormToNumberOfPoints()

        x = pts_pdata.surface.GetBounds()[1] - pts_pdata.surface.GetBounds()[0]
        y = pts_pdata.surface.GetBounds()[3] - pts_pdata.surface.GetBounds()[2]
        z = pts_pdata.surface.GetBounds()[5] - pts_pdata.surface.GetBounds()[4]

        densfilt.SetSampleDimensions(int(x/voxel_size)+1,int(y/voxel_size)+1,int(z/voxel_size)+1)
        densfilt.SetInputData(pts_pdata.surface)
        densfilt.SetModelBounds(pts_pdata.surface.GetBounds())
        densfilt.Update()
        densimage = densfilt.GetOutput()

        if write_density_images:
            writer = vtk.vtkTIFFWriter()
            writer.SetFileName(output_file+'_density_image.tif')
            writer.SetInputData(densimage)
            writer.Update()

        # binarize the image based on density cut off (10% of max density)
        max_dens =  densimage.GetScalarRange()[1]
        thre = vtk.vtkImageThreshold()
        thre.SetInputData(densimage)
        thre.SetInValue(255)
        thre.SetOutValue(0)
        thre.ThresholdByUpper(max_dens*dens_cutoff_threshold)
        thre.Update()
        binim = thre.GetOutput()

        if write_density_images:
            writer = vtk.vtkTIFFWriter()
            writer.SetFileName(output_file+'_binary_image.tif')
            writer.SetInputData(binim)
            writer.Update()
            print('max density is {}'.format(max_dens))
            print('cutoff is {}'.format(dens_cutoff_threshold))
            print('threshold {}'.format(max_dens*dens_cutoff_threshold))


        valid_cluster = []
        for pt in cleaned_pts:
            x = (pt[0]-binim.GetOrigin()[0])/binim.GetSpacing()[0]
            y = (pt[1]-binim.GetOrigin()[1])/binim.GetSpacing()[1]
            z = (pt[2]-binim.GetOrigin()[2])/binim.GetSpacing()[2]
            ind = list(map(int,[x,y,z]))

            if binim.GetScalarComponentAsFloat(ind[0],ind[1],ind[2],0)  >0 : #== M1_region_id+1:
                valid_cluster.append(pt)

        #Landmarks(pts=valid_cluster).write_landmarks\
        #(output_path_landmarks+'{}_{}_{}_{}.landmarksAscii'.format(exp_name,voxel_size,dens_cutoff_threshold,max_dens))


        final_cluster = []
        for pt in valid_cluster:
            dist = distance.euclidean(pt,cluster_center)
            if dist < distance_from_center:
                final_cluster.append(pt)

        cluster_in_bounds = self.get_landmarks_in_bounding_box(Surface(pts=final_cluster).surface.GetBounds(),cleaned_pts)
        if output_file is not None:
            self.write_landmarks(output_file,cluster_in_bounds)

        return cluster_in_bounds,max_dens,dens_cutoff_threshold,max_dens*dens_cutoff_threshold

    def get_landmarks_within_given_surface_using_vtk(self,surface,pts_to_be_used = None,output_file = None,inside_out=True):

        if pts_to_be_used is not None:
            pts_being_used = pts_to_be_used
        else:
            pts_being_used = self.pts

        sel = vtk.vtkSelectEnclosedPoints()
        sel.SetInputData(Surface(pts=pts_being_used).convert_points_to_polydata(pts_being_used))
        sel.SetSurfaceData(surface)
        sel.SetInsideOut(inside_out)
        sel.Update()

        valid_pts = []
        arr =  sel.GetOutput().GetPointData().GetArray("SelectedPoints")
        for i in range((arr.GetNumberOfTuples())):
            if (arr.GetComponent(i,0)) == 1:
                valid_pts.append(pts_being_used[i])

        if output_file is not None:
            self.write_landmarks(output_file,pts_to_be_written=valid_pts)

        return valid_pts

    def get_landmarks_within_given_surface(self,surface,output_file=None):
        obbPia = vtk.vtkOBBTree()
        obbPia.SetDataSet(surface)
        obbPia.BuildLocator()
        inside_pts = []
        for pt in self.pts:
            points1 = vtk.vtkPoints()
            code = obbPia.IntersectWithLine(pt, Vectors().create_pt_along_vector_at_given_distance(10000,pt,[1,0,0]), points1, None)
            points2 = vtk.vtkPoints()
            code = obbPia.IntersectWithLine(pt, Vectors().create_pt_along_vector_at_given_distance(-10000,pt,[1,0,0]), points2, None)

            if points1.GetNumberOfPoints()>0 and points2.GetNumberOfPoints()>0:
                inside_pts.append(pt)
        #Landmarks(pts=inside_pts).write_landmarks(output_root+'inside_pts')
        return inside_pts

    def get_farthest_pt_to_axis(self,axis,pt_list,ignore_negative_direction=False):
        max_dist = 0
        min_edge = [[0,0,0],[1,1,1]]
        min_dist_from_wm = 99999
        min_dist_from_pia = 99999
        pt_on_verctor = []
        farthest_pt = []
        selected_pts = []
        for pt in pt_list:
            dist,dist_from_wm,dist_from_pia, wm_dist_vector = Vectors().get_distance_of_pt_from_line(pt,axis)
            if dist > max_dist:
                if ignore_negative_direction:
                    # allow only if the prijection pt lies on the positive side of the axis
                    # first extend the axis so that we dont fall short
                    cyl = Landmarks().get_cylinder(50,axis)
                    cylinder_srf = Surface(polydata=cyl).create_delunay_surface_3d(return_hull=True,make_cube=False)
                    valid_pts,invalid_pts = Surface().get_landmarks_within_given_surface(cylinder_srf,[wm_dist_vector])
                    if len(valid_pts) > 0:

                        selected_pts.append(valid_pts[0])
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

        return min_edge,max_dist, min_dist_from_wm, min_dist_from_pia, pt_on_vector, farthest_pt, selected_pts

    def get_volume_overlap_from_bounding_box(self,surf_list,landmarks_list):
        outer_cube,inner_cube,vol_ratio = Surface().get_intersection_cube(surf_list)
        aligned_landmarks = Landmarks()
        for landmarks in landmarks_list:
            aligned_landmarks.append_landmarks(landmarks.pts)
        #print(aligned_landmarks.pts)
        union_landmarks = aligned_landmarks.get_landmarks_within_given_surface_using_vtk(outer_cube)
        intersection_landmarks = aligned_landmarks.get_landmarks_within_given_surface_using_vtk(inner_cube)
        if len(union_landmarks)>0:
        #print(union_landmarks)
            return vol_ratio,len(intersection_landmarks)/len(union_landmarks)
        else:
            return vol_ratio,0

    def get_distance_between_centroids(self,data_list,):
        centers = []
        for l in data_list:
            centers.append(np.array(l.pts).mean(axis=0))

        centroid = np.array(centers).mean(axis=0)
        dists = []
        for center in centers:
            #print(center,centroid)
            dists.append(distance.euclidean(center,centroid))
        return dists
