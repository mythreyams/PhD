import sys
import os
sys.path.insert(0,os.path.abspath('..'))
from spatial_graphs.AmiraSpatialGraph import AmiraSpatialGraph,MatchBarrels
from spatial_graphs.Landmarks import Landmarks
from spatial_graphs.Surfaces import Surface
from spatial_graphs.Vectors import Vectors
from spatial_graphs.Alignment import Alignment
from dask import compute,multiprocessing,delayed
import pathlib
import shutil
import glob
import pandas as pd
import vtk
from scipy.spatial import distance
import numpy as np


input_path_spatial_graphs = '/nas1/Data_Mythreya/MotorCortexProject/V9/Registration_Local/Input_Spatial_Graphs/'
surface_resolution = 100
scaling_correction = False
neuron_translated = False


# create output folder structure
HOME = str(pathlib.Path(input_path_spatial_graphs).parent)
if scaling_correction:
    if neuron_translated:
        output_root = HOME + '/Outputs_{}_{}_{}/'.format(surface_resolution,'Scaling_Corrected','Neuron_Translated')
    else:
        output_root = HOME + '/Outputs_{}_{}/'.format(surface_resolution,'Scaling_Corrected')
else:
    output_root = HOME + '/Outputs_{}/'.format(surface_resolution)
#print(output_root)
pathlib.Path((output_root)).mkdir(exist_ok=True)

output_amira_spatial_graphs = output_root + '1_Amira_Spatial_Graphs/'
pathlib.Path((output_amira_spatial_graphs)).mkdir(exist_ok=True)

output_amira_section_spatial_graphs = output_root + '2_Section_Graphs/'
pathlib.Path((output_amira_section_spatial_graphs)).mkdir(exist_ok=True)

output_corrected_spatial_graphs = output_root + '3_Corrected_Spatial_Graphs/'
pathlib.Path((output_corrected_spatial_graphs)).mkdir(exist_ok=True)

output_surfaces = output_root + '4_Surfaces/'
pathlib.Path((output_surfaces)).mkdir(exist_ok=True)

output_neuron_axis = output_root + '5_Neuron_Axis/'
pathlib.Path((output_neuron_axis)).mkdir(exist_ok=True)

output_neuron_depth_stats = output_root + '6_Neuron_Stats/'
pathlib.Path((output_neuron_depth_stats)).mkdir(exist_ok=True)

output_gallery = output_root + '7_Gallery_Neurons/'
pathlib.Path((output_gallery)).mkdir(exist_ok=True)

output_features = output_root + '8_Features/'
pathlib.Path((output_features)).mkdir(exist_ok=True)


def main():

    convert_asc_to_amira()

    correct_spatial_graph()

    create_surfaces()

    generate_axis_field()

    select_neuron_axis()

    register_neuron_local_ref_frame()

    register_neuron_to_vM1_ref_frame()

    register_neuron_to_vS1_ref_frame()

    feature_extraction()








for file in glob.glob(input_path_spatial_graphs+'*'):
    #if os.path.basename(file) == 'RA20150807_2_Cell_B.ASC':
    if not os.path.isdir(file):
        sg = AmiraSpatialGraph(file,file_format='ASC')
        sg.write_spatial_graph(output_amira_spatial_graphs+os.path.basename(file)[:-3]+'am')


# # Correct spatial graph errors if any

# In[ ]:

# split into section graphs
for file in glob.glob(output_amira_spatial_graphs+'/*'):
    #if os.path.basename(file) == 'RA20150807_2_Cell_B.am':
    sg = AmiraSpatialGraph(file,create_section_graphs=True,compute_edge_length=True)
    #print(sg.graph_data.edge_list)
    animal_name = os.path.basename(file)
    pathlib.Path(output_amira_section_spatial_graphs+animal_name).mkdir(exist_ok=True)
    sg.write_all_section_spatial_graphs(output_amira_section_spatial_graphs+animal_name+'/')


# In[ ]:

section_graphs_path = output_amira_section_spatial_graphs
graph_path_3d = output_corrected_spatial_graphs
AXIS_DIRECTION = [1,1,1]
if scaling_correction:
    SCALING_CORRECTION = [2.306/2.35,2.309/2.35,1]
else:
    SCALING_CORRECTION = [1,1,1]
    
for section_folder in sorted(glob.glob(section_graphs_path+'*')):
    print(section_folder)
    sg = AmiraSpatialGraph(section_folder+'/1.am',read_header_only=True,generic_graph=True)
    sg_unscaled = AmiraSpatialGraph(section_folder+'/1.am',read_header_only=True,generic_graph=True)
    found_soma = False
    sections_after_soma = 0
    sections_after_wm = 0
    
    for z in range(len(glob.glob(section_folder+'/*.am'))):
        #print(section_folder+'/{}.am'.format(z+1))
        sg_sec = AmiraSpatialGraph(section_folder+'/{}.am'.format(z+1),axis_directions=AXIS_DIRECTION)
        #print(len(sg_sec.graph_data.edge_list))
        sg_sec_unscaled = AmiraSpatialGraph(section_folder+'/{}.am'.format(z+1),axis_directions=AXIS_DIRECTION,)
        #sg_sec.set_scaling(scaling = SCALING_CORRECTION)
        sg_sec.set_z_coord(z*-100,to_neuron=True)
        sg_sec_unscaled.set_z_coord(z*-100,to_neuron=True)
        
        sg_sec.set_scaling(scaling=SCALING_CORRECTION,to_neuron=True)
        sg_sec.graph_data = sg_sec.combine_subgraphs([sg_sec.pia,sg_sec.wm,                                                       sg_sec.bvs,sg_sec.barrels.all_rows_graphdata,                                                       sg_sec.neuron.all_neurites_subgraphdata])
        
        sg.graph_data = sg.combine_subgraphs([sg.graph_data,sg_sec.graph_data])
        sg_unscaled.graph_data = sg_unscaled.combine_subgraphs([sg_unscaled.graph_data,sg_sec_unscaled.graph_data])
        
#         if len(sg_sec.neuron.soma.vertices)>0:
#             found_soma = True
#         if len(sg_sec.wm.vertices)>0:
#             sections_after_wm = sections_after_wm + 1
#         if found_soma:
#             sections_after_soma = sections_after_soma + 1
        
#         if sections_after_soma > 2 and sections_after_wm > 2:
#             # lets have atleast 2 sections after soma and wm
#             break
    #print(sg.graph_data.edge_pt_coords)
    txmat = [1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1]#get_axis_flip_tx_mat(sg,need_to_flip=True)
    sg.apply_transformation(txmat)
    sg.write_spatial_graph(graph_path_3d+os.path.basename(section_folder))
    
    if scaling_correction and neuron_translated:
        sg_unscaled.apply_transformation(txmat)
        sg_unscaled.write_spatial_graph(graph_path_3d+'_unscaled.am')
        sg = AmiraSpatialGraph(graph_path_3d+os.path.basename(section_folder),)
        sg_unscaled = AmiraSpatialGraph(graph_path_3d+'_unscaled.am')
        os.remove(graph_path_3d+'_unscaled.am')

        # make sure that scaling does not affect soma position 
        # correct soma position is given by the scaled version
        # so tranlate the unscaled neuron to scaled soma position
        target_soma_position = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)
        original_soma_position = np.array(sg_unscaled.neuron.soma.edge_pt_coords).mean(axis=0)
        translation = original_soma_position - target_soma_position
        #print(translation)
        txmat = [1,0,0,0, 0,1,0,0, 0,0,1,0, -translation[0],-translation[1],-translation[2],1]
        sg_unscaled.apply_transformation_to_labelled_subgraph(sg_unscaled.neuron.all_neurites_subgraphdata,np.reshape(txmat,[4,4]))
        sg_unscaled.graph_data = sg_unscaled.combine_subgraphs([sg_unscaled.pia,sg_unscaled.wm,sg_unscaled.barrels.all_rows_graphdata,sg_unscaled.neuron.all_neurites_subgraphdata])
        sg.graph_data = sg.combine_subgraphs([sg.pia,sg.wm,sg.barrels.all_rows_graphdata,sg_unscaled.neuron.all_neurites_subgraphdata])
        sg.write_spatial_graph(graph_path_3d+os.path.basename(section_folder))

    
#     pia = Surface(pts=sg.pia.edge_pt_coords)
#     pia_hull = pia.create_delunay_surface_3d(return_hull=True)
#     #Surface(polydata=pia_hull).write_surface_mesh(graph_path_3d+'_scaled_pia.vtk')
#     pia_unscaled = Surface(pts=sg_unscaled.pia.edge_pt_coords)
#     pia_unscaled_hull = pia_unscaled.create_delunay_surface_3d(return_hull=True)
#     #Surface(polydata=pia_unscaled_hull).write_surface_mesh(graph_path_3d+'_unscaled_pia.vtk')
    
#     # get soma depth from pia  from unscaled version
#     soma_depth = 0
#     intersec_pt = [0,0,0]
#     soma_center = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)
#     #print(soma_center)
#     if len(sg_unscaled.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0:
#         apical_center = np.array(sg_unscaled.neuron.dendrite.apical_dendrite.edge_pt_coords).mean(axis=0)
#         intersec_pt_unscaled,soma_depth_unscaled = Surface(polydata=pia_unscaled_hull).get_vector_intersection_pt(soma_center,apical_center)
#         intersec_pt_scaled,soma_depth_scaled = Surface(polydata=pia_hull).get_vector_intersection_pt(soma_center,apical_center)
#         depth_change = (soma_depth_unscaled - soma_depth_scaled)
#         new_soma_center = pia.create_pt_along_vector_at_given_distance(-depth_change,soma_center,apical_center)
#     else:
#         # if there is no apical then use basal for direction or just leave as it is
#         #basal_center = np.array(sg_unscaled.neuron.dendrite.basal_dendrite.edge_pt_coords).mean(axis=0)
#         #intersec_pt_unscaled,soma_depth_unscaled = Surface(polydata=pia_unscaled_hull).get_vector_intersection_pt(soma_center,basal_center)
#         #intersec_pt_scaled,soma_depth_scaled = Surface(polydata=pia_hull).get_vector_intersection_pt(soma_center,basal_center)
#         #depth_change = (soma_depth_unscaled - soma_depth_scaled)
#         #new_soma_center = pia.create_pt_along_vector_at_given_distance(-depth_change,soma_center,basal_center)
#         new_soma_center = soma_center
#     #Landmarks(pts=[new_soma_center]).write_landmarks(graph_path_3d+'bla.landamarksAscii')
    
#     translation = soma_center - new_soma_center
#     txmat = [1,0,0,0, 0,1,0,0, 0,0,1,0, -translation[0],-translation[1],-translation[2],1]
    
#     sg.apply_transformation_to_labelled_subgraph(sg.neuron.all_neurites_subgraphdata,np.reshape(txmat,[4,4]))
#     sg.graph_data = sg.combine_subgraphs([sg.pia,sg.wm,sg.barrels.all_rows_graphdata,sg.neuron.all_neurites_subgraphdata])
#     sg.write_spatial_graph(graph_path_3d+os.path.basename(section_folder))
    
    


# # Create Pia and WM Surfaces

# In[ ]:

#TODO: fix divide surface... fails sometimes


# In[15]:

@delayed
def create_surface(file,surface_resolution='100'):
    sg = AmiraSpatialGraph(filename=file)
    # flip if required
    #txmat = get_axis_flip_tx_mat(sg)
    #print(txmat)
    #sg.apply_transformation(txmat)
    #sg.write_spatial_graph(root_dst+os.path.basename(file))
    #sg = AmiraSpatialGraph(root_dst+os.path.basename(file))
    #sg_clipped = get_fist_x_z_planes_sg(sg,17)
    #sg_clipped.write_spatial_graph(root_dst+os.path.basename(file)+'_clipped.am')
    os.system('/home/mythreya/project_src/BarrelField3D/DataAnalysis3D/BF3DRecon {} {} {}'.              format(file,output_surfaces+os.path.basename(file)[:-3],str(100)))
    pia = Surface(output_surfaces+os.path.basename(file)[:-3]+'_pia.vtk')
    wm = Surface(output_surfaces+os.path.basename(file)[:-3]+'_WM.vtk')
    
    if pia.surface.GetNumberOfCells()>0 and wm.surface.GetNumberOfCells()>0:
        pia.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_pia_bottom_open_100.vtk')
        wm.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_WM_bottom_open_100.vtk')
        
        os.system('/home/mythreya/project_src/BarrelField3D/DataAnalysis3D/BF3DRecon {} {} {}'.              format(file,output_surfaces+os.path.basename(file)[:-3],str(500)))
        pia_500 = Surface(output_surfaces+os.path.basename(file)[:-3]+'_pia.vtk')
        pia_500.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_pia_bottom_open_500.vtk')
        
        wm_500 = Surface(output_surfaces+os.path.basename(file)[:-3]+'_WM.vtk')
        wm_500_divided = Surface(polydata=wm_500.divide_surface(4,wm_500.surface))
        if wm_500_divided.surface.GetNumberOfCells() > 0:
            wm_500_divided.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_WM_bottom_open_500.vtk')
            flipped_axes = wm_500_divided.create_axis_field(pia,op_sg_name=None,flip_normals=True,return_axes=True,max_allowed_axis_length=4000)
            axes = wm_500_divided.create_axis_field(pia,op_sg_name=None,flip_normals=False,return_axes=True,max_allowed_axis_length=4000)
            sg = AmiraSpatialGraph(generic_graph=True)
            for edge in flipped_axes:
                sg.graph_data.add_edge(edge[0],edge[1])
            for edge in axes:
                sg.graph_data.add_edge(edge[0],edge[1])
            sg.write_spatial_graph(output_surfaces+os.path.basename(file)[:-3]+'_axis_field.am')
            
        else:
            # surface division failed
            flipped_axes = wm.create_axis_field(pia,op_sg_name=None,flip_normals=True,return_axes=True,max_allowed_axis_length=4000)
            axes = wm.create_axis_field(pia,op_sg_name=None,flip_normals=False,return_axes=True,max_allowed_axis_length=4000)
            sg = AmiraSpatialGraph(generic_graph=True)
            for edge in flipped_axes:
                sg.graph_data.add_edge(edge[0],edge[1])
            for edge in axes:
                sg.graph_data.add_edge(edge[0],edge[1])
            sg.write_spatial_graph(output_surfaces+os.path.basename(file)[:-3]+'_axis_field.am')
    else:
        # create a delunay surface as the traditional way has failed
        pia = Surface(pts=sg.pia.edge_pt_coords)
        wm = Surface(pts=sg.wm.edge_pt_coords)
        pia_hull = pia.create_delunay_surface_3d(return_hull=True,output_filename=output_surfaces+os.path.basename(file)[:-3]+'_pia.vtk')
        wm_hull = wm.create_delunay_surface_3d(return_hull=True,output_filename=output_surfaces+os.path.basename(file)[:-3]+'_WM.vtk')
        
        # clip at the bottom
        shutil.copyfile(src=file,dst=output_surfaces+os.path.basename(file)[:-3]+'_barrels.am')
        pia.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_pia_bottom_open.vtk')
        wm.clip_surface_at_given_z(-100,output_filename = output_surfaces+os.path.basename(file)[:-3]+'_WM_bottom_open.vtk')
        wm.create_axis_field(pia,op_sg_name=output_surfaces+os.path.basename(file)[:-3]+'_axis_field.am',flip_normals=False,                                      return_axes=True)


# In[16]:

fn_list_surfs = []
file_names = glob.glob(output_corrected_spatial_graphs+'*.am')
for i in range(len(file_names)):
    #if os.path.basename(file_names[i])=='RA20160323_Cell_D.am':
    fn_list_surfs.append(create_surface(file_names[i]))
    
compute(fn_list_surfs)


# In[ ]:

# see which ones failed
# file_names = glob.glob(output_corrected_spatial_graphs+'*.am')
# for i in range(len(file_names)):
#     file = file_names[i]
#     wm = Surface(output_surfaces+os.path.basename(file)[:-3]+'_WM_bottom_open_500.vtk')
#     if wm.surface.GetNumberOfCells() == 0:
#         print(os.path.basename(file)[:-3])
#         pia = Surface(output_surfaces+os.path.basename(file)[:-3]+'_pia_bottom_open_100.vtk')
#         wm = Surface(output_surfaces+os.path.basename(file)[:-3]+'_WM_bottom_open_100.vtk')
#         wm.create_axis_field(pia,op_sg_name=output_surfaces+os.path.basename(file)[:-3]+'_axis_field.am',flip_normals=False,\
#                                       return_axes=False,max_allowed_axis_length=4000)


# # Stats of input neuron

# In[17]:

def write_edge(edge,filename):
    axis_sg = AmiraSpatialGraph(generic_graph=True)
    axis_sg.graph_data.add_edge(edge[0],edge[1])
    axis_sg.write_spatial_graph(filename)

    return axis_sg
    


# In[18]:

# rotation from a to b
def align_a_to_b(a,b):
    u_a = Vectors().get_unit_vec(a[0],a[1])
    u_b = Vectors().get_unit_vec(b[0],b[1])
    v = np.cross(u_a,u_b)
    s = np.linalg.norm(v)
    c = np.dot(u_a,u_b)
    I = np.array([1,0,0,0,1,0,0,0,1]).reshape([3,3])
    #print(I)
    v_x = np.array([0,-v[2],v[1],v[2],0,-v[0],-v[1],v[0],0]).reshape([3,3])
    R = I + v_x + np.dot(v_x,v_x) * (1/(1+c))
    r_list = []
    for i in range(3):
        for j in range(3):
            r_list.append(R[i,j])
    #print(R)
    return R


# In[19]:

def convert_to_transformation_matrix(rot_mat,translation):
    
    tr_mat = []
    k = 0

    for i in range(3):
        for j in range(3):
            tr_mat.append(rot_mat[i,j])
            k = k + 1
        tr_mat.append(0)
        k = k+1
    tr_mat.append(-translation[0])
    tr_mat.append(-translation[1])
    tr_mat.append(-translation[2])
    tr_mat.append(1)
    return tr_mat


# In[20]:

def get_neuron_stats_and_locally_Register(exp_name):
    df = pd.DataFrame(columns=df_main.columns)
    
    print(i,exp_name)
    sg = AmiraSpatialGraph(output_surfaces+exp_name+'_barrels.am')
    pia_closed = Surface(output_surfaces+exp_name+'_pia.vtk')
    wm_closed = Surface(output_surfaces+exp_name+'_WM.vtk')
    pia = Surface(output_surfaces+exp_name+'_pia_bottom_open_100.vtk')
    wm = Surface(output_surfaces+exp_name+'_WM_bottom_open_100.vtk')

    #axis_field = wm.create_axis_field(pia,op_sg_name=output_surfaces+exp_name+'_axis_field.am',flip_normals=False,\
    #                                  return_axes=True)
    axis_field = AmiraSpatialGraph(output_surfaces+exp_name+'_axis_field.am',generic_graph=True).graph_data.edge_list
    df['Exp_Name'] = [exp_name]
    df['Pia_Area'] = [pia.get_surface_properties(prop='area')]
    df['WM_Area'] = [wm.get_surface_properties(prop='area')]
    df['Volume'] = [pia_closed.get_surface_properties(prop='volume')]
    vec_lens = []
    for edge in axis_field:
        vec_lens.append(Vectors().get_vec_length(edge))
    df['Cortical_Thickness_Mean'] = [np.array(vec_lens).mean(axis=0)]
    df['Cortical_Thickness_Std'] = [np.array(vec_lens).std(axis=0)]
    df['Cortical_Thickness_CoV'] = [100*np.array(vec_lens).std(axis=0)/np.array(vec_lens).mean(axis=0)]

    is_sticking_out = False#pia.is_neuron_sticking_out(sg)
    df['Neuron_Sticking_OutOf_Pia'] = [is_sticking_out]

    # find neuron axis, nearest axis and angle between them
    soma_center = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)
    apical_dend_center = np.array(sg.neuron.dendrite.apical_dendrite.edge_pt_coords).mean(axis=0)
    #print(soma_center,apical_dend_center)
    if len(sg.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0:
        print('Apical dend present in {}'.format(exp_name))
        df['Apical_Dendrite_Present'] = ['True']
        pia_proj_pt_neuron,soma_depth = pia.get_vector_intersection_pt(soma_center,apical_dend_center,extrapolation_len=10000)
        if len(pia_proj_pt_neuron) > 0:
            df['Pia_Intersection_Found'] = True
            neuron_axis_sg = write_edge([soma_center,pia_proj_pt_neuron],output_neuron_axis+exp_name+'_Soma_To_Pia_Axis.am')
            df['Neuron_Depth_from_Pia'] = [soma_depth]
            
        else:
            df['Pia_Intersection_Found'] = False

    else:
        df['Apical_Dendrite_Present'] = ['False']
        df['Pia_Intersection_Found'] = False
    #Landmarks(pts=[soma_center,apical_dend_center]).write_landmarks(output_surfaces+'a')
   
    nearest_axis,min_dist, min_dist_from_wm, min_dist_from_pia, pt_on_vector = Vectors().                    get_nearest_axis_to_pt(soma_center,axis_field,axis_validation_distance=4000)
    
    nearest_axis_sg = write_edge(nearest_axis,output_neuron_axis+exp_name+'_Nearest_Axis.am')
    df['Nearest_Axis_Depth_from_Pia'] = [min_dist_from_pia]
    df['Nearest_Axis_Height_from_WM'] = [min_dist_from_wm]
    df['Nearest_Axis_Length'] = [min_dist_from_pia+min_dist_from_wm]
    if len(sg.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0 and len(pia_proj_pt_neuron) > 0:
        angle = abs(Vectors().get_angle_between_vectors([soma_center,pia_proj_pt_neuron],nearest_axis))
        df['Angle_to_Nearest_Axis'] = [angle]

    if len(sg.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0 and len(pia_proj_pt_neuron)>0: 
        wm_proj_pt_neuron,soma_height = wm.get_vector_intersection_pt(soma_center,apical_dend_center,extrapolation_len=-10000)
        if len(wm_proj_pt_neuron) > 0 and soma_height < 3000:
            df['WM_Present_Below_Neuron'] = ['True']
            write_edge([wm_proj_pt_neuron,pia_proj_pt_neuron],output_neuron_axis+exp_name+'_WM_To_Pia_Axis.am')
            df['Neuron_Height_from_WM'] = [soma_height]
            df['Neuron_Axis_Length'] = [soma_depth+soma_height]

    # Locally transform the nuron on its axis
    if len(sg.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0 and len(pia_proj_pt_neuron) > 0:
        print(soma_center,pia_proj_pt_neuron,apical_dend_center)
        #z_angle = Vectors().get_angle_between_vectors([soma_center,apical_dend_center],[[0,0,0],[0,0,1]])
        rot_mat = align_a_to_b([soma_center,apical_dend_center],[[0,0,0],[0,0,1]])
        translation = pia_proj_pt_neuron
        neuron_axis_sg.apply_rotation(rot_mat,translation=translation)
        neuron_axis_sg.write_spatial_graph(output_gallery+exp_name+'_Rotated_Axis.am')
        #All_Axes.graph_data = All_Axes.combine_subgraphs([All_Axes.graph_data,neuron_axis_sg.graph_data])
    else:
        #z_angle = Vectors().get_angle_between_vectors(nearest_axis,[[0,0,0],[0,0,1]])
        rot_mat = align_a_to_b(nearest_axis,[[0,0,0],[0,0,1]])
        translation = nearest_axis[1]
        nearest_axis_sg.apply_rotation(rot_mat,translation=translation)
        nearest_axis_sg.write_spatial_graph(output_gallery+exp_name+'Rotated_Axis.am')
        #All_Axes.graph_data = All_Axes.combine_subgraphs([All_Axes.graph_data,neuron_axis_sg.graph_data])
    #tr_mat = convert_to_transformation_matrix(rot_mat,translation)
    #print(tr_mat)
    sg.apply_rotation(rot_mat,translation=translation)
    #sg.apply_transformation(tr_mat)
    sg.write_spatial_graph(output_gallery+exp_name+'.am')
    sg_neuron = AmiraSpatialGraph()
    sg_neuron.graph_data = sg_neuron.combine_subgraphs([sg.neuron.all_neurites_subgraphdata])
    sg_neuron.write_spatial_graph(output_gallery+exp_name+'_Neuron.am')
    
    return df,sg_neuron
    


# In[21]:

#@delayed
#def parallel_fun(i):
cols = ['Exp_Name','Pia_Area','WM_Area','Volume','Cortical_Thickness_Mean','Cortical_Thickness_Std','Cortical_Thickness_CoV',        'Neuron_Sticking_OutOf_Pia','Apical_Dendrite_Present','Pia_Intersection_Found','Neuron_Depth_from_Pia',        'Angle_to_Nearest_Axis','Nearest_Axis_Depth_from_Pia','Nearest_Axis_Height_from_WM','Nearest_Axis_Length',        'WM_Present_Below_Neuron','Neuron_Height_from_WM','Neuron_Axis_Length']
df_main = pd.DataFrame(columns=cols)

All_Neurons = AmiraSpatialGraph()
#All_Axes = AmiraSpatialGraph()

for i in range(len(sorted(glob.glob(output_surfaces+'*_pia_bottom_open_100.vtk')))):
    pia_file = sorted(glob.glob(output_surfaces+'*_pia_bottom_open_100.vtk'))[i]
    exp_name = os.path.basename(pia_file)[0:-24]
    #print(exp_name)
    #if exp_name == 'RA20180109_2_Cluster_A_Cell_C_done.am_Segmented_corrected':    
    df,sg_neuron = get_neuron_stats_and_locally_Register(exp_name)

    df_main = df_main.append(df)[df.columns.tolist()]
    All_Neurons.graph_data = All_Neurons.combine_subgraphs([All_Neurons.graph_data,sg_neuron.graph_data])
    
df_main.to_csv(output_root+'Input_Neuron_Stats.csv')
All_Neurons.write_spatial_graph(output_root+'All_Neurons_Rotated.am')
#All_Axes.write_spatial_graph(output_root+'All_Axes_Rotated.am')


# In[ ]:




# # Extract Features

# In[ ]:

def get_width_height_volume(pts):
    surf = Surface(pts = pts)
    bounds = surf.surface.GetBounds()
    width = bounds[1] - bounds[0]
    bredth = bounds[3] - bounds[2]
    height = bounds[5] - bounds[4]
    volume = width * bredth * height
    
    return width,bredth,height,volume
    


# In[ ]:

# Most neuron specific features are extracted using old code
# Here i only update position specific features. i.e. volume,width and height

# read existing feature csv
df = pd.read_csv(output_features+'Features.csv')
for i in df.index:
    filename = sorted(glob.glob(output_gallery+'*_Neuron.am'))[i]
    print(filename)
    sg = AmiraSpatialGraph(filename)
    total_dend_pts = []
    # Apical dend stats
    if len(sg.neuron.dendrite.apical_dendrite.edge_pt_coords) > 0:
        w,b,h,v = get_width_height_volume(sg.neuron.dendrite.apical_dendrite.edge_pt_coords)
        df.at[i,'Apical_Width'] = w
        df.at[i,'Apical_Bredth'] = b
        df.at[i,'Apical_Height'] = h
        df.at[i,'Apical_Volume'] = v
        for pt in sg.neuron.dendrite.apical_dendrite.edge_pt_coords:
            total_dend_pts.append(pt)
    
    # basal dend stats
    w,b,h,v = get_width_height_volume(sg.neuron.dendrite.basal_dendrite.edge_pt_coords)
    df.at[i,'Basal_Width'] = w
    df.at[i,'Basal_Bredth'] = b
    df.at[i,'Basal_Height'] = h
    df.at[i,'Basal_Volume'] = v
    for pt in sg.neuron.dendrite.basal_dendrite.edge_pt_coords:
        total_dend_pts.append(pt)
        
    # Total Dend stats
    w,b,h,v = get_width_height_volume(total_dend_pts)
    df.at[i,'Dendritic_Width'] = w
    df.at[i,'Dendritic_Bredth'] = b
    df.at[i,'Dendritic_Height'] = h
    df.at[i,'Dendritic_Volume'] = v
    
    if len(sg.neuron.axon.edge_pt_coords) > 0:
        w,b,h,v = get_width_height_volume(sg.neuron.axon.edge_pt_coords)
        df.at[i,'Axon_Volume'] = v
    
    #df['Soma_Depth_From_Pia'] = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)[2]
    


# In[ ]:

df


# In[ ]:

df_soma_depth = pd.DataFrame(columns=['Soma_Depth_From_Pia'])
for i in df.index:
    filename = sorted(glob.glob(output_gallery+'*_Neuron.am'))[i]
    print(filename)
    sg = AmiraSpatialGraph(filename)
    df_soma_depth.at[i,'Soma_Depth_From_Pia'] = np.array(sg.neuron.soma.edge_pt_coords).mean(axis=0)[2]
    
    


# In[ ]:

df_final = pd.DataFrame()
df_final = df_final.append(df,ignore_index=True)


# In[ ]:

df_final['Soma_Depth_From_Pia'] = df_soma_depth['Soma_Depth_From_Pia']


# In[ ]:

df_final.to_csv(output_features+'Features_With_Soma_Depth_V2.csv')


# In[ ]:



