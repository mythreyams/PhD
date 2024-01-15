
# coding: utf-8

# In[1]:

import os
import glob
import numpy as np
import pathlib
import SimpleITK as sitk
import scipy.spatial.distance as distance
import shutil


# In[6]:

XY_RES = 0.868
Z_RES = 50

exp_name = 'MG49_lhs/'
input_image_path =  '/nas1/Data_Mythreya/MotorCortexProject/Images_For_NeuN_Count/'#.format(exp_name)
#input_contours_path = '/nas1/Data_Mythreya/MotorCortexProject/V0/0_Inputs/Contours/Rabies/' + exp_name
final_path = '/nas1/Data_Mythreya/MotorCortexProject/NeuN_Counts/{}/'.format(exp_name)
pathlib.Path(final_path).mkdir(exist_ok= True)

output_path = '/usr/tmp_cell_count/'
pathlib.Path(output_path).mkdir(exist_ok= True)
sec_num_start = 18
sec_num_end = 18
#sec_num = 45


# In[91]:

#@delayed
def count_cells(i):
    #print(i)
#     if os.path.exists(final_path+'S_{}/landmarks.landmarkAscii'.format(i)):
#         landmarks = getLandmarkCoords(final_path+'S_{}/landmarks.landmarkAscii'.format(i))
#         if len(landmarks) > 0:
#             # count already exists... so skip this section
#             return
#     print(i)
    
    try:
        pathlib.Path(output_path+'S_{}/'.format(i)).mkdir(exist_ok= True)
        pathlib.Path(final_path+'S_{}/'.format(i)).mkdir(exist_ok= True)
        
        # mask the image so that only cortex is used for counting
        filename = glob.glob(input_image_path+'S{}*.tif'.format(i))
        im = sitk.ReadImage(filename[0])
#         spatial_graph_file = input_contours_path + 'S{:03d}.am'.format(i)
#         sg = AmiraSpatialGraph(spatial_graph_file)
#         wm_surf = Surface(pts=sg.wm.edge_pt_coords)
#         bounds = wm_surf.surface.GetBounds()
#         print(len(sg.wm.edge_pt_coords))
#         #pia,wm = getSectionPiaWMPoints(spatial_graph_file)
#         pia_bin = convertContourCoordsIntoBinaryImage(sg.pia.edge_pt_coords,im,resample_res=0.1)
#         wm_bin = convertContourCoordsIntoBinaryImage(sg.wm.edge_pt_coords,im,resample_res=0.1)

        
#         pia_only = sitk.Mask(im,pia_bin)
#         cortex_im = sitk.Mask(pia_only,sitk.Not(wm_bin))
        sitk.WriteImage(im,output_path+'S_{}/S_{}.tif'.format(i,i))
        with open(output_path+'S_{}'.format(i)+'/arguments_{}.txt'.format(i),'w') as f:
            f.write('<command> CMD_COUNT \n')
            f.write('<ifilename> {}S_{}/S_{}.tif\n'.format(output_path,i,i))
            f.write('<ofilename>    {}\n'.format(output_path+'S_{}/'.format(i)))
            f.write('<cfilename>    {}/CellClusterList.csv\n'.format(output_path+'section_{}/'.format(i)))
            f.write('<scantype> 1\n')
            f.write('<xysampling> 0.868\n')
            f.write('<zsampling> 1\n')
            f.write('<bricksize>  300 m\n')
            f.write('<somaradius> 4\n')
            f.close()
        
        print('/home/mythreya/projects/CellCount/source_code/bugfix/Version20120306_ColocV2_without_clustering/Version20120306_ColocV2_Copy//CellCount         {} > {}'.format(output_path+'S_{}'.format(i)+'/arguments_{}.txt'.format(i),output_path+'S_{}'.format(i)+'/screen_outputs{}.txt'.format(i)))
        os.system('/home/mythreya/projects/CellCount/source_code/bugfix/Version20120306_ColocV2_without_clustering/Version20120306_ColocV2_Copy/CellCount         {} > {} 2>&1'.format(output_path+'S_{}'.format(i)+'/arguments_{}.txt'.format(i),output_path+'S_{}'.format(i)+'/screen_outputs{}.txt'.format(i)))
    
        # copy the final landmarks results into the final result folder
        shutil.copyfile(output_path+'S_{}/screen_outputs{}.txt'.format(i,i), final_path+'S_{}/screen_outputs{}.txt'.format(i,i))
        shutil.copyfile(output_path+'S_{}/landmarks.landmarkAscii'.format(i), final_path+'S_{}/landmarks.landmarkAscii'.format(i))
        shutil.copyfile(output_path+'S_{}/_probs.csv'.format(i), final_path+'S_{}/_probs.csv'.format(i))
        shutil.copyfile(output_path+'S_{}/_histo.csv'.format(i), final_path+'S_{}/_histo.csv'.format(i))
        
        shutil.rmtree(output_path+'S_{}'.format(i))
        
                  
    except RuntimeError:
        print('Image file not found')


# In[108]:

fun_list = []
for i in range(sec_num_start,sec_num_end+1):
    count_cells(i)
    break
#     fun_list.append(count_cells(i))
#     if (len(fun_list)%5)==0:
#         compute(fun_list)
#         fun_list = []


# In[ ]:



