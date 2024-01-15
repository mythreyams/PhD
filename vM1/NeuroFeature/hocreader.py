# FEATURE EXTRACTION
# v2.11 (12/09/2017)
# written by Andrew Johnson
# last modified by Daniel Udvary

import csv
import math
import os
import re
import sys
from datetime import datetime

import functions as fcn
import hocclasses as hc

startTime = datetime.now()
Experiment_ID = None
mtime = os.path.getmtime(sys.argv[0])
last_modified_date = datetime.fromtimestamp(mtime)
multiple_inputs = False
local_orientation = ["None"]

print "****************************************************************************"
print "Using Python version ", sys.version
print "Running python script: ", sys.argv[0]
print "Version: 2.11 Last modified: ", last_modified_date
print "****************************************************************************"
print "INPUTS: "
print "length of input string = ", len(sys.argv)
print "current working directory: ", os.getcwd()

if len(sys.argv) == 4:
    input_hoc_file = str(sys.argv[1])  # Use: python hocreader.py input.hoc -0.0324163,-0.0472919,0.998355 [e.g.]
    # input local orientation at command line
    print "input_hoc_file (4 args) =", input_hoc_file
    upload_path = str(sys.argv[2])
    print "upload_path =", upload_path

    # if transfomration matrix as input u,v,w
    if str(sys.argv[3]).find(',') != -1:
        local_orientation = str(sys.argv[3])
        local_orientation = local_orientation.split(",")  # convert string to list
        print "local orientation = ", local_orientation
    # if path/to/transformation.log as input
    elif str(sys.argv[3])[-4:] == ".log":
        logfile = str(sys.argv[3])
        inlogfile = open(logfile, "r")
        for line in inlogfile:
            if re.search("Local orientation:", line):
                parts = line.split(":")
                parts1 = parts[1]  # Numbers with brackets [1,2,3]
                parts2 = parts1.split("[")[1]  # 1,2,3]
                local_orientation = parts2.split("]")[0]  # 1,2,3
                local_orientation = local_orientation.split(",")  # convert string to list
                break
        print "local orientation =", local_orientation, "extracted from", logfile
    else:
        local_orientation = '0,0,1'  # default local_orientation (none)
        print "WARNING! Could not read second input ", str(sys.argv[3]), "! local orientation set to 0,0,1"
        print "USAGE: python hocreader.py [input.hoc] [upload_path] [transformationMatrix: u,v,w]"
        print "USAGE: python hocreader.py [input.hoc] [upload_path] [transform.log]"
    FILES = [input_hoc_file]
elif len(sys.argv) == 3:
    input_hoc_file = str(sys.argv[1])  # Use: python hocreader.py input.hoc
    upload_path = str(sys.argv[2])
    findlogs = False
    print "input_hoc_file (2 args) =", input_hoc_file
    print "upload_path =", upload_path
    print "local orientation = ", local_orientation
    FILES = [input_hoc_file]
elif len(sys.argv) == 1:
    multiple_inputs = True
    input_hoc_file = fcn.inputHocs()  # Use: hocreader.py (gets input from directory path in function.py)
    # Set findlogs = True to get local orientation from .log files.
    #findlogs = True  # (Will read local orientation from .log files in same path.)
    findlogs = False                     # Otherwise assume no local orientation and no log files present
    print "input_hoc_file (1 arg.) =", input_hoc_file
    FILES = input_hoc_file
    upload_path = "/nas1/Data_Mythreya/MotorCortexProject/V8/Registration_Local/Outputs_200/8_Features/"
    
else:
    print "USAGE: python hocreader.py [input.hoc] [upload_path] [transformationMatrix: u,v,w]"
    print "USAGE: python hocreader.py [input.hoc] [upload_path] [transform.log]"
    print "USAGE: python hocreader.py [input.hoc] [upload_path]"
    sys.exit()

print "****************************************************************************"

# Flags
verbose = 90  # -1   # print voluminously, or not

# Check for local orientation:
if (local_orientation[0] == "None"):
    print "No local orientation detected!"
    barreltrans = False
else:
    # Need to rotate to local vertical.
    # Good reference is: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
    # Used matrices Txz and Tz from there to get rotation matrix below (Tz*Txz)
    #
    # Write rotation matrix in terms of rotation vector (u,v,w):
    #
    #            x_barrel = (u*w*x + w*v*y - uv_squared*z)/sq_root
    #            y_barrel = (-v*x + u*y)/sq_root
    #            z_barrel = u*x + v*y + w*z
    #
    # where (x,y,z) are the "global" coords from the hoc file and
    # (x_barrel, etc) are the rotated coords to the barrel axis (along z_barrel)

    # local orientation vector:
    u = float(local_orientation[0])
    v = float(local_orientation[1])
    w = float(local_orientation[2])
    angle = math.acos(w)

    print "u, v, w, angle = ", u, v, w, angle

    if (u == 0 and v == 0 and angle == 0):
        barreltrans = False
    else:
        barreltrans = True

        # Generate "minimal torsion" rotation matrix based on local orientation
        # rotation matrix = [[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]]
        m = fcn.rotate(angle, u, -v, 0)
        print "Rotation Matrix = ", m

        # Useful quantities for the (old) transformation
        # Nor used anymore
        uv_squared = u ** 2 + v ** 2
        sq_root = math.sqrt(uv_squared)

        print "Local orientation detected & performed"

# Loop over files in FILES list of multiple hoc files:
file_num = 0
for file in FILES:
    file_num += 1
    print "****************************************************************************"
    print "Processing file =", file
    print "File Number =", file_num
    print

    # Check for multiple_inputs (hardwire = true) option, read local_orientation from .log file
    if (multiple_inputs):

        if (findlogs):
            # Find .log file corresponding to the .hoc file:

            # Gerardo's files
            # 500_GP_WR639_cell_1547_SP5C_standardized_registered_D2.hoc
            # 500_GP_WR639_cell_1547_SP5C_standardized_newHB_D2_transform.log
            part1 = file.split("registered")[0]  # 500_GP_WR639_cell_1547_SP5C_standardized_
            logfile = part1 + "newHB_D2_transform.log"  # 500_GP_WR639_cell_1547_SP5C_standardized_newHB_D2_transform.log

            print
            print " ****************************************************************************"
            print " Attempting to find .log file for local orientation extraction: "
            print "     ", logfile
            print " ****************************************************************************"
            print

            # Look for log files to run with individualized local orientations
            # if (open(logfile,"r")):
            try:
                inlogfile = open(logfile, "r")
                barreltrans = True

                print
                print " ****************************************************************************"
                print " Found .log file for local orientation extraction: "
                print "     ", logfile
                print " ****************************************************************************"
                print

                for line in inlogfile:
                    if re.search("Local orientation:", line):
                        parts = line.split(":")
                        parts1 = parts[1];  # Numbers with brackets [1,2,3]
                        parts2 = parts1.split("[")[1];  # 1,2,3]
                        local_orientation = parts2.split("]")[0]  # 1,2,3
                        local_orientation = local_orientation.split(",")  # convert string to list

                        print
                        print "Found Local Orientation Vector from .log file: ", local_orientation
                        print

            # No .log file found, proceed with no local orientation
            except IOError:
                # barreltrans = False # No need to set False, already False
                print
                print " ****************************************************************************"
                print " No .log file found "
                print
                print " ****************************************************************************"
                print

                pass

        if (barreltrans):

            # local orientation vector:
            u = float(local_orientation[0])
            v = float(local_orientation[1])
            w = float(local_orientation[2])
            angle = math.acos(w)

            # Test for case with local orientation given but (0,0,z)
            # (In which case run with no local orientation.)
            if (u == 0 and v == 0):
                barreltrans = False
            else:
                barreltrans = True

                # Generate "minimal torsion" rotation matrix based on local orientation
                # rotation matrix = [[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]]
                m = fcn.rotate(angle, u, -v, 0)

                # Useful quantities for the (old) transformation
                # Not used anymore
                uv_squared = u ** 2 + v ** 2
                sq_root = math.sqrt(uv_squared)

                print "***************** Local orientation found  *****************"
                print "*****************  Performing rotation     *****************"
                print
                # print "*****************uv_squared = ", uv_squared
                # print
                # print "*****************sq_root = ", sq_root
                # print

        # No log files, run with no orintation
        else:
            print
            print "****************************************************************************"
            print
            print "Running with no orientation!!!!!!!!!!!!!!!!!!!!!!!! "
            print
            print "****************************************************************************"
            print

    # Instantiate Neuron as "cell"
    cell = hc.Neuron()
    cell.name = file

    # Initialize Soma
    somaX = 0
    somaY = 0
    somaZ = 0
    soma_x = 0
    soma_y = 0
    soma_z = 0
    ave_somaX = 0
    ave_somaY = 0
    ave_somaZ = 0
    numsomapts = 0
    transform = False
    RuleTscale = False

    # Other initializations
    z_scaling_lookup = {}
    z_scale = False

    m00 = 0
    m01 = 0
    m02 = 0
    m03 = 0

    m10 = 0
    m11 = 0
    m12 = 0
    m13 = 0

    m20 = 0
    m21 = 0
    m22 = 0
    m23 = 0

    m30 = 0
    m31 = 0
    m32 = 0
    m33 = 0

    # Declarations
    hermannhoc = False

    infile = open(file, "r")
    if (verbose == 0):
        print "***********Begin initial loop through file ", file

    # hocfile = infile.readlines()
    numlines = 0

    # Initial run through whole file to determine
    #       * What kind of hoc format
    #       * If there is registration ???
    #       * Transformation matrix ???
    #       * z-scaling lookup ???
    for line in infile:
        numlines += 1
        # theline = line.split()
        # print "line = ",line
        # n = len(line)
        if re.search("CellBuilder-like", line):
            # if re.search("CellBuilder-like", hocfile):
            hermannhoc = True
            if (verbose == 1):
                print
                print "************************************************"
                print "One of Hermann Cuntz's Clone hoc files:"
                print "line = ", line
                print "Setting hermannhoc =", hermannhoc
                print "************************************************"
                print

        # Search for registration factor, TScale
        if re.search("pt3dchange", line):

            RuleTscale = True
            parts = line.split(",")

            p01 = parts[1]
            pI = p01.split("*")
            # print " pI = ", pI

            # For files like ID 65 with different transformation matrix formats
            try:
                TScale = p01 = pI[1]
                TScale = float(TScale)
            except IndexError:
                TScale = None
                RuleTscale = False

        if re.search("transformation_matrix.x", line):

            if (verbose == 0):
                print
                print " Processing Transformation Matrix ..."
                print "                            ... Setting transform = True "
                print

            transform = True

            parts = line.split("=")
            label = parts[0]
            # For files like ID 65 with different transformation matrix formats
            try:
                matrix = float(parts[1])
            except ValueError:
                matrix = None

            # Fuck it (I'm an idiot): Brute force:
            if re.search("\[0\]\[0\]", label):
                m00 = matrix
            if re.search("\[0\]\[1\]", label):
                m01 = matrix
            if re.search("\[0\]\[2\]", label):
                m02 = matrix
                # print "++ m00 = ", m02
            if re.search("\[0\]\[3\]", label):
                m03 = matrix
                # print "++ m03 = ", m03
            if re.search("\[1\]\[0\]", label):
                m10 = matrix
                # print "++ m10 = ", m10
            if re.search("\[1\]\[1\]", label):
                m11 = matrix
                # print "++ m11 = ", m11
            if re.search("\[1\]\[2\]", label):
                m12 = matrix
                # print "++ m12 = ", m12
            if re.search("\[1\]\[3\]", label):
                m13 = matrix
                # print "++ m13 = ", m13
            if re.search("\[2\]\[0\]", label):
                m20 = matrix
                # print "++ m20 = ", m20
            if re.search("\[2\]\[1\]", label):
                m21 = matrix
                # print "++ m21 = ", m21
            if re.search("\[2\]\[2\]", label):
                m22 = matrix
                # print "++ m22 = ", m22
            if re.search("\[2\]\[3\]", label):
                m23 = matrix
                # print "++ m23 = ", m23
            if re.search("\[3\]\[0\]", label):
                m30 = matrix
                # print "++ m30 = ", m30
            if re.search("\[3\]\[1\]", label):
                m31 = matrix
                # print "++ m31 = ", m31
            if re.search("\[3\]\[2\]", label):
                m32 = matrix
                # print "++ m32 = ", m32
            if re.search("\[3\]\[3\]", label):
                m33 = matrix
                # print "++ m33 = ", m33

            if (verbose == 0):
                print "++ m00 = ", m00
                print "++ m01 = ", m01
                print "++ m00 = ", m02
                print "++ m03 = ", m03
                print "++ m10 = ", m10
                print "++ m11 = ", m11
                print "++ m12 = ", m12
                print "++ m13 = ", m13
                print "++ m20 = ", m20
                print "++ m21 = ", m21
                print "++ m22 = ", m22
                print "++ m23 = ", m23
                print "++ m30 = ", m30
                print "++ m31 = ", m31
                print "++ m32 = ", m32
                print "++ m33 = ", m33

        if re.search("z_scaling_lookup.x", line):
            z_scale = True
            part = line.split("=")
            part0 = part[0].split("[")
            part1 = part0[1]
            zDepth = part1.split("]")[0]
            zScalingVal = part[1]
            z_scaling_lookup[zDepth] = zScalingVal

        # Compute average of soma
        if re.search("create", line):
            create_subject = line

        if (re.search("pt3dadd", line) and re.search("soma", create_subject)):

            # Note: Neuroconv exchanges z and y!
            parts = line.split(',')
            # coords in different order for "marcel" or "3D" hoc files:
            # (Obsolete, for invitro data):

            x = float(parts[0].split('(')[1])
            y = float(parts[1])
            z = float(parts[2])
            r = float(parts[3].split(')')[0])

            # transform if local orientation given:
            if (barreltrans):
                rotated_coords = fcn.matrix_multiply(x, y, z, m)

                x_barrel = rotated_coords[0]
                y_barrel = rotated_coords[1]
                z_barrel = rotated_coords[2]

                x = x_barrel
                y = y_barrel
                z = z_barrel

            soma_x = soma_x + x
            soma_y = soma_y + y
            soma_z = soma_z + z

            numsomapts += 1

            ave_somaX = soma_x / numsomapts
            ave_somaY = soma_y / numsomapts
            ave_somaZ = soma_z / numsomapts

    if (verbose == 0):
        print "ave_somaX =", ave_somaX
        print "ave_somaY =", ave_somaY
        print "ave_somaZ =", ave_somaZ

    # No z-scale matrix in hoc file, compute using zLUT table
    if (z_scale == False):
        z_scaling_lookup = fcn.zLUT(z_scaling_lookup)

    ##########################################################################################
    # End for line in infile
    ##########################################################################################

    ###########################################################################################
    # Main loop through file (Separate loops for Neuroconv and Cuntz clone hoc formats        #
    ###########################################################################################

    # set true below of "apical" found as label in hoc file
    basal = False
    apical = False
    axon = False

    infile = open(file, "r")
    if (verbose == 0):
        print
        print
        print
        print "***********Beginning initial loop through file ", file
        print
        print
        print

    seg_count = 0
    pt_count = 0

    #################################################################################
    # Neuroconv hoc files or Non-Hermann-Clone Files...
    #################################################################################
    if (not hermannhoc):
        if (verbose == 0):
            print
            print "************************************************"
            print "Neuroconv or at least non-Hermann hoc file:"
            print "hermannhoc (better be false!) =", hermannhoc
            print "************************************************"
            print

        # import neuroconv # (Hypothetical routine that does all this for non-hermannhoc files)
        axon_bp = 0
        axon_length = 0
        x_axon_all = []
        y_axon_all = []
        z_axon_all = []
        axon_fathers = []
        axon_sons = []

        for line in infile:
            # Search for Soma
            if re.search("Soma with centroid", line):

                parts = line.split(" ")
                if (verbose == 0):
                    print line
                    print "0:", parts[0]
                    print "1:", parts[1]
                    print "2:", parts[2]
                    print "3:", parts[3]
                    print "4:", parts[4]
                    print "5:", parts[5]
                    print "6:", parts[6]
                    print "7:", parts[7]
                    print "8:", parts[8]
                somaX = float(parts[5])
                somaY = float(parts[6])
                somaZ = float(parts[7])

                # print " 3D SOMA x:", somaX ,"y:", somaY ,"z:", somaZ

                # if (transform):
                new_x = m00 * somaX + m01 * somaY + m02 * somaZ + m03
                new_y = m10 * somaX + m11 * somaY + m12 * somaZ + m13
                new_z = m20 * somaX + m21 * somaY + m22 * somaZ + m23

                z_scaling_lookup_offset = -10000
                z_scaling_lookup_delta = 1.000000e-001
                # Also see [***]
                # Needs to be a goddamn string for the dictionary look-up
                # below (*) to not have a missing key error---
                # Don't fucking ask me why.
                new_z_lookupPos = str(int((new_z - z_scaling_lookup_offset) / z_scaling_lookup_delta))

                if (transform):
                    # print "new_z_lookupPos, type = ", new_z_lookupPos, type(new_z_lookupPos)
                    # new_z = z_scaling_lookup.get(new_z_lookupPos, new_z)
                    new_z = float(z_scaling_lookup[new_z_lookupPos])  # (*) Refloat here
                    somaX = new_x  # (I know, idiotic.)
                    somaY = new_y
                    somaZ = new_z
                else:
                    new_z = somaZ

                if (verbose == 0):
                    print
                    print "2. new_x, y, and z = ", new_x, new_y, new_z
                    print

            # New Segment
            if (re.search("create ", line)):

                # From previous segment, record length, print:
                if (pt_count > 0):
                    # print "segment length = ",leg.length()
                    # print "segment surf area = ",leg.length()[1]
                    segmentLength = leg.length()

                    # Calculate axon length on fly:
                    if leg.label == "axon":
                        axon_length = axon_length + segmentLength[0]

                # New Segment:
                secondpart = line.split()[1]
                leg = hc.Segment()
                leg.segmentName = secondpart.split('}')[0]

                # Find and set segment label:
                if re.search("soma", leg.segmentName): leg.label = 'soma'
                if re.search("dend", leg.segmentName): leg.label = 'basal'; basal = True
                if re.search("apical", leg.segmentName): leg.label = 'apical'; apical = True
                if re.search("axon", leg.segmentName): leg.label = 'axon'; axon = True

                if (seg_count == 0):

                    # First "create" line, start segment, sons,
                    # and fathers lists and points dictionary:
                    cell.segments = [leg.segmentName]
                    cell.sons = []
                    cell.father = []
                    cell.points = {}
                    cell.segmentLengths = []
                    cell.label = [leg.label]
                    cell.numpoints = []

                else:

                    # Subsequent "create" lines, append segments
                    # (sons and fathers done in "connect" line below):
                    cell.segments.append(leg.segmentName)
                    cell.segmentLengths.append(segmentLength)
                    cell.label.append(leg.label)
                    cell.numpoints.append(pt_count)

                # reset number of points in segment (tallied in pt3dadd search)
                pt_count = 0

                # Increment number of segments
                seg_count += 1

                # Misc prints
                if (verbose == 2):
                    # if (leg.label == "axon"):
                    # print
                    # print "**************found create line : ", line
                    # print "**************line.split()      = ", line.split()
                    print "**************leg.segmentName   = ", seg_count, ": ", leg.segmentName
                    # print "**************create_subject    = ", create_subject
                    # print

            # Find Fathers and Sons
            if re.search("connect", line):
                parts = line.split(',')
                firstpart = parts[0]
                secondpart = firstpart.split()[1]

                # Segment already created, but not necessarily this one.
                leg.father = parts[1].split('(')[0]

                cell.father.append(leg.father)
                cell.sons.append(leg.segmentName)

                # print "son of ", leg.father, " is ", leg.segmentName
                if re.search("axon", leg.father):
                    axon_fathers.append(leg.father)

            # points in segment
            if re.search("pt3dadd", line):

                # Note: Neuroconv exchanges z and y!
                parts = line.split(',')

                """
                # In vivo data
                # coords in different order for "marcel" or "3D" hoc files:
                #if re.search("marcel" or "_registered_", file):
                if re.search("registered", file):
                    # print "Found marcel or registered file!", file
                    x = float(parts[0].split('(')[1])
                    y = float(parts[1])
                    z = float(parts[2])
                    r = float(parts[3].split(')')[0])

                # Slice (In vitro)
                else:
                    # print "Found nothing!", file
                    x = float(parts[0].split('(')[1])
                    z = float(parts[1])
                    y = float(parts[2])
                    r = float(parts[3].split(')')[0])
                """
                x = float(parts[0].split('(')[1])
                y = float(parts[1])
                z = float(parts[2])
                r = float(parts[3].split(')')[0])

                # transform if local orientation given:
                if (barreltrans):
                    rotated_coords = fcn.matrix_multiply(x, y, z, m)

                    x_barrel = rotated_coords[0]
                    y_barrel = rotated_coords[1]
                    z_barrel = rotated_coords[2]

                    x = x_barrel
                    y = y_barrel
                    z = z_barrel

                # if (transform):
                new_x = m00 * x + m01 * y + m02 * z + m03
                new_y = m10 * x + m11 * y + m12 * z + m13
                new_z = m20 * x + m21 * y + m22 * z + m23

                z_scaling_lookup_offset = -10000
                z_scaling_lookup_delta = 1.000000e-001
                # [***]
                new_z_lookupPos = str(int((new_z - z_scaling_lookup_offset) / z_scaling_lookup_delta))
                new_z = float(z_scaling_lookup[new_z_lookupPos])

                if (transform):
                    x = new_x - somaX
                    y = new_y - somaY
                    z = new_z - somaZ

                if (not transform):
                    x = x - ave_somaX
                    y = y - ave_somaY
                    z = z - ave_somaZ

                polarR = math.sqrt(x ** 2 + y ** 2 + z ** 2)

                if re.search("axon", leg.segmentName):
                    x_axon = x
                    y_axon = y
                    z_axon = z

                    x_axon_all.append(x_axon)
                    y_axon_all.append(y_axon)
                    z_axon_all.append(z_axon)

                # read in a point as a list:
                points = [int(pt_count + 1), float(x), float(y), float(z), float(r), float(polarR)]
                # print "points = ", points[1], points[2], points[3]

                if (pt_count == 0):
                    # Points for a new segment, start list of tuples
                    leg.points = [points]

                else:
                    # Existing segment, append list of tuples
                    leg.points.append(points)

                # Dictionary of segment names and points in class cell:
                cell.points[leg.segmentName] = leg.points

                # Increment number of points in segment:
                pt_count += 1

                if (verbose == 9):
                    print
                    # print "found pt3dadd line     :", line
                    print "**************x        = ", x
                    print "**************y        = ", y
                    print "**************z        = ", z
                    print "**************r        = ", r
                    print "*****leg.segmentName   = ", leg.segmentName
                    # print "*****leg.points        = ", leg.points

        # Last segment (put outside loop instead of searching for "EOF" as below)
        segmentLength = leg.length()

        # Subsequent "create" lines, append segments
        # (sons and fathers done in "connect" line below):
        # cell.segments.append(leg.segmentName)
        cell.segmentLengths.append(segmentLength)
        cell.numpoints.append(pt_count)

        # Find and set segment label:
        if re.search("soma", leg.segmentName): leg.label = 'soma'
        if re.search("dend", leg.segmentName): leg.label = 'basal'
        if re.search("apical", leg.segmentName): leg.label = 'apical'
        if re.search("axon", leg.segmentName): leg.label = 'axon'

        # Calculate axon box volume:
        if (axon):
            axon_box_vol = (max(x_axon_all) - min(x_axon_all)) * (max(y_axon_all) - min(y_axon_all)) * (
                max(z_axon_all) - min(z_axon_all))
        else:
            axon_box_vol = 0

        # Axon branchpoints:
        # Count the number of nonunique fathers in axon_fathers, this is the number of
        # fathers with more than one son = number of branchpoint!
        #
        # Generates a list of tuples the second element is a tally of unique elements.
        # Difference is number of nonunique elements
        #
        # http://stackoverflow.com/questions/9426805/counting-unique-elements-in-a-list-in-python
        #
        # Note: Need python2.6 for enumerate(x,y)!!!!

        # words = ['the', 'counter', 'starts', 'the', 'starts', 'for']
        uniq = set()
        result = []
        # for i, word in enumerate(words, 1):
        for i, segment in enumerate(axon_fathers, 1):
            uniq.add(segment)
            result.append((i, len(uniq)))

        unique = set(axon_fathers)
        num_threes = 0
        for segment in unique:
            if axon_fathers.count(segment) > 2: num_threes += 1

        if (axon):
            # Subtract the number of fathers with 3 or more sons (See below **^**)
            axon_bp = (len(result) - len(uniq)) - num_threes

        else:
            axon_bp = 0

        # Example to show why the number of fathers with three or more sons must be subtracted:
        #
        # result
        # [(1, 1), (2, 2), (3, 3), (4, 3), (5, 3), (6, 4)] num_nonunique = 6 - 4 = 2 (yes!)
        #
        # words = ['the', 'counter', 'starts', 'the', 'starts', 'for', 'the']
        # [(1, 1), (2, 2), (3, 3), (4, 3), (5, 3), (6, 4), (7, 4)] num_nonunique = 7 - 4 = 3 (no!)
        #
        # **^**: Therefore, must count the number of fathers that repeat more
        # than twice and subtract:
        # axon_bp = (len(result) - len(uniq)) - num_two_or_more

        if (verbose == 0):
            print "ave_somaX = ", ave_somaX
            print "ave_somaY = ", ave_somaY
            print "ave_somaZ = ", ave_somaZ
            print "numsomapts = ", numsomapts

            print
            print "End of File Reached!"
            print

            #############################################################################
            # End for line in infile:
            #############################################################################

    #################################################################################
    # Hermann Cuntz hoc files
    #################################################################################
    if (hermannhoc):
        if (verbose == 0):
            print
            print "************************************************"
            print "One of Hermann Cuntz's Clone hoc files:"
            print "Set hermannhoc =", hermannhoc
            print "************************************************"
            print

        connect_count = 0

        #############################################################################
        for line in infile:
            #############################################################################

            # Fathers and Sons (not after Turgenev)
            if re.search("connect", line):

                parts = line.split(',')
                firstpart = parts[0]
                secondpart = parts[1]

                # New Segment
                if (connect_count == 0):
                    cell.sons = [firstpart.split()[1][:-3]]  # Remove "(0)" at the end
                    cell.father = [secondpart.split('(')[0]]
                    connect_count += 1
                else:
                    cell.sons.append(firstpart.split()[1][:-3])  # Remove "(0)" at the end
                    cell.father.append(secondpart.split('(')[0])
                    connect_count += 1

                if (verbose == 0):
                    print
                    print "**************found connect line : ", line
                    print "**************parts              = ", parts
                    print "**************segment            = ", leg.segmentName
                    print "**************sons               =    Don't know yet!"
                    print "**************father             = ", leg.father

            # New Segment
            if re.search("pt3dclear", line):

                # From previous segment, record length, print:
                if (pt_count > 0):
                    segmentLength = leg.length()

                # New Segment:
                firstpart = line.split()[0]
                leg = hc.Segment()
                leg.segmentName = firstpart

                # Find and set segment label:
                if re.search("soma", leg.segmentName): leg.label = 'soma'
                if re.search("basal", leg.segmentName): leg.label = 'basal'
                if re.search("apical", leg.segmentName): leg.label = 'apical'

                if (seg_count == 0):
                    # First "create" line, start segment, sons,
                    # and fathers lists and points dictionary:
                    cell.segments = [leg.segmentName]
                    cell.points = {}
                    cell.segmentLengths = []
                    cell.label = [leg.label]
                    cell.numpoints = []

                else:
                    # Subsequent "create" lines, append segments
                    # (sons and fathers done in "connect" line below):
                    cell.segments.append(leg.segmentName)
                    cell.segmentLengths.append(segmentLength)
                    cell.label.append(leg.label)
                    cell.numpoints.append(pt_count)

                pt_count = 0

                if (verbose == 2):
                    print
                    print "**************found pt3dclear line : ", line
                    print "**************line.split()      = ", line.split()
                    print "**************leg.segmentName   = ", leg.segmentName
                    print

                # Increment number of segments
                seg_count += 1

            # points in segment
            if re.search("pt3dadd", line):

                parts = line.split(',')
                x = float(parts[0].split('(')[1])
                y = float(parts[1])
                z = float(parts[2])
                r = float(parts[3].split(')')[0])

                polarR = math.sqrt(x ** 2 + y ** 2 + z ** 2)

                # read in a point as a list:
                points = [int(pt_count + 1), float(x), float(y), float(z), float(r), float(polarR)]

                if (pt_count == 0):
                    # Points for a new segment, start list
                    leg.points = [points]

                else:
                    # Existing segment, append list
                    leg.points.append(points)

                # Dictionary of segment names and points in class cell:
                cell.points[leg.segmentName] = leg.points

                # Increment number of points in segment:
                pt_count += 1

                if (verbose == 0):
                    print
                    print "found pt3dadd line     :", line
                    print "**************x        = ", x
                    print "**************y        = ", y
                    print "**************z        = ", z
                    print "**************r        = ", r
                    print "*****leg.segmentName   = ", leg.segmentName
                    print "*****leg.points        = ", leg.points

            # if (re.search("proc_basic_shape", line) and  pt_count > 0):
            if (re.search("access ", line) and pt_count > 0):

                segmentLength = leg.length()

                # Subsequent "create" lines, append segments
                # (sons and fathers done in "connect" line below):
                cell.segmentLengths.append(segmentLength)
                cell.numpoints.append(pt_count)

                # Find and set segment label:
                if re.search("soma", leg.segmentName): leg.label = 'soma'
                if re.search("basal", leg.segmentName): leg.label = 'basal'
                if re.search("apical", leg.segmentName): leg.label = 'apical'

                if (verbose == 0):
                    print "End of File Reached."

                    #############################################################################
                    # End "for line in infile:"
                    #############################################################################

    #####################################################################################
    # Summarizing the data         #
    ################################

    # seg_count for whole cell
    cell.segment_counter = seg_count
    # ________________________________________________________________________________________________________________#
    #                           Feature Extractor                                                                    #
    # ________________________________________________________________________________________________________________#
    #                                                                                                                #
    if (verbose == 90):
        print "****************************************************************************"
        print "BEGINNING FEATURE EXTRACTION"
        print "****************************************************************************"

    #####################################################################################
    # Make celldata dictionary (Datastructure)    #
    ###############################################

    ############################################################################
    # Make for loop list to avoid axon, circle and closed segments

    # For doing axons the streamlined way
    gen = (segment for segment in cell.segments if (re.search("axon", segment) == None and \
                                                    re.search("circle", segment) == None and \
                                                    re.search("closed", segment) == None))

    segments = []
    for segment in gen:
        segments.append(segment)

    ############################################################################
    # Make father/son dictionary with 'None' entries for orphans and bachelors

    father = dict(zip(cell.sons, cell.father))
    son = dict(zip(cell.father, cell.sons))

    # axon_fathers = []
    for segment in cell.segments:
        # for segment in gen:

        if (segment not in cell.sons): father[segment] = 'None'
        if (segment not in cell.father): son[segment] = 'None'

    # Also a segment length, label, numpoints etc. dictionary
    segmentlength = dict(zip(cell.segments, cell.segmentLengths))
    segmentlabel = dict(zip(cell.segments, cell.label))
    segmentnumpoints = dict(zip(cell.segments, cell.numpoints))

    ##############################################################################
    # Loop to make celldata dictionary:
    #
    # celldata = {
    #              segment : 0. [sonslist]                    (list of sons of current segment)
    #                        1. father[segment]               (segment's father)
    #                        2. segmentlabel[segment]         ('apical', 'basal', 'axon', or 'soma')
    #                        3. segmentlength[segment][0]     (segment length),
    #                        4. Type                          ('Intersegment', 'Ending', or 'Beginning')
    #                        5. segmentnumpoints[segment]     (number of points in segment)
    #                        6. cell.points[segment]          (coords of points)
    #                        7. segmentlength[segment][1]     (segment surface area)
    #            }
    #
    # Make celldata by looping through every segment in the cell.
    #  * For each segment start a new sonslist
    #  * Go through three branches for each segment.
    #        + In first branch check if son is in keys of father dict AND equals current segment
    #          If so, append sonslist and celldata for Type = 'Intersegment'
    #        + In next branch check if sonslist is empty.
    #          If so, append celldata for Type = 'Ending'
    #        + In next branch check if current segment has no father
    #          If so, append celldata for Type = 'Beginning'

    if (verbose == 90):
        print "Starting celldata construction ..."

    celldata = {}

    for segment in segments:
        # for segment in cell.segments:

        # print "segment in gen 2 =", segment
        # print "segmentLength[segment] = ", segmentLength[segment]

        sonslist = []
        # for son in father.keys():
        for son in cell.segments:

            # Used .strip() to strip off white space
            son = son.strip()

            if (son in father.keys() and father[son].strip() == segment.strip()):
                sonslist.append(son)
                celldata[segment] = \
                    [sonslist, father[segment].strip(), segmentlabel[segment], \
                     segmentlength[segment][0], 'Intersegment', \
                     segmentnumpoints[segment], cell.points[segment], segmentlength[segment][1]]

            elif (sonslist == []):
                celldata[segment] = \
                    [sonslist, father[segment].strip(), segmentlabel[segment], \
                     segmentlength[segment][0], 'Ending', \
                     segmentnumpoints[segment], cell.points[segment], segmentlength[segment][1]]

            elif (father[segment] == 'None'):
                celldata[segment] = \
                    [sonslist, father[segment].strip(), segmentlabel[segment], \
                     segmentlength[segment][0], 'Beginning', \
                     segmentnumpoints[segment], cell.points[segment], segmentlength[segment][1]]

    # print "father.keys() =", father.keys()
    if verbose == 90:
        print "Done celldata dictionary ..."
        print "Starting celldata traversal ..."

    #####################################################################################
    # Traverse celldata, compute lengths     #
    ##########################################

    #####################################################################################
    # Branch Order

    # Pick starter segment with
    # (1) label 'soma',
    # (2) nonempty sonslist, and
    # (3) father segment not containing 'soma'.
    # (This last requirement makes branch order counting in Cuntz hoc files
    # consistent with that of neuroconv hoc files.)

    ##########################################
    # for segment in cell.segments:
    for segment in segments:
        ##########################################

        if (verbose == 0):
            print "celldata[", segment, "][0] (Sons) = ", celldata[segment][0]

        if (celldata[segment][2] == 'soma' and \
                        celldata[segment][0] != [] and \
                        re.search("soma", celldata[segment][1]) == None):

            starter_segment = segment
            if (verbose == 90):
                print "In hocreader:"
                print "Starter segment (root) =", starter_segment
                print "Sons of starter segment (Primagenitors)  =", celldata[starter_segment][0]
                print

            if (starter_segment == []):
                print
                print " ********************** Possible unconnected Soma!          **********************"
                print " ********************** Check HOC file                      **********************"
                print " ********************** Exiting without extracting features **********************"
                print

                sys.exit()

    #####################################################################################
    #####################################################################################
    # initialize branchorder to make it easy to tell when a segment's been visited in
    branchorder = {}
    nonbinary = False
    # for segment in celldata.keys():
    for segment in segments:
        branchorder[segment] = [-1, 0]

    branchorder = fcn.traverse(starter_segment, celldata, branchorder, nonbinary)

    if (verbose == 90):
        print "****************************************************************************"
        print " Done traversing tree                                    "
        print " Using celldata to compute features                      "
        print "****************************************************************************"

    #####################################################################################
    # Use celldata to compute features #
    ####################################

    #####################################################################################
    # Cell totals
    # initialize
    axon_ending = 0
    numaxonpts = 0

    basallength = 0
    apicallength = 0
    dendriticlength = 0  # dendriticlength = basallength + apicallength
    length = 0  # length = dendriticlength + axon_length (not used anywhere)
    surface_area = 0
    bp = 0  # dendritic branching points
    basalbp = 0
    apicalbp = 0
    ending = 0
    basalending = 0
    apicalending = 0
    numpts = 0
    numbasalpts = 0
    numapicalpts = 0

    max_xs = []
    max_ys = []
    max_zs = []

    min_xs = []
    min_ys = []
    min_zs = []

    maxbasal_xs = []
    maxbasal_ys = []
    maxbasal_zs = []

    minbasal_xs = []
    minbasal_ys = []
    minbasal_zs = []

    maxapical_xs = []
    maxapical_ys = []
    maxapical_zs = []

    minapical_xs = []
    minapical_ys = []
    minapical_zs = []

    totalr = 0
    ave_r_segs = []

    x_all = []
    y_all = []
    z_all = []
    r_all = []
    polarR_all = []

    xbp_all = []
    ybp_all = []
    zbp_all = []
    rbp_all = []

    xE_all = []
    yE_all = []
    zE_all = []
    rE_all = []
    polarRE_all = []

    axonx_all = []
    axony_all = []
    axonz_all = []
    axonr_all = []

    apicalx_all = []
    apicaly_all = []
    apicalz_all = []
    apicalr_all = []

    basalx_all = []
    basaly_all = []
    basalz_all = []
    basalr_all = []

    dendx_all = []
    dendy_all = []
    dendz_all = []
    dendr_all = []

    Bp_coords = {}
    E_coords = {}
    bo_S = []
    bo_S_len = []
    bo_I = []
    bo_E = []

    dendx = []
    dendy = []
    dendz = []
    dendr = []

    apicalx = []
    apicaly = []
    apicalz = []
    apicalr = []

    basalx = []
    basaly = []
    basalz = []
    basalr = []

    axonx = []
    axony = []
    axonz = []
    axonr = []

    x = []
    y = []
    z = []
    r = []
    polarR = []
    name = []
    # cellpoints = {}

    primalengths = {}

    #####################################################################################
    # Determine maxBranchName
    #
    # If no apical dendrite specified in hocfile,
    # Find "apical" ("prima"genitor---see traverse.py!) branch
    # either using the longest branch or
    # the branch with maximum z coord (as in php code)
    #
    # (I) First loop through starter segments (primagenitors) and calculate the length of each
    #
    # Set lengthtest true below to test on branch length,
    # false to test on highest z

    # Length of each primagenitor
    lengthprim = 0
    # List of z-values of the primagenitor
    zprima = []

    # if (not hermannhoc):
    if (not apical and not axon):
        for prima in celldata[starter_segment][0]:
            # Exclude axon since it will obviously be the longest!
            if (not re.search("axon", prima)):
                lengthprim = 0
                # print "hermannhoc =", hermannhoc
                # if (hermannhoc): prima = prima[:-3]
                # print "prima = ", prima
                for segment in celldata.keys():
                    # print "   celldata.keys = ", segment
                    # Look for agreement at beginning of string (use re.match):
                    if re.match(prima, segment):
                        # if re.search(prima, segment):
                        # print "       Found segment", segment, "like", prima
                        lengthprim = lengthprim + celldata[segment][3]
                        zp = [row[3] for row in cell.points[segment]]
                        # zprima = zprima + zp
                        zprima.extend(zp)

                maxzprima = max(zprima)

                primalengths[prima] = [lengthprim, maxzprima]
                if (verbose == 0):
                    print "For primagenitor ", prima, " total length = ", primalengths[prima][0], \
                        "and max z =", primalengths[prima][1]
                    # print "primlengths[", prima, "][1]= ", primalengths[prima][1]
            else:
                primalengths[prima] = [0, 0]

        branchlengths = []
        for prima in primalengths.keys():
            print "prima = ", prima
            print "primalengths[prima][0]) = ", primalengths[prima][0]
            branchlengths.append(primalengths[prima][0])

        maxBranch = max(branchlengths)

        # *****************************************************************************
        # Lengthtest for determining the apical dendrite
        #
        # Set true -- to test on branch length,
        #     false -- to test on highest z
        #
        # lengthtest = False # makes ? cells work
        lengthtest = True  # makes Robert's registered cells and human cells work
        # *****************************************************************************

        for prima in primalengths.keys():
            # Test on either max length or max z coord.:
            if (lengthtest):
                if (primalengths[prima][0] == maxBranch):
                    maxBranchName = prima
                    if (verbose == 0):
                        print "Testing on length"
                        print "Longest branch is", maxBranchName, "with length = ", maxBranch
            else:
                if (primalengths[prima][1] == maxzprima):
                    maxBranchName = prima
                    if (verbose == 0):
                        print "testing on max z"
                        print "Branch with largest z coord is", maxBranchName, "with max z =  ", maxzprima

    # else: # hermannhoc is true
    else:  # Apical dendrite specified in hoc file

        maxBranchName = 'apical'
        if (verbose == 0):
            # print " hermannhoc = ", hermannhoc
            print " Using test on string 'apical' to determine apical dendrite...."

    #####################################################################################
    # Loop over segments to build various lists

    ########################################
    for segment in celldata.keys():
        ########################################

        if celldata[segment][2] != 'soma' and \
                (re.search("apical", segment) or re.search("dend", segment)):

            # Make sure not to count segment with only one son as branch points:
            if ((celldata[segment][4] == 'Intersegment') and (len(celldata[segment][0]) > 1)):
                # Increment number of branchpoints
                bp += 1

                # make list of branchorders from branchorder dict to compute BO features:
                bo_I.append(branchorder[segment][0])

                # Find coords of bp (last point of every 'Intersegment'):
                Bp_coords[segment] = cell.points[segment][len(cell.points[segment]) - 1]

            if (celldata[segment][4] == 'Ending'):
                ending += 1

                # make list of branchorders from branchorder dict to compute BO features:
                bo_E.append(branchorder[segment][0])

                # Find coords of ending (use last point of every 'Ending'):
                E_coords[segment] = cell.points[segment][len(cell.points[segment]) - 1]

            numpts = numpts + celldata[segment][5]
            length = length + celldata[segment][3]  # used also for axon below to compute total cell length (not used)
            dendriticlength = dendriticlength + celldata[segment][3]

            # make list of branchorders from branchorder dict to compute BO features:
            bo_S.append(branchorder[segment][0])
            bo_S_len.append(branchorder[segment][1])

            # Box XYZ volume of the whole cell.
            # For handy debug prints set printmaxmin = True
            # x, y, zs for each segment:
            dendx = [row[1] for row in cell.points[segment]]
            dendy = [row[2] for row in cell.points[segment]]
            dendz = [row[3] for row in cell.points[segment]]
            dendr = [row[4] for row in cell.points[segment]]
            polarR = [row[5] for row in cell.points[segment]]

            # Surface area
            surface_area = surface_area + celldata[segment][7]

            x_all.extend(dendx)
            y_all.extend(dendy)
            z_all.extend(dendz)
            r_all.extend(dendr)
            polarR_all.extend(polarR)

            dendx_all.extend(dendx)
            dendy_all.extend(dendy)
            dendz_all.extend(dendz)
            dendr_all.extend(dendr)

            # x, y, zs for each bp:
            if segment in Bp_coords.keys():
                xbp = Bp_coords[segment][1]
                ybp = Bp_coords[segment][2]
                zbp = Bp_coords[segment][3]
                rbp = Bp_coords[segment][4]

                xbp_all.append(xbp)
                ybp_all.append(ybp)
                zbp_all.append(zbp)
                rbp_all.append(rbp)

            # x, y, zs for each E:
            if segment in E_coords.keys():
                xE = E_coords[segment][1]
                yE = E_coords[segment][2]
                zE = E_coords[segment][3]
                rE = E_coords[segment][4]
                polarRE = E_coords[segment][5]

                xE_all.append(xE)
                yE_all.append(yE)
                zE_all.append(zE)
                rE_all.append(rE)
                polarRE_all.append(polarRE)

        ##################################
        # Basal or Dend
        # Pick segment names with "basal" labels or that are NOT the selected apical branch
        if celldata[segment][2] == 'basal':

            basallength = basallength + celldata[segment][3]
            # if celldata[segment][4] == 'Intersegment': basalbp += 1
            # Only increment bp at actual branch points (more than one son)
            if ((celldata[segment][4] == 'Intersegment') and (len(celldata[segment][0]) > 1)): basalbp += 1
            if celldata[segment][4] == 'Ending': basalending += 1
            numbasalpts = numbasalpts + celldata[segment][5]

            # Box volume of basal tree:
            # For handy debug prints set printmaxmin = True
            # x, y, zs for each segment:
            basalx = [row[1] for row in cell.points[segment]]
            basaly = [row[2] for row in cell.points[segment]]
            basalz = [row[3] for row in cell.points[segment]]
            basalr = [row[4] for row in cell.points[segment]]

            basalx_all.extend(basalx)
            basaly_all.extend(basaly)
            basalz_all.extend(basalz)
            basalr_all.extend(r)

        # If no apical dendrite specified, test on branch with most length or largest z:
        if re.search(maxBranchName, segment):

            apicallength = apicallength + celldata[segment][3]
            # if celldata[segment][4] == 'Intersegment': apicalbp += 1
            # Only increment bp at actual branch points (more than one son)
            if ((celldata[segment][4] == 'Intersegment') and (len(celldata[segment][0]) > 1)): apicalbp += 1
            if celldata[segment][4] == 'Ending': apicalending += 1
            numapicalpts = numapicalpts + celldata[segment][5]

            # Box volume of apical tree:
            # For handy debug prints set printmaxmin = True
            # x, y, zs for each segment:
            x = [row[1] for row in cell.points[segment]]
            y = [row[2] for row in cell.points[segment]]
            z = [row[3] for row in cell.points[segment]]
            r = [row[4] for row in cell.points[segment]]

            apicalx_all.extend(x)
            apicaly_all.extend(y)
            apicalz_all.extend(z)
            apicalr_all.extend(r)

        ##################################
        # Axon only
        # Do everything in Just "Ne" or "Dend" and "Basal" if branches so that exceptions
        # aren't needed for axon only case
        # if ((not apical) and (not basal)):
        if (celldata[segment][2] == 'axon'):

            numpts = numpts + celldata[segment][5]
            length = length + celldata[segment][3]

            # For handy debug prints set printmaxmin = True
            # x, y, zs for each axon segment:
            x = [row[1] for row in cell.points[segment]]
            y = [row[2] for row in cell.points[segment]]
            z = [row[3] for row in cell.points[segment]]
            r = [row[4] for row in cell.points[segment]]
            polarR = [row[5] for row in cell.points[segment]]

            # Surface area
            surface_area = surface_area + celldata[segment][7]

            x_all.extend(x)
            y_all.extend(y)
            z_all.extend(z)
            r_all.extend(r)
            polarR_all.extend(polarR)

            # x, y, zs for each bp:
            if segment in Bp_coords.keys():
                xbp = Bp_coords[segment][1]
                ybp = Bp_coords[segment][2]
                zbp = Bp_coords[segment][3]
                rbp = Bp_coords[segment][4]

                xbp_all.append(xbp)
                ybp_all.append(ybp)
                zbp_all.append(zbp)
                rbp_all.append(rbp)

            # x, y, zs for each E:
            if segment in E_coords.keys():
                xE = E_coords[segment][1]
                yE = E_coords[segment][2]
                zE = E_coords[segment][3]
                rE = E_coords[segment][4]
                polarRE = E_coords[segment][5]

                xE_all.append(xE)
                yE_all.append(yE)
                zE_all.append(zE)
                rE_all.append(rE)
                polarRE_all.append(polarRE)

            if celldata[segment][4] == 'Ending': axon_ending += 1
            numaxonpts = numaxonpts + celldata[segment][5]

            axonx_all.extend(x)
            axony_all.extend(y)
            axonz_all.extend(z)
            axonr_all.extend(r)

    ########################################
    # End for segment in celldata.keys():
    ########################################

    #####################################################################################
    # Compute with various lists

    ave_bo_I = fcn.meanstdv(bo_I)[0]
    std_bo_I = fcn.meanstdv(bo_I)[1]

    ##############
    # Th's

    # Create lists x_th1, x_th2, and x_th3, etc.

    # Th3: z > 2*apicalz_max/3
    # Th2: apicalz_max/3 < z < 2*apicalz_max/3
    # Th1: z < apicalz_max/3

    # Set up dictionaries for coords for each of th1, 2, and 3
    #
    # E.g. point_th3[z] = (x, y) for x, y, z where z is defined as above
    # Thus, e.g,. z_th3 = point_th3.keys(1)
    #             y_th3 = point_th3[z_th3][1]
    #             x_th3 = point_th3[z_th3][0]
    #
    # And similarly for th2 and th1

    # Exception handling for case where there are no apical dendrites:
    try:
        apicalz_max = max(apicalz_all)
        apicalz_min = min(apicalz_all)
        apicalz_span = apicalz_max - apicalz_min

    except ValueError:
        apicalz_max = 0
        apicalz_min = 0
        apicalz_span = 0

    # z_th = z_max/3
    # z_th2 = 2*z_th

    point = {}
    point = dict(zip(apicalz_all, zip(apicalx_all, apicaly_all)))

    point_th3 = {}
    point_th2 = {}
    point_th1 = {}

    for z in point.keys():

        # print "point[", z, "] = ", point[z]

        # Th3:
        if z > apicalz_min + 2 * apicalz_span / 3:
            point_th3[z] = point[z]
            # print "point_th3[", z, "] = ", point_th3[z]
            # print " z_th3 = ", z
        # Th2:
        if z > apicalz_min + apicalz_span / 3 and z < apicalz_min + 2 * apicalz_span / 3:
            point_th2[z] = point[z]
            # print "point_th2[", z, "] = ", point_th2[z]
            # print " z_th2 = ", z
        # Th1:
        if z < apicalz_min + apicalz_span / 3:
            point_th1[z] = point[z]
            # print "point_th1[", z, "] = ", point_th1[z]
            # print " z_th1 = ", z

    if (verbose == 0):
        print
        print " **************** apicalz_max = ", apicalz_max
        print " **************** apicalz_min = ", apicalz_min
        print " **************** apicalz_span = ", apicalz_span
        print

    x_th1 = []
    x_th2 = []
    x_th3 = []

    y_th1 = []
    y_th2 = []
    y_th3 = []

    z_th1 = []
    z_th2 = []
    z_th3 = []

    # Exception handling for case where there are no apical dendrites:
    if (apicalz_max != 0) and (len(point_th2) > 0):

        for z in point_th1.keys():
            z_th1.append(z)
            y_th1.append(point_th1[z][1])
            x_th1.append(point_th1[z][0])

        for z in point_th2.keys():
            z_th2.append(z)
            y_th2.append(point_th2[z][1])
            x_th2.append(point_th2[z][0])

        for z in point_th3.keys():
            z_th3.append(z)
            y_th3.append(point_th3[z][1])
            x_th3.append(point_th3[z][0])

    else:

        x_th1 = [0]
        x_th2 = [0]
        x_th3 = [0]

        y_th1 = [0]
        y_th2 = [0]
        y_th3 = [0]

        z_th1 = [0]
        z_th2 = [0]
        z_th3 = [0]

        apicalx_all = [0]
        apicaly_all = [0]
        apicalz_all = [0]
        apicalr_all = [0]

    ##################################
    # Min/maxs/aves
    ##################################

    max_x_th1 = max(x_th1)
    max_x_th2 = max(x_th2)
    max_x_th3 = max(x_th3)

    max_y_th1 = max(y_th1)
    max_y_th2 = max(y_th2)
    max_y_th3 = max(y_th3)

    max_z_th1 = max(z_th1)
    max_z_th2 = max(z_th2)
    max_z_th3 = max(z_th3)

    min_x_th1 = min(x_th1)
    min_x_th2 = min(x_th2)
    min_x_th3 = min(x_th3)

    min_y_th1 = min(y_th1)
    min_y_th2 = min(y_th2)
    min_y_th3 = min(y_th3)

    min_z_th1 = min(z_th1)
    min_z_th2 = min(z_th2)
    min_z_th3 = min(z_th3)

    if (basal or apical):

        # Ne
        max_x = max(x_all)
        max_y = max(y_all)
        max_z = max(z_all)
        max_r = max(r_all)

        min_x = min(x_all)
        min_y = min(y_all)
        min_z = min(z_all)
        min_r = min(r_all)

        ave_x = fcn.meanstdv(x_all)[0]
        ave_y = fcn.meanstdv(y_all)[0]
        ave_z = fcn.meanstdv(z_all)[0]
        ave_r = fcn.meanstdv(r_all)[0]

        std_x = fcn.meanstdv(x_all)[1]
        std_y = fcn.meanstdv(y_all)[1]
        std_z = fcn.meanstdv(z_all)[1]
        std_r = fcn.meanstdv(r_all)[1]

        # Dendritic

        maxdend_x = max(dendx_all)
        maxdend_y = max(dendy_all)
        maxdend_z = max(dendz_all)
        maxdend_r = max(dendr_all)

        mindend_x = min(dendx_all)
        mindend_y = min(dendy_all)
        mindend_z = min(dendz_all)
        mindend_r = min(dendr_all)

    # Zero out features of no dendrites
    else:

        max_x = 0
        max_y = 0
        max_z = 0
        max_r = 0

        min_x = 0
        min_y = 0
        min_z = 0
        min_r = 0

        ave_x = 0
        ave_y = 0
        ave_z = 0
        ave_r = 0

        std_x = 0
        std_y = 0
        std_z = 0
        std_r = 0

        maxdend_x = 0
        maxdend_y = 0
        maxdend_z = 0
        maxdend_r = 0

        mindend_x = 0
        mindend_y = 0
        mindend_z = 0
        mindend_r = 0

    # Apical (Md)
    maxapical_x = max(apicalx_all)
    maxapical_y = max(apicaly_all)
    maxapical_z = max(apicalz_all)
    maxapical_r = max(apicalr_all)

    minapical_x = min(apicalx_all)
    minapical_y = min(apicaly_all)
    minapical_z = min(apicalz_all)
    minapical_r = min(apicalr_all)

    aveapical_x = fcn.meanstdv(apicalx_all)[0]
    aveapical_y = fcn.meanstdv(apicaly_all)[0]
    aveapical_z = fcn.meanstdv(apicalz_all)[0]
    aveapical_r = fcn.meanstdv(apicalr_all)[0]

    stdapical_x = fcn.meanstdv(apicalx_all)[1]
    stdapical_y = fcn.meanstdv(apicaly_all)[1]
    stdapical_z = fcn.meanstdv(apicalz_all)[1]
    stdapical_r = fcn.meanstdv(apicalr_all)[1]

    if (verbose == 0):
        print " min_apicalz (all) = ", minapical_z
        print " max_apicalz (all) = ", maxapical_z

    # Bp
    try:
        max_xbp = max(xbp_all)
        max_ybp = max(ybp_all)
        max_zbp = max(zbp_all)
        max_rbp = max(rbp_all)

        min_xbp = min(xbp_all)
        min_ybp = min(ybp_all)
        min_zbp = min(zbp_all)
        min_rbp = min(rbp_all)

    except ValueError:
        max_xbp = 0
        max_ybp = 0
        max_zbp = 0
        max_rbp = 0

        min_xbp = 0
        min_ybp = 0
        min_zbp = 0
        min_rbp = 0

    ave_xbp = fcn.meanstdv(xbp_all)[0]
    ave_ybp = fcn.meanstdv(ybp_all)[0]
    ave_zbp = fcn.meanstdv(zbp_all)[0]
    ave_rbp = fcn.meanstdv(rbp_all)[0]

    std_xbp = fcn.meanstdv(xbp_all)[1]
    std_ybp = fcn.meanstdv(ybp_all)[1]
    std_zbp = fcn.meanstdv(zbp_all)[1]
    std_rbp = fcn.meanstdv(rbp_all)[1]

    # E
    try:

        max_xE = max(xE_all)
        max_yE = max(yE_all)
        max_zE = max(zE_all)
        max_rE = max(rE_all)
        max_polarRE = max(polarRE_all)

        min_xE = min(xE_all)
        min_yE = min(yE_all)
        min_zE = min(zE_all)
        min_rE = min(rE_all)
        min_polarRE = min(polarRE_all)

        ave_xE = fcn.meanstdv(xE_all)[0]
        ave_yE = fcn.meanstdv(yE_all)[0]
        ave_zE = fcn.meanstdv(zE_all)[0]
        ave_rE = fcn.meanstdv(rE_all)[0]
        ave_polarRE = fcn.meanstdv(polarRE_all)[0]

        std_xE = fcn.meanstdv(xE_all)[1]
        std_yE = fcn.meanstdv(yE_all)[1]
        std_zE = fcn.meanstdv(zE_all)[1]
        std_rE = fcn.meanstdv(rE_all)[1]

        sum_polarRE = sum(polarRE_all)

    except ValueError:

        max_xE = 0
        max_yE = 0
        max_zE = 0
        max_rE = 0
        max_polarRE = 0

        min_xE = 0
        min_yE = 0
        min_zE = 0
        min_rE = 0
        min_polarRE = 0

        ave_xE = 0
        ave_yE = 0
        ave_zE = 0
        ave_rE = 0
        ave_polarRE = 0

        std_xE = 0
        std_yE = 0
        std_zE = 0
        std_rE = 0

        sum_polarRE = 0

    # Basal (Dend)
    if (basal):
        maxbasal_x = max(basalx_all)
        maxbasal_y = max(basaly_all)
        maxbasal_z = max(basalz_all)

        minbasal_x = min(basalx_all)
        minbasal_y = min(basaly_all)
        minbasal_z = min(basalz_all)

        basal_z_span = maxbasal_z - minbasal_z

        avebasal_x = fcn.meanstdv(basalx_all)[0]
        avebasal_y = fcn.meanstdv(basaly_all)[0]

    else:

        maxbasal_x = 0
        maxbasal_y = 0
        maxbasal_z = 0

        minbasal_x = 0
        minbasal_y = 0
        minbasal_z = 0

        basal_z_span = maxbasal_z - minbasal_z

    # if (axon_oldway):
    try:
        maxaxon_x = max(axonx_all)
        maxaxon_y = max(axony_all)
        maxaxon_z = max(axonz_all)

        minaxon_x = min(axonx_all)
        minaxon_y = min(axony_all)
        minaxon_z = min(axonz_all)
    # else:
    except ValueError:
        maxaxon_x = 0
        maxaxon_y = 0
        maxaxon_z = 0

        minaxon_x = 0
        minaxon_y = 0
        minaxon_z = 0

    ##################################
    # Medians
    ##################################

    # Medians
    #
    # The way medians are computed (per legacy php code) is to count the number of points in
    # the lower 1/2, 1/4 and 3/4 and defining Z_Median, Z_1Th_Median, and Z_2Th_Median accordingly
    # (This is done in, e.g., get.S.neuroPart.Y.XZ etc)

    num_half = int(numpts / 2)
    num_th1 = int(numpts / 4)
    num_th2 = int(3 * numpts / 4)

    if (verbose == 0):
        print " length z_all = ", len(z_all)
        print " num_half = ", num_half
        print " num_th1 = ", num_th1
        print " num_th2 = ", num_th2

    z_all.sort()
    zbp_all.sort()
    zE_all.sort()

    # _S_
    z_half_all = []
    z_th1_all = []
    z_th2_all = []
    for i in range(len(z_all)):
        if i <= num_half: z_half_all.append(z_all[i])
        if i <= num_th1: z_th1_all.append(z_all[i])
        if i <= num_th2: z_th2_all.append(z_all[i])

    try:
        # med_z = fcn.meanstdv(z_half_all)[0]
        med_z = sum(z_half_all) / len(z_half_all)
        #  ave_z_th1 = fcn.meanstdv(z_th1_all)[0]
        ave_z_th1 = sum(z_th1_all) / len(z_th1_all)
        # ave_z_th2 = fcn.meanstdv(z_th2_all)[0]
        ave_z_th2 = sum(z_th2_all) / len(z_th2_all)

    except ZeroDivisionError:

        # med_z = fcn.meanstdv(z_half_all)[0]
        med_z = 0
        #  ave_z_th1 = fcn.meanstdv(z_th1_all)[0]
        ave_z_th1 = 0
        # ave_z_th2 = fcn.meanstdv(z_th2_all)[0]
        ave_z_th2 = 0

    # _Bp_
    zbp_half_all = []
    zbp_th1_all = []
    zbp_th2_all = []
    for i in range(len(zbp_all)):
        if i <= num_half: zbp_half_all.append(zbp_all[i])
        if i <= num_th1: zbp_th1_all.append(zbp_all[i])
        if i <= num_th2: zbp_th2_all.append(zbp_all[i])

    # Exception handling for only axon
    try:
        # med_zbp = fcn.meanstdv(zbp_half_all)[0]
        med_zbp = sum(zbp_half_all) / len(zbp_half_all)
        #  ave_zbp_th1 = fcn.meanstdv(zbp_th1_all)[0]
        ave_zbp_th1 = sum(zbp_th1_all) / len(zbp_th1_all)
        # ave_zbp_th2 = fcn.meanstdv(zbp_th2_all)[0]
        ave_zbp_th2 = sum(zbp_th2_all) / len(zbp_th2_all)

    except ZeroDivisionError:
        # med_zbp = fcn.meanstdv(zbp_half_all)[0]
        med_zbp = 0
        #  ave_zbp_th1 = fcn.meanstdv(zbp_th1_all)[0]
        ave_zbp_th1 = 0
        # ave_zbp_th2 = fcn.meanstdv(zbp_th2_all)[0]
        ave_zbp_th2 = 0

    # _E_
    zE_half_all = []
    zE_th1_all = []
    zE_th2_all = []
    for i in range(len(zE_all)):
        if i <= num_half: zE_half_all.append(zE_all[i])
        if i <= num_th1: zE_th1_all.append(zE_all[i])
        if i <= num_th2: zE_th2_all.append(zE_all[i])

    # med_zE = fcn.meanstdv(zE_half_all)[0]
    try:
        med_zE = sum(zE_half_all) / len(zE_half_all)
        #  ave_zE_th1 = fcn.meanstdv(zE_th1_all)[0]
        ave_zE_th1 = sum(zE_th1_all) / len(zE_th1_all)
        # ave_zE_th2 = fcn.meanstdv(zE_th2_all)[0]
        ave_zE_th2 = sum(zE_th2_all) / len(zE_th2_all)

    except ZeroDivisionError:
        med_zE = 0
        ave_zE_th1 = 0
        ave_zE_th2 = 0

    ##################################
    # Box volumes
    ##################################

    box_XYZ_vol = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)
    box_XYZ_vol_dend = (maxdend_x - mindend_x) * (maxdend_y - mindend_y) * (maxdend_z - mindend_z)
    box_XYZ_vol_apical = (maxapical_x - minapical_x) * (maxapical_y - minapical_y) * (maxapical_z - minapical_z)
    box_XYZ_vol_basal = (maxbasal_x - minbasal_x) * (maxbasal_y - minbasal_y) * (maxbasal_z - minbasal_z)
    box_XYZ_vol_axon = (maxaxon_x - minaxon_x) * (maxaxon_y - minaxon_y) * (maxaxon_z - minaxon_z)

    # Th's (all for apical only)
    box_XYZ_vol_th1 = (max_x_th1 - min_x_th1) * (max_y_th1 - min_y_th1) * (max_z_th1 - min_z_th1)
    box_XYZ_vol_th2 = (max_x_th2 - min_x_th2) * (max_y_th2 - min_y_th2) * (max_z_th2 - min_z_th2)
    box_XYZ_vol_th3 = (max_x_th3 - min_x_th3) * (max_y_th3 - min_y_th3) * (max_z_th3 - min_z_th3)

    boxbasal_XYZ_vol = (maxbasal_x - minbasal_x) * (maxbasal_y - minbasal_y) * (maxbasal_z - minbasal_z)
    boxapical_XYZ_vol = (maxapical_x - minapical_x) * (maxapical_y - minapical_y) * (maxapical_z - minapical_z)

    diag_XY = math.sqrt((max_x - min_x) ** 2 + (max_y - min_y) ** 2)
    diag_XY_dend = math.sqrt((maxdend_x - mindend_x) ** 2 + (maxdend_y - mindend_y) ** 2)
    diag_XY_basal = math.sqrt((maxbasal_x - minbasal_x) ** 2 + (maxbasal_y - minbasal_y) ** 2)
    diag_XY_apical = math.sqrt((maxapical_x - minapical_x) ** 2 + (maxapical_y - minapical_y) ** 2)
    # Th's (all for apical only)
    diag_XY_th1 = math.sqrt((max_x_th1 - min_x_th1) ** 2 + (max_y_th1 - min_y_th1) ** 2)
    diag_XY_th2 = math.sqrt((max_x_th2 - min_x_th2) ** 2 + (max_y_th2 - min_y_th2) ** 2)
    diag_XY_th3 = math.sqrt((max_x_th3 - min_x_th3) ** 2 + (max_y_th3 - min_y_th3) ** 2)
    diag_XYZ = math.sqrt((max_x - min_x) ** 2 + (max_y - min_y) ** 2 + (max_z - min_z) ** 2)
    diag_XYZ_th1 = math.sqrt((max_x_th1 - min_x_th1) ** 2 + (max_y_th1 - min_y_th1) ** 2 + (max_z_th1 - min_z_th1) ** 2)
    diag_XYZ_th2 = math.sqrt((max_x_th2 - min_x_th2) ** 2 + (max_y_th2 - min_y_th2) ** 2 + (max_z_th2 - min_z_th2) ** 2)
    diag_XYZ_th3 = math.sqrt((max_x_th3 - min_x_th3) ** 2 + (max_y_th3 - min_y_th3) ** 2 + (max_z_th3 - min_z_th3) ** 2)

    # bp
    box_XYZ_vol_bp = (max_xbp - min_xbp) * (max_ybp - min_ybp) * (max_zbp - min_zbp)

    diag_XY_bp = math.sqrt((max_xbp - min_xbp) ** 2 + (max_ybp - min_ybp) ** 2)
    diag_XYZ_bp = math.sqrt((max_xbp - min_xbp) ** 2 + (max_ybp - min_ybp) ** 2 + (max_zbp - min_zbp) ** 2)

    # E
    box_XYZ_vol_E = (max_xE - min_xE) * (max_yE - min_yE) * (max_zE - min_zE)

    diag_XY_E = math.sqrt((max_xE - min_xE) ** 2 + (max_yE - min_yE) ** 2)
    diag_XYZ_E = math.sqrt((max_xE - min_xE) ** 2 + (max_yE - min_yE) ** 2 + (max_zE - min_zE) ** 2)

    ##################################
    # Surface areas
    ##################################
    surf_area = length * 2 * math.pi * ave_r

    ##################################
    # Centroid
    ##################################

    # Dendritic (Ne)
    ave_centroid_dist_XYZ = math.sqrt(ave_x ** 2 + ave_y ** 2 + ave_z ** 2)
    ave_centroid_dist_XY = math.sqrt(ave_x ** 2 + ave_y ** 2)

    # Apical (Md)
    ave_centroid_dist_XYZ_apical = math.sqrt(aveapical_x ** 2 + aveapical_y ** 2 + aveapical_z ** 2)
    ave_centroid_dist_XY_apical = math.sqrt(aveapical_x ** 2 + aveapical_y ** 2)

    # Basal
    if basal:
        ave_centroid_dist_XY_basal = math.sqrt(avebasal_x ** 2 + avebasal_y ** 2)
    else:
        ave_centroid_dist_XY_basal = 0.0

    # bp
    # Exception handling for only axon
    try:
        xbp_ave = float(sum(xbp_all) / len(xbp_all))
        ybp_ave = float(sum(ybp_all) / len(ybp_all))
        zbp_ave = float(sum(zbp_all) / len(zbp_all))

    except ZeroDivisionError:
        xbp_ave = 0
        ybp_ave = 0
        zbp_ave = 0

    ave_centroid_dist_XYZ_bp = math.sqrt(xbp_ave ** 2 + ybp_ave ** 2 + zbp_ave ** 2)
    ave_centroid_dist_XY_bp = math.sqrt(xbp_ave ** 2 + ybp_ave ** 2)

    # E
    try:
        xE_ave = float(sum(xE_all) / len(xE_all))
        yE_ave = float(sum(yE_all) / len(yE_all))
        zE_ave = float(sum(zE_all) / len(zE_all))

    except ZeroDivisionError:

        xE_ave = 0
        yE_ave = 0
        zE_ave = 0

    ave_centroid_dist_XYZ_E = math.sqrt(xE_ave ** 2 + yE_ave ** 2 + zE_ave ** 2)
    ave_centroid_dist_XY_E = math.sqrt(xE_ave ** 2 + yE_ave ** 2)

    ##################################
    # Cylindrical R
    ##################################

    CyR = [0]
    CyR_bp = [0]
    CyR_E = [0]
    for i in range(len(x_all)):
        CyR.append(math.sqrt(x_all[i] ** 2 + y_all[i] ** 2))

    for i in range(len(xbp_all)):
        CyR_bp.append(math.sqrt(xbp_all[i] ** 2 + ybp_all[i] ** 2))

    for i in range(len(xE_all)):
        CyR_E.append(math.sqrt(xE_all[i] ** 2 + yE_all[i] ** 2))

    CyR_max = max(CyR)
    CyR_ave = fcn.meanstdv(CyR)[0]
    CyR_sd = fcn.meanstdv(CyR)[1]
    CyR_sum = sum(CyR)

    try:
        CyR_max_bp = max(CyR_bp)
    except ValueError:
        CyR_max_bp = 0
    CyR_ave_bp = fcn.meanstdv(CyR_bp)[0]
    CyR_sd_bp = fcn.meanstdv(CyR_bp)[1]
    CyR_sum_bp = sum(CyR_bp)

    CyR_max_E = max(CyR_E)
    CyR_ave_E = fcn.meanstdv(CyR_E)[0]
    CyR_sd_E = fcn.meanstdv(CyR_E)[1]
    CyR_sum_E = sum(CyR_E)

    # Find cell ID from input hoc file name if not specified
    if re.search("/", file):
        # file = file.split("/")[2]
        fileparts = file.split("/")
        # print "fileparts = ", file
        for piece in fileparts:
            if re.search(".hoc", piece):
                hocfilename = piece

        parts = hocfilename.split('.')

        if (not hermannhoc):
            TheID = parts[0]
            ExperimentName = "NeuroConv Hoc File"
        else:
            TheID = parts[1].split('.')[0]
            ExperimentName = "Clone Hoc File"

    # Use firstpart either numeric or otherwise
    else:
        hocfilename = file.split(".")[0]
        TheID = hocfilename

    if verbose == 0:
        print "No ID entered for cell, attempting to extract it from filename...."
        print "For file = ", file
        print "Extracted ID is: ", TheID

    Experiment_ID = TheID

    # Build Features Dictionary out of Names and Values tuples (ordered lists):
    FeatureNames = ()
    FeatureValues = ()
    Features_dict = {}

    # Morphology of Dendrites
    FeatureNames = (
        "Experiment_ID"
        , "Dendritic_Length"  # Ne_S_Sum_Length
        , "Dendritic_Volume"  # Ne_E_Box_XYZ
        , "Dendritic_Width"  # Ne_S_Diagonal_XY
        , "Dendritic_Height"  # Ne_S_Span_Z
        , "Dendritic_Branchpoints"  # Ne_Bp_N
        , "Dendritic_Branchorder_Avg"  # Ne_Bp_Local_Avg_BO
        , "Dendritic_Branchorder_Max"  # Ne_Bp_Local_Max_BO
        , "Dendritic_Endpoint_Avg"  # Ne_E_Avg_PolarR
        , "Dendritic_Endpoint_Max"  # Ne_E_Max_PolarR
        , "Dendritic_Endpoint_Sum"  # Ne_E_Sum_PolarR
        , "Dendritic_Endpoint_Centroid"  # Ne_S_Centroid_Dist_XYZ
        #
        , "Apical_Length"  # Md_S_Sum_Length
        , "Apical_Volume"  # Md_S_Box_XYZ
        , "Apical_Width"  # MdTh3_S_Diagonal_XY
        , "Apical_Height"  # Md_S_Span_Z
        , "Apical_Branchpoints"  # Md_Bp_N
        , "Apical_Endpoint_Centroid"  # Md_S_Centroid_Dist_XYZ
        #
        , "Basal_Length"  # Dend_S_Sum_Length
        , "Basal_Volume"  # Dend_S_Box_XYZ
        , "Basal_Width"  # Dend_S_Diagonal_XY
        , "Basal_Height"  # Dend_S_Span_Z
        , "Basal_Branchpoints"  # Dend_Bp_N
        #
        , "Axon_Branchpoints"
        , "Axon_Length"
        , "Axon_Volume"
    )

    try:
        min_bo_I = min(bo_I)
        max_bo_I = max(bo_I)
    except ValueError:
        min_bo_I = 0
        max_bo_I = 0

    FeatureValues = (Experiment_ID
                     , dendriticlength  # Ne_S_Sum_Length
                     , box_XYZ_vol_dend  # Ne_E_Box_XYZ
                     , diag_XY_dend  # Ne_S_Diagonal_XY
                     , maxdend_z - mindend_z  # Ne_S_Span_Z
                     , bp  # Ne_Bp_N
                     , ave_bo_I - min_bo_I  # Ne_Bp_Local_Avg_BO
                     , max_bo_I - min_bo_I  # Ne_Bp_Local_Max_BO
                     , ave_polarRE  # Ne_E_Avg_PolarR
                     , max_polarRE  # Ne_E_Max_PolarR
                     , sum_polarRE  # Ne_E_Sum_PolarR
                     , ave_centroid_dist_XYZ  # Ne_S_Centroid_Dist_XYZ
                     , apicallength  # Md_S_Sum_Length
                     , box_XYZ_vol_apical  # Md_S_Box_XYZ
                     , diag_XY_apical  # MdTh3_S_Diagonal_XY
                     , apicalz_span  # Md_S_Span_Z
                     , apicalbp  # Md_Bp_N
                     , ave_centroid_dist_XYZ_apical  # Md_S_Centroid_Dist_XYZ
                     , basallength  # Dend_S_Sum_Length
                     , box_XYZ_vol_basal  # Dend_S_Box_XYZ
                     , diag_XY_basal  # Dend_S_Diagonal_XY
                     , basal_z_span  # Dend_S_Span_Z
                     , basalbp  # Dend_Bp_N
                     , axon_bp
                     , axon_length
                     , axon_box_vol
                     )

    Features_dict = dict(zip(FeatureNames, FeatureValues))
    Features_tuple = zip(FeatureNames, FeatureValues)

    # Zero axon features spuriously calculated using "all" variables in axon only case
    if ((not apical) and (not basal)):
        max_polarRE = 0
        ave_polarRE = 0
        sum_polarRE = 0
        box_XYZ_vol_E = 0
        bp = 0
        length = 0
        ave_centroid_dist_XYZ = 0
        diag_XY = 0
        ave_bo_I = 0
        bo_I = [0]

    # Selects which features to print based on morphology
    print "Compartments of Morphology: "
    print "  apical =", apical
    print "  basal =", basal
    print "  axon = ", axon
    for featurename in Features_dict.keys():
        # No Apical and Axon
        # Delete "Basal" and "Apical" features if no Apical dendrite and Axon
        if ((not apical) and axon and ((re.search("Basal", featurename)) or (re.search("Apical", featurename)))):
            Features_dict[featurename] = 0
        # No Apical and No Axon
        # Delete "Basal" and "Apical" features if no Apical dendrite and NO Axon
        if ((not apical) and not axon and ((re.search("Basal", featurename)) or (re.search("Apical", featurename)))):
            Features_dict[featurename] = 0
        # No Axon
        # Delete Axon features if no Axon
        if ((not axon) and (re.search("Axon", featurename))):
            Features_dict[featurename] = 0
        # Only Axon
        # Delete "Dendritic" not deleted in first case
        if (axon and (not apical) and (not basal) and ((re.search("Dendritic", featurename)))):
            print featurename
            # del Features_dict[featurename]
            Features_dict[featurename] = 0

    FeatureNamesActual = ()
    FeatureValuesActual = ()
    for i in range(len(Features_tuple)):
        # Only print features in Features_dict, using Features_tuple to order the list
        if (Features_tuple[i][0] in Features_dict.keys()):
            FeatureNamesActual = FeatureNamesActual + (Features_tuple[i][0],)
            FeatureValuesActual = FeatureValuesActual + (Features_dict[Features_tuple[i][0]],)

    if (verbose == 0):
        print
        print " ****************************************************************************"
        print " ******************* Successful Feature Extraction for cell ID", Experiment_ID
        print " *******************"
        print " ******************* Printing various outputata ..."
        print " ****************************************************************************"
        print

    # Write features to csv file:

    if (multiple_inputs):
        FeaturesName = "Features.csv"
        Local_Orient_Name = "Features_Local_Orientations.csv"
    else:
        FeaturesName = hocfilename.split(".")[0] + ".csv"
        Local_Orient_Name = hocfilename.split(".")[0] + "_local_orient.csv"
    try:
        path_to_csv_file = upload_path + "/" + FeaturesName  # change local Linux path to uploads/features
        path_to_local_orient_file = upload_path + "/" + Local_Orient_Name  # change local Linux path to uploads/features
    except IOError:
        pass

    print
    print "TheID = ", TheID
    print "FeaturesName =", FeaturesName
    print "path/to/csv = " + path_to_csv_file

    if (file_num == 1):
        out = csv.writer(open(path_to_csv_file, "w+"), delimiter=',', quoting=csv.QUOTE_ALL)
        out.writerow(FeatureNamesActual)
        out.writerow(FeatureValuesActual)
        Loc_Orient_Names = ()
        Loc_Orient_Values = ()
        Loc_Orient_Names = ("Experiment_ID", "Local Orientation")
        Loc_Orient_Values = (Experiment_ID, local_orientation)
        out = csv.writer(open(path_to_local_orient_file, "w+"), delimiter=',', quoting=csv.QUOTE_ALL)
        out.writerow(Loc_Orient_Names)
        out.writerow(Loc_Orient_Values)
    else:
        out = csv.writer(open(path_to_csv_file, "a"), delimiter=',', quoting=csv.QUOTE_ALL)
        out.writerow(FeatureValuesActual)
        Loc_Orient_Values = (Experiment_ID, local_orientation)
        out = csv.writer(open(path_to_local_orient_file, "a"), delimiter=',', quoting=csv.QUOTE_ALL)
        out.writerow(Loc_Orient_Values)

    #######################
    # Print Features     #
    ######################

    i = 0
    print "************************** Calculated FEATURES: ************************"
    if (apical): print "                           Apical Dendrites"
    if (basal): print "                           Basal Dendrites"
    if (axon): print "                           Axons"
    print
    for i in range(len(Features_tuple)):
        # Only print features in Features_dict, using Features_tuple to order the list
        if (Features_tuple[i][0] in Features_dict.keys()):
            print "             ", Features_tuple[i][0], "=", Features_dict[Features_tuple[i][0]]
    print
    print "          ... more features:"
    print
    if (axon):
        print "              max/min axon x, y, z =", max(x_axon_all), min(x_axon_all), max(y_axon_all), min(
            y_axon_all), max(z_axon_all), min(z_axon_all)
    print "              Number of threes =", num_threes
    print "              Dendritic_Endpoint_Centroid_XY =", ave_centroid_dist_XY
    print "              Apical_Endpoint_Centroid_XY =", ave_centroid_dist_XY_apical
    print "              Basal_Endpoint_Centroid_XY =", ave_centroid_dist_XY_basal

    #######################################################################################
    # End "for file in FILES:"
    #######################################################################################

    # sys.exit("Got to end of hocreader");
    # print(datetime.now()-startTime)
