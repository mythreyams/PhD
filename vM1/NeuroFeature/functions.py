#######################################################################################

# Flag for branch prints:
printbotrav = False       # print high level diagnostics from traverse
printbo = False          # print everything else
printbo_short = False      # print just segs with no siblings

#from hocreader import printbo


#######################################################################################
#######################################################################################
#######################################################################################
    
def traverse(root, data, branchorder, nonbinary):

    # called from hocreader

    # Driver for traversal through hoc trees, which are binary after a certain point.
    # This function weeds through the "un-binary" part of the tree and calls
    # traverse_binarytree to go through the binary parts.
    #
    # root = segments called "Primagenitors", the last branches of the un-binary part of the
    # tree.  They are the sons of the "root" segment found in hocreader ("starter_segment"s)
    # before calling traverse.
    # It is assumed that from the primagenitor down, the tree is strictly binary.
    #
    # data = celldata dictionary built in hocreader.
    #
    # branchorder = input empty dictionary initialized in hocreader:
    #               branchorder = {segment: 0. bo (branch order)
    #                                       1. bl (branchlength)
    # returns branchorder

    # Loop through the sons of root segment
    # call traverse_binarytree for each son
    # re call this routine if root happens to have more than 2 sons (make sonofroot
    # into "surrogate" father).
        
    import re

    # Set starter segment (root) branch order to zero:

    # In case coming back through via *******[&&&]********
    # (root has more than two sons: nonbinary = True)
    if (nonbinary):
        bo = branchorder[root][0]
        bl = branchorder[root][1]
        branchorder[root] = [bo, bl]
    else:
        branchorder[root] = [0, 0]

    ##########################################
    for sonofroot in data[root][0]:
    ##########################################


        # ignore axons (go to end of loop)
        # print "In traverse: sonofroot = ", sonofroot
        if re.search("axon", sonofroot):
            # print "current son of root is ", sonofroot 
            continue 

        # In case coming back through via *******[&&&]********
        # (root has more than two sons: nonbinary = True)
        if (nonbinary):
            bo = branchorder[root][0]
            bl = branchorder[root][1]
        else:
            bo = 0
            bl = 0
          
        # if "soma" in the label, keep bo = 0 and go back through this fcn
        # with root = sonofroot (Cuntz cells)
        if re.search("soma", sonofroot):

            branchorder[sonofroot] = [bo, bl]

            if (printbotrav):
                print
                print "    *** Found secondary root: ", sonofroot, "recursing through traverse ... ***"
                print

            # if sonofroot has sons, recurse through traverse
            if data[sonofroot][0] != []:
            
                traverse(sonofroot, data, branchorder, nonbinary)

            # "soma" that is an ending (leaf), go to top of loop (next sonofroot)
            else:

                continue

        # if not a soma segment, increment bo (Cuntz cells, some D2 registered cells)
        # and set branchorder for segment:
        if re.search("soma", sonofroot) == None:

            bo += 1
            bl = bl + data[sonofroot][3]
            branchorder[sonofroot] = [bo, bl]
            
            # if sonofroot has sons,
            # traverse binary tree starting with one son of sonofroot (Primagenitor)
            if data[sonofroot][4] != 'Ending':

                if (printbotrav):
                    print
                    print "In traverse nonEnding segment:"
                    print "************* Processing Primagenitor         =", sonofroot
                    print "************* Sons of Primagenitor =", data[sonofroot][0]
                    print

                # Usually there will just be two branches from soma
                if len(data[sonofroot][0]) <= 2:
                    
                    # nonbinary = False
                    if (printbotrav):
                        print
                        print "              Two or fewer branches for ", sonofroot
                        print "              Calling traverse_binarytree with seed ", sonofroot
                        print "              branchorder[", sonofroot, "] =", branchorder[sonofroot] 
                    
                    #________________________________________________________________________________
                    branchorder = traverse_binarytree(sonofroot, None, data, branchorder, sonofroot)
                    #________________________________________________________________________________

                # But in some cases there are more than two branches from soma (e.g. cell 217)
                # In which case recurse through here again
                elif len(data[sonofroot][0]) > 2:

                    nonbinary = True
                    
                    bo = branchorder[root][0]
                    bl = branchorder[root][1]
                    branchorder[root] = [bo, bl]
                    
                    if (printbotrav):
                        print
                        print "              Two or more branches for ", sonofroot, " setting nonbinary = ", nonbinary
                        print "              Calling traverse with root ", sonofroot, "*************"
                        print "              branchorder[", sonofroot, "] =", branchorder[sonofroot] 

                    #________________________________________________________________________________
                    traverse(sonofroot, data, branchorder, nonbinary) # *******[&&&]********
                    #________________________________________________________________________________

            # if an ending (leaf), go to top of loop (next sonofroot):
            else:

                if (printbotrav):
                    print
                    print "In traverse Ending segment:"
                    print "************* Processing Primagenitor         =", sonofroot
                    print "************* Sons of Primagenitor =", data[sonofroot][0], "Duh!"
                    print

                continue                
                
    # print "length of branchorder = ", len(branchorder)
    # print "branchorder = ", branchorder
    ##########################################
    # End: for sonofroot in data[root][0]:
    ##########################################
    return branchorder
        
        
#######################################################################################
#######################################################################################
#######################################################################################

def traverse_binarytree(seed, seedsfather, data, branchorder, primagenitor):
    
    # Move between currentseg and its sons and father
    # keeping track of branchorder [bo, bl] as you go.

    # "up" is towards the soma
    # "down" is away from the soma
    
    if(printbo):
        print "seed = ", seed
        print "called with seedsfather = ", seedsfather
        
    # If called with father, set current seg to seed
    # to be able to recurse through with currentseg
    
    if seedsfather != None:
        
        if (printbo):
            print "Called with seed's father (processing using seed instead of seed's son)"
        
        for son in data[seedsfather][0]:
            if son == seed:
                currentseg = seed
            else:
                siblingseg = son

        # Check to see if only one son (no sibling):
        if len(data[seedsfather][0]) == 1:
            siblingseg = None
                            
    # else called with seedsfather = None: set current seg to 0th son of seed
    else:

        currentseg = data[seed][0][0]

        # Check for only one son (no sibling)
        if len(data[seed][0]) == 1:
            siblingseg = None

        else:
            siblingseg = data[seed][0][1]
        
##         try:
##             siblingseg = data[seed][0][1]

##         # If sibling does not exist:
##         except (IndexError):
##             siblingseg = None
        
        if (printbo):
            print "Called w/o seed's father (processing using seed's son)"
            print " seed       = ", seed
            print " currentseg = ", currentseg
            print " siblingseg = ", siblingseg
        
    # Set current quantities.  bo and bl are set to the father of currentseg
    bo = branchorder[data[currentseg][1]][0]
    bl = branchorder[data[currentseg][1]][1]

    currentsegsons = data[currentseg][0]
    # siblingsegsons = data[siblingseg][0]
    currentsegfather = data[currentseg][1]
    
    currentseglength = data[currentseg][3]
    # siblingseglength = data[siblingseg][3]
    ########################################################################
    # Cases for no siblings (ran into with marcel axons in "3D" hoc files)
    if siblingseg == None:
    ########################################################################

        if (printbo_short):
            print
            print "*************** Segment with no sibling! *******************"
            print "siblingseg =", siblingseg, "for currentseg = ", currentseg
            print
        
        # IN
        # currentseg is 'Intersegment' and no sibling:
        if data[currentseg][4] == 'Intersegment':
            
            if (printbo):
                print
                print " Intersegment-None (IN) Case:"
                print
            
            # Do not increment bo  
            # Only bl for currentseg and sibling:
            # bo += 1
        
            blcs = bl + currentseglength 
            branchorder[currentseg] = [bo, blcs]
        
            # Check to see if both sons are endings:
            for son in currentsegsons:            
            
                # if at least one 'Intersegment', break out of loop and skip else
                if data[son][4] == 'Intersegment':
                    break

                # Else both sons 'Ending': increment leaves and go back up
                else: 

                    if (printbo):
                        print "Sons of intersegment are both endings: currentseg = ", currentseg, data[currentseg][4]
                        print "Currentseg sons = ", son
                        print "Incrementing leaves!"
                        
                        print "Reached two endings (from 'Intersegment' branch)"
                        print "Proceeding back 'up' this blasted tree, which is really down towards the root ..."
                
                    bo = branchorder[currentseg][0] + 1
                    bl = bl + branchorder[currentseg][1]
                    branchorder[son] = [bo, bl]
                
                    #________________________________________________________________________________
                    traverse_binarytree_up(son, data, branchorder, primagenitor)
                    #________________________________________________________________________________

            # At least one son 'Intersegment'.  Find one of them and continue traversing.
            # If fell through from else, neither son "Intersegment' so nothing happens
            for son in currentsegsons:            
                
                if data[son][4] == 'Intersegment':
                    #________________________________________________________________________________
                    traverse_binarytree(son, currentseg, data, branchorder, primagenitor)
                    #________________________________________________________________________________
            
        # EN
        # If currentseg is 'Ending' and no sibling:
        if data[currentseg][4] == 'Ending':

            if (printbo):
                print
                print " Ending-Ending (EE) Case:"
                print
        
            bo += 1
        
            blcs = bl + currentseglength 
            branchorder[currentseg] = [bo, blcs]
            
            if (printbo):
                print "Reached two endings (leaves, Anna Livia!):"
                print "           currentseg = ", currentseg, data[currentseg][4]
                print "           primogenitor = ", primagenitor
                print "Proceeding back 'up' this blasted tree, which is really down towards the root ..."
            
            #________________________________________________________________________________
            traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
            #________________________________________________________________________________
            
            if (printbo):
                print
                print "          Done with Primagenitor (leaves, Anna Livia!)", primagenitor,  ", branch",\
                    data[primagenitor][0][0]       #, ", moving on to other branch", \
                    #data[primagenitor][0][1]
                    # commented out since no sibling in this case (crashing for axons)


                print "Branch order of", currentseg, " = ", branchorder[currentseg]
                
    ########################################################################
    # Cases with siblings (should be all cells except certain "marcel" "3D" cells)
    else: 
    ########################################################################
    
        # currentsegsons = data[currentseg][0]
        siblingsegsons = data[siblingseg][0]
        # currentsegfather = data[currentseg][1]
        
        # currentseglength = data[currentseg][3]
        siblingseglength = data[siblingseg][3]
        
        if (printbo):
            print
            print "          In traverse_binarytree:"
            print "          bo =", bo
            print "          bl =", bl
            # print "          blss =", blss
            # print "          blcs =", blcs
            print "          seed       = ", seed
            print "          currentseg = ", currentseg 
            print "             siblingseg       =", siblingseg 
            print "             currentsegsons   =", currentsegsons 
            print "             currentsegfather =", currentsegfather 
            print "          Type current seg = ", data[currentseg][4]
            print
            
            print "Currseg:", currentseg
            print "In binary tree, leaf's father's father is ", data[currentseg][1]
            print

        # II
        # currentseg is 'Intersegment' and sons are not both Endings, move to next 'Intersegment':
        if data[currentseg][4] == 'Intersegment':

            if (printbo):
                print
                print " Intersegment-Intersegment (II) Case:"
                print
            
            # Increment bo / bl for currentseg and sibling:
            bo += 1
            
            blcs = bl + currentseglength 
            branchorder[currentseg] = [bo, blcs]
            
            blss = bl + siblingseglength 
            branchorder[siblingseg] = [bo, blss]
            
            # Check to see if both sons are endings:
            for son in currentsegsons:            
            
                # if at least one 'Intersegment', break out of loop and skip else
                if data[son][4] == 'Intersegment':
                    break

                # Else both sons 'Ending': increment leaves and go back up
                else: 
                    
                    bo = branchorder[currentseg][0] + 1
                    bl = bl + branchorder[currentseg][1]
                    branchorder[son] = [bo, bl]
                    
                    if (printbo):
                        for segment in branchorder.keys():
                            print "branchorder[", segment, "] =", branchorder[segment][0]
                       
                    if (printbo):
                        print "Sons of intersegment are both endings: currentseg = ", currentseg, data[currentseg][4]
                        print "                                       siblingseg = ", siblingseg, data[siblingseg][4]
                        print "Currentseg sons = ", son
                        print "Incrementing leaves!"
                        
                        print "Reached two endings (from 'Intersegment' branch)"
                        print "Proceeding back 'up' this blasted tree"
                        print "(which is really down towards the root ...)"
                        
                    #________________________________________________________________________________
                    traverse_binarytree_up(son, data, branchorder, primagenitor)
                    #________________________________________________________________________________

            # At least one son 'Intersegment'.  Find it and continue traversing.
            # If fell through from else, neither son "Intersegment' so nothing happens
            for son in currentsegsons:            
                
                if data[son][4] == 'Intersegment':
                    #________________________________________________________________________________
                    traverse_binarytree(son, currentseg, data, branchorder, primagenitor)
                    #________________________________________________________________________________
            
        # EI
        # currentseg is 'Ending', but sibling is 'Intersegment'
        # move on to sons of sibling:
        if data[currentseg][4] == 'Ending' and data[siblingseg][4] == 'Intersegment':
            
            if (printbo):
                print
                print " Ending-Intersegment (EI) Case:"
                print
            
            # Increment bo for currentseg and sibling:
            bo += 1
        
            blcs = bl + currentseglength 
            branchorder[currentseg] = [bo, blcs]
            
            blss = bl + siblingseglength 
            branchorder[siblingseg] = [bo, blss]
            
            #________________________________________________________________________________
            traverse_binarytree(siblingseg, None, data, branchorder, primagenitor)
            #________________________________________________________________________________
                
        # EE
        # If both currentseg and siblingseg are 'Endings':
        # reached the end of the line:
        # go back through tree checking for undone branches
        if data[currentseg][4] == 'Ending' and data[siblingseg][4] == 'Ending':

            if (printbo):
                print
                print " Ending-Ending (EE) Case:"
                print
        
            bo += 1
        
            blcs = bl + currentseglength 
            branchorder[currentseg] = [bo, blcs]
            
            blss = bl + siblingseglength 
            branchorder[siblingseg] = [bo, blss]

                        
            if (printbo):
                print "Reached two endings: currentseg = ", currentseg, data[currentseg][4]
                print "                                     siblingseg = ", siblingseg, data[siblingseg][4]
                print "Proceeding back 'up' this blasted tree"
                print "(which is really down towards the root) ..."
                
            #________________________________________________________________________________
            traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
            #________________________________________________________________________________

            if (printbo):
                print "          Done with Primagenitor ", primagenitor,  ", branch",\
                      data[primagenitor][0][0], ", moving on to other branch", \
                      data[primagenitor][0][1]
                

                print "Branch order of", currentseg, " = ", branchorder[currentseg]

    return branchorder
        
#######################################################################################
#######################################################################################
#######################################################################################

def traverse_binarytree_up(leaf, data, branchorder, primagenitor):

    # global printbo

    # Called from two places in traverse_binarytree

    # Move up from the leaf of a binary tree to father segment
    # and look for an untraversed branch in sibling's sons

    #  set currentseg to father of leaf (data[leaf][1]) (already traversed), set bo / bl:

    currentseg = data[leaf][1]
    bo = branchorder[currentseg][0]
    bl = branchorder[currentseg][1]

    # If at primagenitor, done with whole branch, return:

    if currentseg == primagenitor: 
        if (printbo): print "traverse_binarytree_up: current segment is primagenitor: ", currentseg
        return branchorder
    
    # Find currentseg's sibling from father:

    currentsegfather = data[currentseg][1]

    # print some stuff:
    if (printbo):
        print
        print "                 *** Arriving safely at traverse_binarytree_up! *** "
        print "                      ... Oh, hello Mr. Tyler, going ... Up?"
        print "                      primagenitor = ", primagenitor
        print "                      Leaf = ", leaf
        print "                      currentseg (leaf's father) = ", currentseg
        print "                      currentsegfather = ", currentsegfather
        print "                      branch length of current segment = ", branchorder[currentseg][1]
        print
        
    # Check to see if there are siblings
    if len(data[currentsegfather][0]) > 1:
        
        for son in data[currentsegfather][0]:
            if son != currentseg:
                siblingseg = son

    # currentseg has no sibling, continue up
    else:
                                
        if (printbo):
            print
            print " ********* Aha! No siblings, continuing up ..."
            print "Proceeding back 'up' this blasted tree"
            print "(which is really down towards the root ..."
            print "Calling traverse_binarytree_up with leaf = ", currentseg
            print "                        for primagenitor = ", primagenitor
            print
            
        #________________________________________________________________________________
        traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
        #________________________________________________________________________________
            
        siblingseg = None
        # return


    if siblingseg != None:
        
        if (printbo):
            print
            print "        Running sibling branch ... "
            print "siblingseg = ", siblingseg
            print
        
        # Check if sibling is an ending, if so, go up (we know it's been traversed)

        if data[siblingseg][4] == 'Ending':
            if (printbo):
                print "Looks like we got an ending fellas!"
                print "                           ... Headin back up!"
            #________________________________________________________________________________
            traverse_binarytree_up(siblingseg, data, branchorder, primagenitor)
            #________________________________________________________________________________

        else:

            # Find sons of siblingseg:

            siblingsegsons = data[siblingseg][0]
            if (printbo): print " siblingsegsons = ", siblingsegsons
        
            # Check if sibling seg's sons have been traversed yet.
        
            sonflag = 0
            for son in siblingsegsons:

                if (printbo):
                    print " *********** Sibling son loop: Son         = ", son
                    print " *********** Sibling son loop: BO of Son   = ", branchorder[son]
                    print " *********** Sibling son loop: Type of Son = ", data[son][4]
        
                # if traversed, keep going up
                if branchorder[son][0] != -1:
                    sonflag = 0
                
                # if not traversed, and an ending, increment, and go up
                if branchorder[son][0] == -1 and data[son][4] == 'Ending':
                    sonflag = 1
                    
                # if not traversed, and intersegment, go down
                if branchorder[son][0] == -1 and data[son][4] == 'Intersegment':
                    sonflag = 2

                # Based on above loop, perform the following action:
                
                if (sonflag == 0):
                    if (printbo):
                        print "sonflag = ", sonflag
                        print "Keep going up! current siblingson = ", son
                        print "branchorder[ of son:", son, "][0] = ", branchorder[son][0]
                        print "Leaf = ", leaf
                        print "Leaf's father = ", currentseg 
                        print "Calling traverse_binarytree_up with current leaf's father = ", currentseg
                    #________________________________________________________________________________
                    traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
                    #________________________________________________________________________________
                    # Need break here for "marcel" "3D" hoc files
                    # Not quite sure why
                    break

                if (sonflag == 1):
                    if (printbo):
                        print "sonflag = ", sonflag
                        print "Oops, increment and go back up!  Segment = ", son
                        print "                       Calling with seed = ", currentseg
                    bl = branchorder[siblingseg][1] + data[son][3]
                    branchorder[son] = [bo+1, bl]
                    #________________________________________________________________________________
                    traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
                    #________________________________________________________________________________

                if (sonflag == 2):
                    if (printbo):
                        print "sonflag = ", sonflag
                        print "Oops, go back down!  calling traverse with siblingseg = ", siblingseg
                    #________________________________________________________________________________
                    traverse_binarytree(siblingseg, None, data, branchorder, primagenitor)                  
                    #________________________________________________________________________________
    
##             for son in siblingsegsons:

##                 if (printbo): print " *********** Sibling son loop: Son = ", son
        
##                 # if traversed, keep going up
##                 if branchorder[son][0] != -1:
##                     sonflag = 0
##                     print "sonflag = ", sonflag
##                     if (printbo):
##                         print "Keep going up! current siblingson = ", son
##                         print "branchorder[ of son:", son, "][0] = ", branchorder[son][0]
##                         print "Leaf = ", leaf
##                         print "Leaf's father = ", currentseg 
##                         print "Calling traverse_binarytree_up with current leaf's father = ", currentseg
##                     traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
##                     break
##                     # continue
                
##                 # if not traversed, and an ending, increment, and go up
##                 if branchorder[son][0] == -1 and data[son][4] == 'Ending':
##                     sonflag = 1
##                     print "sonflag = ", sonflag
##                     if (printbo):
##                         print "Oops, increment and go back up!  Segment = ", son
##                         print "                       Calling with seed = ", currentseg
##                     bl = branchorder[siblingseg][1] + data[son][3]
##                     branchorder[son] = [bo+1, bl]
##                     traverse_binarytree_up(currentseg, data, branchorder, primagenitor)
##                     #continue
                    
##                 # if not traversed, and intersegment, go down
##                 if branchorder[son][0] == -1 and data[son][4] == 'Intersegment':
##                     sonflag = 2
##                     print "sonflag = ", sonflag
##                     if (printbo):
##                         print "Oops, go back down!  Segment = ", son
##                     traverse_binarytree(siblingseg, None, data, branchorder, primagenitor)
##                     #continue                    
    
#######################################################################################
#######################################################################################
#######################################################################################

# Calculate mean and standard deviation of data x[]:
#    mean = {\sum_i x_i \over n}
#    std = sqrt(\sum_i (x_i - mean)^2 \over n-1)

def meanstdv(x):

    from math import sqrt
    n, mean, std = len(x), 0, 0

    if n > 1:
        for a in x:
            mean = mean + a
            
        mean = mean / float(n)

        for a in x:
            std = std + (a - mean)**2

        std = sqrt(std / float(n-1))

    else:

        mean = 0
        std = 0

    return mean, std
#######################################################################################
#######################################################################################
#######################################################################################

# Function to determine zScalingVal based on zLUT table in CortexDBA

def zLUT(ScalingDict):

    zDepth, zScalingVal = 0, 0

    for zDepth in range(200000):

        zScalingVal = -9903.83 + zDepth*(11071.9 - (-9903.83))/200000
        # Have to make this a string to get the
        # dictionary lookup to work in hocreader (see [***]
        ScalingDict[str(zDepth)] = zScalingVal
    
    return ScalingDict


#######################################################################################
#######################################################################################
#######################################################################################

# Function to write results into Database tables

def writeCoreFeatures(FeatureDictionary, Database, Username, UserPassword, HostofHosts):

    # FeatureDictionary is a dictionary defined in hocreader with all the relevant features

    import MySQLdb

    # http://www.kitebird.com/articles/pydbapi.html
    # http://mysql-python.sourceforge.net/MySQLdb.html#mysqldb

    #################################################################################
    # Some code to sort a dictionary, kept here for historical interest,
    # remnants from a war waged before I got my head out of my ass and 
    # figured out how to use a ductionary (see below).
    # import operator
    # x = FeatureDictionary
    # sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1))
    # sorted_x will be a list of tuples sorted by the second element in each tuple. 
    # dict(sorted_x) == x
    #
    # E.g.
    #     sorted_x[0][0]=FeatureName[0]
    #     sorted_x[0][1]=FeatureValue[0]
    #
    # print "sorted_x = ", sorted_x
    # print "sorted_x[0][0] = ", sorted_x[0][0]
    # print "sorted_x[0][1] = ", sorted_x[0][0]
    # print 
    #################################################################################

    CellID = FeatureDictionary["Experiment_ID"]

    # Would like to figure out how to have the table name generated dynamically some time
    # table = "corefeatures"

    print "            ********************* In writeCoreFeatures: *************************"
    # print "the DB = ", Database
    # print "user = ", Username
    # print "password = ", UserPassword
    # print "host = ", HostofHosts
    # print "CellID = ", CellID
    # print
    # for featurename in FeatureDictionary.keys():
    #     print featurename, "= ", FeatureDictionary[featurename]

    # Note: Had much trouble with unix socket.
    # Kept getting error:
    # _mysql_exceptions.OperationalError: 
    # (2002, "Can't connect to local MySQL server through socket '/var/lib/mysql/mysql.sock' (2)")
    # To find the unix socket go to a mysql command line and type:
    #     show variables like 'socket';
    
    dB = MySQLdb.connect(
            # host = "10.1.100.47", \
            host = HostofHosts, \
            port = 3306, \
            db = Database, \
            passwd = UserPassword, \
            user = Username,
            # unix_socket = "/opt/lampp/var/mysql/mysql.sock"                       # Linux (Neuromorph-3d) socket
            # unix_socket = "/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock"   # Mac socket
            unix_socket = "/var/lib/mysql/mysql.sock"                               # BB3D socket
        )
        
    try:
        c = dB.cursor()
        # c = dm3.cursor()
        # c = dm3.cursor(MySQLdb.cursors.DictCursor)
        
        c.execute ("SELECT VERSION()")
        row = c.fetchone()
        print
        print "... Proceding to database upload!"
        print "Server MySQL version:", row[0]
        print
        
        # Delete previous entry (if exists)
        # For dm3
        # c.execute("""DELETE FROM `CoreFeatures_python` WHERE `Experiment_ID` = %s""", (CellID,))
        # c.execute("""DELETE FROM `CoreFeatures` WHERE `Experiment_ID` = %s""", (CellID,))
        c.execute("""DELETE FROM `corefeatures` WHERE `Experiment_ID` = %s""", (CellID,))

        # print "FeatureDictionary.keys() =", FeatureDictionary.keys()

        # d = {"spam": "1", etc,etc,"egg": "2"}
        #cols = FeatureDictionary.keys()
        #vals = FeatureDictionary.values()
        #
        #stmt = "INSERT INTO features (%s) VALUES(%s)" % (",".join(cols), ",".join(vals))
        #
        #print "stmt =", stmt
        #
        #c.execute(stmt)

        d = FeatureDictionary
        cols = d.keys()
        vals = d.values()
        sql = "INSERT INTO corefeatures (%s) VALUES (%s)" % (
            ",".join(cols), ",".join(["%s"]*len(vals))
            )
        # sql = 'UPDATE features SET {}'.format(', '.join('{}=%s'.format(k) for k in d.keys()))
        # print sql
        # print vals

        c.execute(sql, vals) 
                                                         
        c.execute("""SELECT 
                    Experiment_ID,\
                    Dendritic_Length, \
                    Dendritic_Volume \
                 FROM corefeatures WHERE Experiment_ID = %s""", (CellID,))
                 # FROM CoreFeatures_python WHERE NeuronID = %s""", (CellID,))


        rows = c.fetchall()
        for row in rows:
            #print "%s, %s, %s, %s" % (row[0], row[1], row[2], row[3])
            print "Number of rows returned : %d" % c.rowcount
        featurewrite = True
            
    except MySQLdb.Error, e:
        print "functions.py (writeCoreFeatures): MySQL Error %d: %s" % (e.args[0], e.args[1])
        featurewrite = False
        
        
    c.close()
    dB.commit()
    dB.close()
    # dm3.commit()
    # dm3.close()

    return featurewrite

#######################################################################################
#######################################################################################
#######################################################################################

# Function to input lots of hocs for hardwire testing

def inputHocs():

    # glob module for retrieving all the files from a directory
    import glob

    inputHocList = ()

    inputHocList = glob.glob("/nas1/Data_Mythreya/MotorCortexProject/V8/Registration_Local/Input_Spatial_Graphs/All_Hocs/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/Sub_Cortical_HocsnLogs/4Andy2/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/Sub_Cortical_HocsnLogs/4Andy/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/Rajeev/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/invivoTracings/PC_Registered_Hocs_from_Robert/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/invivoTracings/PCRegistered_Hocs/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/invivoTracings/D2Registered_BBUpdate_Andy/HocsNLogs/JustHocs/*.hoc")
    # inputHocList = glob.glob("/home/johnsona/invivoTracings/D2Registered_BBUpdate_Andy/HocsNLogs/*.hoc")
    # inputHocList = glob.glob('/home/johnsona/invivoTracings/hocs/*.hoc')
    # inputHocList = glob.glob('/home/johnsona/invivoTracings/hoc_features/*.hoc')
    # inputHocList = glob.glob('/home/johnsona/invivoTracings/BB3D_D2_hocs/unique/*.hoc')
    # inputHocList = glob.glob('/home/johnsona/invivoTracings/test/*.hoc')
    # inputHocList = glob.glob('/disk1/Share/AllCells4Database/NeuroConv_Hoc/*.hoc')
    # inputHocList = glob.glob('/disk1/Share/AllCells4Database/D2_registered_Hoc/*.hoc')

    # inputHocList = ("/disk1/Share/AllCells4Database/D2_registered_Hoc/30_CDK20050414nr6L23_registered_rob_neuron_transform_registered_D2.hoc",
    #                  "/disk1/Share/AllCells4Database/D2_registered_Hoc/04_RMB20040916_nr21_registered_rob_neuron_transform_registered_D2.hoc",
    #                  "/disk1/Share/AllCells4Database/D2_registered_Hoc/06_RMB_nr73_registered_rob_neuron_transform_registered_D2.hoc"
    #                  )


    return inputHocList

#######################################################################################
#######################################################################################
#######################################################################################

# Function to check if string is number

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#######################################################################################
#######################################################################################
#######################################################################################

def rotate(angle, x, y, z):

# Robert's "minimal torsion" rotation algorithm
# Compute a rotation matrix based on a rotation axis and angle
#- Let n be the unit vector for the orientation at the location of the neuron, and z the unit vector pointing in the z-direction
#- compute the rotation axis vector r = n x z , where x is the vector cross-product (order matters!)
#- compute the rotation angle alpha = acos(n * z) , where * is the dot product of the two vectors
#- compute the rotation matrix for rotating points by angle alpha around rotation axis r

    
    import math
    
    # make a normalized quaternion
    w = math.cos(0.5*angle);
    f = math.sin(0.5*angle)/math.sqrt(x*x+y*y+z*z);
    x *= f;
    y *= f;
    z *= f;
    
    # // convert the quaternion to a matrix
    
    ww = w*w;
    wx = w*x;
    wy = w*y;
    wz = w*z;
    
    xx = x*x;
    yy = y*y;
    zz = z*z;
    
    xy = x*y;
    xz = x*z;
    yz = y*z;
    
    s = ww - xx - yy - zz;
    
    m00 = xx*2 + s;
    m10 = (xy + wz)*2;
    m20 = (xz - wy)*2;
    
    m01 = (xy - wz)*2;
    m11 = yy*2 + s;
    m21 = (yz + wx)*2;
    
    m02 = (xz + wy)*2;
    m12 = (yz - wx)*2;
    m22 = zz*2 + s;
    
    # matrix = [[m00, m10, m20], [m01, m11, m21], [m02, m12, m22]]
    matrix = [[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]]
    
    # this->Concatenate(*matrix);
    return matrix

#######################################################################################
#######################################################################################
#######################################################################################

def matrix_multiply(x, y, z, array):

    rotated_x = array[0][0]*x + array[0][1]*y + array[0][2]*z
    rotated_y = array[1][0]*x + array[1][1]*y + array[1][2]*z
    rotated_z = array[2][0]*x + array[2][1]*y + array[2][2]*z

    return rotated_x, rotated_y, rotated_z
