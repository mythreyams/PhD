# Converts output csv file from hocreader.py to csv file that only contains Axonal features, no dendrite features
import sys
import csv

featureNames = ["Experiment_ID","Axon_Branchpoints","Axon_Length","Axon_Volume"]

def main(inputfilenameFeatures,outputfilename):

    strHeader = ()
    strFeatures = ()

    # Read in features from hocreader.py and only write Axon features to csv file
    with open(inputfilenameFeatures,'r') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            for f in featureNames:
                strHeader = strHeader + (f,)
                strFeatures = strFeatures + (row[f],)

    out = csv.writer(open(outputfilename, "w+"), delimiter=',', quoting=csv.QUOTE_ALL)
    out.writerow(strHeader)
    out.writerow(strFeatures)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        inputfilenameFeatures = sys.argv[1]
        outputfilename = sys.argv[2]

        main(inputfilenameFeatures, outputfilename)
    elif len(sys.argv) == 2:
        inputfilenameFeatures = sys.argv[1]
        outputfilename = inputfilenameFeatures

        main(inputfilenameFeatures, outputfilename)
    else:
        print "USAGE: python convertFeaturesToAxonFeatures.py [inputfilenameFeatures] [outputfilename](optional, otherwise outputfilename=inputfilename)"