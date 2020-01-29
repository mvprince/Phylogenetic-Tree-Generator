# Maia Prince 2/17/2019
# Converts Phylip-formatted data to a Newick-formatted evolutionary tree.
# Uses neighbor joining method; includes branch lengths.


# Read in data file
def getData(fileName):
    distances = [] # Temporary list for holding cluster distances
    
    infile = open(fileName, 'r')
    
    # Get number of sequences from 1st line
    s = infile.readline()
    global numSeq
    numSeq = int(s)
    
    # Clean up formatting; store cluster names & distances between clusters
    for line in infile:
        line = line.replace('\n', '')
        lineList = line.split(' ')
        clusterNames.append(lineList[0])
        del lineList[0]
        distances.append(lineList)
    infile.close()
    
    # Store file data in nested dict
    for i in range (0,numSeq):
        clusterDist[clusterNames[i]] = {}
        for j in range(0,numSeq):
            clusterDist[clusterNames[i]][clusterNames[j]] = float(distances[i][j])
# end getData


numSeq = 0 # Number of sequences
clusterDist = {} # Distances between clusters
newick = {} # Newick-formatted strings for each currently existing cluster
clusterNames = [] # All current cluster names

print ('Enter the file name of the Phylip-formatted data (no .txt extension).')
file = input('') + '.txt'
getData(file)

# initialize Newick dict
for name in clusterNames:
    newick[name] = name

# do clustering
while (len(clusterNames) > 2):
    # generate transformed r-values for each cluster to account for different rates of evolution
    transformedDists = {}
    klist = list(clusterDist.keys())
    for k1 in klist:
        r = 0.0
        k2list = list(clusterDist[k1])
        for k2 in k2list:
            r += clusterDist[k1][k2]
        r /= (len(clusterNames)-2)
        transformedDists[k1] = r
    
    # find minimum transition value (amount of evolutionary distance)
    shortestD = float("inf")
    cluster1 = ''
    cluster2 = ''
    klist = list(clusterDist.keys())
    for k1 in klist:
        k2list = list(clusterDist[k1])
        for k2 in k2list:
            # calculate transition value, check if shortest
            tv = clusterDist[k1][k2] - transformedDists[k1] - transformedDists[k2]
            if (tv < shortestD) and (k1 != k2):
                shortestD = tv
                cluster1 = k1
                cluster2 = k2

    # start merging closest clusters
    newClusterName = cluster1 + '-' + cluster2
    
    # do neighbor joining to make new dictionary with key newClusterName
    klist = list(clusterDist.keys())
    clusterDist[newClusterName] = {}
    for k in klist:
        newDist = clusterDist[cluster1][k]
        newDist += clusterDist[cluster2][k]
        newDist -= clusterDist[cluster1][cluster2]
        newDist /= 2
        clusterDist[newClusterName][k] = newDist
    clusterDist[newClusterName][newClusterName] = 0
    
    # update other clusters in nested dictionary to have a distance value for new cluster
    klist = list(clusterDist.keys())
    for k in klist:
        clusterDist[k][newClusterName] = clusterDist[newClusterName][k]
    
    # update Newick dict, including branch lengths
    bLen1 = clusterDist[cluster1][cluster2] + transformedDists[cluster1] - transformedDists[cluster2]
    bLen1 /= 2
    bLen2 = clusterDist[cluster2][cluster1] + transformedDists[cluster2] - transformedDists[cluster1]
    bLen2 /= 2
    newClusterStr = '(' + newick[cluster1] + ':' + str(bLen1) + ','
    newClusterStr += newick[cluster2] + ':' + str(bLen2) + ')'
    newick[newClusterName] = newClusterStr
    
    #Remove old clusters from Newick dict
    del newick[cluster1]
    del newick[cluster2]
    
    # delete all entries involving the two merged clusters in nested dict of distances
    klist = list(clusterDist.keys())
    for k in klist:
        del clusterDist[k][cluster1]
        del clusterDist[k][cluster2]
    del clusterDist[cluster1]
    del clusterDist[cluster2]
    
    # update list of cluster names
    clusterNames.remove(cluster1)
    clusterNames.remove(cluster2)
    clusterNames.append(newClusterName)
# end while; only 2 clusters now remain

print('\nNewick format:')
print('(' + newick[clusterNames[0]] + ',' + newick[clusterNames[1]] + ')' )
 