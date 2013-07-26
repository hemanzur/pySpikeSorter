##########################################################################################
            
def KlustaKwik_call(data, minClust = 2, maxClust = 5):
    ''' data must be an array of observations x dimensions'''
    
    # create a text file with the data. The first line must be the
    # number of dimensions of the data
    f = open('data.fet.1', 'w')
    f.write('%d\n' % data.shape[1])
    for k in data:
        for j in k:
            f.write('%f ' % j)
        f.write('\n')
    f.close()
    
    # call klustakwick with the data 
    import os
    os.system('KlustaKwik data 1 -MinClusters %d -MaxClusters %d' % (minClust, maxClust))

    # read the results
    f = open('data.clu.1','r')
    clusterData = f.readlines()
    f.close()
    clusterData = [int(re.search('[0-9]{1,2}', k).group()) for k in clusterData]

    # the first line is the number of clusters
    nClusters = clusterData[0]
    clusterData.pop(0)
    clusterData = np.array(clusterData)

    # create an array with the indices of each cluster
    clustIndx = []
    for k in range(1, nClusters+1):
        clustIndx.append(np.flatnonzero(clusterData==k))

    return clustIndx

##### AUTOMATIC CLUSTERING WIDGET #####

class AutoClustWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Automatic Clustering')
        vlay = QtGui.QVBoxLayout(self)
        hlay = QtGui.QHBoxLayout()
        self.MinClust = QtGui.QSpinBox()
        self.MinClust.setRange(1,5)
        self.MinClust.setValue(2)
        hlay.addWidget(QtGui.QLabel('Min'))
        hlay.addWidget(self.MinClust)
        self.MaxClust = QtGui.QSpinBox()
        self.MaxClust.setRange(2,10)
        self.MaxClust.setValue(4)
        hlay.addWidget(QtGui.QLabel('Max'))
        hlay.addWidget(self.MaxClust)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        self.ClusteringMethod = QtGui.QComboBox(self)
        self.ClusteringMethod.addItems(['KlustaKwik','KMeans', 'Afinity Propagation',
                                        'Mean-shift', 'Spectral Clustering', 'Hierarchical clustering',
                                        'DBSCAN', 'Gaussian Mixtures'])
        hlay.addWidget(QtGui.QLabel('Method'))
        hlay.addWidget(self.ClusteringMethod)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        AutoClustBtn = QtGui.QPushButton('Preview')
        #AutoClustBtn.clicked.connect(self.AutoClust_proc)
        ClusteringCommitBtn = QtGui.QPushButton('Commit')
        hlay.addWidget(AutoClustBtn)
        hlay.addWidget(ClusteringCommitBtn)
        vlay.addLayout(hlay)

        self.setLayout(vlay)

autoclust = AutoClustWidget()