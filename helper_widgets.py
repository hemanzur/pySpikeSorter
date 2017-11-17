from PyQt4 import QtGui
import guidata
import os
import re
import numpy as np
app = guidata.qapplication()
import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di
from matplotlib.pyplot import cm

colormaps = [k for k in cm.datad.keys() if not k.endswith('_r')]


#==============================================================================
class MergeUnitsWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Merge Units')
        hlay = QtGui.QHBoxLayout()
        self.list1 = QtGui.QListWidget()
        self.list1.setMaximumWidth(100)
        self.list2 = QtGui.QListWidget()
        self.list2.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-next')
        btnRight = QtGui.QPushButton(icon, '')
        btnRight.clicked.connect(self.toRight)
        icon = QtGui.QIcon.fromTheme('go-previous')
        btnLeft = QtGui.QPushButton(icon, '')
        btnLeft.clicked.connect(self.toLeft)
        vlay = QtGui.QVBoxLayout()
        hlay.addWidget(self.list1)
        vlay.addWidget(btnRight)
        vlay.addWidget(btnLeft)
        hlay.addLayout(vlay)
        hlay.addWidget(self.list2)

        vlay = QtGui.QVBoxLayout()
        vlay.addLayout(hlay)

        self.CancelBtn = QtGui.QPushButton('Cancel')
        self.CancelBtn.clicked.connect(self.close)

        self.AcceptBtn = QtGui.QPushButton('Accept')
        self.AcceptBtn.clicked.connect(self.Accept_proc)

        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(self.CancelBtn)
        hlay.addWidget(self.AcceptBtn)

        vlay.addLayout(hlay)
        self.setLayout(vlay)

    #__________________________________________________________________________
    def Accept_proc(self):
        units = [self.list2.item(k).text() for k in range(self.list2.count())]
        self.close()
        return units

    #__________________________________________________________________________
    def toRight(self):
        if self.list1.currentRow() == -1:
            return
        item = self.list1.currentItem()
        self.list1.takeItem(self.list1.currentRow())
        self.list2.addItem(item)

    #__________________________________________________________________________
    def toLeft(self):
        if self.list2.currentRow() == -1:
            return
        item = self.list2.currentItem()
        self.list2.takeItem(self.list2.currentRow())
        self.list1.addItem(item)


#==============================================================================
class MoveUnitsWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Move Units')
        hlay = QtGui.QHBoxLayout()
        self.list = QtGui.QListWidget()
        self.list.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-up')
        btnUp = QtGui.QPushButton(icon, '')
        btnUp.clicked.connect(self.toUp)
        icon = QtGui.QIcon.fromTheme('go-down')
        btnDown = QtGui.QPushButton(icon, '')
        btnDown.clicked.connect(self.toDown)
        vlay = QtGui.QVBoxLayout()
        hlay.addWidget(self.list)
        vlay.addWidget(btnUp)
        vlay.addWidget(btnDown)
        hlay.addLayout(vlay)

        vlay = QtGui.QVBoxLayout()
        vlay.addLayout(hlay)

        self.CancelBtn = QtGui.QPushButton('Cancel')
        self.CancelBtn.clicked.connect(self.close)

        self.AcceptBtn = QtGui.QPushButton('Accept')
        self.AcceptBtn.clicked.connect(self.Accept_proc)

        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(self.CancelBtn)
        hlay.addWidget(self.AcceptBtn)

        vlay.addLayout(hlay)
        self.setLayout(vlay)

    #__________________________________________________________________________
    def Accept_proc(self):
        units = [self.list.item(k).text() for k in range(self.list.count())]
        self.close()
        return units

    #__________________________________________________________________________
    def toUp(self):
        if self.list.currentRow() in [-1, 0]:
            return
        place = self.list.currentRow()
        item = self.list.takeItem(place)
        self.list.insertItem(place - 1, item)
        self.list.setCurrentRow(place - 1)

    #__________________________________________________________________________
    def toDown(self):
        if self.list.currentRow() == -1:
            return
        if self.list.count() - 1 == self.list.currentRow():
            return
        place = self.list.currentRow()
        item = self.list.takeItem(place)
        self.list.insertItem(place + 1, item)
        self.list.setCurrentRow(place + 1)


#==============================================================================
class AutoClustWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Automatic Clustering')
        vlay = QtGui.QVBoxLayout(self)
        hlay = QtGui.QHBoxLayout()
        self.MinClust = QtGui.QSpinBox()
        self.MinClust.setRange(1, 5)
        self.MinClust.setValue(2)
        hlay.addWidget(QtGui.QLabel('Min'))
        hlay.addWidget(self.MinClust)
        self.MaxClust = QtGui.QSpinBox()
        self.MaxClust.setRange(2, 10)
        self.MaxClust.setValue(4)
        hlay.addWidget(QtGui.QLabel('Max'))
        hlay.addWidget(self.MaxClust)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        self.ClusteringMethod = QtGui.QComboBox(self)
        self.ClusteringMethod.addItems(['KlustaKwik', 'KMeans',
                                        'Afinity Propagation',
                                        'Mean-shift', 'Spectral Clustering',
                                        'Hierarchical clustering',
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


#==============================================================================
def KlustaKwik_call(data, minClust=2, maxClust=5):
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
    os.system('KlustaKwik data 1 -MinClusters %d -MaxClusters %d' % (minClust, maxClust))

    # read the results
    f = open('data.clu.1', 'r')
    clusterData = f.readlines()
    f.close()
    clusterData = [int(re.search('[0-9]{1,2}', k).group()) for k in clusterData]

    # the first line is the number of clusters
    nClusters = clusterData[0]
    clusterData.pop(0)
    clusterData = np.array(clusterData)

    # create an array with the indices of each cluster
    clustIndx = []
    for k in range(1, nClusters + 1):
        clustIndx.append(np.flatnonzero(clusterData == k))

    return clustIndx


#==============================================================================
class Settings(dt.DataSet):

    def chDir(self, item, value):
        self.FiguresDir = os.path.split(value)[0] + os.path.sep

    colormaps = [k for k in cm.datad.keys() if not k.endswith('_r')]

    WorkingDir = di.DirectoryItem('Select a Working Dir').set_prop("display", callback=chDir)
    FiguresDir = di.DirectoryItem('Path to save images')
    Figurescolor = di.ColorItem('Fig color', default='black')
    AxesColor = di.ColorItem('Axes color', default='gray').set_pos(col=1)
    LassoColor = di.ColorItem('Lasso color', default='gray')
    DensityCM = di.ChoiceItem('Density Color Map', tuple(colormaps),
                              default = colormaps.index('gist_heat')).set_pos(col=1)


#==============================================================================
class AutocorrOpts(dt.DataSet):
    Mode = di.ChoiceItem('Mode', ['ephys', 'fft', 'time'])
    BinSize = di.IntItem('BinSize', default=1).set_pos(col=1)
    Win0 = di.IntItem('Win 0', default=-150)
    Win1 = di.IntItem('Win 1', default=150).set_pos(col=1)

