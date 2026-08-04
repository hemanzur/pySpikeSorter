from PyQt5 import QtGui, QtWidgets
#import guidata
import os
import re
import numpy as np
#app = guidata.qapplication()
#import guidata.dataset.datatypes as dt
#import guidata.dataset.dataitems as di
from matplotlib.pyplot import cm

colormaps = [k for k in cm.datad.keys() if not k.endswith('_r')]


#==============================================================================
class MergeUnitsWidget(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Merge Units')
        hlay = QtWidgets.QHBoxLayout()
        self.list1 = QtWidgets.QListWidget()
        self.list1.setMaximumWidth(100)
        self.list2 = QtWidgets.QListWidget()
        self.list2.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-next')
        btnRight = QtWidgets.QPushButton(icon, '')
        btnRight.clicked.connect(self.toRight)
        icon = QtGui.QIcon.fromTheme('go-previous')
        btnLeft = QtWidgets.QPushButton(icon, '')
        btnLeft.clicked.connect(self.toLeft)
        vlay = QtWidgets.QVBoxLayout()
        hlay.addWidget(self.list1)
        vlay.addWidget(btnRight)
        vlay.addWidget(btnLeft)
        hlay.addLayout(vlay)
        hlay.addWidget(self.list2)

        vlay = QtWidgets.QVBoxLayout()
        vlay.addLayout(hlay)

        self.CancelBtn = QtWidgets.QPushButton('Cancel')
        self.CancelBtn.clicked.connect(self.close)

        self.AcceptBtn = QtWidgets.QPushButton('Accept')
        self.AcceptBtn.clicked.connect(self.Accept_proc)

        hlay = QtWidgets.QHBoxLayout()
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
class MoveUnitsWidget(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Move Units')
        hlay = QtWidgets.QHBoxLayout()
        self.list = QtWidgets.QListWidget()
        self.list.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-up')
        btnUp = QtWidgets.QPushButton(icon, '')
        btnUp.clicked.connect(self.toUp)
        icon = QtGui.QIcon.fromTheme('go-down')
        btnDown = QtWidgets.QPushButton(icon, '')
        btnDown.clicked.connect(self.toDown)
        vlay = QtWidgets.QVBoxLayout()
        hlay.addWidget(self.list)
        vlay.addWidget(btnUp)
        vlay.addWidget(btnDown)
        hlay.addLayout(vlay)

        vlay = QtWidgets.QVBoxLayout()
        vlay.addLayout(hlay)

        self.CancelBtn = QtWidgets.QPushButton('Cancel')
        self.CancelBtn.clicked.connect(self.close)

        self.AcceptBtn = QtWidgets.QPushButton('Accept')
        self.AcceptBtn.clicked.connect(self.Accept_proc)

        hlay = QtWidgets.QHBoxLayout()
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
class AutoClustWidget(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Automatic Clustering')
        vlay = QtWidgets.QVBoxLayout(self)
        hlay = QtWidgets.QHBoxLayout()
        self.MinClust = QtWidgets.QSpinBox()
        self.MinClust.setRange(1, 5)
        self.MinClust.setValue(2)
        hlay.addWidget(QtWidgets.QLabel('Min'))
        hlay.addWidget(self.MinClust)
        self.MaxClust = QtWidgets.QSpinBox()
        self.MaxClust.setRange(2, 10)
        self.MaxClust.setValue(4)
        hlay.addWidget(QtWidgets.QLabel('Max'))
        hlay.addWidget(self.MaxClust)
        vlay.addLayout(hlay)

        hlay = QtWidgets.QHBoxLayout()
        self.ClusteringMethod = QtWidgets.QComboBox(self)
        self.ClusteringMethod.addItems(['KlustaKwik', 'KMeans',
                                        'Afinity Propagation',
                                        'Mean-shift', 'Spectral Clustering',
                                        'Hierarchical clustering',
                                        'DBSCAN', 'Gaussian Mixtures'])
        hlay.addWidget(QtWidgets.QLabel('Method'))
        hlay.addWidget(self.ClusteringMethod)
        vlay.addLayout(hlay)

        hlay = QtWidgets.QHBoxLayout()
        AutoClustBtn = QtWidgets.QPushButton('Preview')
        #AutoClustBtn.clicked.connect(self.AutoClust_proc)
        ClusteringCommitBtn = QtWidgets.QPushButton('Commit')
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

"""

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

"""

#==============================================================================

class AutocorrOpts(QtWidgets.QWidget):
     def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Autocorrelation Options')
        glay = QtWidgets.QGridLayout(self)

        self.Mode = QtWidgets.QComboBox(self)
        self.Mode.addItems(['ephys', 'fft', 'time'])
        glay.addWidget(QtWidgets.QLabel("Mode"), 0, 0)
        glay.addWidget(self.Mode, 0, 1)

        self.BinSize = QtWidgets.QSpinBox()
        self.BinSize.setRange(1, 20)
        self.BinSize.setValue(1)
        glay.addWidget(QtWidgets.QLabel("Bin Size"), 1, 0)
        glay.addWidget(self.BinSize, 1, 1)

        self.Win0 = QtWidgets.QSpinBox()
        self.Win0.setRange(-400, -50)
        self.Win0.setValue(-150)
        glay.addWidget(QtWidgets.QLabel("Win0"), 2, 0)
        glay.addWidget(self.Win0, 2, 1)

        self.Win1 = QtWidgets.QSpinBox()
        self.Win1.setRange(400, 50)
        self.Win1.setValue(150)
        glay.addWidget(QtWidgets.QLabel("Win1"), 3, 0)
        glay.addWidget(self.Win1, 3, 1)