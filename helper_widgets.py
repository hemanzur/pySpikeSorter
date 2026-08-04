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
        btn_right = QtWidgets.QPushButton(icon, '')
        btn_right.clicked.connect(self.to_right)
        icon = QtGui.QIcon.fromTheme('go-previous')
        btn_left = QtWidgets.QPushButton(icon, '')
        btn_left.clicked.connect(self.to_left)
        vlay = QtWidgets.QVBoxLayout()
        hlay.addWidget(self.list1)
        vlay.addWidget(btn_right)
        vlay.addWidget(btn_left)
        hlay.addLayout(vlay)
        hlay.addWidget(self.list2)

        vlay = QtWidgets.QVBoxLayout()
        vlay.addLayout(hlay)

        self.cancel_btn = QtWidgets.QPushButton('Cancel')
        self.cancel_btn.clicked.connect(self.close)

        self.accept_btn = QtWidgets.QPushButton('Accept')
        self.accept_btn.clicked.connect(self.accept)

        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(self.cancel_btn)
        hlay.addWidget(self.accept_btn)

        vlay.addLayout(hlay)
        self.setLayout(vlay)

    #__________________________________________________________________________
    def accept(self):
        units = [self.list2.item(k).text() for k in range(self.list2.count())]
        self.close()
        return units

    #__________________________________________________________________________
    def to_right(self):
        if self.list1.currentRow() == -1:
            return
        item = self.list1.currentItem()
        self.list1.takeItem(self.list1.currentRow())
        self.list2.addItem(item)

    #__________________________________________________________________________
    def to_left(self):
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
        btn_up = QtWidgets.QPushButton(icon, '')
        btn_up.clicked.connect(self.to_up)
        icon = QtGui.QIcon.fromTheme('go-down')
        btn_down = QtWidgets.QPushButton(icon, '')
        btn_down.clicked.connect(self.to_down)
        vlay = QtWidgets.QVBoxLayout()
        hlay.addWidget(self.list)
        vlay.addWidget(btn_up)
        vlay.addWidget(btn_down)
        hlay.addLayout(vlay)

        vlay = QtWidgets.QVBoxLayout()
        vlay.addLayout(hlay)

        self.cancel_btn = QtWidgets.QPushButton('Cancel')
        self.cancel_btn.clicked.connect(self.close)

        self.accept_btn = QtWidgets.QPushButton('Accept')
        self.accept_btn.clicked.connect(self.accept)

        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(self.cancel_btn)
        hlay.addWidget(self.accept_btn)

        vlay.addLayout(hlay)
        self.setLayout(vlay)

    #__________________________________________________________________________
    def accept(self):
        units = [self.list.item(k).text() for k in range(self.list.count())]
        self.close()
        return units

    #__________________________________________________________________________
    def to_up(self):
        if self.list.currentRow() in [-1, 0]:
            return
        place = self.list.currentRow()
        item = self.list.takeItem(place)
        self.list.insertItem(place - 1, item)
        self.list.setCurrentRow(place - 1)

    #__________________________________________________________________________
    def to_down(self):
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
        self.min_clust = QtWidgets.QSpinBox()
        self.min_clust.setRange(1, 5)
        self.min_clust.setValue(2)
        hlay.addWidget(QtWidgets.QLabel('Min'))
        hlay.addWidget(self.min_clust)
        self.max_clust = QtWidgets.QSpinBox()
        self.max_clust.setRange(2, 10)
        self.max_clust.setValue(4)
        hlay.addWidget(QtWidgets.QLabel('Max'))
        hlay.addWidget(self.max_clust)
        vlay.addLayout(hlay)

        hlay = QtWidgets.QHBoxLayout()
        self.clustering_method = QtWidgets.QComboBox(self)
        self.clustering_method.addItems(['KlustaKwik', 'KMeans',
                                        'Afinity Propagation',
                                        'Mean-shift', 'Spectral Clustering',
                                        'Hierarchical clustering',
                                        'DBSCAN', 'Gaussian Mixtures'])
        hlay.addWidget(QtWidgets.QLabel('Method'))
        hlay.addWidget(self.clustering_method)
        vlay.addLayout(hlay)

        hlay = QtWidgets.QHBoxLayout()
        auto_clust_btn = QtWidgets.QPushButton('Preview')
        #auto_clust_btn.clicked.connect(self.auto_clust)
        clustering_commit_btn = QtWidgets.QPushButton('Commit')
        hlay.addWidget(auto_clust_btn)
        hlay.addWidget(clustering_commit_btn)
        vlay.addLayout(hlay)

        self.setLayout(vlay)


#==============================================================================
def klusta_kwik_call(data, min_clust=2, max_clust=5):
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
    os.system('KlustaKwik data 1 -MinClusters %d -MaxClusters %d' % (min_clust, max_clust))

    # read the results
    f = open('data.clu.1', 'r')
    cluster_data = f.readlines()
    f.close()
    cluster_data = [int(re.search('[0-9]{1,2}', k).group()) for k in cluster_data]

    # the first line is the number of clusters
    n_clusters = cluster_data[0]
    cluster_data.pop(0)
    cluster_data = np.array(cluster_data)

    # create an array with the indices of each cluster
    cluster_indices = []
    for k in range(1, n_clusters + 1):
        cluster_indices.append(np.flatnonzero(cluster_data == k))

    return cluster_indices

"""

#==============================================================================
class settings_dialog(dt.DataSet):

    def ch_dir(self, item, value):
        self.figures_dir = os.path.split(value)[0] + os.path.sep

    colormaps = [k for k in cm.datad.keys() if not k.endswith('_r')]

    working_dir = di.DirectoryItem('Select a Working Dir').set_prop("display", callback=ch_dir)
    figures_dir = di.DirectoryItem('Path to save images')
    figures_color = di.ColorItem('Fig color', default='black')
    axes_color = di.ColorItem('Axes color', default='gray').set_pos(col=1)
    lasso_color = di.ColorItem('Lasso color', default='gray')
    density_cm = di.ChoiceItem('Density Color Map', tuple(colormaps),
                              default = colormaps.index('gist_heat')).set_pos(col=1)


#==============================================================================
class AutocorrOpts(dt.DataSet):
    mode = di.ChoiceItem('mode', ['ephys', 'fft', 'time'])
    bin_size = di.IntItem('bin_size', default=1).set_pos(col=1)
    win_0 = di.IntItem('win 0', default=-150)
    win_1 = di.IntItem('win 1', default=150).set_pos(col=1)

"""

#==============================================================================

class AutocorrOpts(QtWidgets.QWidget):
     def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.setWindowTitle('Autocorrelation Options')
        glay = QtWidgets.QGridLayout(self)

        self.mode = QtWidgets.QComboBox(self)
        self.mode.addItems(['ephys', 'fft', 'time'])
        glay.addWidget(QtWidgets.QLabel("mode"), 0, 0)
        glay.addWidget(self.mode, 0, 1)

        self.bin_size = QtWidgets.QSpinBox()
        self.bin_size.setRange(1, 20)
        self.bin_size.setValue(1)
        glay.addWidget(QtWidgets.QLabel("Bin Size"), 1, 0)
        glay.addWidget(self.bin_size, 1, 1)

        self.win_0 = QtWidgets.QSpinBox()
        self.win_0.setRange(-400, -50)
        self.win_0.setValue(-150)
        glay.addWidget(QtWidgets.QLabel("win_0"), 2, 0)
        glay.addWidget(self.win_0, 2, 1)

        self.win_1 = QtWidgets.QSpinBox()
        self.win_1.setRange(400, 50)
        self.win_1.setValue(150)
        glay.addWidget(QtWidgets.QLabel("win_1"), 3, 0)
        glay.addWidget(self.win_1, 3, 1)
