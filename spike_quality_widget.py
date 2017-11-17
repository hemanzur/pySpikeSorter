# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore
import tables
from matplotlib_widgets import MplWidget
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavToolbar
from matplotlib import rc
import sys
import numpy as np
from m_spike_utils import cross_correlation

rc('font', size = 8)
rc('axes', titlesize = 10, labelsize = 8, labelcolor = 'w')
rc('xtick', labelsize = 8, color = 'w')
rc('ytick', labelsize = 8, color = 'w')

#%%
class Sorting_Quality_Widget(QtGui.QWidget):
    
    def __init__(self, h5file = None):
        QtGui.QWidget.__init__(self)        
        # define a right side control panel
        gLay = QtGui.QGridLayout()
        row = 0
        
        if isinstance(h5file, tables.file.File):
            self.h5file = h5file
            
        elif isinstance(h5file, str):
            self.h5file = str(QtGui.QFileDialog.getOpenFileName(caption='select an h5 file',
                                                                filter='*.h5'))
            if self.h5file:
                self.h5file = tables.openFile(self.h5file, 'r')
        
        elif not h5file:
            self.loadH5FileBtn = QtGui.QPushButton('Load H5File')
            self.loadH5FileBtn.clicked.connect(self.loadH5FileProc)
            gLay.addWidget(self.loadH5FileBtn, row, 0, 1, 2)
            row += 1
            self.setWindowTitle('Spike Sorting Quality Explorer')
                
        self.FirstUnitCombo = QtGui.QComboBox()
        gLay.addWidget(self.FirstUnitCombo, row, 0, 1, 2)
        row += 1
        
        self.selectBtn = QtGui.QPushButton('Select None')
        self.selectBtn.clicked.connect(self.selectProc)
        self.selectBtn.setCheckable(True)
        gLay.addWidget(self.selectBtn, row, 0)
        
        self.plotXCorrBtn = QtGui.QPushButton('Plot xCorr')
        self.plotXCorrBtn.clicked.connect(self.plotXCorr)
        gLay.addWidget(self.plotXCorrBtn, row, 1)
        row += 1
        
        self.UnitsSelector = QtGui.QTableWidget(0, 1)
        self.UnitsSelector.verticalHeader().setVisible(False)
        self.UnitsSelector.horizontalHeader().setVisible(False)
        self.UnitsSelector.setColumnWidth(0, 200)
        gLay.addWidget(self.UnitsSelector, row, 0, 1, 2)
        row += 1
        
        mainLay = QtGui.QHBoxLayout(self)
        mainLay.addLayout(gLay)
        
        # define a left side figure
        vLay = QtGui.QVBoxLayout()
        self.mainFig = MplWidget(self)
        self.mainFig.figure.set_facecolor('k')
        self.ntb = NavToolbar(self.mainFig, self)
        self.ntb.setIconSize(QtCore.QSize(15, 15))
        vLay.addWidget(self.mainFig)
        vLay.addWidget(self.ntb)
    
        mainLay.addLayout(vLay)        
        
        self.show()
        
        self.UnitChecks = []
    
    def loadH5FileProc(self):
        
        if hasattr(self, 'h5file') and\
           isinstance(self.h5file, tables.file.File) and\
           self.h5file.isopen:
               self.h5file.close()
            
        self.h5file = str(QtGui.QFileDialog.getOpenFileName(caption='select an h5 file',
                                                            filter='*.h5',
                                                            directory = '/home/hachi/Desktop/Data/Recording'))
        if self.h5file:
            self.h5file = tables.openFile(self.h5file, 'r')
            self.updateUnitsList()
    
    def selectProc(self):
        if not self.selectBtn.isChecked():
            self.selectBtn.setText('Select None')
            for k in self.UnitChecks:
                k.setChecked(True)
        else:
            self.selectBtn.setText('Select All')
            for k in self.UnitChecks:
                k.setChecked(False)
    
    '''def getUnitsProc(self):
        if not hasattr(self, 'h5file'): return
        if h5file,close(): return
        
        try:
            nodes = self.h5file.listNodes('/Spikes')
        except:
            print 'There is a problem with the H5File'
            
        count = 0
        units = []
        for group in nodes:
            for member in group:
                if member._v_name.find('Unit') != -1:
                    units.append(member)
                    self.UnitsSelector.insertRow(count)
                    count += 1'''

    def updateUnitsList(self):
        
        if not hasattr(self, 'h5file'): return
        
        # clear the FirstUnit Selector
        self.FirstUnitCombo.clear()
        
        # clean the table, kill the checkboxes
        self.UnitsSelector.setRowCount(0)
        for k in self.UnitChecks: k.deleteLater()
            
        try:
            nodes = self.h5file.listNodes('/Spikes')
        except:
            print 'There is a problem with the H5File'
        
        count = 0
        self.UnitChecks = []
        self.unitIDs = []
        for group in nodes:
            for member in group:
                if member._v_name.find('Unit') != -1:
                    self.UnitsSelector.insertRow(count)
                    unitID = group._v_name + ' ' + member._v_name
                    self.UnitChecks.append(QtGui.QCheckBox(unitID))
                    self.UnitsSelector.setCellWidget(count, 0, self.UnitChecks[-1])
                    self.UnitsSelector.setRowHeight(count, 20)
                    self.FirstUnitCombo.addItem(unitID)
                    self.unitIDs.append(unitID)
                    count += 1
    
    def plotXCorr(self):
                
        self.mainFig.figure.clf()
        baseUnit = str(self.FirstUnitCombo.currentText())
        chan = baseUnit[0:8]
        unit = baseUnit[9:]
        
        #get the timestamps for that unit
        baseNode = self.h5file.getNode('/Spikes/'+chan)
        TS = baseNode.TimeStamp.read()
        baseUnitTS = baseNode.__getattr__(unit).Indx.read()
        baseUnitTS = TS[baseUnitTS]
        
        # check wich units to plot
        units2Plot = []
        for k in range(self.UnitsSelector.rowCount()):
            if self.UnitChecks[k].isChecked():
                units2Plot.append( str(self.UnitChecks[k].text()) )
                
        # create a grid of subplots of 8 columns by n rows
        nRows = np.ceil(len(units2Plot)/8.0)
        
        ylim = 0
        axes_list = []
        # iterate over the list of units and plot the crosscorrelation
        for j, k in enumerate(units2Plot):
            axes_list.append(self.mainFig.figure.add_subplot(nRows, 8, j+1))
            
            chan = k[0:8]
            unit = k[9:]
            axes_list[-1].set_title(chan+' '+unit, color = 'w')
            #get the timestamps for that unit
            node = self.h5file.getNode('/Spikes/'+chan)
            TS = node.TimeStamp.read()
            UnitTS = node.__getattr__(unit).Indx.read()
            UnitTS = TS[UnitTS]
            
            r = []
            bin_size = 1
            #r, t = cross_correlation(baseUnitTS, UnitTS, bins = 20, win_lag = [-10, 10])
            for ts in baseUnitTS:
                t = UnitTS - ts
                r.extend(t[(t > -20) & (t < 20)])
            r, t = np.histogram(r, bins = int(40/bin_size) )
            #indx = np.flatnonzero((t>=-200) & (t<=200))
            axes_list[-1].bar(t[:-1], r, edgecolor = 'none', color = 'w')
            #ax.plot(t[indx], r[indx], 'w')
            axes_list[-1].set_xlim(-20, 20)
            ylim = max([ylim, max(r)])
            
            # change the color of the axes to white
            axes_list[-1].tick_params(axis = 'x', colors = 'w')
            axes_list[-1].tick_params(axis = 'y', colors = 'w')
            axes_list[-1].set_axis_bgcolor('none')
            
            for key, spine in axes_list[-1].spines.iteritems():
                spine.set_color('w')

        #for ax in axes_list:
        #    ax.set_ylim(0, ylim)
            
        self.mainFig.figure.tight_layout()
        self.mainFig.figure.canvas.draw()        
                    
if __name__ == '__main__':
    if not QtGui.QApplication.instance():
        app = QtGui.QApplication(sys.argv)
    sqw = Sorting_Quality_Widget()
    #sys.exit(app.exec_())


#%%
from PyQt4 import QtGui
from PyQt4.QtCore import pyqtSignal

class TestWidget(QtGui.QWidget):
    
    Unit = pyqtSignal(str)
    
    def __init__(self):
        QtGui.QWidget.__init__(self)        
        # define a right side control panel
        vlay = QtGui.QVBoxLayout(self)
        self.btn1 = QtGui.QPushButton('Button1')
        vlay.addWidget(self.btn1)
        
        self.btn2 = QtGui.QPushButton('Button2')
        vlay.addWidget(self.btn2)
        
        self.lbl1 = QtGui.QLabel('')
        vlay.addWidget(self.lbl1)
        
        self.lbl2 = QtGui.QLabel('')
        vlay.addWidget(self.lbl2)
        
        #self.btn1.clicked.connect(self.ChangeLbl1)
        #self.btn1.clicked.connect(self.ChangeLbl2)
        
        self.btn1.clicked.connect(lambda: self.Unit.emit('Hola'))
        self.btn2.clicked.connect(lambda: self.Unit.emit('Borrado'))
        self.Unit.connect(self.ChangeLbl1)
        self.Unit.connect(self.ChangeLbl2)
        
        self.show()
        
    def ChangeLbl1(self, s):
        self.lbl1.setText(s+" label1")
        
    def ChangeLbl2(self, s):
        self.lbl2.setText(s+" label2")

test = TestWidget()
#%%        
#Load H5file
h5file = QtGui.QFileDialog.getOpenFileName(directory = RecDir, filter='*h5')
