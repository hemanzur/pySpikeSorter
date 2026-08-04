# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtWidgets
import tables
from matplotlib_widgets import MplWidget
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavToolbar
from matplotlib import rc
import sys
import numpy as np
try:
    from m_spike_utils import cross_correlation
except ImportError:
    cross_correlation = None

rc('font', size = 8)
rc('axes', titlesize = 10, labelsize = 8, labelcolor = 'w')
rc('xtick', labelsize = 8, color = 'w')
rc('ytick', labelsize = 8, color = 'w')

#%%
class SortingQualityWidget(QtWidgets.QWidget):
    
    def __init__(self, h5file = None):
        QtWidgets.QWidget.__init__(self)
        # define a right side control panel
        g_lay = QtWidgets.QGridLayout()
        row = 0
        
        if isinstance(h5file, tables.file.File):
            self.h5file = h5file
            
        elif isinstance(h5file, str):
            self.h5file = QtWidgets.QFileDialog.getOpenFileName(
                caption='select an h5 file', filter='*.h5')[0]
            if self.h5file:
                self.h5file = tables.open_file(self.h5file, 'r')
        
        elif not h5file:
            self.load_h5_file_btn = QtWidgets.QPushButton('Load H5File')
            self.load_h5_file_btn.clicked.connect(self.load_h5_file)
            g_lay.addWidget(self.load_h5_file_btn, row, 0, 1, 2)
            row += 1
            self.setWindowTitle('Spike Sorting Quality Explorer')
                
        self.first_unit_combo = QtWidgets.QComboBox()
        g_lay.addWidget(self.first_unit_combo, row, 0, 1, 2)
        row += 1
        
        self.select_btn = QtWidgets.QPushButton('Select None')
        self.select_btn.clicked.connect(self.select)
        self.select_btn.setCheckable(True)
        g_lay.addWidget(self.select_btn, row, 0)
        
        self.plot_xcorr_btn = QtWidgets.QPushButton('Plot xCorr')
        self.plot_xcorr_btn.clicked.connect(self.plot_xcorr)
        g_lay.addWidget(self.plot_xcorr_btn, row, 1)
        row += 1
        
        self.units_selector = QtWidgets.QTableWidget(0, 1)
        self.units_selector.verticalHeader().setVisible(False)
        self.units_selector.horizontalHeader().setVisible(False)
        self.units_selector.setColumnWidth(0, 200)
        g_lay.addWidget(self.units_selector, row, 0, 1, 2)
        row += 1
        
        main_lay = QtWidgets.QHBoxLayout(self)
        main_lay.addLayout(g_lay)
        
        # define a left side figure
        v_lay = QtWidgets.QVBoxLayout()
        self.main_fig = MplWidget(self)
        self.main_fig.figure.set_facecolor('k')
        self.ntb = NavToolbar(self.main_fig, self)
        self.ntb.setIconSize(QtCore.QSize(15, 15))
        v_lay.addWidget(self.main_fig)
        v_lay.addWidget(self.ntb)
    
        main_lay.addLayout(v_lay)        
        
        self.show()
        
        self.unit_checks = []
    
    def load_h5_file(self):
        
        if hasattr(self, 'h5file') and\
           isinstance(self.h5file, tables.file.File) and\
           self.h5file.isopen:
               self.h5file.close()
            
        self.h5file = QtWidgets.QFileDialog.getOpenFileName(
            caption='select an h5 file',
            filter='*.h5',
            directory='/home/hachi/Desktop/Data/Recording')[0]
        if self.h5file:
            self.h5file = tables.open_file(self.h5file, 'r')
            self.update_units_list()
    
    def select(self):
        if not self.select_btn.isChecked():
            self.select_btn.setText('Select None')
            for k in self.unit_checks:
                k.setChecked(True)
        else:
            self.select_btn.setText('Select All')
            for k in self.unit_checks:
                k.setChecked(False)
    
    '''def getUnitsProc(self):
        if not hasattr(self, 'h5file'): return
        if h5file,close(): return
        
        try:
            nodes = self.h5file.list_nodes('/Spikes')
        except:
            print('There is a problem with the H5File')
            return
            
        count = 0
        units = []
        for group in nodes:
            for member in group:
                if member._v_name.find('Unit') != -1:
                    units.append(member)
                    self.units_selector.insertRow(count)
                    count += 1'''

    def update_units_list(self):
        
        if not hasattr(self, 'h5file'): return
        
        # clear the FirstUnit Selector
        self.first_unit_combo.clear()
        
        # clean the table, kill the checkboxes
        self.units_selector.setRowCount(0)
        for k in self.unit_checks: k.deleteLater()
            
        try:
            nodes = self.h5file.list_nodes('/Spikes')
        except:
            print('There is a problem with the H5File')
            return
        
        count = 0
        self.unit_checks = []
        self.unit_ids = []
        for group in nodes:
            for member in group:
                if member._v_name.find('Unit') != -1:
                    self.units_selector.insertRow(count)
                    unit_id = group._v_name + ' ' + member._v_name
                    self.unit_checks.append(QtWidgets.QCheckBox(unit_id))
                    self.units_selector.setCellWidget(count, 0, self.unit_checks[-1])
                    self.units_selector.setRowHeight(count, 20)
                    self.first_unit_combo.addItem(unit_id)
                    self.unit_ids.append(unit_id)
                    count += 1
    
    def plot_xcorr(self):
                
        self.main_fig.figure.clf()
        base_unit = str(self.first_unit_combo.currentText())
        chan = base_unit[0:8]
        unit = base_unit[9:]
        
        #get the timestamps for that unit
        base_node = self.h5file.get_node('/Spikes/'+chan)
        timestamps = base_node.TimeStamp.read()
        base_unit_timestamps = base_node.__getattr__(unit).Indx.read()
        base_unit_timestamps = timestamps[base_unit_timestamps]
        
        # check wich units to plot
        units_to_plot = []
        for k in range(self.units_selector.rowCount()):
            if self.unit_checks[k].isChecked():
                units_to_plot.append( str(self.unit_checks[k].text()) )
                
        # create a grid of subplots of 8 columns by n rows
        n_rows = np.ceil(len(units_to_plot)/8.0)
        
        y_lim = 0
        axes_list = []
        # iterate over the list of units and plot the crosscorrelation
        for j, k in enumerate(units_to_plot):
            axes_list.append(self.main_fig.figure.add_subplot(int(n_rows), 8, j+1))
            
            chan = k[0:8]
            unit = k[9:]
            axes_list[-1].set_title(chan+' '+unit, color = 'w')
            #get the timestamps for that unit
            node = self.h5file.get_node('/Spikes/'+chan)
            timestamps = node.TimeStamp.read()
            unit_timestamps = node.__getattr__(unit).Indx.read()
            unit_timestamps = timestamps[unit_timestamps]
            
            r = []
            bin_size = 1
            #r, t = cross_correlation(base_unit_timestamps, unit_timestamps, bins = 20, win_lag = [-10, 10])
            for ts in base_unit_timestamps:
                t = unit_timestamps - ts
                r.extend(t[(t > -20) & (t < 20)])
            r, t = np.histogram(r, bins = int(40/bin_size) )
            #indx = np.flatnonzero((t>=-200) & (t<=200))
            axes_list[-1].bar(t[:-1], r, edgecolor = 'none', color = 'w')
            #ax.plot(t[indx], r[indx], 'w')
            axes_list[-1].set_xlim(-20, 20)
            y_lim = max([y_lim, max(r)])
            
            # change the color of the axes to white
            axes_list[-1].tick_params(axis = 'x', colors = 'w')
            axes_list[-1].tick_params(axis = 'y', colors = 'w')
            axes_list[-1].set_facecolor('none')
            
            for key, spine in axes_list[-1].spines.items():
                spine.set_color('w')

        #for ax in axes_list:
        #    ax.set_ylim(0, y_lim)
            
        self.main_fig.figure.tight_layout()
        self.main_fig.figure.canvas.draw()        
                    
if __name__ == '__main__':
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    sqw = SortingQualityWidget()
    #sys.exit(app.exec_())


#%%
from PyQt5.QtCore import pyqtSignal

class TestWidget(QtWidgets.QWidget):
    
    unit_signal = pyqtSignal(str)
    
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        # define a right side control panel
        vlay = QtWidgets.QVBoxLayout(self)
        self.btn1 = QtWidgets.QPushButton('Button1')
        vlay.addWidget(self.btn1)
        
        self.btn2 = QtWidgets.QPushButton('Button2')
        vlay.addWidget(self.btn2)
        
        self.lbl1 = QtWidgets.QLabel('')
        vlay.addWidget(self.lbl1)
        
        self.lbl2 = QtWidgets.QLabel('')
        vlay.addWidget(self.lbl2)
        
        #self.btn1.clicked.connect(self.change_lbl1)
        #self.btn1.clicked.connect(self.change_lbl2)
        
        self.btn1.clicked.connect(lambda: self.unit_signal.emit('Hola'))
        self.btn2.clicked.connect(lambda: self.unit_signal.emit('Borrado'))
        self.unit_signal.connect(self.change_lbl1)
        self.unit_signal.connect(self.change_lbl2)
        
        self.show()
        
    def change_lbl1(self, s):
        self.lbl1.setText(s+" label1")
        
    def change_lbl2(self, s):
        self.lbl2.setText(s+" label2")
