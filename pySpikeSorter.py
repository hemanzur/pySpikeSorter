###!/usr/local/bin/ipython -i

#---------------------------------------------------------------------- IMPORTS

import os
import sys

#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Plain',
#color_scheme='Linux', call_pdb=1)

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtCore import pyqtRemoveInputHook, pyqtRestoreInputHook

def debug_trace():
    pyqtRemoveInputHook()
    from pdb import set_trace
    set_trace()
    

if not QtWidgets.QApplication.instance():
    app = QtWidgets.QApplication([])
else:
    app = QtWidgets.QApplication.instance()

from pyqtgraph import opengl as gl

filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    exec(open(filename).read())


import re
import tables
import numpy as np

# extra widgets import
import matplotlib_widgets
import helper_widgets

from matplotlib import rc
# from matplotlib.mlab import PCA
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

from matplotlib.path import Path
from scipy.spatial import cKDTree
import datetime

import m_BlackrockLib as BL


def pca_scores(data, n_components=3):
    data = np.asarray(data)
    n_components = min(n_components, data.shape[0], data.shape[1])
    if n_components < 1:
        return np.empty((data.shape[0], 0))
    scores = PCA(n_components=n_components).fit_transform(data)
    if scores.shape[1] < 3:
        scores = np.pad(scores, ((0, 0), (0, 3 - scores.shape[1])), mode='constant')
    return scores


def sample_step(total, target):
    return max(1, int(total // target))


def sample_range(total, target):
    return list(range(0, total, sample_step(total, target)))


#==============================================================================
def autocorr(time_stamp, bin_size=20, win=[0, 10000], mode='time',
             time_range=[-200, 200]):

    if not np.any(time_stamp):
        return

    win = np.array(win)
    time_stamp = np.array(time_stamp)
    time_stamp = time_stamp - time_stamp[0]
    timestamps = time_stamp[(time_stamp >= win[0]) & (time_stamp < win[1])]

    if timestamps.size > 1000:
        time_stamp = timestamps

    bin_size = int(bin_size)
    n_bins = int(time_stamp[-1] // bin_size)
    train = np.zeros(n_bins + 1, dtype=np.int16)
    for k in np.floor(time_stamp / bin_size).astype('int'):
        train[k] = train[k] + 1

    if mode == 'time':
        ac = np.correlate(train, train, mode='same')
        x = np.linspace(-time_stamp[-1] / 2, time_stamp[-1] / 2, ac.size)

    elif mode == 'ephys':
        tmp = np.array([])
        for k in time_stamp:
            t = time_stamp - k
            t = t[(t > time_range[0]) & (t < time_range[1])]
            tmp = np.append(tmp, t)
        ac, x = np.histogram(tmp, bins=int(np.diff(time_range) / bin_size))
        x = x[1:]
        ac[np.flatnonzero(x == 0)] = 0

    elif mode == 'fft':
        s = np.fft.fft(train)
        ac = np.abs(np.fft.ifft(s * np.conjugate(s)))
        #ac = ac/(train.size/((win[1]-win[0])/1000))
        mid = ac.size // 2
        ac = np.concatenate([ac[mid:], ac[:mid]])
        x = np.linspace(-time_stamp[-1] / 2, time_stamp[-1] / 2, ac.size)

    return ac, x


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
    if os.system('KlustaKwik data 1 -MinClusters %d -MaxClusters %d'
                 % (min_clust, max_clust)) != 256:
            return

    # wait while klustakwick gets the clusters
    while not os.path.isfile('data.clu.1'):
        continue

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


#==============================================================================
# Spike Sorter Main GUI window

rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

# create instance of imported widgets
if hasattr(helper_widgets, 'settings_dialog'):
    settings = helper_widgets.settings_dialog()
else:
    class _DefaultSettings:
        density_cm = helper_widgets.colormaps.index('gist_heat')

        def edit(self):
            return 0

    settings = _DefaultSettings()
autocorropts = helper_widgets.AutocorrOpts()
autoclust = helper_widgets.AutoClustWidget()


#==============================================================================
class SpikeSorter(QtWidgets.QMainWindow):

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        self.setWindowTitle("pySpikeSorter")
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap('spike_icon.png')))
        self.main_widget = QtWidgets.QWidget(self)
        self.main_layout = QtWidgets.QHBoxLayout(self.main_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(0)

        self.cur_unit = 0
        self.plot_unit_counter = 0
        self.units_list = []
        self.n_units = 0
        self.h5_file_loaded = False
        self.chan_plotted = False
        self.removing_tab = 0
        self.unit_colors = np.array([[1, 0, 0], [0, 0.7, 0], [0, 0.4, 1],
                                    [0.8, 0.6, 0], [0.6, 0, 1], [0, 0.7, 0.7],
                                    [0, 0.5, 1]])
        self.unit_colors = np.tile(self.unit_colors, (10, 1))

        #--------------------------------------------- TOOLBAR ON THE LEFT SIDE

        split1 = QtWidgets.QSplitter(QtCore.Qt.Horizontal, self.main_widget)   # SPLITTER
        self.tools_tab = QtWidgets.QTabWidget()

        tools_tab_1 = QtWidgets.QWidget()
        tools_tab_2 = QtWidgets.QWidget()

        self.tools_tab.addTab(tools_tab_1, 'Main Tools')
        self.tools_tab.addTab(tools_tab_2, 'Chan Tools')
        self.tools_tab.setMaximumWidth(220)

        split1.addWidget(self.tools_tab)

        #---------------------------------------------------------tools_tab No 1
        tools_lay = QtWidgets.QVBoxLayout()

        #-------------------------------------------------------------- FRAME 1
        grp = QtWidgets.QGroupBox('Overview Tools', tools_tab_1)
        vlay = QtWidgets.QVBoxLayout()

        # number of events to overview spin box
        hlay = QtWidgets.QHBoxLayout()
        self.overview_n_events_spin = QtWidgets.QSpinBox()
        self.overview_n_events_spin.setRange(100, 1000)
        self.overview_n_events_spin.setSingleStep(100)
        self.overview_n_events_spin.setValue(500)
        hlay.addWidget(QtWidgets.QLabel('N Events 2 Overview'))
        hlay.addWidget(self.overview_n_events_spin)
        vlay.addLayout(hlay)

        # Y axis limits selector
        hlay = QtWidgets.QHBoxLayout()
        self.overview_y_lims_spin = QtWidgets.QSpinBox()
        self.overview_y_lims_spin.setRange(100, 5000)
        self.overview_y_lims_spin.setSingleStep(100)
        self.overview_y_lims_spin.setValue(2000)
        self.overview_y_lims_spin.editingFinished.connect(self.change_overview_y_lim)
        hlay.addWidget(QtWidgets.QLabel('Overview Axes YLim'))
        hlay.addWidget(self.overview_y_lims_spin)
        vlay.addLayout(hlay)

        #----------------------------------------------------------------------
        btn = QtWidgets.QPushButton('Plot Overview')
        btn.setStyleSheet('QPushButton{background-color: rgba(0,190,0)}')
        btn.clicked.connect(self.load_h5_file)
        vlay.addWidget(btn)

        btn = QtWidgets.QPushButton('Save Overview')
        btn.clicked.connect(self.save_overview_fig)
        vlay.addWidget(btn)

        grp.setLayout(vlay)
        tools_lay.addWidget(grp)

        grp = QtWidgets.QGroupBox('Delete Channel Tools', tools_tab_1)
        vlay = QtWidgets.QVBoxLayout()

        # add mark trash spin and button and link it to a function
        hlay = QtWidgets.QHBoxLayout()
        self.mark_trash_spin = QtWidgets.QSpinBox()
        self.mark_trash_spin.setMinimum(1)
        self.mark_trash_spin.setMaximum(1000000)
        self.mark_trash_spin.setValue(1000)
        hlay.addWidget(QtWidgets.QLabel('Below'))
        hlay.addWidget(self.mark_trash_spin)
        mark_trash_btn = QtWidgets.QPushButton('Mark Trash')
        mark_trash_btn.clicked.connect(self.trash_chans)
        hlay.addWidget(mark_trash_btn)
        vlay.addLayout(hlay)

        # add delete trash chans and link it to a function
        btn = QtWidgets.QPushButton('Delete Trash Chans')
        btn.clicked.connect(self.delete_trash_chans)
        vlay.addWidget(btn)

        grp.setLayout(vlay)
        tools_lay.addWidget(grp)

        #-------------------------------------------------------------- FRAME 2
        grp = QtWidgets.QGroupBox('Channel Plot Options', tools_tab_1)
        vlay = QtWidgets.QVBoxLayout()

        hlay = QtWidgets.QHBoxLayout()
        self.chan_selector = QtWidgets.QComboBox()
        hlay.addWidget(QtWidgets.QLabel('Chan Selector'))
        hlay.addWidget(self.chan_selector)
        vlay.addLayout(hlay)

        plot_chan_btn = QtWidgets.QPushButton('Plot Chan')
        plot_chan_btn.clicked.connect(self.plot_chan)
        vlay.addWidget(plot_chan_btn)

        grp.setLayout(vlay)
        tools_lay.addWidget(grp)

        #------------------------------------------------------------ Group No3

        grp = QtWidgets.QGroupBox('General Tools', tools_tab_1)
        glay = QtWidgets.QGridLayout()

        settings_btn = QtWidgets.QPushButton('Settings')
        settings_btn.clicked.connect(self.settings_dialog)
        glay.addWidget(settings_btn, 0, 0)

        about_btn = QtWidgets.QPushButton('About')
        about_btn.clicked.connect(self.about)
        glay.addWidget(about_btn, 0, 1)

        close_h5_file_btn = QtWidgets.QPushButton('Close H5 File')
        close_h5_file_btn.clicked.connect(self.close_file)
        glay.addWidget(close_h5_file_btn, 1, 0)

        exit_btn = QtWidgets.QPushButton('Exit')
        exit_btn.clicked.connect(self.close_event)
        glay.addWidget(exit_btn, 1, 1)

        convert_file_btn = QtWidgets.QPushButton('Convert File')
        convert_file_btn.clicked.connect(BL.bin_to_h5)
        glay.addWidget(convert_file_btn, 2, 0)

        grp.setLayout(glay)
        tools_lay.addWidget(grp)

        # create an "about" Msg Box
        self.about_msg = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Information,
                                          'about',
                                          u'Spyke Sorter v0.1\nHachi Manzur, 2012')

        tools_lay.addStretch(1)

        tools_tab_1.setLayout(tools_lay)

        #--------------------------------------------------- self.tools_tab No 2
        tools_lay = QtWidgets.QVBoxLayout()

        # group No1
        grp = QtWidgets.QGroupBox('Features Plot Opts', tools_tab_2)
        vlay = QtWidgets.QVBoxLayout()

        # add X and Y features selection combobox
        items = ['PCA1', 'PCA2', 'PCA3', 'Slice1', 'Slice2', 'Time', 'Pk2Pk Amp',
                 'Peak', 'Valley', 'Energy', 'Peak Pt', 'Valley Pt']
        self.x_plot = QtWidgets.QComboBox(grp)
        self.y_plot = QtWidgets.QComboBox(grp)
        self.z_plot = QtWidgets.QComboBox(grp)
        self.x_plot.addItems(items)
        self.y_plot.addItems(items)
        self.z_plot.addItems(items)
        self.z_plot.addItem('Density')
        self.y_plot.setCurrentIndex(1)
        self.z_plot.setCurrentIndex(2)

        # add the X axis combo box
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('X Axis Variable'))
        hlay.addWidget(self.x_plot)
        vlay.addLayout(hlay)

        # add the Y axis combo box
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('Y Axis Variable'))
        hlay.addWidget(self.y_plot)
        vlay.addLayout(hlay)

        # add the Y axis combo box
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('Z Axis Variable'))
        hlay.addWidget(self.z_plot)
        vlay.addLayout(hlay)

        # add a source of what to plot selection combo box
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('What to Plot ?'))
        self.what_to_plot = QtWidgets.QComboBox()
        hlay.addWidget(self.what_to_plot)
        vlay.addLayout(hlay)

        # add two slice selection spin box
        hlay = QtWidgets.QHBoxLayout()
        self.slice_spin_1 = QtWidgets.QSpinBox()
        self.slice_spin_1.setObjectName('Slice1')
        self.slice_spin_1.valueChanged.connect(self.slice_draw)
        hlay.addWidget(QtWidgets.QLabel('Slice 1'))
        hlay.addWidget(self.slice_spin_1)
        self.slice_spin_2 = QtWidgets.QSpinBox()
        self.slice_spin_2.setObjectName('Slice2')
        hlay.addWidget(QtWidgets.QLabel('Slice 2'))
        self.slice_spin_2.valueChanged.connect(self.slice_draw)
        hlay.addWidget(self.slice_spin_2)
        vlay.addLayout(hlay)

        # add a plot density check and a spin box to set the resolution
        hlay = QtWidgets.QHBoxLayout()
        self.plot_density_check = QtWidgets.QCheckBox('Plot Density ?')
        hlay.addWidget(self.plot_density_check)
        self.plot_density_bins = QtWidgets.QSpinBox()
        self.plot_density_bins.setMinimum(50)
        self.plot_density_bins.setMaximum(300)
        self.plot_density_bins.setValue(100)
        hlay.addWidget(self.plot_density_bins)
        vlay.addLayout(hlay)

        # plot only valid wfs check widget
        self.plot_valids_only_check = QtWidgets.QCheckBox('Plot Valids Only')
        self.plot_valids_only_check.setChecked(True)
        vlay.addWidget(self.plot_valids_only_check)

        # label with number of points
        hlay = QtWidgets.QHBoxLayout()
        self.n_pts_label = QtWidgets.QLabel()
        hlay.addWidget(QtWidgets.QLabel('NPoints'))
        hlay.addWidget(self.n_pts_label)
        vlay.addLayout(hlay)

        # number of spikes spin box
        hlay = QtWidgets.QHBoxLayout()
        self.n_pts_spin = QtWidgets.QSpinBox()
        self.n_pts_spin.setRange(10000, 200000)
        self.n_pts_spin.setSingleStep(10000)
        hlay.addWidget(self.n_pts_spin)

        # number of spikes slider
        self.n_pts_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.n_pts_slider.setRange(10000, 200000)
        self.n_pts_slider.setTickInterval(5000)
        self.n_pts_slider.setSingleStep(5000)
        hlay.addWidget(self.n_pts_slider)

        # connect spinner with No-of-spikes slider
        self.n_pts_slider.valueChanged.connect(self.n_pts_spin.setValue)

        # connect slider with No-of-spikes spinner
        self.n_pts_spin.valueChanged.connect(self.n_pts_slider.setValue)

        # set N spikes value
        self.n_pts_slider.setValue(50000)
        vlay.addLayout(hlay)

        hlay = QtWidgets.QHBoxLayout()
        # plot features btn and funcion connection
        self.plot_features_btn = QtWidgets.QPushButton('Plot 2D', grp)
        self.plot_features_btn.clicked.connect(self.plot_features)
        hlay.addWidget(self.plot_features_btn)

        # plot features btn and funcion connection
        self.plot_3d_btn = QtWidgets.QPushButton('Plot 3D', grp)
        self.plot_3d_btn.clicked.connect(self.plot_3d_features)
        hlay.addWidget(self.plot_3d_btn)
        vlay.addLayout(hlay)

        grp.setLayout(vlay)

        tools_lay.addWidget(grp)

        #----------------------------------------------------------- group No 2
        grp = QtWidgets.QGroupBox('Raw waveforms Opts')
        vlay = QtWidgets.QVBoxLayout()

        # number of spikes spin box
        hlay = QtWidgets.QHBoxLayout()
        self.n_spikes_spin = QtWidgets.QSpinBox()
        self.n_spikes_spin.setMaximum(5000)
        self.n_spikes_spin.setMinimum(100)
        self.n_spikes_spin.setSingleStep(100)
        hlay.addWidget(self.n_spikes_spin)

        # number of spikes slider
        self.n_spikes_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.n_spikes_slider.setMaximum(5000)
        self.n_spikes_slider.setMinimum(100)
        self.n_spikes_slider.setSingleStep(100)
        hlay.addWidget(self.n_spikes_slider)

        # connect spinner with No-of-spikes slider
        self.n_spikes_spin.valueChanged.connect(self.n_spikes_slider.setValue)

        # connect slider with No-of-spikes spinner
        self.n_spikes_slider.valueChanged.connect(self.n_spikes_spin.setValue)

        # set N spikes value
        self.n_spikes_slider.setValue(1000)
        vlay.addLayout(hlay)

        # add axes limit spin box
        hlay = QtWidgets.QHBoxLayout()
        self.wave_ax_y_lim_spin = QtWidgets.QSpinBox()
        self.wave_ax_y_lim_spin.setRange(0, 10000)
        self.wave_ax_y_lim_spin.setValue(1000)
        self.wave_ax_y_lim_spin.setSingleStep(100)
        self.wave_ax_y_lim_spin.editingFinished.connect(self.set_wf_plot_lim)
        hlay.addWidget(QtWidgets.QLabel('Axes Y Lim'))
        hlay.addWidget(self.wave_ax_y_lim_spin)
        vlay.addLayout(hlay)

        # create a "plot waveforms" check widget
        self.plot_waveforms_check = QtWidgets.QCheckBox('Plot Raw waveforms ?')
        vlay.addWidget(self.plot_waveforms_check)

        grp.setLayout(vlay)

        tools_lay.addWidget(grp)

        # Automatic clustering box
        #w = QtWidgets.QWidget()
        #autoClustBtn = QtWidgets.QPushButton('Automatic Clustering')
        #autoClustBtn.clicked.connect(autoclust.show)
        #tools_lay.addWidget(autoClustBtn)

        hlay = QtWidgets.QHBoxLayout()
        merge_units_btn = QtWidgets.QPushButton('Merge Units')
        self.merge_units_widget = helper_widgets.MergeUnitsWidget()
        self.merge_units_widget.accept_btn.clicked.connect(self.merge_units)
        merge_units_btn.clicked.connect(self.call_merge_units)
        hlay.addWidget(merge_units_btn)

        move_units_btn = QtWidgets.QPushButton('Move Units')
        self.move_units_widget = helper_widgets.MoveUnitsWidget()
        self.move_units_widget.accept_btn.clicked.connect(self.move_units)
        move_units_btn.clicked.connect(self.call_move_units)
        hlay.addWidget(move_units_btn)

        tools_lay.addLayout(hlay)

        #----------------------------------------------- CHANNEL METAINFO GROUP

        # button to reset a channel
        btn = QtWidgets.QPushButton('Reset Channel')
        btn.clicked.connect(self.reset_chan)
        tools_lay.addWidget(btn)

        # button to reset a channel
        btn = QtWidgets.QPushButton('Autocorr Opts')
        btn.clicked.connect(self.autocorr_opts)
        tools_lay.addWidget(btn)

        tools_lay.addStretch(1)
        tools_tab_2.setLayout(tools_lay)

        #------------------------------------------------ TABBED FIGURES WIDGET

        self.overview_tab_1 = {}
        self.overview_tab_2 = {}

        self.main_fig_tab_widget = QtWidgets.QTabWidget()
        self.main_fig_tab_widget.currentChanged.connect(self.main_fig_tab)
        self.overview_tab_1['main_widget'] = QtWidgets.QWidget(self.main_fig_tab_widget)
        hlay = QtWidgets.QHBoxLayout(self.overview_tab_1['main_widget'])

        self.main_fig_tab_widget.addTab(self.overview_tab_1['main_widget'], 'Channels Overview')

        # overview figure
        self.overview_tab_1['figure'] = matplotlib_widgets.MplWidget()
        self.overview_tab_1['figure'].figure.set_facecolor('k')
        self.overview_tab_1['toolbar'] = matplotlib_widgets.NavToolbar(self.overview_tab_1['figure'], self.overview_tab_1['main_widget'])
        self.overview_tab_1['toolbar'].setIconSize(QtCore.QSize(15, 15))
        vlay = QtWidgets.QVBoxLayout()
        vlay.addWidget(self.overview_tab_1['figure'])
        vlay.addWidget(self.overview_tab_1['toolbar'])
        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(1)
        hlay.addLayout(vlay)
        hlay.setContentsMargins(0, 0, 0, 0)
        hlay.setSpacing(1)

        #------------------------------------------------ OVERVIEW TABLE WIDGET
        self.overview_tab_2['main_widget'] = QtWidgets.QWidget(self.main_fig_tab_widget)
        self.overview_tab_2['overview_table'] = QtWidgets.QTableWidget(0, 6, self.overview_tab_2['main_widget'])
        self.overview_tab_2['overview_table'].setAlternatingRowColors(True)
        self.overview_tab_2['overview_table'].setFont(QtGui.QFont('sans', 8))
        labels = ['Count', 'isTrash', 'MultiUnit?', 'Comments', 'Unsorted', 'Valid']
        self.overview_tab_2['overview_table'].setHorizontalHeaderLabels(labels)
        for k in range(self.overview_tab_2['overview_table'].columnCount()):
            self.overview_tab_2['overview_table'].setColumnWidth(k, 65)
        self.overview_tab_2['overview_table'].setColumnWidth(3, 150)
        self.overview_tab_2['overview_table'].setColumnWidth(2, 75)

        # associate the vertical header click to select the channel
        v_header = self.overview_tab_2['overview_table'].verticalHeader()
        v_header.sectionClicked.connect(self.table_row_changed)

        vlay = QtWidgets.QVBoxLayout(self.overview_tab_2['main_widget'])
        vlay.addWidget(self.overview_tab_2['overview_table'])

        # add a log entry browser
        grp = QtWidgets.QGroupBox('Log Browser')
        grp.setMaximumHeight(100)
        hlay = QtWidgets.QHBoxLayout()
        self.log_combo = QtWidgets.QComboBox()
        self.log_combo.setMinimumWidth(200)
        #self.log_combo.setMinimumHeight(20)
        self.log_combo.currentIndexChanged.connect(self.set_log_text)
        hlay.addWidget(self.log_combo)
        self.log_text_browser = QtWidgets.QTextBrowser()
        #self.log_text_browser.setMaximumHeight(40)
        hlay.addWidget(self.log_text_browser)
        hlay.setContentsMargins(0, 0, 0, 0)
        hlay.setSpacing(1)
        grp.setLayout(hlay)
        vlay.addWidget(grp)

        self.main_fig_tab_widget.addTab(self.overview_tab_2['main_widget'], 'Summary Table')

        #---------------------------------------------------------- CHANNEL TAB
        self.chan_tab = {}
        self.chan_tab['main_widget'] = QtWidgets.QWidget()
        self.main_fig_tab_widget.addTab(self.chan_tab['main_widget'], 'Channel Tab')

        main_h_lay = QtWidgets.QHBoxLayout()

        #------------------------------------------------- RAW WAVEFORMS WIDGET
        # create the mpl widget to plot the raw waveforms
        vlay = QtWidgets.QVBoxLayout()

        # buttons and controls on top of raw waveforms plot
        hlay = QtWidgets.QHBoxLayout()
        self.n_units_spin = QtWidgets.QSpinBox()
        self.n_units_spin.setMaximumHeight(20)
        self.n_units_spin.setMinimum(1)
        self.n_units_spin.setMaximum(10000)
        self.n_units_spin.setValue(1)

        trim_btn = QtWidgets.QPushButton('Trim waveforms')
        trim_btn.clicked.connect(self.activate_trim_waveforms)
        trim_btn.setMaximumHeight(20)

        clean_btn = QtWidgets.QPushButton('Redraw')
        clean_btn.setMaximumHeight(20)

        clean_btn.clicked.connect(self.clean_waves_figure)
        hlay.addStretch(1)
        lbl = QtWidgets.QLabel('waveforms_to_plot:')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.n_units_spin)
        hlay.addWidget(trim_btn)
        hlay.addWidget(clean_btn)
        hlay.addStretch(1)
        vlay.addLayout(hlay)

        #------------------------------------------- waveforms plot and toolbar
        hlay = QtWidgets.QHBoxLayout()
        self.chan_tab['waves_figure'] = matplotlib_widgets.MplWidget()
        self.chan_tab['waves_figure'].figure.set_facecolor('k')
        self.chan_tab['wave_toolbar'] = matplotlib_widgets.NavToolbar(self.chan_tab['waves_figure'],
                                                                    self.chan_tab['main_widget'],
                                                                    coordinates=False)
        self.chan_tab['wave_toolbar'].setIconSize(QtCore.QSize(15, 15))
        self.chan_tab['wave_toolbar'].setOrientation(QtCore.Qt.Vertical)
        self.chan_tab['wave_toolbar'].setMaximumWidth(30)
        hlay.addWidget(self.chan_tab['waves_figure'])
        hlay.addWidget(self.chan_tab['wave_toolbar'])
        hlay.setContentsMargins(0, 0, 0, 0)
        hlay.setSpacing(1)
        vlay.addLayout(hlay)

        #------------------------------------------------------ UNIT TABS WIDGET
        self.chan_tab['unit_tabs_widget'] = QtWidgets.QTabWidget()
        self.chan_tab['unit_tab_bar_widget'] = self.chan_tab['unit_tabs_widget'].tabBar()
        screen = QtWidgets.QApplication.primaryScreen().availableGeometry()
        self.chan_tab['unit_tabs_widget'].setMaximumHeight(int(screen.height() / 4))
        self.chan_tab['unit_figures'] = {}
        self.chan_tab['del_unit_buttons'] = {}
        self.chan_tab['unit_count_label'] = {}
        self.chan_tab['unit_buttons'] = {}
        self.chan_tab['plot_raw_check'] = {}
        self.chan_tab['is_multiunit_check'] = {}
        self.chan_tab['unit_tabs_widget'].currentChanged.connect(self.change_current_unit)
        vlay.addWidget(self.chan_tab['unit_tabs_widget'])

        main_h_lay.addLayout(vlay)

        # configures the waveforms figure
        wavesfig = self.chan_tab['waves_figure'].figure
        ax = wavesfig.add_subplot(111)
        self.trim_waveforms_rect = matplotlib_widgets.MyRectangleSelector(ax,
                                                                        self.trim_waveforms,
                                                                        drawtype='line',
                                                                        useblit=True)
        self.trim_waveforms_rect.set_active(False)
        ax.set_facecolor('k')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        self.sample_waveform, = ax.plot([], color=[.5, .5, .5], linewidth=2)
        self.waveforms, = ax.plot([], animated=True)
        ax.set_ylim(-1000, 1000)
        ax.set_xlim(0, 32)

        # Create Slice plots
        self.slice_1_line = ax.axvline(0, color=[.5, .5, .5])
        self.slice_2_line = ax.axvline(0, color=[.5, .5, .5], linestyle='--')
        ax.grid()
        wavesfig.canvas.mpl_connect('draw_event', self.draw_callback)

        #------------------------------------------------- FEATURES PLOT WIDGET
        main_right_lay = QtWidgets.QVBoxLayout()

        tab = QtWidgets.QTabWidget()

        widget = QtWidgets.QWidget()
        # function buttons on top of the features plot:

        vlay = QtWidgets.QVBoxLayout(widget)

        hlay = QtWidgets.QHBoxLayout()
        hlay.addStretch(1)

        self.add_unit_btn = QtWidgets.QPushButton('Add Unit')
        self.add_unit_btn.setMaximumHeight(20)
        self.add_unit_btn.clicked.connect(self.add_unit)
        hlay.addWidget(self.add_unit_btn)

        # add a "keep" button
        self.keep_btn = QtWidgets.QPushButton('Keep')
        self.keep_btn.setMaximumHeight(20)
        self.keep_btn.setToolTip('Create new unit (only when All waveforms or unsorted are plotted)')
        self.keep_btn.clicked.connect(self.keep)
        hlay.addWidget(self.keep_btn)

        # add an "add region" button
        self.add_region_btn = QtWidgets.QPushButton('Add Region')
        self.add_region_btn.setMaximumHeight(20)
        self.add_region_btn.setToolTip('Add waveforms to the current unit')
        self.add_region_btn.clicked.connect(self.add_region)
        hlay.addWidget(self.add_region_btn)

        # add a "remove region" button
        self.remove_region_btn = QtWidgets.QPushButton('Remove Region')
        self.remove_region_btn.setMaximumHeight(20)
        self.remove_region_btn.clicked.connect(self.remove_region)
        hlay.addWidget(self.remove_region_btn)

        # "set valid waveforms" button
        self.valid_wf_btn = QtWidgets.QPushButton('Set Valid Wfs')
        self.valid_wf_btn.setMaximumHeight(20)
        self.valid_wf_btn.clicked.connect(self.validate_wfs)
        hlay.addWidget(self.valid_wf_btn)

        # "set valid waveforms" button
        self.replot_density_btn = QtWidgets.QPushButton('Replot Density')
        self.replot_density_btn.setMaximumHeight(20)
        self.replot_density_btn.clicked.connect(self.replot_density)
        hlay.addWidget(self.replot_density_btn)

        hlay.addStretch(1)
        vlay.addLayout(hlay)

        # Features figure and toolbar
        self.chan_tab['features_fig'] = matplotlib_widgets.MplWidget()
        self.chan_tab['features_fig'].figure.set_facecolor('k')
        self.chan_tab['features_fig_toolbar'] = matplotlib_widgets.NavToolbar(self.chan_tab['features_fig'].figure.canvas,
                                                                       self.chan_tab['main_widget'])
        self.chan_tab['features_fig_toolbar'].setIconSize(QtCore.QSize(15, 15))
        self.chan_tab['features_fig_toolbar'].setMaximumHeight(30)

        vlay.addWidget(self.chan_tab['features_fig'])
        vlay.addWidget(self.chan_tab['features_fig_toolbar'])

        vlay.setContentsMargins(0, 0, 0, 0)
        vlay.setSpacing(1)

        tab.addTab(widget, '2D')
        main_right_lay.addWidget(tab)

        #-------------------------------------------------------- 3D TAB Widget
        self.widget_3d = gl.GLViewWidget()
        tab.addTab(self.widget_3d, '3D')

        #---------------------------------- Spikes vs time visualization widget
        # add a figure adn axes
        self.time_scroll = {}
        self.time_scroll['figure'] = matplotlib_widgets.MplWidget()
        self.time_scroll['figure'].figure.set_facecolor('k')
        self.time_scroll['draw_fig_cid'] = self.time_scroll['figure'].figure.canvas.mpl_connect('draw_event', self.draw_scroll_fig)
        self.time_scroll['figure'].setMaximumHeight(int(screen.height() / 6))
        self.time_scroll['Ax'] = self.time_scroll['figure'].figure.add_subplot(111)
        self.time_scroll['Ax'].set_facecolor('k')
        self.time_scroll['Ax'].set_ylim(-1500, 1500)
        self.time_scroll['Ax'].set_xticklabels([])
        self.time_scroll['Ax'].set_yticklabels([])
        self.time_scroll['Ax'].set_axis_off()
        self.time_scroll['Plot'], = self.time_scroll['Ax'].plot([], color=[.5, .5, .5])
        self.time_scroll['figure'].figure.tight_layout()
        self.time_scroll['figure'].figure.canvas.draw()

        # add a vertical zoom slider
        self.time_scroll['v_zoom_slider'] = QtWidgets.QSlider(QtCore.Qt.Vertical)
        self.time_scroll['v_zoom_slider'].setMaximumHeight(int(screen.height() / 6))
        self.time_scroll['v_zoom_slider'].setMinimum(100)
        self.time_scroll['v_zoom_slider'].setMaximum(5000)
        self.time_scroll['v_zoom_slider'].setValue(1000)
        self.time_scroll['v_zoom_slider'].valueChanged.connect(self.v_zoom)
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(self.time_scroll['v_zoom_slider'])
        hlay.addWidget(self.time_scroll['figure'])
        main_right_lay.addLayout(hlay)

        # add an horizontal zoom slider
        self.time_scroll['h_zoom_slider'] = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.time_scroll['h_zoom_slider'].setRange(5, 5000)
        self.time_scroll['h_zoom_slider'].setValue(500)
        self.time_scroll['h_zoom_slider'].setSingleStep(5)
        self.time_scroll['h_zoom_slider'].valueChanged.connect(self.h_zoom)
        self.time_scroll['h_zoom_spin'] = QtWidgets.QSpinBox()
        self.time_scroll['h_zoom_spin'].setMinimumWidth(80)
        self.time_scroll['h_zoom_spin'].setMaximumHeight(20)
        self.time_scroll['h_zoom_spin'].setRange(5, 5000)
        self.time_scroll['h_zoom_spin'].setValue(500)
        self.time_scroll['h_zoom_spin'].setSingleStep(10)
        self.time_scroll['h_zoom_spin'].valueChanged.connect(self.time_scroll['h_zoom_slider'].setValue)
        self.time_scroll['h_zoom_slider'].valueChanged.connect(self.time_scroll['h_zoom_spin'].setValue)
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('H Span '))
        hlay.addWidget(self.time_scroll['h_zoom_spin'])
        hlay.addWidget(self.time_scroll['h_zoom_slider'])
        main_right_lay.addLayout(hlay)

        # add a time slider
        self.time_scroll['h_scroll_slider'] = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.time_scroll['h_scroll_slider'].setRange(0, 3000000)
        self.time_scroll['h_scroll_slider'].setSingleStep(int(self.time_scroll['h_zoom_slider'].value() / 10))
        self.time_scroll['h_scroll_slider'].valueChanged.connect(self.h_scroll)
        self.time_scroll['h_spin'] = QtWidgets.QSpinBox()
        self.time_scroll['h_spin'].setRange(0, 3000000)
        self.time_scroll['h_spin'].setMinimumWidth(80)
        self.time_scroll['h_spin'].setMaximumHeight(20)
        self.time_scroll['h_spin'].valueChanged.connect(self.time_scroll['h_scroll_slider'].setValue)
        self.time_scroll['h_scroll_slider'].valueChanged.connect(self.time_scroll['h_spin'].setValue)
        hlay = QtWidgets.QHBoxLayout()
        hlay.addWidget(QtWidgets.QLabel('H Scroll'))
        hlay.addWidget(self.time_scroll['h_spin'])
        hlay.addWidget(self.time_scroll['h_scroll_slider'])
        main_right_lay.addLayout(hlay)

        main_right_lay.setContentsMargins(0, 0, 0, 0)
        main_right_lay.setSpacing(1)
        # add the widget to the main horizontal layout
        main_h_lay.addLayout(main_right_lay)
        main_h_lay.setContentsMargins(1, 1, 1, 1)
        self.chan_tab['main_widget'].setLayout(main_h_lay)

        # create a generic Msg box
        self.msg_box = QtWidgets.QMessageBox()

        # if running in linux set a certain style for the buttons and widgets
        if sys.platform == 'linux2':
            QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Plastique'))

        # add the main tabbed figures widget to the main splitter
        split1.addWidget(self.main_fig_tab_widget)

        # add the splitter to the main layout
        self.main_layout.addWidget(split1)

        # set the layout of the main widget
        #self.main_widget.setLayout(self.main_layout)

        # set the central widget of the application
        self.setCentralWidget(self.main_widget)

        # finally show the object
        self.show()

    #__________________________________________________________________________
    def save_overview_fig(self):
        if self.h5_file_loaded:
            fname = QtWidgets.QFileDialog.getSaveFileName(
                directory=self.h5file.filename[0:-3] + '_sorted.png')[0]
            if fname:
                self.overview_tab_1['figure'].figure.savefig(fname,
                                                           dpi=300,
                                                           facecolor='k')

    #__________________________________________________________________________
    def adjust_plots(self):
        self.time_scroll['figure'].figure.tight_layout()
        self.time_scroll['figure'].figure.canvas.draw()

        self.chan_tab['waves_figure'].figure.tight_layout()
        self.chan_tab['waves_figure'].figure.canvas.draw()

        if len(self.chan_tab['features_fig'].figure.axes) > 0:
            self.chan_tab['features_fig'].figure.tight_layout()
            self.chan_tab['features_fig'].figure.canvas.draw()

        if len(self.overview_tab_1['figure'].figure.axes) > 0:
            self.overview_tab_1['figure'].figure.tight_layout()
            self.overview_tab_1['figure'].figure.canvas.draw()

        for k in self.chan_tab['unit_figures']:
            self.chan_tab['unit_figures'][k].figure.tight_layout()
            self.chan_tab['unit_figures'][k].figure.canvas.draw()

    #__________________________________________________________________________
    def set_wf_plot_lim(self):
        sender = self.sender()
        ax = self.chan_tab['waves_figure'].figure.axes[0]
        lim = sender.value()
        ax.set_ylim(-lim, lim)
        self.chan_tab['waves_figure'].figure.canvas.draw()

    #__________________________________________________________________________
    def h_scroll(self):
        ''' This function gets triggered whenever the user moves the bottom
        scrollbar in the lower right. It helps to explore the raw waveforms'''

        if not self.chan_plotted:
            return

        self.time_scroll['figure'].figure.canvas.restore_region(self.time_scroll['bg'])
        self.chan_tab['waves_figure'].figure.canvas.restore_region(self.chan_tab['waves_fig_bg'])

        v = self.time_scroll['h_scroll_slider'].value()
        h = self.time_scroll['h_zoom_slider'].value()
        indx = np.flatnonzero(np.logical_and(self.cur_ts >= v, self.cur_ts < (v + h)))

        if any(indx):
            # ontain the timestamps corresponding to the indexes
            ts = self.cur_node.TimeStamp[indx]
            # substract the first timestamp to have a 0 based indexing
            ts = ts - v

            # obtain the waveforms to plot
            wf = self.cur_node.Waveforms[indx, :]

            # obtain the length of units to plot
            n = len(indx)

            # create an array of Nones to append
            nones = np.array(n * [None], ndmin=2).transpose()

            # append nones to the waveforms array and reshape it to a vector
            wf = np.append(wf, nones, axis=1).reshape((n * (self.wf_size + 1),))

            # create a vector time, based on the sampling frequency, the
            # the number of points per spike and the timestamp
            ts = np.tile(ts, (self.wf_size, 1)).T + \
                np.tile(np.linspace(0, self.end, self.wf_size), (n, 1))
            ts = np.append(ts, nones, axis=1).reshape((n * (self.wf_size + 1),))

            # set the plot data to the created arrays
            self.time_scroll['Plot'].set_data(ts, wf)

            # set axes limits
            self.time_scroll['Ax'].set_xlim(0, h)
            self.time_scroll['Ax'].draw_artist(self.time_scroll['Plot'])
            self.sample_waveform.set_data(self.waveform_x_axis * n, wf)
            self.chan_tab['waves_figure'].figure.axes[0].draw_artist(self.sample_waveform)

        self.time_scroll['figure'].figure.canvas.blit(self.time_scroll['figure'].figure.bbox)
        self.chan_tab['waves_figure'].figure.canvas.blit(self.chan_tab['waves_figure'].figure.axes[0].bbox)

    #__________________________________________________________________________
    def v_zoom(self):
        v = self.time_scroll['v_zoom_slider'].value()
        self.time_scroll['Ax'].set_ylim(-v, v)
        self.time_scroll['figure'].figure.canvas.restore_region(self.time_scroll['bg'])
        self.time_scroll['Ax'].draw_artist(self.time_scroll['Plot'])
        self.time_scroll['figure'].figure.canvas.blit(self.time_scroll['figure'].figure.bbox)

    #__________________________________________________________________________
    def h_zoom(self):
        v = self.time_scroll['h_zoom_slider'].value()
        self.time_scroll['h_scroll_slider'].setSingleStep(max(1, int(v / 10)))
        self.h_scroll()

    #__________________________________________________________________________
    def draw_scroll_fig(self, event):
        fig = self.time_scroll['figure'].figure
        self.time_scroll['bg'] = fig.canvas.copy_from_bbox(fig.axes[0].bbox)

    #__________________________________________________________________________
    def load_h5_file(self, h5file=None):
        ''' Loads an h5 file that contains all the information about the units:
        waveforms and timestamps '''

        # try to load an h5 file
        #if settings.working_dir:
        #    d = settings.working_dir
        #else:
            #d = ''

        if not h5file:
            h5file = QtWidgets.QFileDialog.getOpenFileName(parent=self,
                                                      caption='Select an H5 File',
                                                      directory="",
                                                      filter='*.h5')[0]

        # in case there is not file selected
        if not h5file:
            return

        # set file loaded var = True
        if hasattr(self, 'h5_file_loaded') and self.h5_file_loaded:
            self.h5file.close()

        # try to open the file
        try:
            self.h5file = tables.open_file(str(h5file), mode='r+')
        except:
            self.msg_box.setIcon(QtWidgets.QMessageBox.Warning)
            self.msg_box.setText('There was a problem opening the H5 file')
            self.msg_box.setWindowTitle('Warning')
            self.msg_box.show()
            return

        # set file loaded var = True
        self.h5_file_loaded = True
        self.file_path = os.path.split(h5file)[0]

        # REPAIR THE H5FILE STRUCTURE

        if self.h5file.__contains__('/Chans'):
            self.h5file.rename_node('/', 'Spikes', name='Chans')

        chan_nodes = self.h5file.list_nodes('/Spikes')

        for k in chan_nodes:
            where = '/Spikes/%s' % k._v_name
            for n in k:
                if re.search('Unit[0-9]{2}(?!_isMultiunit)', n._v_name) and n._c_classid != 'GROUP':
                    unit_name = re.search('Unit[0-9]{2}(?!_isMultiunit)', n._v_name).group()

                    self.h5file.create_group(where=where,
                                             name=unit_name + '_grp')

                    self.h5file.move_node(where=where,
                                         name=unit_name,
                                         newparent='/Spikes/%s/%s' % (k._v_name, unit_name + '_grp'),
                                         newname='Indx')

                    self.h5file.rename_node(where=where,
                                           name=unit_name + '_grp',
                                           newname=unit_name)

                elif re.search('Unit[0-9]{2}_isMultiunit', n._v_name):
                    self.h5file.remove_node(where=where,
                                            name=re.search('Unit[0-9]{2}_isMultiunit', n._v_name).group())

                elif 'tmp' in n._v_name:
                    self.h5file.remove_node(where=where, name=n._v_name, recursive=True)

        # CREATE 'isMultiunit' and 'isBursting' fields
        chan_nodes = self.h5file.list_nodes('/Spikes')

        for k in chan_nodes:
            node = '/Spikes/%s' % k._v_name
            for n in k:
                if 'Unit' in n._v_name and n._c_classid == 'GROUP':
                    parent = node + '/' + n._v_name
                    if not self.h5file.__contains__(parent + '/' + 'isMultiunit'):
                        self.h5file.create_array(parent, 'isMultiunit', False)
                    if not self.h5file.__contains__(parent + '/' + 'isBursting'):
                        self.h5file.create_array(parent, 'isBursting', False)

        # RENAME the "Indexes" field to "Indx"
        chan_nodes = self.h5file.list_nodes('/Spikes')

        for k in chan_nodes:
            for n in k:
                if 'Unit' in n._v_name:
                    node_name = '/Spikes/%s/%s' % (k._v_name, n._v_name)
                    for l in n:
                        if l._v_name == 'Indexes':
                            self.h5file.rename_node(node_name, 'Indx', 'Indexes')

        # save changes to disk
        self.h5file.flush()

        # REPAIR UNIT NAMES #####
        chan_nodes = self.h5file.list_nodes('/Spikes')

        for chan in chan_nodes:
            unit_names = [k for k in chan.__members__ if 'Unit' in k]
            unit_names.sort()
            for j, k in enumerate(unit_names):
                if k != 'Unit%02d' % j:
                    self.h5file.rename_node('/Spikes/%s' % chan._v_name, name=k,
                                           newname='Unit%02d' % j)

        # save changes to disk
        self.h5file.flush()

        # clean the channel figures if something already plotted
        if hasattr(self, 'chan_plotted') and self.chan_plotted:
            self.reset_channel_tab()

        self.plot_overview()

        # clear the Log Browser and load Log info
        self.log_combo.clear()
        self.log_text_browser.clear()
        if self.h5file.__contains__('/Log'):
            nodes = self.h5file.list_nodes('/Log')
            node_names = [k._v_name for k in nodes]
            self.log_combo.addItems(node_names)

        # set window title = to filename
        self.setWindowTitle('Spike Sorter GUI ' + h5file)

    #__________________________________________________________________________
    def plot_overview(self):
        ''' plot an overview of 1000 spikes per channel.
        Also, fills the overview table with the general information about each
        channel'''

        # get the list of nodes inside the "Chans" group
        chan_nodes = self.h5file.list_nodes('/Spikes')

        # get the number of the channels in the file
        self.chans_list = [int(re.search('(?<=Chan_)[0-9]{3}', k._v_name).group()) for k in chan_nodes]

        # get the waveform size (number of points). X is for fast plotting
        self.wf_size = self.h5file.root.Header.WaveformSize.read()
        #x = range(self.wf_size)

        # add items to the channel selector in the toolbar
        self.chan_selector.clear()
        self.chan_selector.addItems(['%d' % k for k in self.chans_list])

        # clean the overview figure
        self.overview_tab_1['figure'].figure.clear()

        # helper to calculate the geometry of the axes
        fig_rows = np.ceil(len(chan_nodes) / 10.0)

        # clear contents of the overview table
        self.overview_tab_2['overview_table'].clearContents()
        c = list(range(self.overview_tab_2['overview_table'].rowCount()))
        c.reverse()
        for k in c:
            self.overview_tab_2['overview_table'].removeRow(k)

        # iterate over the list of channels to add the information to the table
        for j, k in enumerate(chan_nodes):

            # update overveiew table
            self.overview_tab_2['overview_table'].insertRow(j)
            self.overview_tab_2['overview_table'].setRowHeight(j, 20)

            # add an event count
            self.overview_tab_2['overview_table'].setItem(j, 0, QtWidgets.QTableWidgetItem(str(k.TimeStamp.nrows)))

            # add an "isTrash" checkbox
            check = QtWidgets.QCheckBox()
            check.setProperty('Data', self.chans_list[j])
            check.stateChanged.connect(self.set_trash)
            self.overview_tab_2['overview_table'].setCellWidget(j, 1, check)

            # add an "isMultinunit" checkbox
            is_multiunit_check = QtWidgets.QCheckBox()
            is_multiunit_check.setObjectName(k._v_name)
            is_multiunit_check.stateChanged.connect(self.is_multiunit)
            self.overview_tab_2['overview_table'].setCellWidget(j, 2, is_multiunit_check)

            # add information about unsorted units
            if k.__contains__('Unsorted'):
                self.overview_tab_2['overview_table'].setItem(j, 4, QtWidgets.QTableWidgetItem(str(k.Unsorted.nrows)))

            # add information about valif waveforms
            if k.__contains__('ValidWFs'):
                self.overview_tab_2['overview_table'].setItem(j, 5, QtWidgets.QTableWidgetItem(str(k.ValidWFs.nrows)))

            # add info about each unit
            units = [m for m in k.__members__ if re.search('Unit[0-9]{2}', m)]  # obtain unit names
            units.sort()
            if units:  # in case there are units
                for m, n in enumerate(units):
                    if self.overview_tab_2['overview_table'].columnCount() <= (m + 6):
                        self.overview_tab_2['overview_table'].insertColumn(self.overview_tab_2['overview_table'].columnCount())
                        n_cols = self.overview_tab_2['overview_table'].columnCount()
                        self.overview_tab_2['overview_table'].setColumnWidth(n_cols - 1, 65)
                        self.overview_tab_2['overview_table'].setHorizontalHeaderItem(n_cols - 1,
                                                                                   QtWidgets.QTableWidgetItem('Unit%02d' % m))

                    self.overview_tab_2['overview_table'].setItem(j, m + 6,
                                                               QtWidgets.QTableWidgetItem(str(k.__getattr__(n).Indx.nrows)))

            # Create the axes to plot the waveforms
            self.overview_tab_1['figure'].figure.add_subplot(int(fig_rows), 10, int(j + 1))
            self.overview_tab_1['figure'].figure.axes[j].set_yticks([], [])  # eliminate the ticks to have more space
            self.overview_tab_1['figure'].figure.axes[j].set_xticks([], [])  # eliminate the ticks to have more space
            self.overview_tab_1['figure'].figure.axes[j].set_axis_off()
            self.overview_tab_1['figure'].figure.axes[j].set_title('Ch %d' % (self.chans_list[j]),
                                                                 fontsize=10,
                                                                 fontdict={'color': 'w'})

            self.plot_chan_overview(k, axes_to_plot=self.overview_tab_1['figure'].figure.axes[j])

            # check the isTrash widgets and make the axes background yellow
            if k.__contains__('isTrash'):
                if k.isTrash.read():
                    check.setCheckState(2)
                    self.overview_tab_1['figure'].figure.axes[j].set_facecolor([.5, .5, .5])

            if k.__contains__('isMultiunit'):
                if k.isMultiunit.read():
                    is_multiunit_check.setCheckState(2)

        # set the names of the vertical headers
        self.overview_tab_2['overview_table'].setVerticalHeaderLabels(['Ch ' + str(k) for k in self.chans_list])

        # set alternating row colors
        self.overview_tab_2['overview_table'].setAlternatingRowColors(True)

        # connect the clicks on this canvas with the channel select function
        self.overview_tab_1['figure'].figure.canvas.mpl_connect('button_release_event', self.select_channel)

        # tight layout and draw
        self.overview_tab_1['figure'].figure.tight_layout()
        self.overview_tab_1['figure'].figure.canvas.draw()

        # get the sampling frequency
        self.sf = float(self.h5file.root.Header.TimeStamp_Res.read())

        self.step = self.wf_size + 1

        # set boolean variable

    #__________________________________________________________________________
    def plot_chan_overview(self, node, axes_to_plot):
        '''Helper function that plots the unsorted as well as the sorted events
        in a given axes on the overview figure'''

        # get the number of events to plot
        n_events = self.overview_n_events_spin.value()
        waveforms = node.Waveforms.read()

        # iterate over the members of a node
        for k in node:
            if not re.search('unsorted|Unit[0-9]{2}', k._v_name):
                continue
            
            #print(k._v_children().keys())
            
            # read the indices first:
            if 'Unit' in k._v_name and 'Indx' in k:
                if k.Indx.nrows >= n_events:
                    indx = k.Indx.read(start=0, stop=k.Indx.nrows, step=int(k.Indx.nrows / n_events))
                else:
                    indx = k.Indx.read()
            elif 'Unsorted' in k._v_name:
                if k.nrows >= n_events:
                    indx = k.read(start=0, stop=k.nrows, step=int(k.nrows / n_events))
                else:
                    indx = k.read()

            # obtain the waveforms
            wf = waveforms[indx, :]
            if not wf.any():
                continue

            # faster plotting strategy:

            # obtain the length of units to plot
            n = len(indx)

            # create an array of Nones to append
            nones = np.array(n * [None], ndmin=2).T

            # append nones to the waveforms array and reshape it to a vector
            wf = np.append(wf, nones, axis=1).reshape((n * (self.wf_size + 1),))

            # create the x axis
            ts = list(range(self.wf_size))
            ts.append(None)

            # choose the color and the zorder according to the type of unit
            if k._v_name == 'Unsorted':
                color = 'w'
                zorder = 1
            elif 'Unit' in k._v_name:
                # get the unit number
                zorder = int(re.search('[0-9]{2}', k._v_name).group())
                color = self.unit_colors[zorder, :]
                zorder = 100 - zorder

            # get the list of plots in the particular axes
            l = [l for l in axes_to_plot.lines if str(l.get_label()) == k._v_name]

            # if a plot with a label equal to the name of the unit exist, the update the data
            if len(l) > 0:
                l[0].set_data(ts * n, wf)
            # if not create one
            else:
                axes_to_plot.plot(ts * n, wf, color=color, rasterized=True,
                               alpha=0.5, label=k._v_name, zorder=zorder)

        # set the limits of the axes
        axes_to_plot.set_ylim(-self.overview_y_lims_spin.value(), self.overview_y_lims_spin.value())

        # add a small text box with the event count
        bbox_props = dict(boxstyle='round', fc='0.75', ec='0.25', alpha=0.8)
        axes_to_plot.text(0.5, 0, 'Count: %d' % node.TimeStamp.nrows, transform=axes_to_plot.transAxes,
                       color='k', bbox=bbox_props, size=10, ha='center')

    #__________________________________________________________________________
    def change_overview_y_lim(self):

        if not self.h5_file_loaded:
            return

        lim = self.overview_y_lims_spin.value()
        for k in self.overview_tab_1['figure'].figure.axes:
            k.set_ylim(-lim, lim)

        self.overview_tab_1['figure'].figure.canvas.draw()

    #__________________________________________________________________________
    def plot_chan(self):

        # exit if ther is no H5 file loaded
        if not self.h5_file_loaded:
            return

        # clean the channels tab
        self.reset_channel_tab()

        # reset Units list
        self.units_list = []

        #pdb.set_trace()

        # load waveforms for a specific channel
        self.cur_chan = int(self.chan_selector.currentText())
        #n_spikes = self.n_spikes_slider.value()
        self.cur_node_name = '/Spikes/Chan_%03d' % self.cur_chan
        self.cur_node = self.h5file.get_node(self.cur_node_name)
        self.cur_waveforms = self.cur_node.Waveforms.read()
        self.cur_ts = self.cur_node.TimeStamp.read()
        self.time_scroll['h_scroll_slider'].setMaximum(int(self.cur_ts[-1]))
        self.time_scroll['h_spin'].setMaximum(int(self.cur_ts[-1]))
        self.unit_nodes = [k for k in self.h5file.list_nodes(self.cur_node_name) if re.search('Unit[0-9]{2}', k._v_name)]

        # get the indices of the unsorted. If there are no, create one
        if not self.cur_node.__contains__('Unsorted'):
            self.unsorted = np.arange(len(self.cur_ts))
            self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)
        else:
            self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()

        #set the unit names in the combo box
        self.what_to_plot.clear()
        self.what_to_plot.addItems(['All waveforms', 'Sorted', 'Unsorted'])
        if self.unit_nodes:
            self.what_to_plot.addItems([k._v_name for k in self.unit_nodes])

        # set the axis limits to apropriately view the unit
        v = self.wave_ax_y_lim_spin.value()
        self.chan_tab['waves_figure'].figure.axes[0].set_ylim(-v, v)

        # get the waveform size for this specific waveform
        self.wf_size = self.h5file.root.Header.WaveformSize.read()

        # cheack whether to plot raw waveforms
        '''
        if self.plot_waveforms_check.checkState() == 2:
            if self.cur_waveforms.shape[1] < 10000:
                self.waveforms_to_plot = self.cur_waveforms
            else:
                indx = np.int32(np.linspace(0,self.cur_waveforms.shape[0]-1,10000))
                self.waveforms_to_plot = self.cur_waveforms[indx,:]

            for k in range(self.waveforms_to_plot.shape[0]):
                self.waveforms.set_data(range(self.wf_size), self.waveforms_to_plot[k,:])
                self.chan_tab['waves_figure'].figure.axes[0].draw_artist(self.waveforms)
                self.chan_tab['waves_figure'].figure.canvas.blit(self.chan_tab['waves_figure'].figure.axes[0].bbox)'''

        # grab background from the waveforms figure to make animations
        self.chan_tab['waves_fig_bg'] = self.chan_tab['waves_figure'].figure.canvas.copy_from_bbox(self.chan_tab['waves_figure'].figure.axes[0].bbox)
        self.main_fig_tab_widget.setTabText(2, 'Chan %02d' % self.cur_chan)

        # calculate pca
        pc = pca_scores(self.cur_waveforms)

        # put data in a KDTree to easily calculate distance with the cursor
        self.xy_data = cKDTree(pc[:, 0:2], 1000)
        self.chan_tab['pca'] = pc

        # set the internal variable to true
        self.chan_plotted = True

        # copy the number of events in the channel into a label to see if the user
        # wants to decimate or plot them all
        self.n_pts_label.setText(str(self.cur_ts.size))

        # read the plotting parameters in the "Chan Tools" tab to plot
        # the selected feature
        self.plot_features()

        if self.chan_tab['unit_tabs_widget'].count() > 0:
            self.chan_tab['unit_tabs_widget'].setCurrentIndex(0)
            self.cur_unit_name = str(self.chan_tab['unit_tabs_widget'].tabText(0))
            self.cur_unit = int(re.search('(?<=Unit)[0-9]{2}', self.cur_unit_name).group())

        self.waveform_x_axis = list(range(self.wf_size))
        self.waveform_x_axis.append(None)
        self.end = 1000 * self.wf_size / self.sf

        # save h5file changes to disk
        self.h5file.flush()

    #__________________________________________________________________________
    def set_trash(self):
        sender = self.sender()
        try:
            chan = sender.property('Data').toPyObject()
        except:
            chan = sender.property('Data')

        node_name = '/Spikes/Chan_%03d' % chan

        indx = self.chans_list.index(chan)

        if self.h5file.get_node(node_name).__contains__('isTrash'):
            self.h5file.remove_node(node_name, 'isTrash')

        if sender.checkState() in [1, 2]:
            self.h5file.create_array(node_name, 'isTrash', True)
            self.overview_tab_1['figure'].figure.axes[indx].set_facecolor('y')
        elif sender.checkState() == 0:
            self.h5file.create_array(node_name, 'isTrash', False)
            self.overview_tab_1['figure'].figure.axes[indx].set_facecolor('w')

        # save changes to disk
        self.h5file.flush()

    #__________________________________________________________________________
    def is_multiunit(self):
        sender = self.sender()
        node_name = '/Spikes/%s' % sender.objectName()

        if self.h5file.get_node(node_name).__contains__('isMultiunit'):
            self.h5file.remove_node(node_name, 'isMultiunit')

        if sender.checkState() in [1, 2]:
            self.h5file.create_array(node_name, 'isMultiunit', True)
        elif sender.checkState() == 0:
            self.h5file.create_array(node_name, 'isMultiunit', False)

        # save changes to disk
        self.h5file.flush()

    #__________________________________________________________________________
    def trash_chans(self):
        '''Utility function to mark the channels with fewer than a defined number
        of units'''

        # check whether an h5file has been loaded
        if not self.h5_file_loaded:
            return

        # obtain parameters
        n = self.mark_trash_spin.value()
        chans = self.h5file.list_nodes('/Spikes')

        # iterate over nodes in h5file; if fewer than n mark as trash
        for l, k in enumerate(chans):
            j = int(re.search('(?<=Chan_)[0-9]{3}', k._v_name).group())
            if k.TimeStamp.nrows < n:
                self.overview_tab_1['figure'].figure.axes[l].set_facecolor('y')
                self.overview_tab_2['overview_table'].cellWidget(l, 1).setChecked(True)
                if self.h5file.get_node('/Spikes', 'Chan_%03d' % j).__contains__('isTrash'):
                    self.h5file.remove_node('/Spikes/Chan_%03d' % j, 'isTrash')
                self.h5file.create_array('/Spikes/Chan_%03d' % j, 'isTrash', True)

            else:
                self.overview_tab_1['figure'].figure.axes[l].set_facecolor('w')
                self.overview_tab_2['overview_table'].cellWidget(l, 1).setChecked(False)
                if self.h5file.get_node('/Spikes', 'Chan_%03d' % j).__contains__('isTrash'):
                    self.h5file.remove_node('/Spikes/Chan_%03d' % j, 'isTrash')
                self.h5file.create_array('/Spikes/Chan_%03d' % j, 'isTrash', False)

        # save changes to disk
        self.h5file.flush()

        #update the overview
        self.overview_tab_1['figure'].figure.canvas.draw()

    #__________________________________________________________________________
    def delete_trash_chans(self):
        # check whether an h5file has been loaded
        if not self.h5_file_loaded:
            return

        chans = self.h5file.list_nodes('/Spikes')
        chans.reverse()
        n = list(range(len(chans)))
        n.reverse()

        delchans = []
        for j, k in zip(n, chans):
            state = self.overview_tab_2['overview_table'].cellWidget(j, 1).checkState()

            if state == 2:
                delchans.append(k._v_name)
                self.h5file.remove_node('/Spikes', k._v_name, recursive=True)

        if len(delchans) > 0:
            self.add_log('Deleted channels: ' + str(delchans))

        self.plot_overview()

    #__________________________________________________________________________
    def reset_chan(self):
        ''' check whether an h5file has been loaded '''
        if not self.h5_file_loaded or not self.chan_plotted:
            return

        for k in self.h5file.list_nodes(self.cur_node_name):
            if k._v_name not in ['Waveforms', 'TimeStamp', 'isTrash']:
                self.h5file.remove_node(self.cur_node_name, k._v_name, recursive=True)

        self.plot_chan()
        self.add_log('%s resetted' % self.cur_node_name)

    #__________________________________________________________________________
    def add_log(self, message):
        ''' add log to to keep a history of changes to the file '''

        if not self.h5_file_loaded:
            return

        if not self.h5file.__contains__('/Log'):
            self.h5file.create_group('/', 'Log', 'History of changes')
        name = 'Entry_%s_%s_%s_%s_%s_%s' % datetime.datetime.now().timetuple()[0:6]
        self.h5file.create_array('/Log', name, message)

        # save changes to disk
        self.h5file.flush()

        #add the item to the log browser
        self.log_combo.addItem(name)

    #__________________________________________________________________________
    def set_log_text(self):
        if self.log_combo.currentIndex() == -1:
            return
        node = str(self.log_combo.currentText())
        if node:
            log = self.h5file.get_node('/Log', node).read()
            #self.log_text_browser.setText(log)

    #__________________________________________________________________________
    def close_file(self):
        ''' close the h5 file'''
        if not self.h5_file_loaded:
            return
        self.h5file.flush()
        self.h5file.close()
        self.h5_file_loaded = False
        print('h5 File closed')

    #__________________________________________________________________________
    def select_channel(self, event):
        ''' selects a channel when axes are clicked'''
        if event.in_axes:
            chan = int(re.search('(?<=Ch )[0-9]{1,3}', event.in_axes.get_title()).group())
            c = [int(self.chan_selector.itemText(k)) for k in range(self.chan_selector.count())].index(chan)
            self.chan_selector.setCurrentIndex(c)

    #__________________________________________________________________________
    def table_row_changed(self, sel):
        self.chan_selector.setCurrentIndex(sel)

    #__________________________________________________________________________
    def settings_dialog(self):
        ''' edit paths'''
        if settings.edit() == 1:
            self.working_dir = settings.working_dir

    #__________________________________________________________________________
    def autocorr_opts(self):
        autocorropts.show()

    #__________________________________________________________________________
    def about(self):
        ''' opens a small dialog with information about the software'''
        self.about_msg.show()

    #__________________________________________________________________________
    def nearest_point(self, event):
        ''' when right button clicked over the features window, calculates the closest
        point and plots its corresponding waveform'''

        if event.button == 3 and event.in_axes and self.chan_tab['features_fig_toolbar'].mode == '':
            featuresax = self.chan_tab['features_fig'].figure.axes[0]
            wavesax = self.chan_tab['waves_figure'].figure.axes[0]

            if self.plot_unit_counter >= self.n_units_spin.value():
                self.chan_tab['waves_figure'].figure.canvas.restore_region(self.chan_tab['waves_fig_bg'])
                self.plot_unit_counter = 0

            for k in self.chan_tab['features_fig_bg']:
                self.chan_tab['features_fig'].figure.canvas.restore_region(k)

            _, res = self.xy_data.query([event.xdata, event.ydata], 1)
            self.cursor.set_data(self.xy_data.data[res, 0], self.xy_data.data[res, 1])
            self.sample_waveform.set_data(list(range(self.wf_size)),
                                         self.cur_waveforms[self.data_index[res], :])

            featuresax.draw_artist(self.cursor)
            wavesax.draw_artist(self.sample_waveform)
            self.chan_tab['features_fig'].figure.canvas.blit(featuresax.bbox)
            self.chan_tab['waves_figure'].figure.canvas.blit(wavesax.bbox)
            self.plot_unit_counter += 1

    #__________________________________________________________________________
    def draw_callback(self, event):
        ''' any draw callback triggers the capture of the figure background for
        using it in the animations '''

        if not self.chan_plotted:
            return

        if event.canvas == self.chan_tab['features_fig'].figure.canvas:
            bg = []
            for k in self.chan_tab['features_fig'].figure.axes:
                bg.append(self.chan_tab['features_fig'].figure.canvas.copy_from_bbox(k.bbox))
            self.chan_tab['features_fig_bg'] = bg
        elif event.canvas == self.chan_tab['waves_figure'].figure.canvas:
            self.chan_tab['waves_fig_bg'] = self.chan_tab['waves_figure'].figure.canvas.copy_from_bbox(self.chan_tab['waves_figure'].figure.axes[0].bbox)

    #__________________________________________________________________________
    def plot_features(self):
        ''' determines what 2 plot based on the user choices'''

        # obtain labels and return if are the same
        xlabel = self.x_plot.currentText()
        ylabel = self.y_plot.currentText()
        if xlabel == ylabel:
            return

        curchan = int(self.chan_selector.currentText())

        if self.plot_valids_only_check.checkState() == 2 and \
           self.cur_node.__contains__('ValidWFs'):
            print('you selected to plot only the valid Wfs')

        what_to_plot = str(self.what_to_plot.currentText())  # string containing what to plot
        self.cur_node_name = '/Spikes/Chan_%03d' % curchan
        #unit_nodes = [k for k in self.h5file.list_nodes(self.cur_node_name) if re.search('Unit[0-9]{2}', k._v_name)]

        if what_to_plot in ['All waveforms', 'Sorted']:
            self.data_index = list(range(self.cur_ts.size))
            pc = self.chan_tab['pca']
        elif what_to_plot == 'Unsorted':
            self.data_index = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
            pc = pca_scores(self.cur_waveforms[self.data_index, :])
        elif re.search('Unit', what_to_plot):
            self.cur_unit_name = what_to_plot
            self.data_index = self.h5file.get_node(self.cur_node_name, what_to_plot).Indx.read()
            pc = pca_scores(self.cur_waveforms[self.data_index, :])

        # save what is the feature
        self.cur_feature_plot = what_to_plot

        # get the choice for the x axis
        if xlabel == 'PCA1':
            x = pc[:, 0]

        elif xlabel == 'PCA2':
            x = pc[:, 1]

        elif xlabel == 'PCA3':
            x = pc[:, 2]

        elif xlabel == 'Slice1':
            x = self.cur_waveforms[self.data_index, self.slice_spin_1.value()]
            x = x / 100.0

        elif xlabel == 'Slice2':
            x = self.cur_waveforms[self.data_index, self.slice_spin_2.value()]
            x = x / 100.0

        elif xlabel == 'Energy':
            x = np.sum(np.power(self.cur_waveforms[self.data_index, :], 2), axis=1)
            x = x / 1000000.0

        elif xlabel == 'Peak':
            x = self.cur_waveforms[self.data_index, :].max(axis=1)
            x = x / 100.0

        elif xlabel == 'Valley':
            x = self.cur_waveforms[self.data_index, :].min(axis=1)
            x = x / 100.0

        elif xlabel == 'Peak Pt':
            x = self.cur_waveforms[self.data_index, :].argmax(axis=1)

        elif xlabel == 'Valley Pt':
            x = self.cur_waveforms[self.data_index, :].argmin(axis=1)

        elif xlabel == 'Pk2Pk Amp':
            x = self.cur_waveforms[self.data_index, :].max(axis=1) - self.cur_waveforms[self.data_index, :].min(axis=1)
            x = x / 100.0

        elif xlabel == 'Time':
            x = self.cur_ts[self.data_index]
            x = x / 60000.0

        # get the choice for the y axis
        if ylabel == 'PCA1':
            y = pc[:, 0]

        elif ylabel == 'PCA2':
            y = pc[:, 1]

        elif ylabel == 'PCA3':
            y = pc[:, 2]

        elif ylabel == 'Slice1':
            y = self.cur_waveforms[self.data_index, self.slice_spin_1.value()]
            y = y / 100.0

        elif ylabel == 'Slice2':
            y = self.cur_waveforms[self.data_index, self.slice_spin_2.value()]
            y = y / 100.0

        elif ylabel == 'Energy':
            y = np.sum(np.power(self.cur_waveforms[self.data_index, :], 2), axis=1)
            y = y / 1000000.0

        elif ylabel == 'Peak':
            y = self.cur_waveforms[self.data_index, :].max(axis=1)
            y = y / 100.0

        elif ylabel == 'Valley':
            y = self.cur_waveforms[self.data_index, :].min(axis=1)
            y = y / 100.0

        elif ylabel == 'Peak Pt':
            y = self.cur_waveforms[self.data_index, :].argmax(axis=1)

        elif ylabel == 'Valley Pt':
            y = self.cur_waveforms[self.data_index, :].argmin(axis=1)

        elif ylabel == 'Pk2Pk Amp':
            y = self.cur_waveforms[self.data_index, :].max(axis=1) - self.cur_waveforms[self.data_index, :].min(axis=1)
            y = y / 100.0

        elif ylabel == 'Time':
            y = self.cur_ts[self.data_index]
            y = y / 60000.0

        n_axes = len(self.chan_tab['features_fig'].figure.axes)

        #n_spikes = self.n_spikes_slider.value()

        title = '%s: %s vs %s' % (what_to_plot, xlabel, ylabel)
        # obtain the axis limits if we are plotting the same variables
        same_limits = False
        if n_axes > 0 and \
           self.chan_tab['features_fig'].figure.axes[0].get_title() == title:
            same_limits = True
            x_lim = self.chan_tab['features_fig'].figure.axes[0].get_x_lim()
            y_lim = self.chan_tab['features_fig'].figure.axes[0].get_y_lim()

        # plot only on one axes
        if self.plot_density_check.checkState() == 0:
            if n_axes == 0:
                ax1 = self.chan_tab['features_fig'].figure.add_subplot(111)
                ax1.set_facecolor('k')

            elif n_axes == 1:
                ax1 = self.chan_tab['features_fig'].figure.axes[0]
                ax1.cla()
                ax1.set_facecolor('k')

            elif n_axes >= 2:
                self.chan_tab['features_fig'].figure.clear()
                ax1 = self.chan_tab['features_fig'].figure.add_subplot(111)
                ax1.set_facecolor('k')

        # create 2 subplots to host the density
        elif self.plot_density_check.checkState() == 2:
            if n_axes == 0:
                ax1 = self.chan_tab['features_fig'].figure.add_subplot(121)
                ax2 = self.chan_tab['features_fig'].figure.add_subplot(122, sharex=ax1, sharey=ax1)

            elif n_axes == 1:
                self.chan_tab['features_fig'].figure.clear()
                ax1 = self.chan_tab['features_fig'].figure.add_subplot(121)
                ax2 = self.chan_tab['features_fig'].figure.add_subplot(122, sharex=ax1, sharey=ax1)

            elif n_axes == 2:
                ax1 = self.chan_tab['features_fig'].figure.axes[0]
                ax2 = self.chan_tab['features_fig'].figure.axes[1]
                ax1.cla()
                ax2.cla()

            ax2.set_facecolor('k')
            # create and plot a 2d histogram

        # setup the axes
        ax1.set_title(title, fontdict={'color': 'w'})
        ax1.tick_params(color=[.5, .5, .5])
        for k in ax1.spines.values():
            k.set_edgecolor([.5, .5, .5])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_facecolor('k')

        self.cursor, = ax1.plot([], 's', mfc='none', ms=6, mec='r',
                                animated=True, label='sample')

        # iterate over the members of that channel
        if what_to_plot == 'All waveforms':
            nodes = self.h5file.list_nodes(self.cur_node_name)
            for leaf in nodes:
                if leaf._v_name == 'Unsorted':
                    # select only some indices to plot
                    if leaf.nrows > self.n_pts_spin.value():
                        indx = leaf.read(0, leaf.nrows, sample_step(leaf.nrows, self.n_pts_spin.value()))
                    else:
                        indx = leaf.read()

                    # plot unsorted
                    ax1.plot(x[indx], y[indx], ',',
                             color=[.5, .5, .5], label='data_unsorted')

                unit = re.search('(?<=Unit)[0-9]{2}', leaf._v_name)
                
                if unit:
                    unit_index = int(unit.group())
                    unit_color = tuple(self.unit_colors[unit_index])
                    # select some units to plot
                    if leaf.Indx.nrows > self.n_pts_spin.value():
                        indx = leaf.Indx.read(0, leaf.Indx.nrows, sample_step(leaf.Indx.nrows, self.n_pts_spin.value()))
                    else:
                        indx = leaf.Indx.read()

                    ax1.plot(x[indx], y[indx], ',', label='data_' + leaf._v_name,
                             rasterized=True,
                             color=unit_color,
                             mec=unit_color)

                    # add unit to the tab widget
                    self.units_table_add_unit(leaf._v_name)

        elif what_to_plot == 'Sorted':
            nodes = self.h5file.list_nodes(self.cur_node_name)
            for leaf in nodes:
                unit = re.search('(?<=Unit)[0-9]{2}', leaf._v_name)
                
                if unit:
                    unit_index = int(unit.group())
                    unit_color = tuple(self.unit_colors[unit_index])
                    # select some units to plot
                    if leaf.Indx.nrows > self.n_pts_spin.value():
                        indx = leaf.Indx.read(0, leaf.Indx.nrows, sample_step(leaf.Indx.nrows, self.n_pts_spin.value()))
                    else:
                        indx = leaf.Indx.read()

                    ax1.plot(x[indx], y[indx], ',', label='data_' + leaf._v_name,
                             rasterized=True,
                             color=unit_color,
                             mec=unit_color,
                             zorder=10)

                    # add unit to the tab widget
                    self.units_table_add_unit(leaf._v_name)

        # to plot the unsorted channels
        elif what_to_plot == 'Unsorted':
            lx = len(x)
            # select some units to plot
            if lx > self.n_pts_spin.value():
                indx = sample_range(lx, self.n_pts_spin.value())
            else:
                indx = list(range(lx))
            ax1.plot(x[indx], y[indx], ',', color=[.5, .5, .5],
                     label='data_unsorted',
                     rasterized=True,
                     zorder=10)

        # plot a specific unit
        
        elif re.search('Unit', what_to_plot):
            unit = re.search('(?<=Unit)[0-9]{0,2}', what_to_plot).group()
            unit_index = int(unit)
            unit_color = tuple(self.unit_colors[unit_index])
            lx = len(x)
            # select some units to plot
            if lx > self.n_pts_spin.value():
                indx = sample_range(lx, self.n_pts_spin.value())
            else:
                indx = list(range(lx))

            ax1.plot(x[indx], y[indx], ',', label='data_' + what_to_plot,
                     rasterized=True,
                     color=unit_color,
                     mec=unit_color,
                     zorder=10)

            # add unit to the tab widget
            self.units_table_add_unit(what_to_plot)

        if same_limits:
            ax1.set_ylim(y_lim)
            ax1.set_xlim(x_lim)
        else:
            ax1.relim()
            ax1.autoscale_view(True, True, True)

        # vertical and horizontal lines @ x and y = 0
        ax1.axvline(0, color=[.5, .5, .5], zorder=0)
        ax1.axhline(0, color=[.5, .5, .5], zorder=0)

        # 60 minute line
        if xlabel == 'Time':
            ax1.axvline(60, color='gray', zorder=0)

        # create KDTree objet from the selected data for fast search
        self.xy_data = cKDTree(np.array([x, y]).T)

        # connect figure to the motion notify function
        if not ax1.callbacks.callbacks or not hasattr(self, 'ax_zoom_cid'):
            self.ax_zoom_cid = ax1.callbacks.connect('y_lim_changed', self.axes_zoom)

        # connect figure to the motion notify function
        if not hasattr(self, 'motion_cid'):
            self.motion_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('motion_notify_event',
                                                                                   self.nearest_point)

        # connect figure to the draw function
        if not hasattr(self, 'draw_cid'):
            self.draw_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('draw_event',
                                                                                 self.draw_callback)

        # plot density if checked
        if self.plot_density_check.checkState() == 2:
            self.replot_density()

        # set tight layout and redraw figure
        self.chan_tab['features_fig'].figure.tight_layout()
        self.chan_tab['features_fig'].figure.canvas.draw()

    #__________________________________________________________________________
    def plot_3d_features(self):

        # obtain labels and return if are the same
        xlabel = self.x_plot.currentText()
        ylabel = self.y_plot.currentText()
        zlabel = self.z_plot.currentText()
        if xlabel == ylabel or xlabel == zlabel or ylabel == zlabel:
            return

        curchan = int(self.chan_selector.currentText())

        if self.plot_valids_only_check.checkState() and self.cur_node.__contains__('ValidWFs'):
            print('you selected to plot only the valid Wfs')

        what_to_plot = str(self.what_to_plot.currentText())  # string containing what to plot
        self.cur_node_name = '/Spikes/Chan_%03d' % curchan

        if what_to_plot == 'All waveforms':
            self.data_index = list(range(self.cur_ts.size))
            pc = self.chan_tab['pca']
        elif what_to_plot == 'Unsorted':
            self.data_index = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
            pc = pca_scores(self.cur_waveforms[self.data_index, :])
        elif re.search('Unit', what_to_plot):
            self.cur_unit_name = what_to_plot
            self.data_index = self.h5file.get_node(self.cur_node_name, what_to_plot).Indx.read()
            pc = pca_scores(self.cur_waveforms[self.data_index, :])
        elif what_to_plot == 'Sorted':
            return

        # save what is the feature
        self.cur_feature_plot = what_to_plot

        # get the choice for the x axis
        if xlabel == 'PCA1':
            x = pc[:, 0]

        elif xlabel == 'PCA2':
            x = pc[:, 1]

        elif xlabel == 'PCA3':
            x = pc[:, 2]

        elif xlabel == 'Slice1':
            x = self.cur_waveforms[self.data_index, self.slice_spin_1.value()] / 100.0
            x = x / 100.0

        elif xlabel == 'Slice2':
            x = self.cur_waveforms[self.data_index, self.slice_spin_2.value()] / 100.0
            x = x / 100.0

        elif xlabel == 'Energy':
            x = np.sum(np.power(self.cur_waveforms[self.data_index, :], 2), axis=1)
            x = x / 1000000.0

        elif xlabel == 'Peak':
            x = self.cur_waveforms[self.data_index, :].max(axis=1)
            x = x / 100.0

        elif xlabel == 'Valley':
            x = self.cur_waveforms[self.data_index, :].min(axis=1)
            x = x / 100.0

        elif xlabel == 'Pk2Pk Amp':
            x = self.cur_waveforms[self.data_index, :].max(axis=1) - self.cur_waveforms[self.data_index, :].min(axis=1)
            x = x / 100.0

        elif xlabel == 'Time':
            x = self.cur_ts[self.data_index]
            x = x / 60000.0

        # get the choice for the y axis
        if ylabel == 'PCA1':
            y = pc[:, 0]

        elif ylabel == 'PCA2':
            y = pc[:, 1]

        elif ylabel == 'PCA3':
            y = pc[:, 2]

        elif ylabel == 'Slice1':
            y = self.cur_waveforms[self.data_index, self.slice_spin_1.value()] / 100.0
            y = y / 100.0

        elif ylabel == 'Slice2':
            y = self.cur_waveforms[self.data_index, self.slice_spin_2.value()] / 100.0
            y = y / 100.0

        elif ylabel == 'Energy':
            y = np.sum(np.power(self.cur_waveforms[self.data_index, :], 2), axis=1)
            y = y / 1000000.0

        elif ylabel == 'Peak':
            y = self.cur_waveforms[self.data_index, :].max(axis=1)
            y = y / 100.0

        elif ylabel == 'Valley':
            y = self.cur_waveforms[self.data_index, :].min(axis=1)
            y = y / 100.0

        elif ylabel == 'Pk2Pk Amp':
            y = self.cur_waveforms[self.data_index, :].max(axis=1) - self.cur_waveforms[self.data_index, :].min(axis=1)
            y = y / 100.0

        elif ylabel == 'Time':
            y = self.cur_ts[self.data_index]
            y = y / 60000.0

        # get the choice for the z axis
        if zlabel == 'PCA1':
            z = pc[:, 0]

        elif zlabel == 'PCA2':
            z = pc[:, 1]

        elif zlabel == 'PCA3':
            z = pc[:, 2]

        elif zlabel == 'Slice1':
            z = self.cur_waveforms[self.data_index, self.slice_spin_1.value()] / 100.0
            z = z / 100.0

        elif zlabel == 'Slice2':
            z = self.cur_waveforms[self.data_index, self.slice_spin_2.value()] / 100.0
            z = z / 100.0

        elif zlabel == 'Energy':
            z = np.sum(np.power(self.cur_waveforms[self.data_index, :], 2), axis=1)
            z = z / 1000000.0

        elif zlabel == 'Peak':
            z = self.cur_waveforms[self.data_index, :].max(axis=1)
            z = z / 100.0

        elif zlabel == 'Valley':
            z = self.cur_waveforms[self.data_index, :].min(axis=1)
            z = z / 100.0

        elif zlabel == 'Pk2Pk Amp':
            z = self.cur_waveforms[self.data_index, :].max(axis=1) - self.cur_waveforms[self.data_index, :].min(axis=1)
            z = z / 100.0

        elif zlabel == 'Time':
            z = self.cur_ts[self.data_index]
            z = z / 60000.0

        if what_to_plot == 'All waveforms' and self.cur_node.__contains__('ValidWFs'):
            valid = self.cur_node.ValidWFs.read()
            x = x[valid]
            y = y[valid]
            z = z[valid]

        items = self.widget_3d.items
        for i in items:
            self.widget_3d.removeItem(i)

        grid = gl.GLGridItem()
        self.widget_3d.addItem(grid)

        if zlabel != 'Density':
            handle = gl.GLScatterPlotItem(pos=np.array([x, y, z]).T,
                                          size=np.ones(x.size),
                                          color=(1.0, 0.0, 0.0, 1.0),
                                          pxmode=True)
            self.widget_3d.addItem(handle)

        else:
            #pass

            # obtain axes and first axes limits
            ax1 = self.chan_tab['features_fig'].figure.axes[0]
            #x_lim = ax1.get_x_lim()
            #y_lim = ax1.get_y_lim()

            # search for the unsorted or the units plots to obatin data
            x_points = []
            y_points = []
            for k in ax1.get_children():
                if re.search('unsorted|Unit', str(k.get_label())):
                    data = k.get_data()
                    x_points.extend(data[0])
                    y_points.extend(data[1])
            xy_points = np.array([x_points, y_points]).T

            # check wich points are inside the axes
            verts = ax1.viewLim.corners()
            verts[2] = ax1.viewLim.corners()[3]
            verts[3] = ax1.viewLim.corners()[2]

            in_poly = Path(verts).contains_points(xy_points)

            # create a 2d histogram of the data and scale it logaritmically
            h, _, _ = np.histogram2d(xy_points[in_poly, 0],
                                     xy_points[in_poly, 1],
                                     bins=self.plot_density_bins.value(),
                                     density=False)
            h[h <= 0] = 1
            h = np.log10(h)
            x, y = h.shape
            x, y = np.arange(x), np.arange(y)
            handle = gl.GLSurfacePlotItem(x, y, z=10 * h / h.max(),
                                          shader='heightColor')
            handle.translate(-x.size / 2.0, -y.size / 2.0, 0.0)
            handle.scale(1, 1, 2)
            self.widget_3d.addItem(handle)

    #__________________________________________________________________________
    def validate_wfs(self):
        ''' obtains the coordinates of the current feature axis, and uses it to
        determine wich points lay inside it.
        It also saves the data to the h5file'''

        # exits if there is no h5 file loaded or channel plotted
        if not self.h5_file_loaded or not self.chan_plotted:
            return

        # get axes handle and limits
        ax = self.chan_tab['features_fig'].figure.axes[0]
        x_lim = ax.get_x_lim()
        y_lim = ax.get_y_lim()

        # obtain coordinates of the current axes and uses that to build a polygon
        xyverts = [[x_lim[0], y_lim[0]], [x_lim[0], y_lim[1]], [x_lim[1], y_lim[1]], [x_lim[1], y_lim[0]]]

        # obtain the indices of the waveforms inside the polygon
        p = Path(xyverts).contains_points(self.xy_data.data)

        # in case no points were inside the axes
        if len(p) == 0:
            self.msg_box.setIcon(QtWidgets.QMessageBox.Warning)
            self.msg_box.setText('There were no selected points')
            self.msg_box.setWindowTitle('Warning')
            self.msg_box.show()
            return

        self.valid_wfs = np.flatnonzero(p)
        self.invalid_wfs = np.flatnonzero(~p)

        # remove the 'ValidWFs' field if it already exists
        if self.h5file.get_node(self.cur_node_name).__contains__('ValidWFs'):
            self.h5file.remove_node(self.cur_node_name, 'ValidWFs')

        # remove the 'InvalidWFs' field if it already exists
        if self.h5file.get_node(self.cur_node_name).__contains__('InvalidWFs'):
            self.h5file.remove_node(self.cur_node_name, 'InvalidWFs')

        # save the valid_wfs indices to the h5file
        self.h5file.create_array(self.cur_node_name, 'ValidWFs', self.valid_wfs)

        # save the invalid_wfs indices to the h5file
        self.h5file.create_array(self.cur_node_name, 'InvalidWFs', self.invalid_wfs)

        # save changes to disk
        self.h5file.flush()

        # update the information on the overview table
        row = self.chan_selector.currentIndex()
        item = QtWidgets.QTableWidgetItem(str(self.valid_wfs.size))
        self.overview_tab_2['overview_table'].takeItem(row, 5)
        self.overview_tab_2['overview_table'].setItem(row, 5, item)

    #__________________________________________________________________________
    def replot_density(self):
        ''' replot density using all the resolution only to the visible points'''

        # check whether the number of axes in the figure
        if len(self.chan_tab['features_fig'].figure.axes) != 2:
            return

        # obtain axes and first axes limits
        ax1 = self.chan_tab['features_fig'].figure.axes[0]
        ax2 = self.chan_tab['features_fig'].figure.axes[1]
        x_lim = ax1.get_x_lim()
        y_lim = ax1.get_y_lim()

        # search for the unsorted or the units plots to obatin data
        x_points = []
        y_points = []
        for k in ax1.get_children():
            if re.search('unsorted|Unit', str(k.get_label())):
                data = k.get_data()
                x_points.extend(data[0])
                y_points.extend(data[1])
        xy_points = np.array([x_points, y_points]).T

        # check wich points are inside the axes
        verts = ax1.viewLim.corners()
        verts[2] = ax1.viewLim.corners()[3]
        verts[3] = ax1.viewLim.corners()[2]

        in_poly = Path(verts).contains_points(xy_points)

        # create a 2d histogram of the data and scale it logaritmically
        h, xd, yd = np.histogram2d(xy_points[in_poly, 0], xy_points[in_poly, 1],
                                   bins=self.plot_density_bins.value(),
                                   density=False)
        h[h <= 0] = 1
        h = np.log10(h)

        # clean axes No2 and plot the 2d histogram
        ax2.cla()
        cmap = helper_widgets.colormaps[settings.density_cm]
        ax2.pcolormesh(xd, yd, h.transpose(), cmap=cmap)

        # set axis limits
        ax2.set_xlim(x_lim)
        ax2.set_ylim(y_lim)

        # remove tick labels
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])

        # create vertical and horizontal lines at 0
        ax2.axvline(0, color=[.5, .5, .5])
        ax2.axhline(0, color=[.5, .5, .5])

        # redraw the figure
        self.chan_tab['features_fig'].figure.canvas.draw()

    #__________________________________________________________________________
    def axes_zoom(self, ax):

        x_points = []
        y_points = []
        for k in ax.get_children():
            if re.search('unsorted|Unit', str(k.get_label())):
                data = k.get_data()
                x_points.extend(data[0])
                y_points.extend(data[1])
        #xy_points = np.array([x_points,y_points]).transpose()
        # check wich points are inside the axes
        verts = ax.viewLim.corners()
        verts[2] = ax.viewLim.corners()[3]
        verts[3] = ax.viewLim.corners()[2]

#        in_poly = points_inside_poly(xy_points, verts)
#        w = self.cur_waveforms[in_poly,:]
#        self.chan_tab['waves_figure'].figure.axes[0].set_ylim(w.min(), w.max())
#        self.chan_tab['waves_figure'].figure.canvas.draw()

        if len(self.chan_tab['features_fig'].figure.axes) == 2:
            self.replot_density()

    #__________________________________________________________________________
    def auto_clust(self):
        if not self.h5_file_loaded or not self.chan_plotted:
            return

        if self.xy_data.data.shape[1] > 2:
            data = self.xy_data.data[:, 0:2]
        else:
            data = self.xy_data.data

        cluster_indices = klusta_kwik_call(data, self.min_clust.value(), self.max_clust.value())

        fig = self.chan_tab['features_fig'].figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.set_facecolor('k')

        for k in range(len(cluster_indices)):
            ax.plot(data[cluster_indices[k], 0], data[cluster_indices[k], 1], '.',
                    label='clust %d' % k)

        ax.legend(fancybox=True, mode='expand', ncol=max(1, int(len(cluster_indices) / 2)),
                  loc=9, prop={'size': 10})
        ax.grid(color='grey')
        fig.canvas.draw()
        self.sender().parentWidget().close()

    #__________________________________________________________________________
    def trim_waveforms(self, eclick, erelease):

        # first check whether there's any waveform plotted
        # if it is visible, and if it is the current unit
        for k in self.chan_tab['waves_figure'].figure.axes[0].lines:
            if 'Unit' in k.get_label() and k.get_visible() and \
               k.get_label() == self.cur_unit_name:
                break
        else:
            print("No units found in the plot ...")
            self.trim_waveforms_rect.set_active(False)
            return

        # get the indices
        indx = self.h5file.get_node('/Spikes/Chan_%03d/%s' % (self.cur_chan, self.cur_unit_name), 'Indx').read()
        data = self.cur_waveforms[indx, :]

        # get line equation y = mx + n
        x1 = eclick.xdata
        x2 = erelease.xdata
        y1 = eclick.ydata
        y2 = erelease.ydata

        # return if is a point and not a line
        if x1 == x2:
            self.trim_waveforms_rect.set_active(False)
            return

        m = (y2 - y1) / (x2 - x1)
        n = y1 - m * x1

        # get the y value of nearest integer x:
        x = np.array([x1, x2])
        x.sort()
        x_data = list(range(self.wf_size))
        indx1 = np.flatnonzero(x_data > x[0]).min()
        indx2 = np.flatnonzero(x_data < x[1]).max()
        y = np.array([m * x_data[k] + n for k in range(indx1, indx2)])
        #print x, y

        # get the data bounded by the indices
        data2 = data[:, indx1:indx2]
        #print data2.shape, y.shape
        t = data2 - y
        #print t
        t = np.array(t)

        # get the indices that intersect the line
        intersect = []
        for j, k in enumerate(t):
            if not (np.all(k < 0) or np.all(k > 0)):
                intersect.append(j)

        # update the node containing the unit indices
        self.h5file.remove_node(self.cur_node_name + '/' + self.cur_unit_name, 'Indx')
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'Indx', np.delete(indx, intersect))

        # add the remaining points to the unsorted indexes
        self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
        self.unsorted = np.append(self.unsorted, indx[intersect])
        self.unsorted.sort()

        # update the unsorted in the h5file
        self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # save changes to disk
        self.h5file.flush()

        # update the information in the overview table
        row = self.chan_selector.currentIndex()
        self.overview_tab_2['overview_table'].takeItem(row, self.cur_unit + 6)
        lbl = QtWidgets.QTableWidgetItem(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))
        self.overview_tab_2['overview_table'].setItem(row, self.cur_unit + 6, lbl)

        # update the information on the unit label
        self.chan_tab['unit_count_label'][self.cur_unit_name].setText(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))

        # replot the features
        self.plot_features()

        # replot waveforms
        self.plot_unit_waveforms()

        # replot the unit avg waveform, histogram and autocorrelation
        self.plot_unit_figure()

        eclick.in_axes.figure.canvas.draw()

        self.trim_waveforms_rect.set_active(False)

    #__________________________________________________________________________
    def activate_trim_waveforms(self):
        self.trim_waveforms_rect.set_active(True)

    #__________________________________________________________________________
    def plot_unit_waveforms(self):

        # get unit name and number
        unit_name = str(self.chan_tab['unit_tabs_widget'].tabText(self.chan_tab['unit_tabs_widget'].currentIndex()))
        unit_no = int(re.search('(?<=Unit)[0-9]{2}', unit_name).group())

        # get axes handle and children labels
        fig = self.chan_tab['waves_figure'].figure
        ax = fig.axes[0]
        children_labels = [str(k.get_label()) for k in ax.get_children()]

        # get the number of spikes to plot
        n_spikes = self.n_spikes_spin.value()

        node = self.cur_node.__getattr__(self.cur_unit_name)
        nrows = node.Indx.nrows

        if nrows > n_spikes:
            unit_indices = node.Indx.read(start=0, stop=nrows, step=sample_step(nrows, n_spikes))
        else:
            unit_indices = node.Indx.read()

        # obtain the length of units to plot
        n = len(unit_indices)

        # create an array of Nones to append
        nones = np.array(n * [None], ndmin=2).T

        # create the x indexes
        ts = np.tile(np.arange(self.wf_size), (n, 1))
        ts = np.append(ts, nones, axis=1).reshape((n * (self.wf_size + 1),))

        # get the waveforms, append nones, and reshape it to a vector
        wf = self.cur_node.Waveforms[unit_indices, :]
        wf = np.append(wf, nones, axis=1).reshape((n * (self.wf_size + 1),))

        # create the plot if it doesn't exists
        if unit_name not in children_labels:
            ax.plot(ts, wf, color=self.unit_colors[unit_no, :], alpha=0.7,
                    label=unit_name)

        # if exists update the data
        elif unit_name in children_labels:
            for k in self.chan_tab['waves_figure'].figure.axes[0].get_lines():
                if k.get_label() == self.cur_unit_name:
                    break

            k.set_data(ts, wf)
            k.set_visible(True)

        fig.canvas.draw()

    #__________________________________________________________________________
    def add_unit(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.chan_plotted:
            return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.chan_tab['features_fig'].figure.axes[0].get_title())
        if not re.search('waveforms|unsorted', title):
            return

        # return if a tool is selected in the toolbar
        if self.chan_tab['features_fig_toolbar'].mode != '':
            return

        # create a new lasso instance
        self.lasso_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.lasso_add_unit)

    #__________________________________________________________________________
    def keep(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.chan_plotted:
            return

        # check if a unit is plotted
        title = str(self.chan_tab['features_fig'].figure.axes[0].get_title())
        if not re.search('Unit', title):
            return

        self.what_to_plot.count()

        # return if a tool is selected in the toolbar
        if self.chan_tab['features_fig_toolbar'].mode != '':
            return

        # create a new lasso instance
        self.lasso_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.lasso_keep)

    #__________________________________________________________________________
    def add_region(self):
        ''' starts a lasso instance to draw a line around a ROI'''

        # check whether there is a channel ploted
        if not self.chan_plotted:
            return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.chan_tab['features_fig'].figure.axes[0].get_title())
        if not re.search('waveforms|unsorted', title):
            return

        # return if a tool is selected in the toolbar
        if self.chan_tab['features_fig_toolbar'].mode != '':
            return

        # create a new lasso instance
        self.lasso_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.lasso_add_region)

    #__________________________________________________________________________
    def remove_region(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.chan_plotted:
            return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.chan_tab['features_fig'].figure.axes[0].get_title())
        if not re.search('Unit', title):
            return

        # return if a tool is selected in the toolbar
        if self.chan_tab['features_fig_toolbar'].mode != '':
            return

        # create a new lasso instance
        self.lasso_cid = self.chan_tab['features_fig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.lasso_remove_region)

    #__________________________________________________________________________
    def lasso_add_unit(self, event):
        if self.chan_tab['features_fig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        if event.in_axes is None or event.button != 1:
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        # create a lasso instance
        self.lasso = matplotlib_widgets.MyLasso(event.in_axes, (event.xdata, event.ydata),
                                                self.lasso_callback_add_unit,
                                                color='gray', lw=1)
        self.chan_tab['features_fig'].figure.canvas.widgetlock(self.lasso)

    #__________________________________________________________________________
    def lasso_keep(self, event):
        if self.chan_tab['features_fig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        if event.in_axes is None or event.button != 1:
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        self.keep_btn.setCheckable(True)
        self.keep_btn.setChecked(True)
        self.lasso = matplotlib_widgets.MyLasso(event.in_axes, (event.xdata, event.ydata),
                                                self.lasso_callback_keep,
                                                color='gray', lw=1)
        self.chan_tab['features_fig'].figure.canvas.widgetlock(self.lasso)

    #__________________________________________________________________________
    def lasso_add_region(self, event):
        if self.chan_tab['features_fig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        if event.in_axes is None or event.button != 1:
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        self.lasso = matplotlib_widgets.MyLasso(event.in_axes, (event.xdata, event.ydata),
                                                self.lasso_callback_add_region,
                                                color='gray', lw=1)
        self.chan_tab['features_fig'].figure.canvas.widgetlock(self.lasso)

    #__________________________________________________________________________
    def lasso_remove_region(self, event):
        if self.chan_tab['features_fig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        if event.in_axes is None or event.button != 1:
            if hasattr(self, 'lasso_cid'):
                self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
                del self.lasso_cid
            return

        self.lasso = matplotlib_widgets.MyLasso(event.in_axes, (event.xdata, event.ydata),
                                                self.lasso_callback_remove_region,
                                                color='gray', lw=1)
        self.chan_tab['features_fig'].figure.canvas.widgetlock(self.lasso)

    #__________________________________________________________________________
    def lasso_callback_add_unit(self, verts):

        # disconnect Lasso callback
        self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
        self.chan_tab['features_fig'].figure.canvas.draw_idle()
        del self.lasso_cid

        # release widget lock
        self.chan_tab['features_fig'].figure.canvas.widgetlock.release(self.lasso)

        # delete lasso
        del self.lasso

        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n, int(n / 25))]

        #pdb.set_trace()

        # get the axes handle
        ax = self.chan_tab['features_fig'].figure.axes[0]

        if re.search('Waveforms', ax.get_title()):
            # test which points are inside the lasso
            xy_points = self.xy_data.data[self.unsorted, :]
        elif re.search('Unsorted', ax.get_title()):
            xy_points = self.xy_data.data

        p = Path(self.verts).contains_points(xy_points)

        # in case there were no points selected
        if len(p) == 0:
            self.msg_box.setIcon(QtWidgets.QMessageBox.Warning)
            self.msg_box.setText('There were no selected points')
            self.msg_box.setWindowTitle('Warning')
            self.msg_box.show()
            return

        # set the unit name
        self.n_units = len(self.units_list)
        self.cur_unit_name = 'Unit%02d' % self.n_units

        # look for the unsorted plot handle in the axes
        for k in self.chan_tab['features_fig'].figure.axes[0].get_children():
            if re.search('Unsorted', str(k.get_label())):
                break

        # obtain the unsorted points
        unsorted_data = xy_points[~p, :]
        len_unsorted = len(unsorted_data)

        # select some indices to plot
        if len_unsorted > self.n_pts_spin.value():
            indx = sample_range(len_unsorted, self.n_pts_spin.value())
        else:
            indx = list(range(len_unsorted))

        # replot the unsorted without the corresponding points to the new unit
        k.set_data(unsorted_data[:, 0][indx], unsorted_data[:, 1][indx])
        ax.draw_artist(k)

        # select some indices to plot
        unit_data = xy_points[p, :]
        len_unit = len(unit_data)

        if len_unit > self.n_pts_spin.value():
            indx = sample_range(len_unit, self.n_pts_spin.value())
        else:
            indx = list(range(len_unit))

        ax.plot(unit_data[:, 0][indx], unit_data[:, 1][indx],
                linestyle = '',
                marker = ',',
                mfc = list(self.unit_colors[self.n_units]),
                mec = list(self.unit_colors[self.n_units]),
                label = 'data_' + self.cur_unit_name)

        self.n_units += 1

        # if unit name not in combo box add it
        if self.cur_unit_name not in [str(self.what_to_plot.itemText(k)) for k in range(self.what_to_plot.count())]:
            self.what_to_plot.addItem(self.cur_unit_name)

        # add the indexes of the current unit to the h5file
        if self.h5file.get_node(self.cur_node_name).__contains__(self.cur_unit_name):
            self.h5file.remove_node(self.cur_node_name, self.cur_unit_name, recursive=True)
        self.h5file.create_group(self.cur_node_name, self.cur_unit_name)
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'Indx', self.unsorted[p])
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'isMultiunit', False)
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'isBursting', False)

        # update the list of unsorted indexes
        self.unsorted = self.unsorted[~p]

        # update  the indexes of the unsorted units
        if self.h5file.get_node(self.cur_node_name).__contains__('Unsorted'):
            self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # save changes to disk
        self.h5file.flush()

        # add log
        #self.add_log('%s %s added' % (self.cur_node_name, self.cur_unit_name))

        # add unit to the units tab widget
        self.units_table_add_unit(self.cur_unit_name)

        self.chan_tab['unit_figures'][self.cur_unit_name].figure.tight_layout()
        self.chan_tab['unit_figures'][self.cur_unit_name].figure.canvas.draw()

        # update the overview figure
        for k in self.overview_tab_1['figure'].figure.axes:
            if str(self.cur_chan) in k.get_title():
                break

        self.plot_chan_overview(self.cur_node, axes_to_plot=k)
        for l in k.lines:
            k.draw_artist(l)

    #__________________________________________________________________________
    def lasso_callback_keep(self, verts):
        # disconnect Lasso callback from figure
        self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
        self.chan_tab['features_fig'].figure.canvas.draw_idle()
        del self.lasso_cid

        # release the lock from the lasso
        self.chan_tab['features_fig'].figure.canvas.widgetlock.release(self.lasso)

        # erase lasso
        del self.lasso

        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n, int(n / 25))]

        # test which points lay inside the polygon
        p = Path(self.verts).contains_points(self.xy_data.data)

        # change to not checked
        self.keep_btn.setChecked(False)
        self.keep_btn.setCheckable(False)

        # check how many points were selected
        if len(np.flatnonzero(p)) <= self.wf_size:
            print("Didn't doo anything: Too few points selected")
            return

        # create a KDTree object for efficient neighbor search
        self.xy_data = cKDTree(self.xy_data.data[p, :])

        # get unitname and number from the axes title
        ax = self.chan_tab['features_fig'].figure.axes[0]
        self.cur_unit_name = re.search('Unit[0-9]{2}', ax.get_title()).group()
        self.cur_unit = int(re.search('(?<=Unit)[0-9]{2}', ax.get_title()).group())

        # update plot:
        for k in ax.get_children():
            if re.search(str(k.get_label), self.cur_unit_name):
                k.set_data(self.xy_data.data[:, 0], self.xy_data.data[:, 1])
                ax.draw_artist(k)
                break

        # return if no points selected
        if len(p) < 1:
            return

        node_name = self.cur_node_name + '/' + self.cur_unit_name

        # obtain the unit data
        unit_points = self.h5file.get_node(node_name, 'Indx').read()

        # update the node containing the unit indices
        self.h5file.remove_node(node_name, 'Indx')
        self.h5file.create_array(node_name, 'Indx', unit_points[p])

        # add the remaining points to the unsorted indexes
        self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
        self.unsorted = np.append(self.unsorted, unit_points[~p])
        self.unsorted.sort()

        # update the unsorted in the h5file
        self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # save changes to disk
        self.h5file.flush()

        # replot the unit avg waveform, histogram and autocorrelation
        self.plot_unit_figure()

        # replot the waveforms
        self.plot_unit_waveforms()

        # replot the features
        self.plot_features()

        # update the information in the overview table
        row = self.chan_selector.currentIndex()
        self.overview_tab_2['overview_table'].takeItem(row, self.cur_unit + 6)
        lbl = QtWidgets.QTableWidgetItem(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))
        self.overview_tab_2['overview_table'].setItem(row, self.cur_unit + 6, lbl)

        # update the information on the unit label
        self.chan_tab['unit_count_label'][self.cur_unit_name].setText(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))

    #__________________________________________________________________________
    def lasso_callback_add_region(self, verts):

        # disconnect the lasso from the canvas and redraw the figure
        self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
        self.chan_tab['features_fig'].figure.canvas.draw_idle()
        del self.lasso_cid

        # release widget lock
        self.chan_tab['features_fig'].figure.canvas.widgetlock.release(self.lasso)

        # delete lasso handle
        del self.lasso

        # get the vertices of the polygon to the object and downsample them if too many
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[sample_range(n, 25)]

        # check whether there is any unit
        if not hasattr(self, 'cur_unit_name') or not self.cur_node.__contains__(self.cur_unit_name):
            return

        # get the axes handle
        ax = self.chan_tab['features_fig'].figure.axes[0]

        # get the unsorted
        self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()

        # check what is plotted on the axes
        if re.search('Waveforms', str(ax.get_title())):
            # test which points are inside the lasso
            p = Path(self.verts).contains_points(self.xy_data.data[self.unsorted, :])
            self.xy_data = cKDTree(self.xy_data.data[self.unsorted, :][p])

        elif re.search('Unsorted', str(ax.get_title())):
            # test which points are inside the lasso
            p = Path(self.verts).contains_points(self.xy_data.data)
            self.xy_data = cKDTree(self.xy_data.data[p, :])

        # update plot:
        for k in ax.get_children():
            if re.search(str(k.get_label), self.cur_unit_name):
                k.set_data(self.xy_data.data[:, 0], self.xy_data.data[:, 1])
                ax.draw_artist(k)
                break

        # if more than 0 selected points
        if len(p) > self.wf_size:
            indx = self.unsorted[p]
        else:
            print("Didn't add any unit: Too few points selected")
            return

        # update the unsorted
        self.unsorted = self.unsorted[~p]
        self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # update the plots
        for k in ax.get_children():
            if re.search('Unsorted', str(k.get_label())):
                pass
            elif re.search(self.cur_unit_name, str(k.get_label())):
                pass

        # update the unit information in the file
        unit = self.h5file.get_node(self.cur_node_name + '/' + self.cur_unit_name, 'Indx').read()
        self.h5file.remove_node(self.cur_node_name + '/' + self.cur_unit_name, 'Indx')
        # append the new indexes to the waveform and sort
        unit = np.append(unit, indx)
        unit.sort()
        # create a new array in the h5file to hold the updated unit information
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'Indx', unit)

        # save changes to disk
        self.h5file.flush()

        # update the information in the overview table
        row = self.chan_selector.currentIndex()
        self.overview_tab_2['overview_table'].takeItem(row, self.cur_unit + 6)
        lbl = QtWidgets.QTableWidgetItem(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))
        self.overview_tab_2['overview_table'].setItem(row, self.cur_unit + 6, lbl)

        # update the information on the unit label
        self.chan_tab['unit_count_label'][self.cur_unit_name].setText(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))

        # replot the unit avg waveform, histogram and autocorrelation
        self.plot_unit_figure()

        # replot the features
        self.plot_features()

    #__________________________________________________________________________
    def lasso_callback_remove_region(self, verts):

        # disconnect Lasso callback from figure
        self.chan_tab['features_fig'].figure.canvas.mpl_disconnect(self.lasso_cid)
        self.chan_tab['features_fig'].figure.canvas.draw_idle()
        del self.lasso_cid

        # release the lock from the lasso
        self.chan_tab['features_fig'].figure.canvas.widgetlock.release(self.lasso)

        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[sample_range(n, 25)]

        # test which points lay inside the polygon
        p = Path(self.verts).contains_points(self.xy_data.data)

        # return if no points selected
        if len(p) < 1:
            return

        # get unitname and number from the axes title
        ax = self.chan_tab['features_fig'].figure.axes[0]
        self.cur_unit_name = re.search('Unit[0-9]{2}', ax.get_title()).group()
        self.cur_unit = int(re.search('(?<=Unit)[0-9]{2}', ax.get_title()).group())

        # obtain the unit data
        unit_points = self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.read()

        # update the node containing the unit indexes
        self.h5file.remove_node(self.cur_node_name + '/' + self.cur_unit_name, 'Indx')
        self.h5file.create_array(self.cur_node_name + '/' + self.cur_unit_name, 'Indx', unit_points[~p])

        # add the remaining points to the unsorted indexes
        self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
        self.unsorted = np.append(self.unsorted, unit_points[p])
        self.unsorted.sort()

        # update the unsorted in the h5file
        self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # save changes to disk
        self.h5file.flush()

        # update the information in the overview table
        row = self.chan_selector.currentIndex()
        self.overview_tab_2['overview_table'].takeItem(row, self.cur_unit + 6)
        lbl = QtWidgets.QTableWidgetItem(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))
        self.overview_tab_2['overview_table'].setItem(row, self.cur_unit + 6, lbl)

        # update the information on the unit label
        self.chan_tab['unit_count_label'][self.cur_unit_name].setText(str(self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.nrows))

        # replot the features
        self.plot_features()

        # replot the waveforms
        self.plot_unit_waveforms()

        # replot the unit avg waveform, histogram and autocorrelation
        self.plot_unit_figure()

        #plot_chan_overview(self.cur_node, axes_to_plot)

        # erase lasso
        del self.lasso

    #__________________________________________________________________________
    def units_table_add_unit(self, unit_name):
        ''' creates a new tab per each new unit'''

        # check whether that tab already exists
        for k in list(range(self.chan_tab['unit_tabs_widget'].count())):
            if unit_name == self.chan_tab['unit_tabs_widget'].tabText(k):
                return

        self.cur_unit_name = unit_name

        # create a widget and a layout
        widget = QtWidgets.QWidget()
        vlay = QtWidgets.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setContentsMargins(0, 0, 0, 0)

        # get a unit number
        unit_no = int(re.search('(?<=Unit)[0-9]{2}', unit_name).group())

        # add the unit number to a list
        self.units_list.append(unit_no)

        # create a btn to change unit color
        hlay = QtWidgets.QHBoxLayout()
        hlay.setContentsMargins(0, 0, 0, 0)
        hlay.addStretch(1)
        self.chan_tab['unit_buttons'][unit_name] = QtWidgets.QPushButton('Unit %02d' % unit_no)
        self.chan_tab['unit_buttons'][unit_name].setMaximumHeight(20)
        self.chan_tab['unit_buttons'][unit_name].clicked.connect(self.change_unit_color)
        self.chan_tab['unit_buttons'][unit_name].setStyleSheet('QPushButton {background: rgb%s}' % str(tuple(np.int16(255 * self.unit_colors[unit_no]))))
        hlay.addWidget(self.chan_tab['unit_buttons'][unit_name])
        hlay.addStretch(1)

        # plot-raw check button
        self.chan_tab['plot_raw_check'][unit_name] = QtWidgets.QCheckBox()
        self.chan_tab['plot_raw_check'][unit_name].setObjectName(str(unit_no))
        self.chan_tab['plot_raw_check'][unit_name].setChecked(False)
        self.chan_tab['plot_raw_check'][unit_name].setMaximumHeight(20)
        self.chan_tab['plot_raw_check'][unit_name].stateChanged.connect(self.set_waveform_visible)
        lbl = QtWidgets.QLabel('Plot Raw ?')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.chan_tab['plot_raw_check'][unit_name])
        hlay.addStretch(1)

        # is Multiunit check button
        self.chan_tab['is_multiunit_check'][unit_name] = QtWidgets.QCheckBox()
        self.chan_tab['is_multiunit_check'][unit_name].setObjectName(str(unit_no))
        self.chan_tab['is_multiunit_check'][unit_name].setChecked(False)
        self.chan_tab['is_multiunit_check'][unit_name].setMaximumHeight(20)
        self.chan_tab['is_multiunit_check'][unit_name].stateChanged.connect(self.set_is_multiunit)
        lbl = QtWidgets.QLabel('isMultiunit ?')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.chan_tab['is_multiunit_check'][unit_name])
        hlay.addStretch(1)

        # set the checkstate of the 'isMultiunit' check according to what is saved in the h5file
        if self.h5file.get_node('/Spikes/Chan_%03d/Unit%02d' % (self.cur_chan, unit_no)).__contains__('isMultiunit'):
            is_multiunit_value = self.h5file.get_node('/Spikes/Chan_%03d/Unit%02d' % (self.cur_chan, unit_no),
                                              'isMultiunit').read()
            if is_multiunit_value:
                self.chan_tab['is_multiunit_check'][unit_name].setChecked(True)
            else:
                self.chan_tab['is_multiunit_check'][unit_name].setChecked(False)
        else:
            self.h5file.create_array('/Spikes/Chan_%03d/Unit%02d' % (self.cur_chan, unit_no), 'isMultiunit', False)

        # add a label with the waveform count
        lbl = QtWidgets.QLabel('Count')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        self.chan_tab['unit_count_label'][unit_name] = QtWidgets.QLabel('%d' % self.h5file.get_node(self.cur_node_name, unit_name).Indx.nrows)
        self.chan_tab['unit_count_label'][unit_name].setMaximumHeight(20)
        hlay.addWidget(self.chan_tab['unit_count_label'][unit_name])
        hlay.addStretch(1)

        # add delete-unit button
        self.chan_tab['del_unit_buttons'][unit_name] = QtWidgets.QPushButton('Del Unit')
        self.chan_tab['del_unit_buttons'][unit_name].setObjectName(unit_name)
        self.chan_tab['del_unit_buttons'][unit_name].setMaximumHeight(20)
        self.chan_tab['del_unit_buttons'][unit_name].clicked.connect(self.del_unit)
        hlay.addWidget(self.chan_tab['del_unit_buttons'][unit_name])
        hlay.addStretch(1)

        vlay.addLayout(hlay)

        # add the figure widget
        self.chan_tab['unit_figures'][unit_name] = matplotlib_widgets.MplWidget()
        self.chan_tab['unit_figures'][unit_name].setObjectName(unit_name)  # set the name of the object
        self.chan_tab['unit_figures'][unit_name].figure.set_facecolor('k')
        n = matplotlib_widgets.NavToolbar(self.chan_tab['unit_figures'][unit_name].figure.canvas,
                                          widget, coordinates=False)
        n.setIconSize(QtCore.QSize(12, 12))
        n.setOrientation(QtCore.Qt.Vertical)

        vlay.addWidget(self.chan_tab['unit_figures'][unit_name])
        #vlay.addWidget(n)

        hlay = QtWidgets.QHBoxLayout()
        hlay.setSpacing(0)
        hlay.setContentsMargins(2, 2, 2, 2)
        hlay.addLayout(vlay)
        hlay.addWidget(n)

        widget.setLayout(hlay)

        # Plot the data
        self.plot_unit_figure()

        #if unit_name == 'Unit00':
        #    self.chan_tab['unit_tabs_widget'].removeTab(0)

        self.chan_tab['unit_tabs_widget'].addTab(widget, unit_name)
        indx = self.chan_tab['unit_tabs_widget'].count() - 1

        color = QtGui.QColor(*np.int32(self.unit_colors[indx] * 255))
        self.chan_tab['unit_tabs_widget'].tabBar().setTabTextColor(indx, color)

        # update the information in the overview table
        row = self.chan_selector.currentIndex()

        if self.overview_tab_2['overview_table'].columnCount() <= (unit_no + 6):
            self.overview_tab_2['overview_table'].insertColumn(self.overview_tab_2['overview_table'].columnCount())
            n_cols = self.overview_tab_2['overview_table'].columnCount()
            self.overview_tab_2['overview_table'].setColumnWidth(n_cols - 1, 65)
            self.overview_tab_2['overview_table'].setHorizontalHeaderItem(n_cols - 1,
                                                                       QtWidgets.QTableWidgetItem('Unit%02d' % unit_no))

        self.overview_tab_2['overview_table'].takeItem(row, unit_no + 6)
        lbl = QtWidgets.QTableWidgetItem(str(self.h5file.get_node(self.cur_node_name, unit_name).Indx.nrows))
        self.overview_tab_2['overview_table'].setItem(row, unit_no + 6, lbl)

        # update the unsorted number in the overview table
        self.overview_tab_2['overview_table'].takeItem(row, 4)
        lbl = QtWidgets.QTableWidgetItem(str(len(self.unsorted)))
        self.overview_tab_2['overview_table'].setItem(row, 4, lbl)

    #__________________________________________________________________________
    def plot_unit_figure(self):

        # get a unit name and number
        unit_no = int(re.search('(?<=Unit)[0-9]{2}', self.cur_unit_name).group())

        # find the figure that has a particular name
        fig = self.chan_tab['unit_figures'][self.cur_unit_name].figure

        # check whether we have to create axes
        if len(fig.axes) == 0:
            ax0 = fig.add_subplot(131)
            ax1 = fig.add_subplot(132)
            ax2 = fig.add_subplot(133)
        else:
            ax0 = fig.axes[0]
            ax0.cla()
            ax1 = fig.axes[1]
            ax1.cla()
            ax2 = fig.axes[2]
            ax2.cla()

        # set the axis background color
        ax0.set_facecolor('k')
        ax1.set_facecolor('k')
        ax2.set_facecolor('k')

        # PLOT AVERAGE WAVEFORM #####
        x = list(range(self.wf_size))
        p = self.h5file.get_node(self.cur_node_name, self.cur_unit_name).Indx.read()
        m = self.cur_waveforms[p, :].mean(axis=0)
        s = self.cur_waveforms[p, :].std(axis=0)
        mn = self.cur_waveforms[p, :].min(axis=0)
        mx = self.cur_waveforms[p, :].max(axis=0)

        # plot average waveform
        ax0.plot(x, m, color=self.unit_colors[unit_no], lw=2, label=self.cur_unit_name)

        #plot shaded area of 3 standard devoations around it
        ax0.fill_between(x, m + 3 * s, m - 3 * s, color=self.unit_colors[unit_no],
                         alpha=0.5, label=self.cur_unit_name)

        #plot maximum and minimum boundaries
        ax0.fill_between(x, mx, mn, color=self.unit_colors[unit_no], alpha=0.35,
                         label=self.cur_unit_name)
        ax0.set_xlim(0, self.wf_size - 1)
        ax0.set_yticklabels([])
        ax0.grid(color=[.5, .5, .5])
        ax0.tick_params(color=[.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax0.spines.values():
            k.set_edgecolor([.5, .5, .5])

        # PLOT ISI HISTOGRAM #####
        delta_ts = np.diff(self.cur_ts[p])
        delta_ts = delta_ts[delta_ts < 100]

        ld = len(delta_ts)
        if ld > 1000:
            indx = list(range(0, 1000))
        else:
            indx = list(range(ld))

        if len(delta_ts[indx]) > 0:
            ax1.hist(delta_ts[indx], bins=100, range=[0, 100], ec='none',
                     color=self.unit_colors[unit_no], label=self.cur_unit_name)
            ax1.tick_params(color=[.5, .5, .5], labelcolor=[.5, .5, .5])

        for k in ax1.spines.values():
            k.set_edgecolor([.5, .5, .5])

        wf_width = self.wf_size * 1000 / self.sf
        try:
            collision = 100 * np.flatnonzero(delta_ts < 1.5 * wf_width).size / float(delta_ts.size)
            # put a "percentage of collision" label
            ax1.text(0.5, 0.01, u'Collision = %0.2f %%' % collision,
                     transform=ax1.transAxes, color='w', size=10, ha='center')
        except:
            pass

        ax1.set_xlim(0, 100)

        # PLOT AUTOCORRELATION #####

        ts = self.cur_ts[p]

        time = 25000

        ts = ts[np.flatnonzero(ts < time)]

        ts11 = np.tile(ts, (ts.size, 1))
        ts22 = np.tile(ts, (ts.size, 1)).transpose()
        x = ts11 - ts22
        ac, lags = np.histogram(x.flatten(), bins=100, range=(-500, 500))
        ac[np.flatnonzero(lags == 0)] = 0.0
        ax2.bar(lags[0:-1], ac, width=np.diff(lags)[0], edgecolor='none',
                color=self.unit_colors[unit_no])

        '''if ts.size > 1000: ts = ts[0:1000]

        ac, x = autocorr(ts, bin_size = 20, win = [0,10000],
                         mode = 'fft', time_range = [-150, 150])
        ac[ac.argmax()] = 0
        ax2.plot(x, ac, color = self.unit_colors[unit_no], lw = 2)'''

        ax2.set_xlim(-500, 500)
        #ax2.set_ylim(0, ac.max())
        ax2.tick_params(color=[.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax2.spines.values():
            k.set_edgecolor([.5, .5, .5])
        ax2.set_yticklabels([])

        self.chan_tab['unit_figures'][self.cur_unit_name].figure.tight_layout()
        self.chan_tab['unit_figures'][self.cur_unit_name].figure.canvas.draw()

    #__________________________________________________________________________
    def set_waveform_visible(self):
        ''' makes the raw waveform of each unit visible or invisible'''

        sender = self.sender()
        state = sender.checkState()
        name = int(sender.objectName())

        # get unit name and number
        unit_name = str(self.chan_tab['unit_tabs_widget'].tabText(self.chan_tab['unit_tabs_widget'].currentIndex()))
        unit_no = int(re.search('(?<=Unit)[0-9]{2}', unit_name).group())

        # get axes handle and children labels
        ax = self.chan_tab['waves_figure'].figure.axes[0]
        children_labels = [str(k.get_label()) for k in ax.get_children()]

        # get the node to read from
        node = self.h5file.get_node(self.cur_node_name + '/' + unit_name, 'Indx')

        # get the number of spikes to plot
        n_spikes = self.n_spikes_spin.value()

        if state == 2:  # if checked
            nrows = node.nrows
            if nrows > n_spikes:
                unit_indices = node.read(start=0, stop=nrows, step=sample_step(nrows, n_spikes))
            else:
                unit_indices = node.read()

            # obtain the length of units to plot
            n = len(unit_indices)

            # create an array of Nones to append
            nones = np.array(n * [None], ndmin=2).T

            # create the x indexes
            ts = np.tile(np.arange(self.wf_size), (n, 1))
            ts = np.append(ts, nones, axis=1).reshape((n * (self.wf_size + 1),))

            # get the waveforms, append nones, and reshape it to a vector
            wf = self.cur_node.Waveforms[unit_indices, :]
            wf = np.append(wf, nones, axis=1).reshape((n * (self.wf_size + 1),))

            # create the plot if it doesn't exists
            if unit_name not in children_labels:
                ax.plot(ts, wf, color=self.unit_colors[unit_no, :],
                        alpha=0.7, label=unit_name)

            # if exists update the data
            elif unit_name in children_labels:
                for k in self.chan_tab['waves_figure'].figure.axes[0].get_children():
                    if k.get_label() == 'Unit%02d' % name:
                        break

                k.set_data(ts, wf)
                k.set_visible(True)

        elif state == 0:  # if unchecked
            for k in ax.get_children():
                if str(k.get_label()) == unit_name:
                    k.set_visible(False)

        # set axes limit
        lim = self.wave_ax_y_lim_spin.value()
        ax.set_ylim(-lim, lim)

        # finally redraw the figure
        self.chan_tab['waves_figure'].figure.canvas.draw()

    #__________________________________________________________________________
    def set_is_multiunit(self):

        sender = self.sender()
        #state  = sender.checkState()
        unit_no = int(sender.objectName())

        # current node name
        node_name = '/Spikes/Chan_%03d/Unit%02d' % (self.cur_chan, unit_no)

        # eliminate 'isMultiunit' f already exists
        if self.h5file.get_node(node_name).__contains__('isMultiunit'):
            self.h5file.remove_node(node_name, 'isMultiunit')

        # create a new "isMultiunt" array to hold the value of the cehckbox
        self.h5file.create_array(node_name, 'isMultiunit', bool(sender.checkState()))

        # save changes to disk
        self.h5file.flush()

    #__________________________________________________________________________
    def exchange_unit_name(self, initial, final):
        tb = self.chan_tab['unit_tab_bar_widget']

        # get the names of the changed tabs
        old_name_bg_tab = tb.tabText(final)
        new_name_bg_tab = tb.tabText(initial)

        # change the name of the background tab
        tb.setTabText(final, new_name_bg_tab)

        # change the name of the front tab the oldname of the unit
        tb.setTabText(tb.currentIndex(), old_name_bg_tab)

        # PROPAGATE CHANGES TO THE H5FILE #####

        # first change the background moved unit name to "tmpUnitData"
        self.h5file.rename_node(where='/Spikes/Chan_%03d' % self.cur_chan,
                               name=old_name_bg_tab, newname='tmpUnitData',
                               overwrite=True)

        # second change the front moved unit name to the old name of the background unit
        self.h5file.rename_node(where='/Spikes/Chan_%03d' % self.cur_chan,
                               name=new_name_bg_tab, newname=old_name_bg_tab,
                               overwrite=True)

        # third change the background moved unit name to its new name
        self.h5file.rename_node(where='/Spikes/Chan_%03d' % self.cur_chan,
                               name='tmpUnitData', newname=new_name_bg_tab,
                               overwrite=True)

        # CHANGE THE NAME OF THE FIGURES #####

        # first change the figure name of the background unit to "tmpFigName"
        for k in self.chan_tab['unit_figures']:
            if k.objectName() == old_name_bg_tab:
                k.setObjectName('tmpFigName')
                break

        # second, change the front tab figure name to the old background tab name
        for k in self.chan_tab['unit_figures']:
            if k.objectName() == new_name_bg_tab:
                k.setObjectName(old_name_bg_tab)
                break

        # third, change the figname of the background tab to the new one
        for k in self.chan_tab['unit_figures']:
            if k.objectName() == 'tmpFigName':
                k.setObjectName(new_name_bg_tab)
                break

        # CHANGE UNIT COLOR #####
        self.change_unit_color(unit_name=new_name_bg_tab,
                                  color=tuple(np.append(np.int32(self.unit_colors[final] * 255), 255)))

        self.change_unit_color(unit_name=old_name_bg_tab,
                                  color=tuple(np.append(np.int32(self.unit_colors[initial] * 255), 255)))

    #__________________________________________________________________________
    def repair_unit_names(self):

        unit_names = [k for k in self.cur_node.__members__ if 'Unit' in k]

        for j, k in enumerate(unit_names):
            if k != 'Unit%02d' % j:
                self.h5file.rename_node(self.cur_node_name, k, 'Unit%02d' % j)

        self.h5file.flush()

    #__________________________________________________________________________
    def call_merge_units(self):
        if not self.h5_file_loaded:
            return
        self.merge_units_widget.list1.clear()
        self.merge_units_widget.list2.clear()

        units_list = [k for k in self.cur_node.__members__ if 'Unit' in k]
        units_list.sort()
        self.merge_units_widget.list1.addItems(units_list)
        self.merge_units_widget.show()

    #__________________________________________________________________________
    def merge_units(self):

        # get the list of units to merge
        units_to_merge = [str(self.merge_units_widget.list2.item(k).text()) for k in range(self.merge_units_widget.list2.count())]
        # sort the names
        units_to_merge.sort()
        # if fewer than 2 return
        if len(units_to_merge) < 2:
            return

        # store the unit indexes in a list, sort them, and trnasform it into an array
        new_unit = []
        for k in units_to_merge:
            new_unit.extend(self.cur_node.__getattr__(k).Indx.read())
        new_unit.sort()
        new_unit = np.array(new_unit)

        # remove all the listed units from the h5file
        for k in units_to_merge:
            self.h5file.remove_node(self.cur_node_name, k, recursive=True)

        # create a group with the name of the first unit in the list, and
        # add all the indices of that
        self.h5file.create_group(self.cur_node_name, units_to_merge[0])
        self.h5file.create_array(self.cur_node_name + '/' + units_to_merge[0], 'Indx', new_unit)
        self.h5file.create_array(self.cur_node_name + '/' + units_to_merge[0], 'isMultiunit', False)
        self.h5file.create_array(self.cur_node_name + '/' + units_to_merge[0], 'isBursting', False)

        # save changes to disk
        self.h5file.flush()

        # add log
        self.add_log('%s %s merged' % (self.cur_node_name, str(units_to_merge)))

        # REMOVE ALL THE GRAPHICAL ELEMENTS #####
        # get the axes to remove from
        ax = self.chan_tab['waves_figure'].figure.axes[0]

        for k in units_to_merge[1:]:
            # remove the tabs
            for tab_index in range(self.chan_tab['unit_tabs_widget'].count()):
                if str(self.chan_tab['unit_tabs_widget'].tabText(tab_index)) == k:
                    self.chan_tab['unit_tabs_widget'].removeTab(tab_index)

            # remove unit figure
            self.chan_tab['unit_figures'][k].figure.clear()
            self.chan_tab['unit_figures'][k].close()
            self.chan_tab['unit_figures'].pop(k, 0)

            # removes the unitname from the what2 plot list
            for n in range(self.what_to_plot.count()):
                if self.what_to_plot.itemText(n) == k:
                    self.what_to_plot.removeItem(n)

            # eliminate the raw waveforms from the plot
            for line in ax.lines:
                if k in line.get_label():
                    line.remove()

            # remove the unit from the list
            unit_no = int(re.search('[0-9]{2}', k).group())
            self.units_list.remove(unit_no)

        # update the information in the overview table
        #self.overview_tab_2['overview_table'].takeItem(self.chans_list.index(self.cur_chan),
        #                                            unit_no+4)

        # redraw the waveforms figure
        self.chan_tab['waves_figure'].figure.canvas.draw()

        # replot features
        self.plot_features()

        # add the merged unit to the table
        self.units_table_add_unit(units_to_merge[0])

    #__________________________________________________________________________
    def call_move_units(self):
        if not self.h5_file_loaded:
            return
        self.move_units_widget.list.clear()
        units_list = [k for k in self.cur_node.__members__ if 'Unit' in k]
        units_list.sort()
        self.move_units_widget.list.addItems(units_list)
        self.move_units_widget.show()

    #__________________________________________________________________________
    def move_units(self):

        # first get the needed changes
        old = []
        new = []
        for k in range(self.move_units_widget.list.count()):
            if 'Unit%02d' % k != str(self.move_units_widget.list.item(k).text()):
                old.append(str(self.move_units_widget.list.item(k).text()))
                new.append('Unit%02d' % k)

        # in case no changes are needed
        if len(old) == 0:
            return

        # RENAME ALL THE UNITS AND GRAPHICAL ELEMENTS TO "_tmp" #####
        for k in self.cur_node.__members__:
            if 'Unit' in k:
                # rename the nodes
                self.h5file.rename_node(self.cur_node_name, name=k, newname=k + '_tmp')

                for key in ['unit_figures', 'unit_count_label', 'del_unit_buttons',
                            'plot_raw_check', 'unit_buttons', 'is_multiunit_check']:
                    self.chan_tab[key][k + '_tmp'] = self.chan_tab[key][k]
                    self.chan_tab[key].pop(k, 0)  # remove

        for j, k in zip(old, new):
            self.change_unit_name(j + '_tmp', k)

        # move everything back
        for k in self.cur_node.__members__:
            if '_tmp' in k:
                if k.replace('_tmp', '') in self.cur_node.__members__:
                    self.h5file.remove_node(self.cur_node_name, name=k)
                    for key in ['unit_figures', 'unit_count_label', 'del_unit_buttons',
                                'plot_raw_check', 'unit_buttons', 'is_multiunit_check']:
                        self.chan_tab[key].deleteLater()
                        self.chan_tab[key].pop(k, 0)
                else:
                    self.h5file.rename_node(self.cur_node_name, name=k, newname=k.replace('_tmp', ''))
                    for key in ['unit_figures', 'unit_count_label', 'del_unit_buttons',
                                'plot_raw_check', 'unit_buttons', 'is_multiunit_check']:
                        self.chan_tab[key][k.replace('_tmp', '')] = self.chan_tab[key][k]
                        self.chan_tab[key].pop(k, 0)

        # save changes to disk
        self.h5file.flush()

    #__________________________________________________________________________
    def change_unit_name(self, old_name, new_name):

        # rename node
        self.h5file.rename_node(self.cur_node_name, name=old_name, newname=new_name,
                               overwrite=True)

        # get the unit numbers from the names
        old_unit_no = int(re.search('[0-9]{2}', old_name).group())
        new_unit_no = int(re.search('[0-9]{2}', new_name).group())

        # move the tab and change its name
        self.chan_tab['unit_tab_bar_widget'].setTabText(new_unit_no, new_name)
        self.chan_tab['unit_tab_bar_widget'].moveTab(old_unit_no, new_unit_no)

        for key in ['unit_figures', 'unit_count_label', 'del_unit_buttons',
                    'plot_raw_check', 'unit_buttons', 'is_multiunit_check']:
            self.chan_tab[key][new_name] = self.chan_tab[key][old_name]
            self.chan_tab[key][new_name].setObjectName(new_name)
            self.chan_tab[key].pop(old_name, 0)

        # change color of the unit
        self.change_unit_color(new_name, color=255 * self.unit_colors[new_unit_no])

    #__________________________________________________________________________
    def clean_waves_figure(self):
        self.chan_tab['waves_figure'].figure.canvas.draw()

    #__________________________________________________________________________
    def units_table_add_row(self):
        self.cur_unit = self.chan_tab['unit_tabs_widget'].currentIndex()
        self.cur_unit_name = self.chan_tab['unit_tabs_widget'].tabText(self.cur_unit)

    #__________________________________________________________________________
    def del_unit(self):

        if not self.h5_file_loaded or not self.chan_plotted:
            return

        # get sender
        sender = self.sender()

        # get unit name and number
        unit_name = str(sender.objectName())
        unit_no = int(re.search('(?<=Unit)[0-9]{2}', unit_name).group())

        # remove the unit from the list
        self.units_list.remove(unit_no)

        # get the indexes of the unit
        indx = self.h5file.get_node(self.cur_node_name, unit_name).Indx.read()

        # get unsorted, append the indexes from the unit and update that
        # to the h5file
        self.unsorted = self.h5file.get_node(self.cur_node_name, 'Unsorted').read()
        self.unsorted = np.append(self.unsorted, indx)
        self.unsorted.sort()
        self.h5file.remove_node(self.cur_node_name, 'Unsorted')
        self.h5file.remove_node(self.cur_node_name, unit_name, recursive=True)
        self.h5file.create_array(self.cur_node_name, 'Unsorted', self.unsorted)

        # add log
        self.add_log('%s %s deleted' % (self.cur_node_name, unit_name))

        # remove the tab
        for tab_index in range(self.chan_tab['unit_tabs_widget'].count()):
            if str(self.chan_tab['unit_tabs_widget'].tabText(tab_index)) == unit_name:
                break
        self.chan_tab['unit_tabs_widget'].removeTab(tab_index)

        # close and remove unit figure
        plt.close(self.chan_tab['unit_figures'][unit_name].figure)
        self.chan_tab['unit_figures'].pop(unit_name, 0)

        # removes the unitname from the what2 plot list
        for n in range(self.what_to_plot.count()):
            if self.what_to_plot.itemText(n) == unit_name:
                self.what_to_plot.removeItem(n)

        # update the information in the overview table
        self.overview_tab_2['overview_table'].takeItem(self.chans_list.index(self.cur_chan),
                                                    unit_no + 4)

        # eliminate the raw waveforms from the plot
        ax = self.chan_tab['waves_figure'].figure.axes[0]
        for line in ax.lines:
            if unit_name in line.get_label():
                line.remove()
                break

        # redraw the waveforms figure
        self.chan_tab['waves_figure'].figure.canvas.draw()

        # replot features
        self.plot_features()

    #__________________________________________________________________________
    def change_unit_color(self, unit_name=None, color=None):
        ''' Change unit color utility function
        inputs:
            unit_name : string containing the unit name
            color    : must be a four element RGB tuple from 0 to 255, for example,
                       the output of getRgB() output from a Qt Color instance.
                       The fourth element is the alpha (usually = to 255)'''

        if unit_name in [None, False]:
            sender = self.sender()
            unit_name = str(sender.text()).replace(' ', '')

        unit_no = int(re.search('[0-9]{1,3}', unit_name).group())

        if not np.any(color):
            c = QtWidgets.QColorDialog()
            color = c.getColor(sender.palette().color(1))
            if not color.isValid():
                return

        if isinstance(color, QtGui.QColor):
            qt_color = color
        else:
            qt_color = QtGui.QColor(*color)

        mpl_color = np.array(qt_color.getRgb()[0:3]) / 255.0

        if isinstance(self.sender(), QtWidgets.QPushButton) and \
           'Unit' in self.sender().text():
            self.sender().setStyleSheet('QPushButton {background: rgb%s}' % str(qt_color.getRgb()[0:3]))

        self.unit_colors[unit_no, :] = mpl_color

        # get the figure with a name equal to the current unit
        ax = self.chan_tab['unit_figures'][unit_name].figure.axes

        # iterate over axes to change colors
        for k in ax:
            for j in k.lines:
                j.set_color(mpl_color)
            for j in k.collections:
                j.set_color(mpl_color)
            for j in k.patches:
                j.set_color(mpl_color)

        # search a figure with a specific name
        self.chan_tab['unit_figures'][unit_name].figure.canvas.draw()

        # change the color of the raw waveforms
        for k in self.chan_tab['waves_figure'].figure.axes[0].lines:
            if re.search('Unit%02d' % unit_no, str(k.get_label())):
                k.set_color(mpl_color)
        self.chan_tab['waves_figure'].figure.canvas.draw()

        # change the color in the features plot
        for k in self.chan_tab['features_fig'].figure.axes[0].lines:
            if re.search('Unit%02d' % unit_no, str(k.get_label())):
                k.set_color(mpl_color)
        self.chan_tab['features_fig'].figure.canvas.draw()

        self.chan_tab['unit_tabs_widget'].tabBar().setTabTextColor(unit_no, qt_color)

    #__________________________________________________________________________
    def reset_channel_tab(self):
        ''' reset the units tab'''

        self.n_units = 0

        # clear the unit figures
        for k in self.chan_tab['unit_figures']:
            plt.close(self.chan_tab['unit_figures'][k].figure)
        self.chan_tab['unit_figures'] = {}

        # clean the button dictionaries
        for key in ['del_unit_buttons', 'unit_count_label', 'unit_buttons', 'plot_raw_check', 'is_multiunit_check']:
            for k in self.chan_tab[key].keys():
                self.chan_tab[key][k].deleteLater()
            self.chan_tab[key] = {}

        # Reset waves_figure canvas
        ax = self.chan_tab['waves_figure'].figure.axes[0]
        ax.cla()
        self.sample_waveform, = ax.plot([], color=[.5, .5, .5], lw=2,
                                       animated=True)
        ax.set_ylim(-1000, 1000)
        ax.set_xlim(-2, self.wf_size + 1)
        ax.tick_params(color=[.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax.spines.values():
            k.set_edgecolor([.5, .5, .5])
        self.slice_1_line = ax.axvline(0, color=[.5, .5, .5])
        self.slice_2_line = ax.axvline(0, color=[.5, .5, .5], linestyle='--')
        ax.grid(color=[.5, .5, .5])
        self.chan_tab['waves_figure'].figure.tight_layout()
        self.chan_tab['waves_figure'].figure.canvas.draw()

        # clean the 3d widget
        for k in self.widget_3d.items:
            self.widget_3d.removeItem(k)
        #self.Fig3d.clf()

        # set the current indexes of the X and Y variable-selecting comboboxes
        self.x_plot.setCurrentIndex(0)
        self.y_plot.setCurrentIndex(1)
        self.z_plot.setCurrentIndex(2)

        # reset Units list
        self.units_list = []

        # reset the units tabbed widget
        tabs = list(range(self.chan_tab['unit_tabs_widget'].count()))
        tabs.reverse()
        if len(tabs) > 0:
            for k in tabs:
                self.chan_tab['unit_tabs_widget'].removeTab(k)

        # reset the time scroll widget and axes
        self.time_scroll['v_zoom_slider'].setValue(1000)
        self.time_scroll['h_zoom_slider'].setValue(500)
        self.time_scroll['h_scroll_slider'].setValue(0)
        self.time_scroll['figure'].figure.axes[0].cla()
        self.time_scroll['Ax'].set_xticklabels([])
        self.time_scroll['Ax'].set_yticklabels([])
        self.time_scroll['figure'].figure.canvas.draw()

        # reset label
        self.n_pts_label.setText('')

        # reset the features figure:
        self.chan_tab['features_fig'].figure.clf()
        self.chan_tab['features_fig'].figure.canvas.draw()

        # delete KDTree object
        if hasattr(self, 'xy_data'):
            del self.xy_data

        if hasattr(self, 'cur_waveforms'):
            del self.cur_waveforms

        if hasattr(self, 'cur_ts'):
            del self.cur_ts

        # remove the pca from the dictionarys
        self.chan_tab.pop('pca', 0)

        # reset the channel tab name
        self.main_fig_tab_widget.setTabText(2, 'Channel Tab')

    #__________________________________________________________________________
    def slice_draw(self):

        sender = self.sender()
        fig = self.chan_tab['waves_figure'].figure
        ax = fig.axes[0]

        if sender.objectName() == 'Slice1':
            self.slice_1_line.set_xdata([sender.value()])

        elif sender.objectName() == 'Slice2':
            self.slice_2_line.set_xdata([sender.value()])

        ax.draw_artist(ax.patch)

        for k in ax.get_lines():
            ax.draw_artist(k)

        for k in ax.get_xgridlines():
            ax.draw_artist(k)

        for k in ax.get_ygridlines():
            ax.draw_artist(k)

        ax.draw_artist(ax.spines['top'])
        ax.draw_artist(ax.spines['left'])
        fig.canvas.update()
        fig.canvas.flush_events()
        #self.chan_tab['waves_fig_bg'] = fig.canvas.copy_from_bbox(ax.bbox)

    #__________________________________________________________________________
    def change_current_unit(self):
        '''set the current unit'''
        self.cur_unit = self.chan_tab['unit_tabs_widget'].currentIndex()
        self.cur_unit_name = str(self.chan_tab['unit_tabs_widget'].tabText(self.cur_unit))
        for k in range(self.what_to_plot.count()):
            if str(self.what_to_plot.itemText(k)) == self.cur_unit_name:
                self.what_to_plot.setCurrentIndex(k)
                break

    #__________________________________________________________________________
    def main_fig_tab(self):
        '''Change the toolbar tab acording to the selected view'''
        cur_tab = self.main_fig_tab_widget.currentIndex()
        cur_tab_name = str(self.main_fig_tab_widget.tabText(cur_tab))
        if cur_tab_name == 'Channels Overview' or cur_tab_name == 'Summary Table':
            self.tools_tab.setCurrentIndex(0)
        elif re.search('Chan [0-9]{1,2}', cur_tab_name):
            self.tools_tab.setCurrentIndex(1)

    #__________________________________________________________________________
    def close_event(self, *event):
        '''Close the h5 file before killing the window.'''
        if self.h5_file_loaded:
            self.h5file.close()
        self.deleteLater()

    # Qt calls this exact method name for window close events.
    def closeEvent(self, *event):
        self.close_event(*event)


#==============================================================================
#if __name__ == '__main__':
#if not QtWidgets.QApplication.instance():
#    app = QtWidgets.QApplication(sys.argv)
#else:
#    app = QtWidgets.QApplication.instance()
spikesorter = SpikeSorter()
sys.exit(app.exec_())
