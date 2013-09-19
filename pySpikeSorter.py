## Initialization Routine

'''
Required Packages:
    - Numpy
    - Matplotlib
    - Mayavi2
    - Scipy
    - Spatial (KDTree)
    - PyQt4
    - guidata
    - pytables
'''

################ IMPORTS #################################################

import sys, os, re, tables
import numpy as np

# extra widgets import
import mayavi_widgets
import matplotlib_widgets
import helper_widgets

from m_Spike_Utils import autocorr
from PyQt4 import QtGui, QtCore

from matplotlib import rc
from matplotlib.mlab import PCA
from matplotlib import pyplot as plt

from matplotlib.path import Path
from scipy.spatial import cKDTree
import datetime

########## UTILITY FUNCTIONS #######################################################
        
## Spike Sorter Main GUI Window

rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

# create instance of imported widgets
settings     = helper_widgets.Settings()
autocorropts = helper_widgets.AutocorrOpts()
autoclust    = helper_widgets.AutoClustWidget()
    
class pySpikeSorter(QtGui.QMainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.setWindowTitle("Spike Sorter GUI")
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap('spike_icon.png')))
        self.main_widget = QtGui.QWidget(self)
        self.MainLayout  = QtGui.QHBoxLayout(self.main_widget)
        self.MainLayout.setMargin(0)
        self.MainLayout.setSpacing(0)
        
        ###########  SET SEVERAL VARIABLES ##############################################

        self.CurUnit      = 0
        self.PlotUnitCounter =0
        self.UnitsList    = []
        self.NUnits       = 0
        self.H5FileLoaded = False
        self.ChanPlotted  = False
        self.RemovingTab  = 0
        self.UnitColors   = np.array([[  1,   0,   0],
                                      [  0, 0.7,   0],
                                      [  0, 0.4,   1],
                                      [0.8, 0.6,   0],
                                      [0.6,   0,   1],
                                      [  0, 0.7, 0.7],
                                      [  0, 0.5,   1]])
        self.UnitColors = np.tile(self.UnitColors, (10,1))

        ########### TOOLBAR ON THE LEFT SIDE ################################

        split1 = QtGui.QSplitter(QtCore.Qt.Horizontal, self.main_widget)   ## SPLITTER
        self.ToolsTab = QtGui.QTabWidget()

        ToolsTab1 = QtGui.QWidget()
        ToolsTab2 = QtGui.QWidget()

        self.ToolsTab.addTab(ToolsTab1, 'Main Tools')
        self.ToolsTab.addTab(ToolsTab2, 'Chan Tools')
        self.ToolsTab.setMaximumWidth(220)

        split1.addWidget(self.ToolsTab)

        ########### self.ToolsTab No 1 #############################################
        toolslay = QtGui.QVBoxLayout()

        ### FRAME 1

        grp = QtGui.QGroupBox('Overview Tools', ToolsTab1)
        vlay = QtGui.QVBoxLayout()

        # number of events to overview spin box
        hlay = QtGui.QHBoxLayout()
        self.OverviewNEventsSpin = QtGui.QSpinBox()
        self.OverviewNEventsSpin.setRange(100, 1000)
        self.OverviewNEventsSpin.setSingleStep(100)
        self.OverviewNEventsSpin.setValue(500)
        hlay.addWidget(QtGui.QLabel('N Events 2 Overview'))
        hlay.addWidget(self.OverviewNEventsSpin)
        vlay.addLayout(hlay)

        # Y axis limits selector
        hlay = QtGui.QHBoxLayout()
        self.OverviewYLimsSpin = QtGui.QSpinBox()
        self.OverviewYLimsSpin.setRange(100, 5000)
        self.OverviewYLimsSpin.setSingleStep(100)
        self.OverviewYLimsSpin.setValue(2000)
        self.OverviewYLimsSpin.editingFinished.connect(self.ChangeOverviewYLim_Proc)
        hlay.addWidget(QtGui.QLabel('Overview Axes YLim'))
        hlay.addWidget(self.OverviewYLimsSpin)
        vlay.addLayout(hlay)

        
        btn = QtGui.QPushButton('Plot Overview')
	btn.setStyleSheet('QPushButton{background-color: rgba(0,190,0)}')
        btn.clicked.connect(self.LoadH5File)
        vlay.addWidget(btn)

        btn = QtGui.QPushButton('Save Overview')
        btn.clicked.connect(self.SaveOverviewFig_proc)
        vlay.addWidget(btn)

        grp.setLayout(vlay)
        toolslay.addWidget(grp)

        grp = QtGui.QGroupBox('Delete Channel Tools', ToolsTab1)
        vlay = QtGui.QVBoxLayout()
        
        # add mark trash spin and button and link it to a function
        hlay = QtGui.QHBoxLayout()
        self.MarkTrashSpin = QtGui.QSpinBox()
        self.MarkTrashSpin.setMinimum(1)
        self.MarkTrashSpin.setMaximum(1000000)
        self.MarkTrashSpin.setValue(1000)
        hlay.addWidget(QtGui.QLabel('Below'))
        hlay.addWidget(self.MarkTrashSpin)
        MarkTrashBtn  = QtGui.QPushButton('Mark Trash')
        MarkTrashBtn.clicked.connect(self.TrashChans_proc)
        hlay.addWidget(MarkTrashBtn)
        vlay.addLayout(hlay)

        # add delete trash chans and link it to a function
        btn = QtGui.QPushButton('Delete Trash Chans')
        btn.clicked.connect(self.DeleteTrashChans_proc)
        vlay.addWidget(btn)
        
        grp.setLayout(vlay)
        toolslay.addWidget(grp)

        ### FRAME 2
        grp = QtGui.QGroupBox('Channel Plot Options', ToolsTab1)
        vlay = QtGui.QVBoxLayout()

        hlay = QtGui.QHBoxLayout()
        self.ChanSelector = QtGui.QComboBox()
        hlay.addWidget(QtGui.QLabel('Chan Selector'))
        hlay.addWidget(self.ChanSelector)
        vlay.addLayout(hlay)

        PlotChanBtn = QtGui.QPushButton('Plot Chan')
        PlotChanBtn.clicked.connect(self.PlotChanProc)
        vlay.addWidget(PlotChanBtn)

        grp.setLayout(vlay)
        toolslay.addWidget(grp)

        ######### Group No3 #####################

        grp = QtGui.QGroupBox('General Tools', ToolsTab1)
        glay = QtGui.QGridLayout()
        
	'''
        convertFileBtn = QtGui.QPushButton('Convert File')
        convertFileBtn.clicked.connect(self.ConvertNevFile)
        glay.addWidget(convertFileBtn, 0,0)

        bin2H5Btn = QtGui.QPushButton('Bin2H5')
        bin2H5Btn.clicked.connect(bin2h5.bin2h5)
        glay.addWidget(bin2H5Btn, 0, 1)
	'''	
	
        setSettigsBtn = QtGui.QPushButton('Settings')
        setSettigsBtn.clicked.connect(self.Settings)
        glay.addWidget(setSettigsBtn, 0, 0)

        aboutBtn = QtGui.QPushButton('About')
        aboutBtn.clicked.connect(self.About)
        glay.addWidget(aboutBtn, 0, 1)

        closeH5FileBtn = QtGui.QPushButton('Close H5 File')
        closeH5FileBtn.clicked.connect(self.CloseFile)
        glay.addWidget(closeH5FileBtn, 1, 0)
        
        exitBtn = QtGui.QPushButton('Exit')
        exitBtn.clicked.connect(self.closeEvent)
        glay.addWidget(exitBtn, 1, 1)

        grp.setLayout(glay)
        toolslay.addWidget(grp)

        # create an "About" Msg Box
        self.AboutMsg = QtGui.QMessageBox(QtGui.QMessageBox.Information,
                                          'About',
                                          u'Spyke Sorter v0.1\nHachi Manzur, 2012')

        ''' Code to have a menu bar (commented to gain space)
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')

        openFileAction = QtGui.QAction('&Open File', self)
        openFileAction.setStatusTip('Open a File')
        openFileAction.setShortcut('Ctrl+o')
        openFileAction.triggered.connect(self.SelFile)
        fileMenu.addAction(openFileAction)

        bin2H5Action = QtGui.QAction('&bin2H5', self)
        bin2H5Action.setStatusTip('Transofmr binary files into H5')
        bin2H5Action.triggered.connect(bin2h5)
        fileMenu.addAction(bin2H5Action)
        
        closeFileAction = QtGui.QAction('&Close H5 File', self)
        closeFileAction.setStatusTip('Close an H5 File')
        closeFileAction.setShortcut('Ctrl+x')
        closeFileAction.triggered.connect(self.CloseFile)
        fileMenu.addAction(closeFileAction)

        settingsAction = QtGui.QAction('&Settings', self)
        settingsAction.setStatusTip('Settings')
        settingsAction.triggered.connect(self.Settings)
        fileMenu.addAction(settingsAction)

        exitAction = QtGui.QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        fileMenu.addAction(exitAction)

        helpMenu = menubar.addMenu('&?')
        aboutAction = QtGui.QAction('&About', self)
        aboutAction.setStatusTip('About')
        aboutAction.triggered.connect(self.About)
        helpMenu.addAction(aboutAction)'''
      
        toolslay.addStretch(1)

        ToolsTab1.setLayout(toolslay)

        ########### self.ToolsTab No 2 #############################################

        toolslay = QtGui.QVBoxLayout()

        # group No1
        grp = QtGui.QGroupBox('Features Plot Opts', ToolsTab2)
        vlay = QtGui.QVBoxLayout()

        # add X and Y features selection combobox
        items=['PCA1','PCA2','PCA3',
               'Slice1','Slice2','Time','Pk2Pk Amp',
               'Peak','Valley','Energy',
               'Peak Pt','Valley Pt']
        self.XPlot = QtGui.QComboBox(grp)
        self.YPlot = QtGui.QComboBox(grp)
        self.ZPlot = QtGui.QComboBox(grp)
        self.XPlot.addItems(items)
        self.YPlot.addItems(items)
        self.ZPlot.addItems(items)
        self.ZPlot.addItem('Density')
        self.YPlot.setCurrentIndex(1)
        self.ZPlot.setCurrentIndex(2)

        # add the X axis combo box
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('X Axis Variable'))
        hlay.addWidget(self.XPlot)
        vlay.addLayout(hlay)

        # add the Y axis combo box
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('Y Axis Variable'))
        hlay.addWidget(self.YPlot)
        vlay.addLayout(hlay)

        # add the Y axis combo box
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('Z Axis Variable'))
        hlay.addWidget(self.ZPlot)
        vlay.addLayout(hlay)

        # add a source of what to plot selection combo box
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('What to Plot ?'))
        self.What2Plot = QtGui.QComboBox()
        hlay.addWidget(self.What2Plot)
        vlay.addLayout(hlay)

        # add two slice selection spin box
        hlay = QtGui.QHBoxLayout()
        self.SliceSpBx1 = QtGui.QSpinBox()
        self.SliceSpBx1.setObjectName('Slice1')
        self.SliceSpBx1.valueChanged.connect(self.SliceDraw)
        hlay.addWidget(QtGui.QLabel('Slice 1'))
        hlay.addWidget(self.SliceSpBx1)
        self.SliceSpBx2 = QtGui.QSpinBox()
        self.SliceSpBx2.setObjectName('Slice2')
        hlay.addWidget(QtGui.QLabel('Slice 2'))
        self.SliceSpBx2.valueChanged.connect(self.SliceDraw)
        hlay.addWidget(self.SliceSpBx2)
        vlay.addLayout(hlay)

        # add a plot density check and a spin box to set the resolution
        hlay = QtGui.QHBoxLayout()
        self.PlotDensityCheck = QtGui.QCheckBox('Plot Density ?')
        hlay.addWidget(self.PlotDensityCheck)
        self.PlotDensityBins = QtGui.QSpinBox()
        self.PlotDensityBins.setMinimum(50)
        self.PlotDensityBins.setMaximum(300)
        self.PlotDensityBins.setValue(50)
        hlay.addWidget(self.PlotDensityBins)
        vlay.addLayout(hlay)

        # plot only valid Wfs check widget
        self.PlotValidsOnlyCheck = QtGui.QCheckBox('Plot Valids Only')
        self.PlotValidsOnlyCheck.setChecked(True)
        vlay.addWidget(self.PlotValidsOnlyCheck)

        # label with number of points
        hlay = QtGui.QHBoxLayout()
        self.nPtsLabel = QtGui.QLabel()
        hlay.addWidget(QtGui.QLabel('NPoints'))
        hlay.addWidget(self.nPtsLabel)
        vlay.addLayout(hlay)

        # number of spikes spin box
        hlay = QtGui.QHBoxLayout()
        self.nPtsSpin = QtGui.QSpinBox()
        self.nPtsSpin.setRange(10000, 200000)
        self.nPtsSpin.setSingleStep(10000)
        hlay.addWidget(self.nPtsSpin)

        # number of spikes slider
        self.nPtsSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.nPtsSlider.setRange(10000, 200000)
        self.nPtsSlider.setTickInterval(5000)
        self.nPtsSlider.setSingleStep(5000)
        hlay.addWidget(self.nPtsSlider)

        # connect spinner with No-of-spikes slider
        self.nPtsSlider.valueChanged.connect(self.nPtsSpin.setValue)

        # connect slider with No-of-spikes spinner
        self.nPtsSpin.valueChanged.connect(self.nPtsSlider.setValue)
        
        # set N spikes value
        self.nPtsSlider.setValue(50000)
        vlay.addLayout(hlay)

        hlay = QtGui.QHBoxLayout()
        # plot features btn and funcion connection
        self.PlotFeaturesBtn = QtGui.QPushButton('Plot 2D', grp)
        self.PlotFeaturesBtn.clicked.connect(self.PlotFeatures)
        hlay.addWidget(self.PlotFeaturesBtn)

        # plot features btn and funcion connection
        self.Plot3DBtn = QtGui.QPushButton('Plot 3D', grp)
        self.Plot3DBtn.clicked.connect(self.Plot3DFeatures)
        hlay.addWidget(self.Plot3DBtn)

        vlay.addLayout(hlay)

        grp.setLayout(vlay)
        
        toolslay.addWidget(grp)

        ### group No 2 ###
        grp  = QtGui.QGroupBox('Raw Waveforms Opts')
        vlay = QtGui.QVBoxLayout()

        # number of spikes spin box
        hlay = QtGui.QHBoxLayout()
        self.NSpikesSpin = QtGui.QSpinBox()
        self.NSpikesSpin.setMaximum(5000)
        self.NSpikesSpin.setMinimum(100)
        self.NSpikesSpin.setSingleStep(100)
        hlay.addWidget(self.NSpikesSpin)

        # number of spikes slider
        self.NSpikesSlider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.NSpikesSlider.setMaximum(5000)
        self.NSpikesSlider.setMinimum(100)
        self.NSpikesSlider.setSingleStep(100)
        hlay.addWidget(self.NSpikesSlider)

        # connect spinner with No-of-spikes slider
        self.NSpikesSpin.valueChanged.connect(self.NSpikesSlider.setValue)
        
        # connect slider with No-of-spikes spinner
        self.NSpikesSlider.valueChanged.connect(self.NSpikesSpin.setValue)
        
        # set N spikes value
        self.NSpikesSlider.setValue(1000)
        vlay.addLayout(hlay)

        # add axes limit spin box
        hlay = QtGui.QHBoxLayout()
        self.WaveAxYLim_Spin = QtGui.QSpinBox()
        self.WaveAxYLim_Spin.setRange(0, 10000)
        self.WaveAxYLim_Spin.setValue(1000)
        self.WaveAxYLim_Spin.setSingleStep(100)
        self.WaveAxYLim_Spin.editingFinished.connect(self.SetWfPlotLim_proc)
        hlay.addWidget(QtGui.QLabel('Axes Y Lim'))
        hlay.addWidget(self.WaveAxYLim_Spin)
        vlay.addLayout(hlay)

        # create a "plot waveforms" check widget
        self.PlotWaveformsCheck = QtGui.QCheckBox('Plot Raw Waveforms ?')
        vlay.addWidget(self.PlotWaveformsCheck)

        grp.setLayout(vlay)

        toolslay.addWidget(grp)

        # Automatic clustering box
        #w = QtGui.QWidget()
        autoClustBtn = QtGui.QPushButton('Automatic Clustering')
        autoClustBtn.clicked.connect(autoclust.show)
        toolslay.addWidget(autoClustBtn)

        hlay = QtGui.QHBoxLayout()
        mergeUnitsBtn = QtGui.QPushButton('Merge Units')
        self.MergeUnitsWidget = helper_widgets.MergeUnitsWidget()
        self.MergeUnitsWidget.AcceptBtn.clicked.connect(self.MergeUnits_proc)
        mergeUnitsBtn.clicked.connect(self.CallMergeUnits_proc)
        hlay.addWidget(mergeUnitsBtn)
        
        moveUnitsBtn = QtGui.QPushButton('Move Units')
        self.MoveUnitsWidget = helper_widgets.MoveUnitsWidget()
        self.MoveUnitsWidget.AcceptBtn.clicked.connect(self.MoveUnits_proc)
        moveUnitsBtn.clicked.connect(self.CallMoveUnits_proc)
        hlay.addWidget(moveUnitsBtn)
        
        toolslay.addLayout(hlay)
        
        ##### CHANNEL METAINFO GROUP #####

        '''grp = QtGui.QGroupBox('Channel Metainfo')
        vlay = QtGui.QVBoxLayout()
        grp.setLayout(vlay)
        toolslay.addWidget(grp)
        # button to tight_layout() on resize of the main window
        #btn = QtGui.QPushButton('Adjust')
        #btn.clicked.connect(self.AdjustPlots_proc)
        #toolslay.addWidget(btn)'''

        # button to reset a channel
        btn = QtGui.QPushButton('Reset Channel')
        btn.clicked.connect(self.ResetChan_proc)
        toolslay.addWidget(btn)

        # button to reset a channel
        btn = QtGui.QPushButton('Autocorr Opts')
        btn.clicked.connect(self.AutocorrOpts)
        toolslay.addWidget(btn)

        # replot units
        #btn = QtGui.QPushButton('Replot Unit')
        #btn.clicked.connect(self.PlotUnitFigure_proc)
        #toolslay.addWidget(btn)
        
        toolslay.addStretch(1)
        ToolsTab2.setLayout(toolslay)

        ########### TABBED FIGURES WIDGET #################################

        self.OverviewTab1 = {}
        self.OverviewTab2 = {}

        self.MainFigTab = QtGui.QTabWidget()
        self.MainFigTab.currentChanged.connect(self.MainFigTabProc)
        self.OverviewTab1['MainWidget']  = QtGui.QWidget(self.MainFigTab)
        hlay = QtGui.QHBoxLayout(self.OverviewTab1['MainWidget'])

        self.MainFigTab.addTab(self.OverviewTab1['MainWidget'],'Channels Overview')

        # overview figure
        self.OverviewTab1['Figure']  = matplotlib_widgets.MplWidget()
        self.OverviewTab1['Figure'].figure.set_facecolor('k')
        self.OverviewTab1['Toolbar'] = matplotlib_widgets.NavToolbar(self.OverviewTab1['Figure'], self.OverviewTab1['MainWidget'])
        self.OverviewTab1['Toolbar'].setIconSize(QtCore.QSize(15,15))
        vlay = QtGui.QVBoxLayout()
        vlay.addWidget(self.OverviewTab1['Figure'])
        vlay.addWidget(self.OverviewTab1['Toolbar'])
        vlay.setMargin(0)
        vlay.setSpacing(1)
        hlay.addLayout(vlay)
        hlay.setMargin(0)
        hlay.setSpacing(1)
        
        ########### OVERVIEW TABLE WIDGET #################################

        self.OverviewTab2['MainWidget']    = QtGui.QWidget(self.MainFigTab)
        self.OverviewTab2['OverviewTable'] = QtGui.QTableWidget(0,6,self.OverviewTab2['MainWidget'])
        self.OverviewTab2['OverviewTable'].setAlternatingRowColors(True)
        self.OverviewTab2['OverviewTable'].setFont(QtGui.QFont('sans',8))
        labels = ['Count','isTrash','MultiUnit?','Comments','Unsorted','Valid']
        self.OverviewTab2['OverviewTable'].setHorizontalHeaderLabels(labels)
        for k in range(self.OverviewTab2['OverviewTable'].columnCount()):
            self.OverviewTab2['OverviewTable'].setColumnWidth(k,65)
        self.OverviewTab2['OverviewTable'].setColumnWidth(3,150)
        self.OverviewTab2['OverviewTable'].setColumnWidth(2,75)

        # associate the vertical header click to select the channel
        vHeader = self.OverviewTab2['OverviewTable'].verticalHeader()
        vHeader.sectionClicked.connect(self.TableRowChanged_proc)

        vlay = QtGui.QVBoxLayout(self.OverviewTab2['MainWidget'])
        vlay.addWidget(self.OverviewTab2['OverviewTable'])

        # add a log entry browser
        grp = QtGui.QGroupBox('Log Browser')
        grp.setMaximumHeight(100)
        hlay = QtGui.QHBoxLayout()
        self.LogCombo = QtGui.QComboBox()
        self.LogCombo.setMinimumWidth(200)
        #self.LogCombo.setMinimumHeight(20)
        self.LogCombo.currentIndexChanged.connect(self.SetLogText_proc)
        hlay.addWidget(self.LogCombo)
        self.LogTextBrowser = QtGui.QTextBrowser()
        #self.LogTextBrowser.setMaximumHeight(40)
        hlay.addWidget(self.LogTextBrowser)
        hlay.setMargin(0)
        hlay.setSpacing(1)
        grp.setLayout(hlay)
        vlay.addWidget(grp)
        
        self.OverviewTab2['MainWidget'].setLayout(hlay)
        self.MainFigTab.addTab(self.OverviewTab2['MainWidget'],'Summary Table')
        
        ##### CHANNEL TAB #####

        self.ChanTab = {}
        curfigtab = self.MainFigTab.currentIndex()-1
        self.ChanTab['MainWidget'] = QtGui.QWidget()
        self.MainFigTab.addTab(self.ChanTab['MainWidget'],'Channel Tab')
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal, self.ChanTab['MainWidget'])

        mainHLay = QtGui.QHBoxLayout()

        ##### RAW WAVEFORMS WIDGET #####
        # create the mpl widget to plot the raw waveforms
        vlay = QtGui.QVBoxLayout()

        # buttons and controls on top of raw waveforms plot
        hlay = QtGui.QHBoxLayout()
        self.NUnitsSpin = QtGui.QSpinBox()
        self.NUnitsSpin.setMaximumHeight(20)
        self.NUnitsSpin.setMinimum(1)
        self.NUnitsSpin.setMaximum(10000)
        self.NUnitsSpin.setValue(1)
                
        TrimBtn  = QtGui.QPushButton('Trim Waveforms')
        TrimBtn.clicked.connect(self.ActivateTrimWaveforms_proc)
        TrimBtn.setMaximumHeight(20)
        
        CleanBtn = QtGui.QPushButton('Redraw')
        CleanBtn.setMaximumHeight(20)
        
        CleanBtn.clicked.connect(self.CleanWavesFigure_proc)
        hlay.addStretch(1)
        lbl = QtGui.QLabel('Waveforms2Plot:')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.NUnitsSpin)
        hlay.addWidget(TrimBtn)
        hlay.addWidget(CleanBtn)
        hlay.addStretch(1)
        vlay.addLayout(hlay)

        # waveforms plot and toolbar
        hlay = QtGui.QHBoxLayout()
        self.ChanTab['WavesFigure'] = matplotlib_widgets.MplWidget()
        self.ChanTab['WavesFigure'].figure.set_facecolor('k')
        self.ChanTab['WaveToolbar'] = matplotlib_widgets.NavToolbar(self.ChanTab['WavesFigure'], self.ChanTab['MainWidget'], coordinates=False)
        self.ChanTab['WaveToolbar'].setIconSize(QtCore.QSize(15,15))
        self.ChanTab['WaveToolbar'].setOrientation(QtCore.Qt.Vertical)
        self.ChanTab['WaveToolbar'].setMaximumWidth(30)
##        self.ChanTab['WaveToolbar'].setFloatable(True)    
##        self.ChanTab['WaveToolbar'].setMovable(True)
        hlay.addWidget(self.ChanTab['WavesFigure'])
        hlay.addWidget(self.ChanTab['WaveToolbar'])
        hlay.setMargin(0)
        hlay.setSpacing(1)
        vlay.addLayout(hlay)

        ###### UNIT TABS WIDGET ######
        
        self.ChanTab['UnitTabsWidget'] = QtGui.QTabWidget()
        #self.ChanTab['UnitTabsWidget'].setMovable(True)
        self.ChanTab['UnitTabBarWidget'] = self.ChanTab['UnitTabsWidget'].tabBar()
        #self.ChanTab['UnitTabBarWidget'].tabMoved.connect(self.ExchangeUnitName_proc)
        self.ChanTab['UnitTabsWidget'].setMaximumHeight(QtGui.QApplication.desktop().availableGeometry().height()/4)
        self.ChanTab['UnitFigures']      = {}
        self.ChanTab['DelUnitBtns']      = {}
        self.ChanTab['UnitCountLabel']   = {}
        self.ChanTab['UnitBtns']         = {}
        self.ChanTab['PlotRawCheck']     = {}
        self.ChanTab['isMultiunitCheck'] = {}
        self.ChanTab['UnitTabsWidget'].currentChanged.connect(self.ChangeCurrentUnit_proc)
        vlay.addWidget(self.ChanTab['UnitTabsWidget'])

        mainHLay.addLayout(vlay)
        
        # configures the waveforms figure
        wavesfig = self.ChanTab['WavesFigure'].figure
        ax = wavesfig.add_subplot(111)
        self.trimWaveformsRect = matplotlib_widgets.MyRectangleSelector(ax, self.TrimWaveforms_proc, drawtype='line', useblit=True)
        self.trimWaveformsRect.set_active(False)
        ax.set_axis_bgcolor('k')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        self.SampleWaveform, = ax.plot([], color=[.5,.5,.5], linewidth=2)
        self.Waveforms, = ax.plot([], animated=True)
        ax.set_ylim(-1000, 1000)
        ax.set_xlim(0,32)

        # Create Slice plots
        self.Slice1Ln = ax.axvline(0, color=[.5, .5, .5])
        self.Slice2Ln = ax.axvline(0, color=[.5, .5, .5], linestyle = '--')
        #ax.set_title('Raw Waveforms')
        ax.grid()
        wavesfig.canvas.mpl_connect('draw_event',self.draw_callback)

        ###### FEATURES PLOT WIDGET #################

        mainRightLay = QtGui.QVBoxLayout()
        
        tab = QtGui.QTabWidget()

        widget = QtGui.QWidget()
        # function buttons on top of the features plot:
        
        vlay = QtGui.QVBoxLayout(widget)

        hlay = QtGui.QHBoxLayout()
        hlay.addStretch(1)

        self.AddUnitBtn = QtGui.QPushButton('Add Unit')
        self.AddUnitBtn.setMaximumHeight(20)
        self.AddUnitBtn.clicked.connect(self.AddUnit_proc)
        hlay.addWidget(self.AddUnitBtn)

        # add a "keep" button
        self.KeepBtn = QtGui.QPushButton('Keep')
        self.KeepBtn.setMaximumHeight(20)
        self.KeepBtn.setToolTip('Create new unit (only when All waveforms or Unsorted are plotted)')
        self.KeepBtn.clicked.connect(self.Keep_proc)
        hlay.addWidget(self.KeepBtn)

        # add an "add region" button
        self.AddRegionBtn = QtGui.QPushButton('Add Region')
        self.AddRegionBtn.setMaximumHeight(20)
        self.AddRegionBtn.setToolTip('Add waveforms to the current unit')
        self.AddRegionBtn.clicked.connect(self.AddRegion_proc)
        hlay.addWidget(self.AddRegionBtn)

        # add a "remove region" button
        self.RemoveRegionBtn = QtGui.QPushButton('Remove Region')
        self.RemoveRegionBtn.setMaximumHeight(20)
        self.RemoveRegionBtn.clicked.connect(self.RemoveRegion_proc)
        hlay.addWidget(self.RemoveRegionBtn)

        # "set valid waveforms" button
        self.ValidWFBtn = QtGui.QPushButton('Set Valid WFs')
        self.ValidWFBtn.setMaximumHeight(20)
        self.ValidWFBtn.clicked.connect(self.ValidateWFs_proc)
        hlay.addWidget(self.ValidWFBtn)

        # "set valid waveforms" button
        self.ReplotDensityBtn = QtGui.QPushButton('Replot Density')
        self.ReplotDensityBtn.setMaximumHeight(20)
        self.ReplotDensityBtn.clicked.connect(self.ReplotDensity_proc)
        hlay.addWidget(self.ReplotDensityBtn)

        hlay.addStretch(1)
        vlay.addLayout(hlay)

        # Features figure and toolbar
        self.ChanTab['FeaturesFig']    = matplotlib_widgets.MplWidget()
        self.ChanTab['FeaturesFig'].figure.set_facecolor('k')
        self.ChanTab['FeaturesFigNtb'] = matplotlib_widgets.NavToolbar(self.ChanTab['FeaturesFig'].figure.canvas,
                                                            self.ChanTab['MainWidget'])
        self.ChanTab['FeaturesFigNtb'].setIconSize(QtCore.QSize(15,15))
        self.ChanTab['FeaturesFigNtb'].setMaximumHeight(30)

        vlay.addWidget(self.ChanTab['FeaturesFig'])
        vlay.addWidget(self.ChanTab['FeaturesFigNtb'])

        vlay.setMargin(0)
        vlay.setSpacing(1)
        
        widget.setLayout(vlay)
        tab.addTab(widget,'2D')
        mainRightLay.addWidget(tab)

        ###### 3D TAB Widget  #######################################
        self.Widget3d = mayavi_widgets.MayaviQWidget()
        tab.addTab(self.Widget3d,'3D')
        self.Fig3d = self.Widget3d.visualization.scene.mlab

        ######## Spikes vs time visualization widget #################

        # add a figure adn axes
        self.TimeScroll = {}
        self.TimeScroll['Figure']  = matplotlib_widgets.MplWidget()
        self.TimeScroll['Figure'].figure.set_facecolor('k')
        self.TimeScroll['DrawFigCID']  = self.TimeScroll['Figure'].figure.canvas.mpl_connect('draw_event', self.DrawScrollFig_Func)
        self.TimeScroll['Figure'].setMaximumHeight(QtGui.QApplication.desktop().availableGeometry().height()/6)
        self.TimeScroll['Ax'] = self.TimeScroll['Figure'].figure.add_subplot(111)
        self.TimeScroll['Ax'].set_axis_bgcolor('k')
        self.TimeScroll['Ax'].set_ylim(-1500,1500)
        self.TimeScroll['Ax'].set_xticklabels([])
        self.TimeScroll['Ax'].set_yticklabels([])
        self.TimeScroll['Ax'].set_axis_off()
        self.TimeScroll['Plot'], = self.TimeScroll['Ax'].plot([],color=[.5, .5, .5])
        self.TimeScroll['Figure'].figure.tight_layout()
        self.TimeScroll['Figure'].figure.canvas.draw()

        # add a vertical zoom slider
        self.TimeScroll['VZoom'] = QtGui.QSlider(QtCore.Qt.Vertical)
        self.TimeScroll['VZoom'].setMaximumHeight(QtGui.QApplication.desktop().availableGeometry().height()/6)
        self.TimeScroll['VZoom'].setMinimum(100)
        self.TimeScroll['VZoom'].setMaximum(5000)
        self.TimeScroll['VZoom'].setValue(1000)
        self.TimeScroll['VZoom'].valueChanged.connect(self.VZoom_Func)
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(self.TimeScroll['VZoom'])
        hlay.addWidget(self.TimeScroll['Figure'])
        mainRightLay.addLayout(hlay)

        # add an horizontal zoom slider
        self.TimeScroll['HZoom']  = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.TimeScroll['HZoom'].setRange(5, 5000)
        self.TimeScroll['HZoom'].setValue(500)
        self.TimeScroll['HZoom'].setSingleStep(5)
        self.TimeScroll['HZoom'].valueChanged.connect(self.HZoom_Func)
        self.TimeScroll['HZoomSpin'] = QtGui.QSpinBox()
        self.TimeScroll['HZoomSpin'].setMinimumWidth(80)
        self.TimeScroll['HZoomSpin'].setMaximumHeight(20)
        self.TimeScroll['HZoomSpin'].setRange(5, 5000)
        self.TimeScroll['HZoomSpin'].setValue(500)
        self.TimeScroll['HZoomSpin'].setSingleStep(10)
        self.TimeScroll['HZoomSpin'].valueChanged.connect(self.TimeScroll['HZoom'].setValue)
        self.TimeScroll['HZoom'].valueChanged.connect(self.TimeScroll['HZoomSpin'].setValue)
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('H Span '))
        hlay.addWidget(self.TimeScroll['HZoomSpin'])
        hlay.addWidget(self.TimeScroll['HZoom'])
        mainRightLay.addLayout(hlay)

        # add a time slider
        self.TimeScroll['HScroll'] = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.TimeScroll['HScroll'].setRange(0, 3000000)
        self.TimeScroll['HScroll'].setSingleStep(self.TimeScroll['HZoom'].value()/10)
        self.TimeScroll['HScroll'].valueChanged.connect(self.HScroll_Func)
        self.TimeScroll['HSpin'] = QtGui.QSpinBox()
        self.TimeScroll['HSpin'].setRange(0, 3000000)
        self.TimeScroll['HSpin'].setMinimumWidth(80)
        self.TimeScroll['HSpin'].setMaximumHeight(20)
        self.TimeScroll['HSpin'].valueChanged.connect(self.TimeScroll['HScroll'].setValue)
        self.TimeScroll['HScroll'].valueChanged.connect(self.TimeScroll['HSpin'].setValue)
        hlay = QtGui.QHBoxLayout()
        hlay.addWidget(QtGui.QLabel('H Scroll'))
        hlay.addWidget(self.TimeScroll['HSpin'])
        hlay.addWidget(self.TimeScroll['HScroll'])
        mainRightLay.addLayout(hlay)

        mainRightLay.setMargin(0)
        mainRightLay.setSpacing(1)
        # add the widget to the main horizontal layout
        mainHLay.addLayout(mainRightLay)
        mainHLay.setMargin(1)
        self.ChanTab['MainWidget'].setLayout(mainHLay)

        # create a generic Msg box
        self.MsgBox = QtGui.QMessageBox()
               
        # if running in linux set a certain style for the buttons and widgets
        if sys.platform == 'linux2':
            QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Plastique'))

        # add the main tabbed figures widget to the main splitter
        split1.addWidget(self.MainFigTab)

        # add the splitter to the main layout
        self.MainLayout.addWidget(split1)

        # set the layout of the main widget
        self.main_widget.setLayout(self.MainLayout)

        # set the central widget of the application
        self.setCentralWidget(self.main_widget)

    ########################################################################################################
    ################### METHODS ############################################################################
    ########################################################################################################

    def SaveOverviewFig_proc(self):
        if self.H5FileLoaded:
            fname = str(QtGui.QFileDialog.getSaveFileName(directory = self.h5file.filename[0:-3]+'_sorted.jpg'))
            if fname:
                self.OverviewTab1['Figure'].figure.savefig(fname, dpi = 300, facecolor = 'k')
        
    ########################################################################################################
        
    
    def resizeEvent(self, *event):
        # reimplementation of the resize event to apply tight_layout to
        # all the figures
        pass

    ########################################################################################################
    
    def AdjustPlots_proc(self):
        self.TimeScroll['Figure'].figure.tight_layout()
        self.TimeScroll['Figure'].figure.canvas.draw()
        
        self.ChanTab['WavesFigure'].figure.tight_layout()
        self.ChanTab['WavesFigure'].figure.canvas.draw()

        if len(self.ChanTab['FeaturesFig'].figure.axes)>0:
            self.ChanTab['FeaturesFig'].figure.tight_layout()
            self.ChanTab['FeaturesFig'].figure.canvas.draw()

        if len(self.OverviewTab1['Figure'].figure.axes)>0:
            self.OverviewTab1['Figure'].figure.tight_layout()
            self.OverviewTab1['Figure'].figure.canvas.draw()

        for k in self.ChanTab['UnitFigures']:
            self.ChanTab['UnitFigures'][k].figure.tight_layout()
            self.ChanTab['UnitFigures'][k].figure.canvas.draw()

    ########################################################################################################
        
    def SetWfPlotLim_proc(self):
        sender = self.sender()
        ax = self.ChanTab['WavesFigure'].figure.axes[0]
        curlim = ax.get_ylim()
        lim = sender.value()           
        ax.set_ylim(-lim, lim)
        self.ChanTab['WavesFigure'].figure.canvas.draw()
        
        
    ########################################################################################################
            
    def HScroll_Func(self):
        ''' This function gets triggered whenever the user moves the bottom
        scrollbar in the lower right. It helps to explore the raw waveforms'''
        
        if not self.ChanPlotted: return
        self.TimeScroll['Figure'].figure.canvas.restore_region(self.TimeScroll['bg'])
        self.ChanTab['WavesFigure'].figure.canvas.restore_region(self.ChanTab['WavesFigBG'])

        v = self.TimeScroll['HScroll'].value()
        h = self.TimeScroll['HZoom'].value()
        indx = np.flatnonzero(np.logical_and(self.CurTs>=v, self.CurTs<(v+h)))

        if any(indx):
            # ontain the timestamps corresponding to the indexes
            Ts = self.CurNode.TimeStamp[indx]
            # substract the first timestamp to have a 0 based indexing
            Ts = Ts-v

            # obtain the waveforms to plot
            Wf = self.CurNode.Waveforms[indx,:]

            # obtain the length of units to plot
            n = len(indx)
            
            # create an array of Nones to append
            nones = np.array(n*[None], ndmin=2).transpose()

            # append nones to the waveforms array and reshape it to a vector
            Wf = np.append(Wf, nones, axis=1).reshape((n*(self.WfSize+1),))

            # create a vector time based on the sampling frequency, the
            # the number of points per spike and the timestamp
            Ts = np.tile(Ts, (self.WfSize,1)).transpose() + \
                 np.tile(np.linspace(0, self.End, self.WfSize),(n,1))
            Ts = np.append(Ts, nones, axis=1).reshape((n*(self.WfSize+1),))
            
            # set the plot data to the created arrays
            self.TimeScroll['Plot'].set_data(Ts, Wf)

            # set axes limits
            self.TimeScroll['Ax'].set_xlim(0, h)
            self.TimeScroll['Ax'].draw_artist(self.TimeScroll['Plot'])
            self.SampleWaveform.set_data(self.WaveformXax*n, Wf)
            self.ChanTab['WavesFigure'].figure.axes[0].draw_artist(self.SampleWaveform)
        
        self.TimeScroll['Figure'].figure.canvas.blit(self.TimeScroll['Figure'].figure.bbox)
        self.ChanTab['WavesFigure'].figure.canvas.blit(self.ChanTab['WavesFigure'].figure.axes[0].bbox)

    ########################################################################################################
        
    def VZoom_Func(self):
        v = self.TimeScroll['VZoom'].value()
        self.TimeScroll['Ax'].set_ylim(-v,v)
        self.TimeScroll['Figure'].figure.canvas.restore_region(self.TimeScroll['bg'])
        self.TimeScroll['Ax'].draw_artist(self.TimeScroll['Plot'])
        self.TimeScroll['Figure'].figure.canvas.blit(self.TimeScroll['Figure'].figure.bbox)

    ########################################################################################################
        
    def HZoom_Func(self):
        v = self.TimeScroll['HZoom'].value()
        self.TimeScroll['HScroll'].setSingleStep(v/10)
        self.HScroll_Func()

    ########################################################################################################
        
    def DrawScrollFig_Func(self, event):
        self.TimeScroll['bg'] = self.TimeScroll['Figure'].figure.canvas.copy_from_bbox(self.TimeScroll['Figure'].figure.axes[0].bbox)

    ########################################################################################################
        
    def LoadH5File(self):
        ''' Loads an h5 file that contains all the information about the units:
        waveforms and timestamps '''
        
        # try to load an h5 file
        if settings.WorkingDir: d = settings.WorkingDir
        else: d = ''
        
        f = QtGui.QFileDialog.getOpenFileName(parent = self,
                                              caption = 'Select an H5 File',
                                              directory = d,
                                              filter = '*.h5')

        # in case there is not file selected
        if not f: return

        # set file loaded var = True
        if hasattr(self, 'H5FileLoaded') and self.H5FileLoaded == True:
            self.h5file.close()
            
        # try to open the file
        try:
            self.h5file = tables.openFile(str(f), mode = 'r+')
        except:
            self.MsgBox.setIcon(QtGui.QMessageBox.Warning)
            self.MsgBox.setText('There was a problem opening the H5 file')
            self.MsgBox.setwindowTitle('Warning')
            self.MsgBox.show()
            return

        # set file loaded var = True
        self.H5FileLoaded = True
        self.FilePath = os.path.split(f)[0]
        
        ### REPAIR THE H5FILE STRUCTURE

        if self.h5file.__contains__('/Chans'):
            self.h5file.renameNode('/', 'Spikes', name = 'Chans')

        chanNodes = self.h5file.listNodes('/Spikes')

        for k in chanNodes:
            where = '/Spikes/%s' % k._v_name
            for n in k:                
                if re.search('Unit[0-9]{2}(?!_isMultiunit)', n._v_name) and n._c_classId != 'GROUP':
                    unitName = re.search('Unit[0-9]{2}(?!_isMultiunit)', n._v_name).group()
                    
                    self.h5file.createGroup(where = where, name = unitName+'_grp')
                    self.h5file.moveNode(where = where, name = unitName,
                                    newparent = '/Spikes/%s/%s' % (k._v_name, unitName+'_grp'),
                                    newname = 'Indx')
                    self.h5file.renameNode(where = where, name = unitName+'_grp', newname = unitName)

                elif re.search('Unit[0-9]{2}_isMultiunit', n._v_name):
                    self.h5file.removeNode(where = where, name = re.search('Unit[0-9]{2}_isMultiunit', n._v_name).group())

                elif n._v_name.find('tmp') != -1:
                    self.h5file.removeNode(where = where, name = n._v_name, recursive = True)

        ### CREATE 'isMultiunit' and 'isBursting' fields
        chanNodes = self.h5file.listNodes('/Spikes')

        for k in chanNodes:
            node = '/Spikes/%s' % k._v_name
            for n in k:
                if n._v_name.find('Unit') != -1 and n._c_classId == 'GROUP':
                    parent = node + '/' + n._v_name
                    if not self.h5file.__contains__(parent+'/'+'isMultiunit'):
                        self.h5file.createArray(parent, 'isMultiunit', False)
                    if not self.h5file.__contains__(parent+'/'+'isBursting'):
                        self.h5file.createArray(parent, 'isBursting', False)

        ### RENAME the "Indexes" field to "Indx"
        chanNodes = self.h5file.listNodes('/Spikes')

        for k in chanNodes:
            for n in k:
                if n._v_name.find('Unit') != -1:
                    nodeName = '/Spikes/%s/%s' % (k._v_name, n._v_name)
                    for l in n:
                        if l._v_name == 'Indexes':
                            self.h5file.renameNode(nodeName, 'Indx', 'Indexes')

        # save changes to disk              
        self.h5file.flush()

        ##### REPAIR UNIT NAMES #####
        chanNodes = self.h5file.listNodes('/Spikes')

        for chan in chanNodes:
            unitNames = [k for k in chan.__members__ if k.find('Unit')!=-1]
            unitNames.sort()
            for j,k in enumerate(unitNames):
                if k != 'Unit%02d' % j:
                    self.h5file.renameNode('/Spikes/%s' % chan._v_name , name = k, newname = 'Unit%02d' % j)

        # save changes to disk
        self.h5file.flush()
        ###############################################
        
        # clean the channel figures if something already plotted
        if hasattr(self, 'ChanPlotted') and self.ChanPlotted == True:
            self.ResetChannelTab_proc()
            
        self.PlotOverview()

        # clear the Log Browser and load Log info
        self.LogCombo.clear()
        self.LogTextBrowser.clear()
        if self.h5file.__contains__('/Log'):
            nodes = self.h5file.listNodes('/Log')
            nodeNames = [k._v_name for k in nodes]
            self.LogCombo.addItems(nodeNames)
            
        # set window title = to filename
        self.setWindowTitle('Spike Sorter GUI '+f)

    ########################################################################################################
        
    def PlotOverview(self):
        ''' plot an overview of 1000 spikes per channel.
        Also, fills the overview table with the general information about each
        channel'''
        
        # get the list of nodes inside the "Chans" group
        chanNodes = self.h5file.listNodes('/Spikes')
        
        # get the number of the channels in the file
        self.ChansList = [int(re.search('(?<=Chan_)[0-9]{3}',k._v_name).group()) for k in chanNodes]
                
        # get the waveform size (number of points). X is for fast plotting
        self.WfSize = self.h5file.root.Header.WaveformSize.read()
        x = range(self.WfSize)
        
        # add items to the channel selector in the toolbar
        self.ChanSelector.clear()
        self.ChanSelector.addItems(['%d' % k for k in self.ChansList])

        # clean the overview figure
        self.OverviewTab1['Figure'].figure.clear()
        
        # helper to calculate the geometry of the axes
        figrows = np.ceil(len(chanNodes) / 10.0)

        # clear contents of the overview table
        self.OverviewTab2['OverviewTable'].clearContents()
        c = range(self.OverviewTab2['OverviewTable'].rowCount())
        c.reverse()
        for k in c: self.OverviewTab2['OverviewTable'].removeRow(k)
        
        
        # iterate over the list of channels to add the information to the table
        for j, k in enumerate(chanNodes):

            # update overveiew table
            self.OverviewTab2['OverviewTable'].insertRow(j)
            self.OverviewTab2['OverviewTable'].setRowHeight(j,20)

            # add an event count
            self.OverviewTab2['OverviewTable'].setItem(j,0, QtGui.QTableWidgetItem(str(k.TimeStamp.nrows)))

            # add an "isTrash" checkbox
            check = QtGui.QCheckBox()
            check.setProperty('Data', self.ChansList[j])
            check.stateChanged.connect(self.setTrash_proc)
            self.OverviewTab2['OverviewTable'].setCellWidget(j,1,check)

            # add an "isMultinunit" checkbox
            isMultiunitCheck = QtGui.QCheckBox()
            isMultiunitCheck.setObjectName(k._v_name)
            isMultiunitCheck.stateChanged.connect(self.isMultiunit_proc)
            self.OverviewTab2['OverviewTable'].setCellWidget(j,2,isMultiunitCheck)
            
            # add information about unsorted units
            if k.__contains__('Unsorted'):
                self.OverviewTab2['OverviewTable'].setItem(j, 4, QtGui.QTableWidgetItem(str(k.Unsorted.nrows)))

            # add information about valif waveforms
            if k.__contains__('ValidWFs'):
                self.OverviewTab2['OverviewTable'].setItem(j, 5, QtGui.QTableWidgetItem(str(k.ValidWFs.nrows)))

            # add info about each unit                
            units = [m for m in k.__members__ if re.search('Unit[0-9]{2}', m)] # obtain unit names
            units.sort()
            if units: # in case there are units
                for m,n in enumerate(units):
                    if self.OverviewTab2['OverviewTable'].columnCount() <= (m+6):
                        self.OverviewTab2['OverviewTable'].insertColumn(self.OverviewTab2['OverviewTable'].columnCount())
                        nCols = self.OverviewTab2['OverviewTable'].columnCount()
                        self.OverviewTab2['OverviewTable'].setColumnWidth(nCols-1, 65)
                        self.OverviewTab2['OverviewTable'].setHorizontalHeaderItem(nCols-1,
                                                                                   QtGui.QTableWidgetItem('Unit%02d' % m))
                                                                       
                    self.OverviewTab2['OverviewTable'].setItem(j, m+6,
                                                               QtGui.QTableWidgetItem(str(k.__getattr__(n).Indx.nrows)))

            # Create the axes to plot the waveforms
            self.OverviewTab1['Figure'].figure.add_subplot(figrows, 10, j+1)
            self.OverviewTab1['Figure'].figure.axes[j].set_yticks([],[]) # eliminate the ticks to have more space
            self.OverviewTab1['Figure'].figure.axes[j].set_xticks([],[]) # eliminate the ticks to have more space
            self.OverviewTab1['Figure'].figure.axes[j].set_axis_off()
            self.OverviewTab1['Figure'].figure.axes[j].set_title('Ch %d' % (self.ChansList[j]),
                                                                 fontsize = 10,
                                                                 fontdict={'color':'w'})

            self.PlotChanOverview_proc(k, axes2Plot = self.OverviewTab1['Figure'].figure.axes[j])
            
            # check the isTrash widgets and make the axes background yellow
            if k.__contains__('isTrash'):
                if k.isTrash.read() == True:
                    check.setCheckState(2)
                    self.OverviewTab1['Figure'].figure.axes[j].set_axis_bgcolor([.5, .5, .5])

            if k.__contains__('isMultiunit'):
                if k.isMultiunit.read()==True:
                    isMultiunitCheck.setCheckState(2)
        
        # set the names of the vertical headers
        self.OverviewTab2['OverviewTable'].setVerticalHeaderLabels(['Ch '+str(k) for k in self.ChansList])

        # set alternating row colors
        self.OverviewTab2['OverviewTable'].setAlternatingRowColors(True)
        
        # connect the clicks on this canvas with the channel select function
        self.OverviewTab1['Figure'].figure.canvas.mpl_connect('button_release_event', self.SelChannel)
        
        # tight layout and draw
        self.OverviewTab1['Figure'].figure.tight_layout()
        self.OverviewTab1['Figure'].figure.canvas.draw()
        
        # get the sampling frequency
        self.Sf     = float(self.h5file.root.Header.TimeStamp_Res.read())
        
        self.Step   = self.WfSize+1

        # set boolean variable

    ########################################################################################################

    def PlotChanOverview_proc(self, node, axes2Plot):
        '''Helper function that plots the unsorted as well as the sorted events
        in a given axes on the overview figure'''

        # get the number of events to plot
        nEvents = self.OverviewNEventsSpin.value()
        Waveforms = node.Waveforms.read()
        
        # iterate over the members of a node
        for k in node:
            if not re.search('Unsorted|Unit[0-9]{2}', k._v_name): continue

            # read the indices first:
            if k._v_name.find('Unit') != -1:
                if k.Indx.nrows >= nEvents:
                    indx = k.Indx.read(start = 0, stop = k.Indx.nrows, step = k.Indx.nrows/nEvents)
                else:
                    indx = k.Indx.read()
            elif k._v_name.find('Unsorted') != -1:
                if k.nrows >= nEvents:
                    indx = k.read(start = 0, stop = k.nrows, step = k.nrows/nEvents)
                else:
                    indx = k.read()

            # obtain the waveforms
            Wf = Waveforms[indx,:]
            if not Wf.any(): continue
            
            #### faster plotting strategy:
            
            # obtain the length of units to plot
            n = len(indx)
            
            # create an array of Nones to append
            nones = np.array(n*[None], ndmin=2).transpose()

            # append nones to the waveforms array and reshape it to a vector
            Wf = np.append(Wf, nones, axis=1).reshape((n*(self.WfSize+1),))

            # create the x axis
            Ts = range(self.WfSize)
            Ts.append(None)

            # choose the color and the zorder according to the type of unit
            if k._v_name == 'Unsorted':
                color = 'w'
                zorder = 1
            elif k._v_name.find('Unit')!=-1:
                # get the unit number
                zorder = int(re.search('[0-9]{2}', k._v_name).group())
                color = self.UnitColors[zorder,:]
                zorder = 100-zorder

            # get the list of plots in the particular axes
            l = [l for l in axes2Plot.lines if str(l.get_label()) == k._v_name]

            # if a plot with a label equal to the name of the unit exist, the update the data
            if len(l)>0:
                l[0].set_data(Ts*n, Wf)
            # if not create one
            else:
                axes2Plot.plot(Ts*n, Wf, color = color, rasterized =  True,
                               alpha = 0.5, label = k._v_name, zorder = zorder)

        # set the limits of the axes
        axes2Plot.set_ylim(-self.OverviewYLimsSpin.value(), self.OverviewYLimsSpin.value())

        # add a small text box with the event count
        bbox_props = dict(boxstyle='round', fc='0.75', ec='0.25', alpha=0.8)
        axes2Plot.text(0.5, 0, 'Count: %d' % node.TimeStamp.nrows, transform = axes2Plot.transAxes,
                       color = 'k', bbox = bbox_props, size = 10, ha = 'center')

    ##################################################################################################

    def ChangeOverviewYLim_Proc(self):

        if not self.H5FileLoaded: return
        
        lim = self.OverviewYLimsSpin.value()
        for k in self.OverviewTab1['Figure'].figure.axes:
            k.set_ylim(-lim, lim)

        self.OverviewTab1['Figure'].figure.canvas.draw()
        
    ##################################################################################################
            
    def PlotChanProc(self):

        # exit if ther is no H5 file loaded
        if not self.H5FileLoaded: return
        
        # clean the channels tab
        self.ResetChannelTab_proc()

        # reset Units list
        self.UnitsList = []
        
        #pdb.set_trace()
        
        # load waveforms for a specific channel
        self.CurChan = int(self.ChanSelector.currentText())
        nspikes = self.NSpikesSlider.value()
        self.CurNodeName = '/Spikes/Chan_%03d' % self.CurChan
        self.CurNode     = self.h5file.getNode(self.CurNodeName)
        self.CurWaveforms = self.CurNode.Waveforms.read()
        self.CurTs        = self.CurNode.TimeStamp.read()
        self.TimeScroll['HScroll'].setMaximum(int(self.CurTs[-1]))
        self.TimeScroll['HSpin'].setMaximum(int(self.CurTs[-1]))
        self.unitNodes = [k for k in self.h5file.listNodes(self.CurNodeName) if re.search('Unit[0-9]{2}',k._v_name)]

        # get the indices of the unsorted. If there are no, create one
        if not self.CurNode.__contains__('Unsorted'):
            self.Unsorted = np.arange(len(self.CurTs))
            self.h5file.createArray(self.CurNodeName,'Unsorted', self.Unsorted)
        else:
            self.Unsorted = self.h5file.getNode(self.CurNodeName,'Unsorted').read()

        #set the unit names in the combo box
        self.What2Plot.clear()
        self.What2Plot.addItems(['All Waveforms', 'Sorted', 'Unsorted'])
        if self.unitNodes:
            self.What2Plot.addItems([k._v_name for k in self.unitNodes])

        # set the axis limits to apropriately view the unit
        v = self.WaveAxYLim_Spin.value()
        self.ChanTab['WavesFigure'].figure.axes[0].set_ylim(-v, v)
            
        # get the waveform size for this specific waveform
        self.WfSize = self.h5file.root.Header.WaveformSize.read()

        # cheack whether to plot raw waveforms
        '''
        if self.PlotWaveformsCheck.checkState() == 2:
            if self.CurWaveforms.shape[1] < 10000:
                self.Waveforms2Plot = self.CurWaveforms
            else:
                indx = np.int32(np.linspace(0,self.CurWaveforms.shape[0]-1,10000))
                self.Waveforms2Plot = self.CurWaveforms[indx,:]

            for k in range(self.Waveforms2Plot.shape[0]):
                self.Waveforms.set_data(range(self.WfSize), self.Waveforms2Plot[k,:])
                self.ChanTab['WavesFigure'].figure.axes[0].draw_artist(self.Waveforms)
                self.ChanTab['WavesFigure'].figure.canvas.blit(self.ChanTab['WavesFigure'].figure.axes[0].bbox)'''

        # grab background from the Waveforms Figure to make animations
        self.ChanTab['WavesFigBG'] = self.ChanTab['WavesFigure'].figure.canvas.copy_from_bbox(self.ChanTab['WavesFigure'].figure.axes[0].bbox)
        self.MainFigTab.setTabText(2,'Chan %02d' % self.CurChan)

        # calculate PCA
        pc = PCA(self.CurWaveforms)

        # put data in a KDTree to easily calculate distance with the cursor
        self.XYData = cKDTree(pc.Y[:,0:2],1000)
        self.ChanTab['PCA'] = pc.Y

        # set the internal variable to true
        self.ChanPlotted = True

        # copy the number of events in the channel into a label to see if the user
        # wants to decimate or plot them all
        self.nPtsLabel.setText(str(self.CurTs.size))

        # read the plotting parameters in the "Chan Tools" tab to plot
        # the selected feature
        self.PlotFeatures()

        if self.ChanTab['UnitTabsWidget'].count()>0:
            self.ChanTab['UnitTabsWidget'].setCurrentIndex(0)
            self.CurUnitName = str(self.ChanTab['UnitTabsWidget'].tabText(0))
            self.CurUnit = int(re.search('(?<=Unit)[0-9]{2}',self.CurUnitName).group())

        self.WaveformXax = range(self.WfSize)
        self.WaveformXax.append(None)
        self.End = 1000*self.WfSize/self.Sf

        # save h5file changes to disk
        self.h5file.flush()

    ########################################################################################################
        
    def setTrash_proc(self):
        sender = self.sender()
        chan = sender.property('Data')
        nodeName = '/Spikes/Chan_%03d' % chan

        indx = self.ChansList.index(chan)
        
        if self.h5file.getNode(nodeName).__contains__('isTrash'):
            self.h5file.removeNode(nodeName, 'isTrash')
                
        if sender.checkState() in [1,2]:
            self.h5file.createArray(nodeName, 'isTrash', True)
            self.OverviewTab1['Figure'].figure.axes[indx].set_axis_bgcolor('y')
        elif sender.checkState() == 0:
            self.h5file.createArray(nodeName, 'isTrash', False)
            self.OverviewTab1['Figure'].figure.axes[indx].set_axis_bgcolor('w')

        # save changes to disk
        self.h5file.flush()

    ########################################################################################################
        
    def isMultiunit_proc(self):
        sender = self.sender()
        nodeName = '/Spikes/%s' % sender.objectName()

        if self.h5file.getNode(nodeName).__contains__('isMultiunit'):
            self.h5file.removeNode(nodeName, 'isMultiunit')

        if sender.checkState() in [1,2]:
            self.h5file.createArray(nodeName, 'isMultiunit', True)
        elif sender.checkState() == 0:
            self.h5file.createArray(nodeName, 'isMultiunit', False)
            
        # save changes to disk
        self.h5file.flush()

    ########################################################################################################
        
    def TrashChans_proc(self):
        '''Utility function to mark the channels with fewer than a defined number
        of units'''
        
        # check whether an h5file has been loaded
        if not self.H5FileLoaded: return

        # obtain parameters
        n     = self.MarkTrashSpin.value()
        chans = self.h5file.listNodes('/Spikes')

        # iterate over nodes in h5file; if fewer than n mark as trash
        for l,k in enumerate(chans):
            j = int(re.search('(?<=Chan_)[0-9]{3}',k._v_name).group())
            if k.TimeStamp.nrows < n:
                self.OverviewTab1['Figure'].figure.axes[l].set_axis_bgcolor('y')
                self.OverviewTab2['OverviewTable'].cellWidget(l, 1).setChecked(True)
                if self.h5file.getNode('/Spikes','Chan_%03d' % j).__contains__('isTrash'):
                    self.h5file.removeNode('/Spikes/Chan_%03d' % j, 'isTrash')
                self.h5file.createArray('/Spikes/Chan_%03d' % j, 'isTrash', True)
                
            else:
                self.OverviewTab1['Figure'].figure.axes[l].set_axis_bgcolor('w')
                self.OverviewTab2['OverviewTable'].cellWidget(l, 1).setChecked(False)
                if self.h5file.getNode('/Spikes','Chan_%03d' % j).__contains__('isTrash'):
                    self.h5file.removeNode('/Spikes/Chan_%03d' % j, 'isTrash')
                self.h5file.createArray('/Spikes/Chan_%03d' % j, 'isTrash', False)

        # save changes to disk
        self.h5file.flush()

        #update the overview
        self.OverviewTab1['Figure'].figure.canvas.draw()

    ########################################################################################################        

    def DeleteTrashChans_proc(self):
        # check whether an h5file has been loaded
        if not self.H5FileLoaded: return
        
        chans = self.h5file.listNodes('/Spikes')
        chans.reverse()
        n = range(len(chans))
        n.reverse()

        delchans = []
        for j,k in zip(n, chans):
            state = self.OverviewTab2['OverviewTable'].cellWidget(j,1).checkState()
            
            if state == 2:
                delchans.append(k._v_name)
                self.h5file.removeNode('/Spikes', k._v_name, recursive=True)
                    #self.OverviewTab1['Figure'].figure.delaxes(self.OverviewTab1['Figure'].figure.axes[j])

        if len(delchans)>0:
            self.AddLog('Deleted channels: '+str(delchans))

        self.PlotOverview()

    ########################################################################################################
        
    def ResetChan_proc(self):
        ''' check whether an h5file has been loaded '''
        if not self.H5FileLoaded or not self.ChanPlotted: return
        
        for k in self.h5file.listNodes(self.CurNodeName):
            if k._v_name not in ['Waveforms','TimeStamp','isTrash']:
                self.h5file.removeNode(self.CurNodeName, k._v_name, recursive = True)

        self.PlotChanProc()
        self.AddLog('%s resetted' % self.CurNodeName)

    ########################################################################################################
        
    def AddLog(self, message):
        ''' add log to to keep a history of changes to the file '''

        if not self.H5FileLoaded: return
        
        if not self.h5file.__contains__('/Log'):
            self.h5file.createGroup('/', 'Log', 'History of changes')
        name = 'Entry_%s_%s_%s_%s_%s_%s' % datetime.datetime.now().timetuple()[0:6]
        self.h5file.createArray('/Log', name, message)

        # save changes to disk
        self.h5file.flush()

        #add the item to the log browser
        self.LogCombo.addItem(name)
    
    ########################################################################################################
        
    def SetLogText_proc(self):
        if self.LogCombo.currentIndex == -1: return
        node = str(self.LogCombo.currentText())
        if node:
            log = self.h5file.getNode('/Log',node).read()
            self.LogTextBrowser.setText(log)

    ########################################################################################################

    def CloseFile(self):
        ''' close the h5 file'''
        if self.H5FileLoaded == False:  return
        self.h5file.flush()
        self.h5file.close()
        self.H5FileLoaded = False
        print 'h5 File closed'

    ########################################################################################################
        
    def SelChannel(self, event):
        ''' selects a channel when axes are clicked'''
        if event.inaxes:
            chan = int(re.search('(?<=Ch )[0-9]{1,3}',event.inaxes.get_title()).group())
            c = [int(self.ChanSelector.itemText(k)) for k in range(self.ChanSelector.count())].index(chan)
            self.ChanSelector.setCurrentIndex(c)

    ########################################################################################################
            
    def TableRowChanged_proc(self, sel):
        self.ChanSelector.setCurrentIndex(sel)

    ########################################################################################################
        
    def Settings(self):
        ''' edit paths'''
        if settings.edit() == 1:
            self.WorkingDir = settings.WorkingDir

    ########################################################################################################
            
    def AutocorrOpts(self):
        if autocorropts.edit() == 1:
            pass

    ########################################################################################################
        
    def About(self):
        ''' opens a small dialog with information about the software'''
        self.AboutMsg.show()

    ########################################################################################################
        
    def NearestPoint(self, event):
        ''' when right button clicked over the features window, calculates the closest
        point and plots its corresponding waveform'''
        
        if event.button==3 and event.inaxes and self.ChanTab['FeaturesFigNtb'].mode=='':
            featuresax = self.ChanTab['FeaturesFig'].figure.axes[0]
            wavesax    = self.ChanTab['WavesFigure'].figure.axes[0]

            if self.PlotUnitCounter >= self.NUnitsSpin.value():
                self.ChanTab['WavesFigure'].figure.canvas.restore_region(self.ChanTab['WavesFigBG'])
                self.PlotUnitCounter = 0

            for k in self.ChanTab['FeaturesFigBG']:
                self.ChanTab['FeaturesFig'].figure.canvas.restore_region(k)

            _,res = self.XYData.query([event.xdata, event.ydata], 1)
            self.cursor.set_data(self.XYData.data[res,0], self.XYData.data[res,1])
            self.SampleWaveform.set_data(range(self.WfSize), self.CurWaveforms[self.dataIndx[res],: ])

            featuresax.draw_artist(self.cursor)
            wavesax.draw_artist(self.SampleWaveform)
            self.ChanTab['FeaturesFig'].figure.canvas.blit(featuresax.bbox)
            self.ChanTab['WavesFigure'].figure.canvas.blit(wavesax.bbox)
            self.PlotUnitCounter+=1

    ########################################################################################################
            
    def draw_callback(self, event):
        ''' any draw callback triggers the capture of the figure background for
        using it in the animations '''
        
        if not self.ChanPlotted: return

        if event.canvas == self.ChanTab['FeaturesFig'].figure.canvas:
            bg = []
            for k in self.ChanTab['FeaturesFig'].figure.axes:
                bg.append(self.ChanTab['FeaturesFig'].figure.canvas.copy_from_bbox(k.bbox))
            self.ChanTab['FeaturesFigBG'] = bg
        elif event.canvas == self.ChanTab['WavesFigure'].figure.canvas:
            self.ChanTab['WavesFigBG'] = self.ChanTab['WavesFigure'].figure.canvas.copy_from_bbox(self.ChanTab['WavesFigure'].figure.axes[0].bbox)

    ########################################################################################################
            
    def PlotFeatures(self):
        ''' determines what 2 plot based on the user choices'''

        # obtain labels and return if are the same
        xlabel = self.XPlot.currentText()
        ylabel = self.YPlot.currentText()
        if xlabel == ylabel: return
        
        curchan = int(self.ChanSelector.currentText())
        
        if self.PlotValidsOnlyCheck.checkState()==2 and \
           self.CurNode.__contains__('ValidWFs'):
            print 'you selected to plot only the valid WFs'

        What2Plot = str(self.What2Plot.currentText()) # string containing what to plot
        self.CurNodeName = '/Spikes/Chan_%03d' % curchan
        unitNodes = [k for k in self.h5file.listNodes(self.CurNodeName) if re.search('Unit[0-9]{2}', k._v_name)]
        

        if What2Plot in ['All Waveforms', 'Sorted']:
            self.dataIndx = range(self.CurTs.size)
            pc  = self.ChanTab['PCA']
        elif What2Plot == 'Unsorted':
            self.dataIndx = self.h5file.getNode(self.CurNodeName,'Unsorted').read()
            pc = PCA(self.CurWaveforms[self.dataIndx,:])
            pc = pc.Y
        elif re.search('Unit', What2Plot):
            self.CurUnitName = What2Plot
            self.dataIndx = self.h5file.getNode(self.CurNodeName, What2Plot).Indx.read()
            pc = PCA(self.CurWaveforms[self.dataIndx,:])
            pc = pc.Y

        # save what is the feature
        self.CurFeaturePlot = What2Plot

        # get the choice for the x axis
        if xlabel == 'PCA1':
            x = pc[:,0]

        elif xlabel == 'PCA2':
            x = pc[:,1]

        elif xlabel == 'PCA3':
            x = pc[:,2]

        elif xlabel == 'Slice1':
            x = self.CurWaveforms[self.dataIndx, self.SliceSpBx1.value()]

        elif xlabel == 'Slice2':
            x = self.CurWaveforms[self.dataIndx, self.SliceSpBx2.value()]

        elif xlabel == 'Energy':
            x = np.sum(np.power(self.CurWaveforms[self.dataIndx,:], 2), axis=1)
            x = x/1000000.0

        elif xlabel == 'Peak':
            x = self.CurWaveforms[self.dataIndx,:].max(axis=1)
            x = x/100.0

        elif xlabel == 'Valley':
            x = self.CurWaveforms[self.dataIndx,:].min(axis=1)
            x = x/100.0

        elif xlabel == 'Peak Pt':
            x = self.CurWaveforms[self.dataIndx,:].argmax(axis=1)

        elif xlabel == 'Valley Pt':
            x = self.CurWaveforms[self.dataIndx,:].argmin(axis=1)
            
        elif xlabel == 'Pk2Pk Amp':
            x = self.CurWaveforms[self.dataIndx,:].max(axis=1)-self.CurWaveforms[self.dataIndx,:].min(axis=1)
            x = x/100.0

        elif xlabel == 'Time':
            x = self.CurTs[self.dataIndx]
            x = x/60000.0


        # get the choice for the y axis
        if ylabel == 'PCA1':
            y = pc[:,0]

        elif ylabel == 'PCA2':
            y = pc[:,1]

        elif ylabel == 'PCA3':
            y = pc[:,2]

        elif ylabel == 'Slice1':
            y = self.CurWaveforms[self.dataIndx, self.SliceSpBx1.value()]

        elif ylabel == 'Slice2':
            y = self.CurWaveforms[self.dataIndx, self.SliceSpBx2.value()]

        elif ylabel == 'Energy':
            y = np.sum(np.power(self.CurWaveforms[self.dataIndx,:], 2), axis=1)
            y = y/1000000.0

        elif ylabel == 'Peak':
            y = self.CurWaveforms[self.dataIndx,:].max(axis=1)
            y = y/100.0

        elif ylabel == 'Valley':
            y = self.CurWaveforms[self.dataIndx,:].min(axis=1)
            y = y/100.0

        elif ylabel == 'Peak Pt':
            y = self.CurWaveforms[self.dataIndx,:].argmax(axis=1)

        elif ylabel == 'Valley Pt':
            y = self.CurWaveforms[self.dataIndx,:].argmin(axis=1)

        elif ylabel == 'Pk2Pk Amp':
            y = self.CurWaveforms[self.dataIndx,:].max(axis=1)-self.CurWaveforms[self.dataIndx,:].min(axis=1)
            y = y/100.0

        elif ylabel == 'Time':
            y = self.CurTs[self.dataIndx]
            y = y/60000.0

        naxes = len(self.ChanTab['FeaturesFig'].figure.axes)

        nspikes = self.NSpikesSlider.value()

        title = '%s: %s vs %s' % (What2Plot, xlabel, ylabel)
        # obtain the axis limits if we are plotting the same variables
        same_limits = False
        if naxes > 0 and \
           self.ChanTab['FeaturesFig'].figure.axes[0].get_title() == title:
            same_limits = True
            xlim = self.ChanTab['FeaturesFig'].figure.axes[0].get_xlim()
            ylim = self.ChanTab['FeaturesFig'].figure.axes[0].get_ylim()
        
        # plot only on one axes
        if self.PlotDensityCheck.checkState() == 0:
            if naxes == 0:
                ax1 = self.ChanTab['FeaturesFig'].figure.add_subplot(111)
                ax1.set_axis_bgcolor('k')
            elif naxes == 1:
                ax1 = self.ChanTab['FeaturesFig'].figure.axes[0]
                ax1.cla()
                ax1.set_axis_bgcolor('k')
            elif naxes >= 2:
                self.ChanTab['FeaturesFig'].figure.clear()
                ax1 = self.ChanTab['FeaturesFig'].figure.add_subplot(111)
                ax1.set_axis_bgcolor('k')
                
        # create 2 subplots to host the density
        elif self.PlotDensityCheck.checkState() == 2:
            if naxes == 0:
                ax1 = self.ChanTab['FeaturesFig'].figure.add_subplot(121)
                ax2 = self.ChanTab['FeaturesFig'].figure.add_subplot(122, sharex=ax1, sharey=ax1)
            elif naxes == 1:
                self.ChanTab['FeaturesFig'].figure.clear()
                ax1 = self.ChanTab['FeaturesFig'].figure.add_subplot(121)
                ax2 = self.ChanTab['FeaturesFig'].figure.add_subplot(122, sharex=ax1, sharey=ax1)
            elif naxes == 2:
                ax1 = self.ChanTab['FeaturesFig'].figure.axes[0]
                ax2 = self.ChanTab['FeaturesFig'].figure.axes[1]
                ax1.cla()
                ax2.cla()
                
            ax2.set_axis_bgcolor('k')
            # create and plot a 2d histogram
            
            
        # setup the axes
        ax1.set_title(title, fontdict={'color':'w'})
        ax1.tick_params(color = [.5, .5, .5])
        for k in ax1.spines.values(): k.set_edgecolor([.5, .5, .5])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_axis_bgcolor('k')
        
        self.cursor, = ax1.plot([], 's', mfc='none', ms=6, mec='r',
                                animated=True, label='sample')

        # iterate over the members of that channel
        if What2Plot == 'All Waveforms':
            nodes = self.h5file.listNodes(self.CurNodeName)
            for leaf in nodes:
                if leaf._v_name == 'Unsorted':
                    # select only some indices to plot
                    if  leaf.nrows > self.nPtsSpin.value():
                        indx = leaf.read(0, leaf.nrows, leaf.nrows/self.nPtsSpin.value())
                    else:
                        indx = leaf.read()

                    # plot unsorted
                    ax1.plot(x[indx,:], y[indx,:], ',',
                             color = [.5, .5, .5], label = 'data_Unsorted')
                    
                unit = re.search('(?<=Unit)[0-9]{2}', leaf._v_name)
                if unit:
                    # select some units to plot
                    if leaf.Indx.nrows > self.nPtsSpin.value():
                        indx = leaf.Indx.read(0, leaf.Indx.nrows, leaf.Indx.nrows/self.nPtsSpin.value())
                    else:
                        indx = leaf.Indx.read()
                        
                    ax1.plot(x[indx,:], y[indx,:], ',', label = 'data_'+leaf._v_name,
                             rasterized = True,
                             color = self.UnitColors[int(unit.group()),:],
                             mec = self.UnitColors[int(unit.group()),:])

                    # add unit to the tab widget
                    self.UnitsTable_AddUnit(leaf._v_name)

        elif What2Plot == 'Sorted':
            nodes = self.h5file.listNodes(self.CurNodeName)
            for leaf in nodes:
                unit = re.search('(?<=Unit)[0-9]{2}', leaf._v_name)
                if unit:
                    # select some units to plot
                    if leaf.Indx.nrows > self.nPtsSpin.value():
                        indx = leaf.Indx.read(0, leaf.Indx.nrows, leaf.Indx.nrows/self.nPtsSpin.value())
                    else:
                        indx = leaf.Indx.read()
                        
                    ax1.plot(x[indx,:], y[indx,:], ',', label = 'data_'+leaf._v_name,
                             rasterized = True,
                             color = self.UnitColors[int(unit.group()),:],
                             mec = self.UnitColors[int(unit.group()),:])

                    # add unit to the tab widget
                    self.UnitsTable_AddUnit(leaf._v_name)

        # to plot the unsorted channels
        elif What2Plot == 'Unsorted':
            lx = len(x)
            # select some units to plot
            if lx > self.nPtsSpin.value():
                indx = range(0, lx, lx/self.nPtsSpin.value())
            else:
                indx = range(lx)
            ax1.plot(x[indx], y[indx], ',', color = [.5, .5, .5]
                     , label = 'data_Unsorted',
                     rasterized = True)

        # plot a specific unit
        elif re.search('Unit', What2Plot):
            unit = re.search('(?<=Unit)[0-9]{0,2}', What2Plot).group()
            lx = len(x)
            # select some units to plot
            if lx > self.nPtsSpin.value():
                indx = range(0, lx, lx/self.nPtsSpin.value())
            else:
                indx = range(lx)
                
            ax1.plot(x[indx], y[indx], ',', label = 'data_'+What2Plot,
                     rasterized = True,
                     color = self.UnitColors[int(unit),:],
                     mec = self.UnitColors[int(unit),:])

            # add unit to the tab widget
            self.UnitsTable_AddUnit(What2Plot)

        if same_limits == True:
            ax1.set_ylim(ylim)
            ax1.set_xlim(xlim)
        else:
            ax1.relim()
            ax1.autoscale_view(True,True,True)

        # vertical and horizontal lines @ x and y = 0
        ax1.axvline(0, color=[.5, .5, .5])
        ax1.axhline(0, color=[.5, .5, .5])
        
        # create KDTree objet from the selected data for fast search 
        self.XYData = cKDTree(np.array([x,y]).transpose())

        # connect figure to the motion notify function
        if not ax1.callbacks.callbacks or not hasattr(self, 'axZoomCID'):
            self.axZoomCID = ax1.callbacks.connect('ylim_changed', self.AxesZoom_proc)

        # connect figure to the motion notify function
        if not hasattr(self, 'motionCID'):
            self.motionCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('motion_notify_event',
                                                                                   self.NearestPoint)

        # connect figure to the draw function
        if not hasattr(self, 'drawCID'):
            self.drawCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('draw_event',
                                                                                 self.draw_callback)

        # plot density if checked
        if self.PlotDensityCheck.checkState() == 2:
            self.ReplotDensity_proc()
        
        # set tight layout and redraw figure
        self.ChanTab['FeaturesFig'].figure.tight_layout()
        self.ChanTab['FeaturesFig'].figure.canvas.draw()

    ########################################################################################################
        
    def Plot3DFeatures(self):

        # obtain labels and return if are the same
        xlabel = self.XPlot.currentText()
        ylabel = self.YPlot.currentText()
        zlabel = self.ZPlot.currentText()
        if xlabel==ylabel or xlabel==zlabel or ylabel==zlabel: return
        
        curchan = int(self.ChanSelector.currentText())
        
        if self.PlotValidsOnlyCheck.checkState()==2 and \
           self.CurNode.__contains__('ValidWFs'):
            print 'you selected to plot only the valid WFs'

        What2Plot = str(self.What2Plot.currentText()) # string containing what to plot
        self.CurNodeName = '/Spikes/Chan_%03d' % curchan
                

        if What2Plot == 'All Waveforms':
            self.dataIndx = range(self.CurTs.size)
            pc  = self.ChanTab['PCA']
        elif What2Plot == 'Unsorted':
            self.dataIndx = self.h5file.getNode(self.CurNodeName,'Unsorted').read()
            pc = PCA(self.CurWaveforms[self.dataIndx,:])
            pc = pc.Y
        elif re.search('Unit',What2Plot):
            self.CurUnitName = What2Plot
            self.dataIndx = self.h5file.getNode(self.CurNodeName, What2Plot).Indx.read()
            pc = PCA(self.CurWaveforms[self.dataIndx,:])
            pc = pc.Y
        elif What2Plot == 'Sorted':
            return

        # save what is the feature
        self.CurFeaturePlot = What2Plot

        # get the choice for the x axis
        if xlabel == 'PCA1':
            x = pc[:,0]

        elif xlabel == 'PCA2':
            x = pc[:,1]

        elif xlabel == 'PCA3':
            x = pc[:,2]

        elif xlabel == 'Slice1':
            x = self.CurWaveforms[self.dataIndx, self.SliceSpBx1.value()]

        elif xlabel == 'Slice2':
            x = self.CurWaveforms[self.dataIndx, self.SliceSpBx2.value()]

        elif xlabel == 'Energy':
            x = np.sum(np.power(self.CurWaveforms[self.dataIndx,:], 2), axis=1)
            x = x/1000000.0

        elif xlabel == 'Peak':
            x = self.CurWaveforms[self.dataIndx,:].max(axis=1)
            x = x/100.0

        elif xlabel == 'Valley':
            x = self.CurWaveforms[self.dataIndx,:].min(axis=1)
            x = x/100.0

        elif xlabel == 'Pk2Pk Amp':
            x = self.CurWaveforms[self.dataIndx,:].max(axis=1)-self.CurWaveforms[self.dataIndx,:].min(axis=1)
            x = x/100.0

        elif xlabel == 'Time':
            x = self.CurTs[self.dataIndx]
            x = x/60000.0


        # get the choice for the y axis
        if ylabel == 'PCA1':
            y = pc[:,0]

        elif ylabel == 'PCA2':
            y = pc[:,1]

        elif ylabel == 'PCA3':
            y = pc[:,2]

        elif ylabel == 'Slice1':
            y = self.CurWaveforms[self.dataIndx, self.SliceSpBx1.value()]

        elif ylabel == 'Slice2':
            y = self.CurWaveforms[self.dataIndx, self.SliceSpBx2.value()]

        elif ylabel == 'Energy':
            y = np.sum(np.power(self.CurWaveforms[self.dataIndx,:], 2), axis=1)
            y = y/1000000.0

        elif ylabel == 'Peak':
            y = self.CurWaveforms[self.dataIndx,:].max(axis=1)
            y = y/100.0

        elif ylabel == 'Valley':
            y = self.CurWaveforms[self.dataIndx,:].min(axis=1)
            y = y/100.0

        elif ylabel == 'Pk2Pk Amp':
            y = self.CurWaveforms[self.dataIndx,:].max(axis=1)-self.CurWaveforms[self.dataIndx,:].min(axis=1)
            y = y/100.0

        elif ylabel == 'Time':
            y = self.CurTs[self.dataIndx]
            y = y/60000.0

        # get the choice for the z axis
        if zlabel == 'PCA1':
            z = pc[:,0]

        elif zlabel == 'PCA2':
            z = pc[:,1]

        elif zlabel == 'PCA3':
            z = pc[:,2]

        elif zlabel == 'Slice1':
            z = self.CurWaveforms[self.dataIndx, self.SliceSpBx1.value()]

        elif zlabel == 'Slice2':
            z = self.CurWaveforms[self.dataIndx, self.SliceSpBx2.value()]

        elif zlabel == 'Energy':
            z = np.sum(np.power(self.CurWaveforms[self.dataIndx,:], 2), axis=1)
            z = y/1000000.0

        elif zlabel == 'Peak':
            z = self.CurWaveforms[self.dataIndx,:].max(axis=1)
            z = y/100.0

        elif zlabel == 'Valley':
            z = self.CurWaveforms[self.dataIndx,:].min(axis=1)
            z = y/100.0

        elif zlabel == 'Pk2Pk Amp':
            z = self.CurWaveforms[self.dataIndx,:].max(axis=1)-self.CurWaveforms[self.dataIndx,:].min(axis=1)
            z = y/100.0

        elif zlabel == 'Time':
            z = self.CurTs[self.dataIndx]
            z = y/60000.0

        if What2Plot == 'All Waveforms' and self.CurNode.__contains__('ValidWFs'):
            valid = self.CurNode.ValidWFs.read()
            x = x[valid]
            y = y[valid]
            z = z[valid]


        if zlabel != 'Density':
            self.Fig3d.clf()
            self.Fig3d.points3d(x,y,z, mode = 'point')
            self.Fig3d.axes()
        else:
            # obtain axes and first axes limits
            ax1 = self.ChanTab['FeaturesFig'].figure.axes[0]
            xlim = ax1.get_xlim()
            ylim = ax1.get_ylim()
            
            # search for the unsorted or the units plots to obatin data
            xpts = []; ypts = []
            for k in ax1.get_children():
                if re.search('Unsorted|Unit', str(k.get_label())):
                    data = k.get_data()
                    xpts.extend(data[0])
                    ypts.extend(data[1])
            xypoints = np.array([xpts,ypts]).transpose()

            # check wich points are inside the axes
            verts = ax1.viewLim.corners()
            verts[2] = ax1.viewLim.corners()[3]
            verts[3] = ax1.viewLim.corners()[2]

            
            inpoly = Path(verts).contains_points(xypoints)

            # create a 2d histogram of the data and scale it logaritmically
            h,xd,yd = np.histogram2d(xypoints[inpoly,0], xypoints[inpoly,1],
                                     bins = self.PlotDensityBins.value(), normed = False)
            h[h<=0] = 1
            h = np.log10(h)
            self.Fig3d.clf()
            self.Fig3d.barchart(h*10)
            self.Fig3d.axes()
        
    ########################################################################################################
        
    def ValidateWFs_proc(self):
        ''' obtains the coordinates of the current feature axis, and uses it to
        determine wich points lay inside it.
        It also saves the data to the h5file'''

        # exits if there is no h5 file loaded or channel plotted
        if not self.H5FileLoaded or not self.ChanPlotted: return

        # get axes handle and limits
        ax = self.ChanTab['FeaturesFig'].figure.axes[0]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # obtain coordinates of the current axes and uses that to build a polygon
        xyverts = [[xlim[0], ylim[0]], [xlim[0], ylim[1]], [xlim[1], ylim[1]], [xlim[1], ylim[0]]]

        # obtain the indices of the waveforms inside the polygon
        p = Path(xyverts).contains_points(self.XYData.data)

        # in case no points were inside the axes
        if len(p)==0:
            self.MsgBox.setIcon(QtGui.QMessageBox.Warning)
            self.MsgBox.setText('There were no selected points')
            self.MsgBox.setwindowTitle('Warning')
            self.MsgBox.show()
            return
        
        self.ValidWFs   = np.flatnonzero(p)
        self.InvalidWFs = np.flatnonzero(~p)

        # remove the 'ValidWFs' field if it already exists
        if self.h5file.getNode(self.CurNodeName).__contains__('ValidWFs'):
            self.h5file.removeNode(self.CurNodeName, 'ValidWFs')

        # remove the 'InvalidWFs' field if it already exists
        if self.h5file.getNode(self.CurNodeName).__contains__('InvalidWFs'):
            self.h5file.removeNode(self.CurNodeName, 'InvalidWFs')
            
        # save the ValidWFs indices to the h5file
        self.h5file.createArray(self.CurNodeName, 'ValidWFs', self.ValidWFs)

        # save the InvalidWFs indices to the h5file
        self.h5file.createArray(self.CurNodeName, 'InvalidWFs', self.InvalidWFs)

        # save changes to disk
        self.h5file.flush()

        # update the information on the overview table
        row = self.ChanSelector.currentIndex()
        item = QtGui.QTableWidgetItem(str(self.ValidWFs.size))
        self.OverviewTab2['OverviewTable'].takeItem(row, 5)
        self.OverviewTab2['OverviewTable'].setItem(row, 5, item)

    ########################################################################################################
        
    def ReplotDensity_proc(self):
        ''' replot density using all the resolution only to the visible points'''

        # check whether the number of axes in the figure
        if len(self.ChanTab['FeaturesFig'].figure.axes)!=2: return

        # obtain axes and first axes limits
        ax1 = self.ChanTab['FeaturesFig'].figure.axes[0]
        ax2 = self.ChanTab['FeaturesFig'].figure.axes[1]
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()
        
        # search for the unsorted or the units plots to obatin data
        xpts = []; ypts = []
        for k in ax1.get_children():
            if re.search('Unsorted|Unit', str(k.get_label())):
                data = k.get_data()
                xpts.extend(data[0])
                ypts.extend(data[1])
        xypoints = np.array([xpts,ypts]).transpose()

        # check wich points are inside the axes
        verts = ax1.viewLim.corners()
        verts[2] = ax1.viewLim.corners()[3]
        verts[3] = ax1.viewLim.corners()[2]
        
        inpoly = Path(verts).contains_points(xypoints)

        # create a 2d histogram of the data and scale it logaritmically
        h,xd,yd = np.histogram2d(xypoints[inpoly,0], xypoints[inpoly,1],
                                 bins = self.PlotDensityBins.value(), normed = False)
        h[h<=0] = 1
        h = np.log10(h)
        
        # clean axes No2 and plot the 2d histogram
        ax2.cla()
        ax2.pcolor(xd, yd, h.transpose(), cmap =  helper_widgets.colormaps[settings.DensityCM])

        # set axis limits
        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim)

        # remove tick labels
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        
        # create vertical and horizontal lines at 0
        ax2.axvline(0,color=[.5, .5, .5])
        ax2.axhline(0,color=[.5, .5, .5])

        # redraw the figure
        self.ChanTab['FeaturesFig'].figure.canvas.draw()

    ########################################################################################################
        
    def AxesZoom_proc(self, ax):

        xpts = []; ypts = []
        for k in ax.get_children():
            if re.search('Unsorted|Unit', str(k.get_label())):
                data = k.get_data()
                xpts.extend(data[0])
                ypts.extend(data[1])
        xypoints = np.array([xpts,ypts]).transpose()
        # check wich points are inside the axes
        verts = ax.viewLim.corners()
        verts[2] = ax.viewLim.corners()[3]
        verts[3] = ax.viewLim.corners()[2]
        
##        inpoly = points_inside_poly(xypoints, verts)
##        w = self.CurWaveforms[inpoly,:]
##        self.ChanTab['WavesFigure'].figure.axes[0].set_ylim(w.min(), w.max())
##        self.ChanTab['WavesFigure'].figure.canvas.draw()

        if len(self.ChanTab['FeaturesFig'].figure.axes) == 2:
            self.ReplotDensity_proc()

    ########################################################################################################
            
    def AutoClust_proc(self):
        if not self.H5FileLoaded or not self.ChanPlotted: return

        if self.XYData.data.shape[1] > 2:
            data= self.XYData.data[:,0:2]
        else:
            data = self.XYData.data
            
        clustIndx = KlustaKwik_call(data, self.MinClust.value(), self.MaxClust.value())
        
        fig = self.ChanTab['FeaturesFig'].figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.set_axis_bgcolor('k')
      
        for k in range(len(clustIndx)):
            ax.plot(data[clustIndx[k],0], data[clustIndx[k],1], '.', label = 'clust %d' % k)

        ax.legend(fancybox=True, mode='expand', ncol=len(clustIndx)/2, loc=9, prop={'size':10})
        ax.grid(color='grey')
        fig.canvas.draw()
        sender().parentWidget().close()
        
    ########################################################################################################

    def TrimWaveforms_proc(self, eclick, erelease):

        # first check whether there's any waveform plotted
        # if it is visible, and if it is the current unit
        unitFound = False
        for k in self.ChanTab['WavesFigure'].figure.axes[0].get_children():
            if isinstance(k.get_label(), str)\
               and re.search('Unit[0-9]{2}', k.get_label()) and k.get_visible():
                if self.CurUnit == int(re.search('(?<=Unit)[0-9]{2}', k.get_label()).group()):
                    unitFound = True
                    break

        if not unitFound: return

        # get the indices
        indx = self.h5file.getNode('/Spikes/Chan_%03d/%s' % (self.CurChan, self.CurUnitName), 'Indx').read()
        data = self.CurWaveforms[indx,:]

        # get line equation y = mx + n
        x1 = eclick.xdata
        x2 = erelease.xdata
        y1 = eclick.ydata
        y2 = erelease.ydata

        # return if is a point and not a line
        if x1 == x2: return
        
        m = (y2 - y1)/(x2 - x1)
        n = y1 - m*x1

        # get the y value of nearest integer x:
        x = np.array([x1, x2])
        x.sort()
        xData = range(self.WfSize)
        indx1 = np.flatnonzero(xData > x[0]).min()
        indx2 = np.flatnonzero(xData < x[1]).max()
        y = np.array([m*xData[k]+n for k in range(indx1, indx2)])
        #print x, y

        # get the data bounded by the indices 
        data2 = data[:,indx1:indx2]
        #print data2.shape, y.shape
        t = data2-y
        #print t
        t = np.array(t)

        # get the indices that intersect the line
        intersect = []
        for j,k in enumerate(t):
            if (np.all(k<0) or np.all(k>0)) == False:
                intersect.append(j)

        print intersect

        # update the node containing the unit indexes
        self.h5file.removeNode(self.CurNodeName + '/' + self.CurUnitName, 'Indx')
        self.h5file.createArray(self.CurNodeName + '/' + self.CurUnitName, 'Indx', np.delete(indx, intersect))

        # add the remaining points to the unsorted indexes
        self.Unsorted = self.h5file.getNode(self.CurNodeName, 'Unsorted').read()
        self.Unsorted = np.append(self.Unsorted, indx[intersect])
        self.Unsorted.sort()
        
        # update the unsorted in the h5file
        self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # save changes to disk
        self.h5file.flush()
        
        # update the information in the overview table
        row = self.ChanSelector.currentIndex()
        self.OverviewTab2['OverviewTable'].takeItem(row, self.CurUnit+6)
        lbl = QtGui.QTableWidgetItem(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))
        self.OverviewTab2['OverviewTable'].setItem(row, self.CurUnit+6, lbl)

        # update the information on the unit label
        self.ChanTab['UnitCountLabel'][self.CurUnitName].setText(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))

        # replot the features
        self.PlotFeatures()
        
        # replot the unit avg waveform, histogram and autocorrelation
        self.PlotUnitFigure_proc()
        
        '''
        for k in intersect:
            deleted.append(int(re.search('(?<=data_)[0-9]{3}', handles[k].get_label()).group()))
            handles[k].remove()'''
        eclick.inaxes.figure.canvas.draw()
        
        self.trimWaveformsRect.set_active(False)

    def ActivateTrimWaveforms_proc(self):
        self.trimWaveformsRect.set_active(True)

    ########################################################################################################
        
    def AddUnit_proc(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.ChanPlotted: return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.ChanTab['FeaturesFig'].figure.axes[0].get_title())
        if not re.search('Waveforms|Unsorted', title): return

        # return if a tool is selected in the toolbar
        if self.ChanTab['FeaturesFigNtb'].mode!='': return

        # create a new lasso instance
        self.LassoCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.LassoAddUnit_proc)
        
    ########################################################################################################
        
    def Keep_proc(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.ChanPlotted: return

        # check if a unit is plotted
        title = str(self.ChanTab['FeaturesFig'].figure.axes[0].get_title())
        if not re.search('Unit', title): return

        self.What2Plot.count()
        
        # return if a tool is selected in the toolbar
        if self.ChanTab['FeaturesFigNtb'].mode!='': return

        # create a new lasso instance
        self.LassoCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.LassoKeep_proc)
        
    ########################################################################################################
        
    def AddRegion_proc(self):
        ''' starts a lasso instance to draw a line around a ROI'''

        # check whether there is a channel ploted
        if not self.ChanPlotted: return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.ChanTab['FeaturesFig'].figure.axes[0].get_title())
        if not re.search('Waveforms|Unsorted', title): return

        # return if a tool is selected in the toolbar
        if self.ChanTab['FeaturesFigNtb'].mode!='': return

        # create a new lasso instance
        self.LassoCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.LassoAddRegion_proc)

    ########################################################################################################

    def RemoveRegion_proc(self):
        ''' starts a lasso instance to draw a line around a ROI'''
        # check whether there is a channel ploted
        if not self.ChanPlotted: return

        # check if what is plotted is all waveforms or unsorted
        title = str(self.ChanTab['FeaturesFig'].figure.axes[0].get_title())
        if not re.search('Unit', title): return

        # return if a tool is selected in the toolbar
        if self.ChanTab['FeaturesFigNtb'].mode!='': return

        # create a new lasso instance
        self.LassoCID = self.ChanTab['FeaturesFig'].figure.canvas.mpl_connect('button_press_event',
                                                                              self.LassoRemoveRegion_proc)
        
    ########################################################################################################
        
    def LassoAddUnit_proc(self, event):
        if self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        if event.inaxes is None or event.button !=1:
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return

        # create a lasso instance
        self.lasso = matplotlib_widgets.MyLasso(event.inaxes, (event.xdata, event.ydata),
                                     self.LassoCallback_AddUnit, color = 'gray', lw=1)
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock(self.lasso)

    ########################################################################################################
        
    def LassoKeep_proc(self, event):
        if self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        if event.inaxes is None or event.button !=1:
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        self.lasso = matplotlib_widgets.MyLasso(event.inaxes, (event.xdata, event.ydata),
                                     self.LassoCallback_Keep, color = 'gray', lw=1)
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock(self.lasso)

    ########################################################################################################
        
    def LassoAddRegion_proc(self, event):
        if self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        if event.inaxes is None or event.button !=1:
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        self.lasso = matplotlib_widgets.MyLasso(event.inaxes, (event.xdata, event.ydata),
                                     self.LassoCallback_AddRegion, color = 'gray', lw=1)
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock(self.lasso)

    ########################################################################################################
        
    def LassoRemoveRegion_proc(self, event):
        if self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.locked():
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
        
        if event.inaxes is None or event.button !=1:
            if hasattr(self, 'LassoCID'):
                self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
                del self.LassoCID
            return
       
        self.lasso = matplotlib_widgets.MyLasso(event.inaxes, (event.xdata, event.ydata),
                                     self.LassoCallback_RemoveRegion, color = 'gray', lw=1)
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock(self.lasso)
        
    ########################################################################################################
        
    def LassoCallback_AddUnit(self, verts):

        # disconnect Lasso callback
        self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
        self.ChanTab['FeaturesFig'].figure.canvas.draw_idle()
        del self.LassoCID
        
        # release widget lock
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.release(self.lasso)

        # delete lasso
        del self.lasso
        
        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n, n/25)]

        #pdb.set_trace()
        
        # get the axes handle
        ax = self.ChanTab['FeaturesFig'].figure.axes[0]

        if re.search('Waveforms', ax.get_title()):
            # test which points are inside the lasso
            xypoints = self.XYData.data[self.Unsorted,:]
        elif re.search('Unsorted', ax.get_title()):
            xypoints = self.XYData.data

        p = Path(self.verts).contains_points(xypoints)
        
        # in case there were no points selected
        if len(p)==0:
            self.MsgBox.setIcon(QtGui.QMessageBox.Warning)
            self.MsgBox.setText('There were no selected points')
            self.MsgBox.setwindowTitle('Warning')
            self.MsgBox.show()
            return

        # set the unit name
        self.NUnits = len(self.UnitsList)
        self.CurUnitName = 'Unit%02d' % self.NUnits

        # look for the unsorted plot handle in the axes
        for k in self.ChanTab['FeaturesFig'].figure.axes[0].get_children():
            if re.search('Unsorted', str(k.get_label())):
                break

        # obtain the unsorted points
        unsortedData = xypoints[~p,:]
        lunsort = len(unsortedData)

        # select some indices to plot
        if lunsort > self.nPtsSpin.value():
            indx = range(0, lunsort, lunsort/self.nPtsSpin.value())
        else:
            indx = range(lunsort)
        
        # replot the unsorted without the corresponding points to the new unit
        k.set_data(unsortedData[:,0][indx], unsortedData[:,1][indx])
        ax.draw_artist(k)

        # select some indices to plot
        unitData = xypoints[p,:]
        lunit = len(unitData)

        if lunit > self.nPtsSpin.value():
            indx = range(0, lunit, lunit/self.nPtsSpin.value())
        else:
            indx = range(lunit)
            
        ax.plot(unitData[:,0][indx], unitData[:,1][indx],
                linestyle = '',
                marker = ',',
                mfc = self.UnitColors[self.NUnits],
                mec = self.UnitColors[self.NUnits],
                label = 'data_'+self.CurUnitName)

        self.NUnits+=1

        # if unit name not in combo box add it
        if self.CurUnitName not in [str(self.What2Plot.itemText(k)) for k in range(self.What2Plot.count())]:
            self.What2Plot.addItem(self.CurUnitName)

        # add the indexes of the current unit to the h5file
        if self.h5file.getNode(self.CurNodeName).__contains__(self.CurUnitName):
            self.h5file.removeNode(self.CurNodeName, self.CurUnitName, recursive = True)
        self.h5file.createGroup(self.CurNodeName, self.CurUnitName)
        self.h5file.createArray(self.CurNodeName  + '/' + self.CurUnitName, 'Indx', self.Unsorted[p])
        self.h5file.createArray(self.CurNodeName  + '/' + self.CurUnitName, 'isMultiunit', False)
        self.h5file.createArray(self.CurNodeName  + '/' + self.CurUnitName, 'isBursting', False)

        # update the list of unsorted indexes
        self.Unsorted = self.Unsorted[~p]
        
        # update  the indexes of the unsorted units
        if self.h5file.getNode(self.CurNodeName).__contains__('Unsorted'):
            self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # save changes to disk
        self.h5file.flush()

        # add log
        self.AddLog('%s %s added' % (self.CurNodeName, self.CurUnitName))
        
        # add unit to the units tab widget
        self.UnitsTable_AddUnit(self.CurUnitName)

        self.ChanTab['UnitFigures'][self.CurUnitName].figure.tight_layout()
        self.ChanTab['UnitFigures'][self.CurUnitName].figure.canvas.draw()
        
        # update the overview figure
        for j,k in enumerate(self.OverviewTab1['Figure'].figure.axes):
            if k.get_title().find(str(self.CurChan))!=-1:
                break
        self.PlotChanOverview_proc(self.CurNode, axes2Plot = k)
        for l in k.lines:
            k.draw_artist(l)

    ########################################################################################################
        
    def LassoCallback_Keep(self, verts):
        # disconnect Lasso callback from figure
        self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
        self.ChanTab['FeaturesFig'].figure.canvas.draw_idle()
        del self.LassoCID
        
        # release the lock from the lasso
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.release(self.lasso)

        # erase lasso
        del self.lasso
        
        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n ,n/25)]
        
        # test which points lay inside the polygon
        p = Path(self.verts).contains_points(self.XYData.data)
        self.XYData = cKDTree(self.XYData.data[p,:])

        # get unitname and number from the axes title
        ax = self.ChanTab['FeaturesFig'].figure.axes[0]
        self.CurUnitName = re.search('Unit[0-9]{2}', ax.get_title()).group()
        self.CurUnit = int(re.search('(?<=Unit)[0-9]{2}', ax.get_title()).group())
        
        # update plot:
        for k in ax.get_children():
            if re.search(str(k.get_label), self.CurUnitName):
                k.set_data(self.XYData.data[:,0], self.XYData.data[:,1])
                ax.draw_artist(k)
                break
        
        # return if no points selected
        if len(p) < 1: return

        nodeName = self.CurNodeName + '/' + self.CurUnitName
        
        # obtain the unit data
        unitPts = self.h5file.getNode(nodeName, 'Indx').read()

        # update the node containing the unit indexes
        self.h5file.removeNode(nodeName, 'Indx')
        self.h5file.createArray(nodeName, 'Indx', unitPts[p])

        # add the remaining points to the unsorted indexes
        self.Unsorted = self.h5file.getNode(self.CurNodeName, 'Unsorted').read()
        self.Unsorted = np.append(self.Unsorted, unitPts[~p])
        self.Unsorted.sort()
        
        # update the unsorted in the h5file
        self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # save changes to disk
        self.h5file.flush()

        # replot the unit avg waveform, histogram and autocorrelation
        self.PlotUnitFigure_proc()

        # replot the features
        self.PlotFeatures()

        # update the information in the overview table
        row = self.ChanSelector.currentIndex()
        self.OverviewTab2['OverviewTable'].takeItem(row, self.CurUnit+6)
        lbl = QtGui.QTableWidgetItem(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))
        self.OverviewTab2['OverviewTable'].setItem(row, self.CurUnit+6, lbl)

        # update the information on the unit label
        self.ChanTab['UnitCountLabel'][self.CurUnitName].setText(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))

    ########################################################################################################
        
    def LassoCallback_AddRegion(self, verts):

        # disconnect the lasso from the canvas and redraw the figure
        self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
        self.ChanTab['FeaturesFig'].figure.canvas.draw_idle()
        del self.LassoCID

        # release widget lock
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.release(self.lasso)

        # delete lasso handle
        del self.lasso
        
        # get the vertices of the polygon to the object and downsample them if too many
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n, n/25)]
        
        # check whether there is any unit
        if not hasattr(self, 'CurUnitName') or not self.CurNode.__contains__(self.CurUnitName):
            return

        # get the axes handle
        ax = self.ChanTab['FeaturesFig'].figure.axes[0]

        # get the unsorted
        self.Unsorted = self.h5file.getNode(self.CurNodeName, 'Unsorted').read()
    
        # check what is plotted on the axes
        if re.search('Waveforms', str(ax.get_title())):
            # test which points are inside the lasso
            p = Path(self.verts).contains_points(self.XYData.data[self.Unsorted,:])
            self.XYData = cKDTree(self.XYData.data[self.Unsorted,:][p])
            
        elif re.search('Unsorted', str(ax.get_title())):
            # test which points are inside the lasso
            p = Path(self.verts).contains_points(self.XYData.data)
            self.XYData = cKDTree(self.XYData.data[p,:])

        # update plot:
        for k in ax.get_children():
            if re.search(str(k.get_label), self.CurUnitName):
                k.set_data(self.XYData.data[:,0], self.XYData.data[:,1])
                ax.draw_artist(k)
                break
        
        # if more than 0 selected points
        if len(p) > 0:
            indx = self.Unsorted[p]
        else:
            return

        # update the unsorted
        self.Unsorted = self.Unsorted[~p]
        self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # update the plots
        for k in ax.get_children():
            if re.search('Unsorted', str(k.get_label())):
                pass
            elif re.search(self.CurUnitName, str(k.get_label())):
                pass
    
        # update the unit information in the file
        unit = self.h5file.getNode(self.CurNodeName + '/' + self.CurUnitName, 'Indx').read()
        self.h5file.removeNode(self.CurNodeName + '/' + self.CurUnitName, 'Indx')
        # append the new indexes to the waveform and sort
        unit = np.append(unit, indx)
        unit.sort()
        # create a new array in the h5file to hold the updated unit information
        self.h5file.createArray(self.CurNodeName+ '/' + self.CurUnitName, 'Indx', unit)

        # save changes to disk
        self.h5file.flush()

        # update the information in the overview table
        row = self.ChanSelector.currentIndex()
        self.OverviewTab2['OverviewTable'].takeItem(row, self.CurUnit+6)
        lbl = QtGui.QTableWidgetItem(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))
        self.OverviewTab2['OverviewTable'].setItem(row, self.CurUnit+6, lbl)

        # update the information on the unit label
        self.ChanTab['UnitCountLabel'][self.CurUnitName].setText(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))
        
        # replot the unit avg waveform, histogram and autocorrelation
        self.PlotUnitFigure_proc()

        # replot the features
        self.PlotFeatures()
        
    ########################################################################################################
        
    def LassoCallback_RemoveRegion(self, verts):

        # disconnect Lasso callback from figure
        self.ChanTab['FeaturesFig'].figure.canvas.mpl_disconnect(self.LassoCID)
        self.ChanTab['FeaturesFig'].figure.canvas.draw_idle()
        del self.LassoCID
        
        # release the lock from the lasso
        self.ChanTab['FeaturesFig'].figure.canvas.widgetlock.release(self.lasso)

        # copy the vertices of the polygon to the object and downsample them
        n = len(verts)
        self.verts = np.array(verts)
        if n > 25:
            self.verts = self.verts[range(0, n ,n/25)]
        
        # test which points lay inside the polygon
        p = Path(self.verts).contains_points(self.XYData.data)

        # return if no points selected
        if len(p) < 1: return

        # get unitname and number from the axes title
        ax = self.ChanTab['FeaturesFig'].figure.axes[0]
        self.CurUnitName = re.search('Unit[0-9]{2}', ax.get_title()).group()
        self.CurUnit = int(re.search('(?<=Unit)[0-9]{2}', ax.get_title()).group())
        
        # obtain the unit data
        unitPts = self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.read()

        # update the node containing the unit indexes
        self.h5file.removeNode(self.CurNodeName + '/' + self.CurUnitName, 'Indx')
        self.h5file.createArray(self.CurNodeName + '/' + self.CurUnitName, 'Indx', unitPts[~p])

        # add the remaining points to the unsorted indexes
        self.Unsorted = self.h5file.getNode(self.CurNodeName, 'Unsorted').read()
        self.Unsorted = np.append(self.Unsorted, unitPts[p])
        self.Unsorted.sort()
        
        # update the unsorted in the h5file
        self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # save changes to disk
        self.h5file.flush()
        
        # update the information in the overview table
        row = self.ChanSelector.currentIndex()
        self.OverviewTab2['OverviewTable'].takeItem(row, self.CurUnit+6)
        lbl = QtGui.QTableWidgetItem(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))
        self.OverviewTab2['OverviewTable'].setItem(row, self.CurUnit+6, lbl)

        # update the information on the unit label
        self.ChanTab['UnitCountLabel'][self.CurUnitName].setText(str(self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.nrows))

        # replot the features
        self.PlotFeatures()
        
        # replot the unit avg waveform, histogram and autocorrelation
        self.PlotUnitFigure_proc()

        #PlotChanOverview_proc(self.CurNode, axes2Plot)

        # erase lasso
        del self.lasso
        
    ########################################################################################################
        
    def UnitsTable_AddUnit(self, unitName):
        ''' creates a new tab per each new unit'''

        # check whether that tab already exists
        for k in range(self.ChanTab['UnitTabsWidget'].count()):
            if unitName == self.ChanTab['UnitTabsWidget'].tabText(k):
                return

        self.CurUnitName = unitName
        
        # create a widget and a layout
        widget = QtGui.QWidget()
        vlay   = QtGui.QVBoxLayout()
        vlay.setSpacing(2)
        vlay.setMargin(0)

        # get a unit number
        unitNo = int(re.search('(?<=Unit)[0-9]{2}', unitName).group())

        # add the unit number to a list
        self.UnitsList.append(unitNo)

        # create a btn to change unit color
        hlay = QtGui.QHBoxLayout()
        hlay.setMargin(0)
        hlay.addStretch(1)
        self.ChanTab['UnitBtns'][unitName] = QtGui.QPushButton('Unit %02d' % unitNo)
        self.ChanTab['UnitBtns'][unitName].setMaximumHeight(20)
        self.ChanTab['UnitBtns'][unitName].clicked.connect(self.ChangeUnitColor_proc)
        self.ChanTab['UnitBtns'][unitName].setStyleSheet('QPushButton {background: rgb%s}' % str(tuple(np.int16(255*self.UnitColors[unitNo]))))
        hlay.addWidget(self.ChanTab['UnitBtns'][unitName])
        hlay.addStretch(1)

        # plot-raw check button
        self.ChanTab['PlotRawCheck'][unitName] = QtGui.QCheckBox()
        self.ChanTab['PlotRawCheck'][unitName].setObjectName(str(unitNo))
        self.ChanTab['PlotRawCheck'][unitName].setChecked(False)
        self.ChanTab['PlotRawCheck'][unitName].setMaximumHeight(20)
        self.ChanTab['PlotRawCheck'][unitName].stateChanged.connect(self.SetWaveformVisible_proc)
        lbl = QtGui.QLabel('Plot Raw ?')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.ChanTab['PlotRawCheck'][unitName])
        hlay.addStretch(1)

        # is Multiunit check button
        self.ChanTab['isMultiunitCheck'][unitName] = QtGui.QCheckBox()
        self.ChanTab['isMultiunitCheck'][unitName].setObjectName(str(unitNo))
        self.ChanTab['isMultiunitCheck'][unitName].setChecked(False)
        self.ChanTab['isMultiunitCheck'][unitName].setMaximumHeight(20)
        self.ChanTab['isMultiunitCheck'][unitName].stateChanged.connect(self.SetisMultiunit_proc)
        lbl = QtGui.QLabel('isMultiunit ?')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        hlay.addWidget(self.ChanTab['isMultiunitCheck'][unitName])
        hlay.addStretch(1)

        # set the checkstate of the 'isMultiunit' check according to what is saved in the h5file
        if self.h5file.getNode('/Spikes/Chan_%03d/Unit%02d' % (self.CurChan, unitNo)).__contains__('isMultiunit'):
            isMultiunit  = self.h5file.getNode('/Spikes/Chan_%03d/Unit%02d' % (self.CurChan, unitNo),
                                               'isMultiunit').read()
            if isMultiunit: self.ChanTab['isMultiunitCheck'][unitName].setChecked(True)
            else: self.ChanTab['isMultiunitCheck'][unitName].setChecked(False)
        else:         
            self.h5file.createArray('/Spikes/Chan_%03d/Unit%02d' % (self.CurChan, unitNo), 'isMultiunit', False)

        # add a label with the waveform count
        lbl = QtGui.QLabel('Count')
        lbl.setMaximumHeight(20)
        hlay.addWidget(lbl)
        self.ChanTab['UnitCountLabel'][unitName] = QtGui.QLabel('%d' % self.h5file.getNode(self.CurNodeName, unitName).Indx.nrows)
        self.ChanTab['UnitCountLabel'][unitName].setMaximumHeight(20)
        hlay.addWidget(self.ChanTab['UnitCountLabel'][unitName])
        hlay.addStretch(1)

        # add delete-unit button
        self.ChanTab['DelUnitBtns'][unitName] = QtGui.QPushButton('Del Unit')
        self.ChanTab['DelUnitBtns'][unitName].setObjectName(unitName)
        self.ChanTab['DelUnitBtns'][unitName].setMaximumHeight(20)
        self.ChanTab['DelUnitBtns'][unitName].clicked.connect(self.DelUnit_proc)
        hlay.addWidget(self.ChanTab['DelUnitBtns'][unitName])
        hlay.addStretch(1)
        
        vlay.addLayout(hlay)

        # add the figure widget
        self.ChanTab['UnitFigures'][unitName] = matplotlib_widgets.MplWidget()
        self.ChanTab['UnitFigures'][unitName].setObjectName(unitName) # set the name of the object
        self.ChanTab['UnitFigures'][unitName].figure.set_facecolor('k')
        n = matplotlib_widgets.NavToolbar(self.ChanTab['UnitFigures'][unitName].figure.canvas, widget, coordinates = False)
        n.setIconSize(QtCore.QSize(12,12))
        n.setOrientation(QtCore.Qt.Vertical)

        vlay.addWidget(self.ChanTab['UnitFigures'][unitName])
        #vlay.addWidget(n)

        hlay = QtGui.QHBoxLayout()
        hlay.setSpacing(0)
        hlay.setMargin(2)
        hlay.addLayout(vlay)
        hlay.addWidget(n)
        
        widget.setLayout(hlay)

        # Plot the data
        self.PlotUnitFigure_proc()
        
        #if unitName == 'Unit00':
        #    self.ChanTab['UnitTabsWidget'].removeTab(0)

        self.ChanTab['UnitTabsWidget'].addTab(widget, unitName)
        indx = self.ChanTab['UnitTabsWidget'].count()-1
        self.ChanTab['UnitTabsWidget'].tabBar().setTabTextColor(indx,
                                                                QtGui.QColor(int(self.UnitColors[indx,0]*255),
                                                                             int(self.UnitColors[indx,1]*255),
                                                                             int(self.UnitColors[indx,2]*255)))


        # update the information in the overview table
        row = self.ChanSelector.currentIndex()

        if self.OverviewTab2['OverviewTable'].columnCount() <= (unitNo+6):
            self.OverviewTab2['OverviewTable'].insertColumn(self.OverviewTab2['OverviewTable'].columnCount())
            nCols = self.OverviewTab2['OverviewTable'].columnCount()
            self.OverviewTab2['OverviewTable'].setColumnWidth(nCols-1, 65)
            self.OverviewTab2['OverviewTable'].setHorizontalHeaderItem(nCols-1,
                                                                       QtGui.QTableWidgetItem('Unit%02d' % unitNo))
        
        self.OverviewTab2['OverviewTable'].takeItem(row, unitNo+6)
        lbl = QtGui.QTableWidgetItem(str(self.h5file.getNode(self.CurNodeName, unitName).Indx.nrows))
        self.OverviewTab2['OverviewTable'].setItem(row, unitNo+6, lbl)

        # update the unsorted number in the overview table
        self.OverviewTab2['OverviewTable'].takeItem(row, 4)
        lbl = QtGui.QTableWidgetItem(str(len(self.Unsorted)))
        self.OverviewTab2['OverviewTable'].setItem(row, 4, lbl)

    ########################################################################################################
        
    def PlotUnitFigure_proc(self):
        
        # get a unit name and number
        unitNo = int(re.search('(?<=Unit)[0-9]{2}', self.CurUnitName).group())

        # find the figure that has a particular name
        fig = self.ChanTab['UnitFigures'][self.CurUnitName].figure
        
        # check whether we have to create axes
        if len(fig.axes) == 0:
            ax0 = fig.add_subplot(131)
            ax1 = fig.add_subplot(132)
            ax2 = fig.add_subplot(133)
        else:
            ax0 = fig.axes[0]; ax0.cla()
            ax1 = fig.axes[1]; ax1.cla()
            ax2 = fig.axes[2]; ax2.cla()

        # set the axis background color
        ax0.set_axis_bgcolor('k')
        ax1.set_axis_bgcolor('k')
        ax2.set_axis_bgcolor('k')
        
        ##### PLOT AVERAGE WAVEFORM #####
        x = range(self.WfSize)
        p = self.h5file.getNode(self.CurNodeName, self.CurUnitName).Indx.read()
        m = self.CurWaveforms[p,:].mean(axis=0)
        s = self.CurWaveforms[p,:].std(axis=0)
        mn = self.CurWaveforms[p,:].min(axis=0)
        mx = self.CurWaveforms[p,:].max(axis=0)

        # plot average waveform
        ax0.plot(x, m  , color = self.UnitColors[unitNo], lw=2, label = self.CurUnitName)

        #plot shaded area of 3 standard devoations around it 
        ax0.fill_between(x, m+3*s, m-3*s, color = self.UnitColors[unitNo], alpha=0.5, label = self.CurUnitName)

        #plot maximum and minimum boundaries
        ax0.fill_between(x, mx, mn, color = self.UnitColors[unitNo], alpha=0.35, label = self.CurUnitName)
        ax0.set_xlim(0, self.WfSize-1)
        ax0.set_yticklabels([])
        ax0.grid(color = [.5, .5, .5])
        ax0.tick_params(color = [.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax0.spines.values(): k.set_edgecolor([.5, .5, .5])
        
        ##### PLOT ISI HISTOGRAM #####
        dts = np.diff(self.CurTs[p])
        dts = dts[dts<100]

        ld = len(dts)
        if ld > 1000:
            indx = range(0, 1000)
        else:
            indx = range(ld)

        ax1.hist(dts[indx], bins = 100, range = [0,100], ec='none',
                 color=self.UnitColors[unitNo], label = self.CurUnitName)
        ax1.tick_params(color = [.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax1.spines.values(): k.set_edgecolor([.5, .5, .5])
        WfWidth = self.WfSize*1000/self.Sf
        try:
            collision = 100*np.flatnonzero(dts < 1.5*WfWidth ).size/np.float(dts.size)
            # put a "percentage of collision" label
            ax1.text(0.5, 0.01, u'Collision = %0.2f %%' % collision,
                     transform = ax1.transAxes, color = 'w', size=10, ha = 'center')
        except:
            pass
        ax1.set_xlim(0,100)

        ##### PLOT AUTOCORRELATION #####
        
        ts = self.CurTs[p]
        if ts.size > 1000: ts = ts[0:1000]
            
        ac, x = autocorr(ts, binSize = 10, Win = [0,10000],
                         mode = 'fft', Range = [-150, 150])
        ax2.plot(x, ac, color = self.UnitColors[unitNo], lw = 2)
        ax2.set_xlim(-200, 200)
        ax2.tick_params(color = [.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax2.spines.values(): k.set_edgecolor([.5, .5, .5])
        ax2.set_yticklabels([])

        self.ChanTab['UnitFigures'][self.CurUnitName].figure.tight_layout()
        self.ChanTab['UnitFigures'][self.CurUnitName].figure.canvas.draw()

    ########################################################################################################
        
    def SetWaveformVisible_proc(self):
        ''' makes the raw waveform of each unit visible or invisible'''
       
        sender = self.sender()
        state  = sender.checkState()
        name   = int(sender.objectName())

        # get unit name and number
        unitName = str(self.ChanTab['UnitTabsWidget'].tabText(self.ChanTab['UnitTabsWidget'].currentIndex()))
        unitNo   = int(re.search('(?<=Unit)[0-9]{2}', unitName).group())

        # get axes handle and children labels
        ax = self.ChanTab['WavesFigure'].figure.axes[0]
        childrenLabels = [str(k.get_label()) for k in ax.get_children()]

        # get the node to read from
        node = self.h5file.getNode(self.CurNodeName +'/'+ unitName, 'Indx')

        # get the number of spikes to plot
        nspikes = self.NSpikesSpin.value()
        
        if state == 2: # if checked
            nrows    = node .nrows
            if nrows > nspikes:
                unitIndx = node.read(start = 0, stop = nrows, step = nrows/nspikes)
            else:
                unitIndx = node.read()
            
            # obtain the length of units to plot
            n = len(unitIndx)

            # create an array of Nones to append
            nones = np.array(n*[None], ndmin=2).transpose()

            # create the x indexes
            Ts = np.tile(np.arange(self.WfSize),(n,1))
            Ts = np.append(Ts, nones, axis = 1).reshape((n*(self.WfSize+1),))
            
            # get the waveforms, append nones, and reshape it to a vector
            Wf = self.CurNode.Waveforms[unitIndx,:]
            Wf = np.append(Wf, nones, axis=1).reshape((n*(self.WfSize+1),))

            # create the plot if it doesn't exists
            if unitName not in childrenLabels:
                ax.plot(Ts, Wf,
                        color = self.UnitColors[unitNo,:],
                        alpha = 0.7,
                        label = unitName)

            # if exists update the data
            elif unitName in childrenLabels:
                for k in self.ChanTab['WavesFigure'].figure.axes[0].get_children():
                    if k.get_label()=='Unit%02d' % name:
                        break
                    
                k.set_data(Ts, Wf)
                k.set_visible(True)

        elif state == 0: # if unchecked
            for k in ax.get_children():
                if str(k.get_label()) == unitName:
                    k.set_visible(False)

        # set axes limit
        lim = self.WaveAxYLim_Spin.value()
        ax.set_ylim(-lim, lim)
        
        # finally redraw the figure
        self.ChanTab['WavesFigure'].figure.canvas.draw()

    ########################################################################################################
        
    def SetisMultiunit_proc(self):

        sender = self.sender()
        state  = sender.checkState()
        unitNo = int(sender.objectName())

        # current node name
        nodeName = '/Spikes/Chan_%03d/Unit%02d' % (self.CurChan, unitNo)

        # eliminate 'isMultiunit' f already exists
        if self.h5file.getNode(nodeName).__contains__('isMultiunit'):
            self.h5file.removeNode(nodeName, 'isMultiunit')

        # get the value of the checkbox and save it to "val"
        if sender.checkState() == 2: val = True            
        else: val = False

        # create a new "isMultiunt" array to hold the value of the cehckbox
        self.h5file.createArray(nodeName, 'isMultiunit', val)

        # save changes to disk
        self.h5file.flush()

    ########################################################################################################
            
    def ExchangeUnitName_proc(self, initial, final):
        tb = self.ChanTab['UnitTabBarWidget']

        # get the names of the changed tabs
        oldNameBgTab    = tb.tabText(final)
        newNameBgTab    = tb.tabText(initial)

        # change the name of the background tab
        tb.setTabText(final, newNameBgTab)

        # change the name of the front tab the oldname of the unit
        tb.setTabText(tb.currentIndex(), oldNameBgTab)

        ###### PROPAGATE CHANGES TO THE H5FILE #####

        # first change the background moved unit name to "tmpUnitData"
        self.h5file.renameNode(where = '/Spikes/Chan_%03d' % self.CurChan, name=oldNameBgTab,
                               newname = 'tmpUnitData', overwrite=True)

        # second change the front moved unit name to the old name of the background unit
        self.h5file.renameNode(where = '/Spikes/Chan_%03d' % self.CurChan, name=newNameBgTab,
                               newname = oldNameBgTab, overwrite=True)

        # third change the background moved unit name to its new name
        self.h5file.renameNode(where = '/Spikes/Chan_%03d' % self.CurChan, name='tmpUnitData',
                               newname = newNameBgTab, overwrite=True)

        ##### CHANGE THE NAME OF THE FIGURES #####

        # first change the figure name of the background unit to "tmpFigName"
        for k in self.ChanTab['UnitFigures']:
            if k.objectName() == oldNameBgTab:
                k.setObjectName('tmpFigName')
                break

        # second, change the front tab figure name to the old background tab name
        for k in self.ChanTab['UnitFigures']:
            if k.objectName() == newNameBgTab:
                k.setObjectName(oldNameBgTab)
                break

        # third, change the figname of the background tab to the new one
        for k in self.ChanTab['UnitFigures']:
            if k.objectName() == 'tmpFigName':
                k.setObjectName(newNameBgTab)
                break

        ##### CHANGE UNIT COLOR #####
        self.ChangeUnitColor_proc(unitName = newNameBgTab,
                                  color = tuple(np.append(np.int32(self.UnitColors[final]*255),255)))
        self.ChangeUnitColor_proc(unitName = oldNameBgTab,
                                  color = tuple(np.append(np.int32(self.UnitColors[initial]*255),255)))
        
            
    ########################################################################################################

    def RepairUnitNames_proc(self):

        unitNames = [k for k in self.CurNode.__members__ if k.find(Unit)!=-1]

        for j,k in enumerate(unitNames):
            if k != 'Unit%02d' % j:
                self.h5file.renameNode(self.CurChan, k, 'Unit%02d' % j)

        self.h5file.flush()

    ########################################################################################################
        
    def CallMergeUnits_proc(self):
        if not self.H5FileLoaded == True: return
        self.MergeUnitsWidget.list1.clear()
        self.MergeUnitsWidget.list2.clear()

        unitsList = [k for k in self.CurNode.__members__ if k.find('Unit')!=-1]
        unitsList.sort()
        self.MergeUnitsWidget.list1.addItems(unitsList)
        self.MergeUnitsWidget.show()
        
    ########################################################################################################
        
    def MergeUnits_proc(self):

        # get the list of units to merge
        units2Merge = [self.MergeUnitsWidget.list2.item(k).text() for k in range(self.MergeUnitsWidget.list2.count())]
        # sort the names
        units2Merge.sort()
        # if fewer than 2 return
        if len(units2Merge) < 2: return

        # store the unit indexes in a list, sort them, and trnasform it into an array
        newUnit = []
        for k in units2Merge:
            newUnit.extend(self.CurNode.__getattr__(k).Indx.read())
        newUnit.sort()
        newUnit = np.array(newUnit)

        # remove all the listed units from the h5file
        for k in units2Merge:
            self.h5file.removeNode(self.CurNodeName, k, recursive = True)

        # create a group with the name of the first unit in the list, and
        # add all the indices of that
        self.h5file.createGroup(self.CurNodeName, units2Merge[0])
        self.h5file.createArray(self.CurNodeName + '/' + units2Merge[0], 'Indx', newUnit)
        self.h5file.createArray(self.CurNodeName + '/' + units2Merge[0], 'isMultiunit', False)
        self.h5file.createArray(self.CurNodeName + '/' + units2Merge[0], 'isBursting', False)
        
        # save changes to disk
        self.h5file.flush()

        # add log
        self.AddLog('%s %s merged' % (self.CurNodeName, str(units2Merge)))

        ##### REMOVE ALL THE GRAPHICAL ELEMENTS ##### 
        # get the axes to remove from
        ax = self.ChanTab['WavesFigure'].figure.axes[0]
        
        for k in units2Merge[1:]:
            # remove the tabs
            for tabIndx in range(self.ChanTab['UnitTabsWidget'].count()):
                if str(self.ChanTab['UnitTabsWidget'].tabText(tabIndx)) == k:
                    self.ChanTab['UnitTabsWidget'].removeTab(tabIndx)

            # remove unit figure
            self.ChanTab['UnitFigures'][k].figure.clear()
            self.ChanTab['UnitFigures'][k].close()
            self.ChanTab['UnitFigures'].pop(k, 0)
                
            # removes the unitname from the what2 plot list
            for n in range(self.What2Plot.count()):
                if self.What2Plot.itemText(n) == k:
                    self.What2Plot.removeItem(n)

            # eliminate the raw waveforms from the plot
            for a in ax.get_children():
                if str(a.get_label()).find(k)!=-1:
                    a.remove()

            # remove the unit from the list
            unitNo = int(re.search('[0-9]{2}', k).group())
            self.UnitsList.remove(unitNo)
            
        # update the information in the overview table
        #self.OverviewTab2['OverviewTable'].takeItem(self.ChansList.index(self.CurChan),
        #                                            unitNo+4)
   
        # redraw the waveforms figure
        self.ChanTab['WavesFigure'].figure.canvas.draw()            

        # replot features
        self.PlotFeatures()
        
        # add the merged unit to the table                        
        self.UnitsTable_AddUnit(units2Merge[0])        

    ########################################################################################################

    def CallMoveUnits_proc(self):
        if not self.H5FileLoaded == True: return
        self.MoveUnitsWidget.list.clear()
        unitsList = [k for k in self.CurNode.__members__ if k.find('Unit')!=-1]
        unitsList.sort()
        self.MoveUnitsWidget.list.addItems(unitsList)
        self.MoveUnitsWidget.show()
        
    ########################################################################################################

    def MoveUnits_proc(self):

        # first get the needed changes
        old = []; new = []
        for k in range(self.MoveUnitsWidget.list.count()):
            if 'Unit%02d' % k != str(self.MoveUnitsWidget.list.item(k).text()):
                old.append(str(self.MoveUnitsWidget.list.item(k).text()))
                new.append('Unit%02d' % k)

        # in case no changes are needed
        if len(old) == 0: return
        
        ##### RENAME ALL THE UNITS AND GRAPHICAL ELEMENTS TO "_tmp" #####
        for k in self.CurNode.__members__:
            if k.find('Unit') != -1:
                # rename the nodes
                self.h5file.renameNode(self.CurNodeName, name = k, newname = k+'_tmp')

                for key in ['UnitFigures', 'UnitCountLabel', 'DelUnitBtns',
                            'PlotRawCheck', 'UnitBtns', 'isMultiunitCheck']:
                    self.ChanTab[key][k+'_tmp'] = self.ChanTab[key][k]
                    self.ChanTab[key].pop(k, 0) # remove

        for j,k in zip(old, new):
            self.ChangeUnitName_proc(j+'_tmp', k)

        # move everything back
        for k in self.CurNode.__members__:
            if k.find('_tmp') != -1:
                if k.replace('_tmp','') in self.CurNode.__members__:
                    self.h5file.removeNode(self.CurNodeName, name = k)
                    for key in ['UnitFigures', 'UnitCountLabel', 'DelUnitBtns',
                                'PlotRawCheck', 'UnitBtns', 'isMultiunitCheck']:
                        self.ChanTab[key].deleteLater()
                        self.ChanTab[key].pop(k, 0)
                else:
                    self.h5file.renameNode(self.CurNodeName, name = k, newname = k.replace('_tmp',''))
                    for key in ['UnitFigures', 'UnitCountLabel', 'DelUnitBtns',
                                'PlotRawCheck', 'UnitBtns', 'isMultiunitCheck']:
                        self.ChanTab[key][k.replace('_tmp','')] = self.ChanTab[key][k]
                        self.ChanTab[key].pop(k, 0)
                    
        # save changes to disk
        self.h5file.flush()
            
    ########################################################################################################
    
    def ChangeUnitName_proc(self, oldName, newName):

        # rename node                    
        self.h5file.renameNode(self.CurNodeName, name = oldName, newname = newName, overwrite = True)
                    
        # get the unit numbers from the names
        oldUnitNo = int(re.search('[0-9]{2}', oldName).group())
        newUnitNo = int(re.search('[0-9]{2}', newName).group())

        # move the tab and change its name
        self.ChanTab['UnitTabBarWidget'].setTabText(newUnitNo, newName)
        self.ChanTab['UnitTabBarWidget'].moveTab(oldUnitNo, newUnitNo)

        for key in ['UnitFigures', 'UnitCountLabel', 'DelUnitBtns',
                    'PlotRawCheck', 'UnitBtns', 'isMultiunitCheck']:
            self.ChanTab[key][newName] = self.ChanTab[key][oldName]
            self.ChanTab[key][newName].setObjectName(newName)
            self.ChanTab[key].pop(oldName,0)

        # change color of the unit
        self.ChangeUnitColor_proc(newName, color = 255*self.UnitColors[newUnitNo])
        

    ########################################################################################################
    
    def CleanWavesFigure_proc(self):
        self.ChanTab['WavesFigure'].figure.canvas.draw()

    ########################################################################################################
        
    def UnitsTable_AddRow(self):
        self.CurUnit     = self.ChanTab['UnitTabsWidget'].currentIndex()
        self.CurUnitName = self.ChanTab['UnitTabsWidget'].tabText(self.CurUnit)

    ########################################################################################################
        
    def DelUnit_proc(self):
        if not self.H5FileLoaded or not self.ChanPlotted: return

        # get sender
        sender = self.sender()

        # get unit name and number
        unitName = str(sender.objectName())
        unitNo = int(re.search('(?<=Unit)[0-9]{2}', unitName).group())

        # remove the unit from the list
        self.UnitsList.remove(unitNo)
        
        # get the indexes of the unit
        indx = self.h5file.getNode(self.CurNodeName, unitName).Indx.read()

        # get unsorted, append the indexes from the unit and update that
        # to the h5file
        self.Unsorted = self.h5file.getNode(self.CurNodeName, 'Unsorted').read()
        self.Unsorted = np.append(self.Unsorted, indx)
        self.Unsorted.sort()
        self.h5file.removeNode(self.CurNodeName, 'Unsorted')
        self.h5file.removeNode(self.CurNodeName, unitName, recursive = True)
        self.h5file.createArray(self.CurNodeName, 'Unsorted', self.Unsorted)

        # add log
        self.AddLog('%s %s deleted' % (self.CurNodeName, unitName))

        # remove the tab
        for tabIndx in range(self.ChanTab['UnitTabsWidget'].count()):
            if str(self.ChanTab['UnitTabsWidget'].tabText(tabIndx)) == unitName:
                break
        self.ChanTab['UnitTabsWidget'].removeTab(tabIndx)

        # close and remove unit figure
        plt.close(self.ChanTab['UnitFigures'][unitName].figure)
        self.ChanTab['UnitFigures'].pop(unitName, 0)
                
        # removes the unitname from the what2 plot list
        for n in range(self.What2Plot.count()):
            if self.What2Plot.itemText(n) == unitName:
                self.What2Plot.removeItem(n)

        # update the information in the overview table
        self.OverviewTab2['OverviewTable'].takeItem(self.ChansList.index(self.CurChan),
                                                    unitNo+4)

        # eliminate the raw waveforms from the plot
        ax = self.ChanTab['WavesFigure'].figure.axes[0]
        for k in ax.get_children():
            if str(k.get_label()).find(unitName)!=-1:
                k.remove()
                break

        # redraw the waveforms figure
        self.ChanTab['WavesFigure'].figure.canvas.draw()            

        # replot features
        self.PlotFeatures()
            

    ########################################################################################################
            
    def ChangeUnitColor_proc(self, unitName = None, color = None):
        ''' Change unit color utility function
        inputs:
            unitName : string containing the unit name
            color    : must be a four element RGB tuple from 0 to 255, for example,
                       the output of getRgB() output from a Qt Color instance.
                       The fourth element is the alpha (usually = to 255)'''
        
        if unitName in [None, False]:
            sender = self.sender()
            unitName = str(sender.text()).replace(' ','')

        unitNo = int(re.search('[0-9]{1,3}', unitName).group())

        if not np.any(color):
            c = QtGui.QColorDialog()
            color = c.getColor(sender.palette().color(1))
            if not color.isValid(): return
            
        if isinstance(color, QtGui.QColor):
            qtColor  = color
        else:
            qtColor  = QtGui.QColor(color[0], color[1], color[2])
    
        mplColor = np.array(qtColor.getRgb()[0:3])/255.0
        
        if type(self.sender()) == QtGui.QPushButton and str(self.sender().text()).find('Unit') != -1:
            self.sender().setStyleSheet('QPushButton {background: rgb%s}' % str(qtColor.getRgb()[0:3]))

        self.UnitColors[unitNo,:] = mplColor

        # get the figure with a name equal to the current unit
        ax = self.ChanTab['UnitFigures'][unitName].figure.axes
        
        # iterate over axes to change colors
        for k in ax:
            for j in k.lines:
                j.set_color(mplColor)
            for j in k.collections:
                j.set_color(mplColor)
            for j in k.patches:
                j.set_color(mplColor)

        # search a figure with a specific name
        self.ChanTab['UnitFigures'][unitName].figure.canvas.draw()

        # change the color of the raw waveforms
        for k in self.ChanTab['WavesFigure'].figure.axes[0].lines:
            if re.search('Unit%02d' % unitNo, str(k.get_label())):
                k.set_color(mplColor)
        self.ChanTab['WavesFigure'].figure.canvas.draw()

        # change the color in the features plot
        for k in self.ChanTab['FeaturesFig'].figure.axes[0].lines:
            if re.search('Unit%02d' % unitNo, str(k.get_label())):
                k.set_color(mplColor)
        self.ChanTab['FeaturesFig'].figure.canvas.draw()

        self.ChanTab['UnitTabsWidget'].tabBar().setTabTextColor(unitNo, qtColor)
        
    ########################################################################################################
        
    def ResetChannelTab_proc(self):
        ''' reset the units tab'''

        self.NUnits = 0
                
        # clear the unit figures
        for k in self.ChanTab['UnitFigures']:
            plt.close(self.ChanTab['UnitFigures'][k].figure)
        self.ChanTab['UnitFigures'] = {}

        # clean the button dictionaries
        for key in ['DelUnitBtns', 'UnitCountLabel', 'UnitBtns', 'PlotRawCheck', 'isMultiunitCheck']:
            for k in self.ChanTab[key].keys():
                self.ChanTab[key][k].deleteLater()
            self.ChanTab[key] = {}
        
        # Reset WavesFigure canvas
        ax = self.ChanTab['WavesFigure'].figure.axes[0]
        ax.cla()
        self.SampleWaveform, = ax.plot([], color= [.5,.5,.5], lw=2, animated=True)
        ax.set_ylim(-1000, 1000)
        ax.set_xlim(0, self.WfSize)
        ax.tick_params(color = [.5, .5, .5], labelcolor=[.5, .5, .5])
        for k in ax.spines.values(): k.set_edgecolor([.5, .5, .5])        
        self.Slice1Ln = ax.axvline(0, color=[.5, .5, .5])
        self.Slice2Ln = ax.axvline(0, color=[.5, .5, .5], linestyle = '--')
        ax.grid(color = [.5, .5, .5])
        self.ChanTab['WavesFigure'].figure.tight_layout()
        self.ChanTab['WavesFigure'].figure.canvas.draw()

        # clean the 3d widget
        self.Fig3d.clf()
        
        # set the current indexes of the X and Y variable-selecting comboboxes
        self.XPlot.setCurrentIndex(0)
        self.YPlot.setCurrentIndex(1)
        self.ZPlot.setCurrentIndex(2)

        # reset Units list
        self.UnitsList = []

        # reset the units tabbed widget
        tabs = range(self.ChanTab['UnitTabsWidget'].count())
        tabs.reverse()
        if len(tabs)>0:
            for k in tabs:
                self.ChanTab['UnitTabsWidget'].removeTab(k)
            
        # reset the time scroll widget and axes
        self.TimeScroll['VZoom'].setValue(1000)
        self.TimeScroll['HZoom'].setValue(500)
        self.TimeScroll['HScroll'].setValue(0)
        self.TimeScroll['Figure'].figure.axes[0].cla()
        self.TimeScroll['Ax'].set_xticklabels([])
        self.TimeScroll['Ax'].set_yticklabels([])
        self.TimeScroll['Figure'].figure.canvas.draw()

        # reset label
        self.nPtsLabel.setText('')

        # reset the features figure:
        self.ChanTab['FeaturesFig'].figure.clf()
        self.ChanTab['FeaturesFig'].figure.canvas.draw()

        # delete KDTree object
        if hasattr(self, 'XYData'):       del self.XYData
        if hasattr(self, 'CurWaveforms'): del self.CurWaveforms
        if hasattr(self, 'CurTs'):        del self.CurTs

        # remove the PCA from the dictionarys
        self.ChanTab.pop('PCA',0)

        # reset the channel tab name
        self.MainFigTab.setTabText(2,'Channel Tab')

    ########################################################################################################
        
    def SliceDraw(self):
        sender = self.sender()
        if sender.objectName() == 'Slice1':
            self.Slice1Ln.set_xdata(sender.value())
            self.ChanTab['WavesFigure'].figure.axes[0].draw_artist(self.Slice1Ln)
        elif sender.objectName() == 'Slice2':
            self.Slice2Ln.set_xdata(sender.value())
            self.ChanTab['WavesFigure'].figure.axes[0].draw_artist(self.Slice2Ln)

    ########################################################################################################
        
    def ChangeCurrentUnit_proc(self):
        '''set the current unit'''
        self.CurUnit     = self.ChanTab['UnitTabsWidget'].currentIndex()
        self.CurUnitName = str(self.ChanTab['UnitTabsWidget'].tabText(self.CurUnit))
        for k in range(self.What2Plot.count()):
            if str(self.What2Plot.itemText(k))==self.CurUnitName:
                self.What2Plot.setCurrentIndex(k)
                break

    ########################################################################################################
            
    def MainFigTabProc(self):
        '''Change the toolbar tab acording to the selected view'''
        curtab = self.MainFigTab.currentIndex()
        curtabname = str(self.MainFigTab.tabText(curtab))
        if curtabname == 'Channels Overview' or curtabname == 'Summary Table':
            self.ToolsTab.setCurrentIndex(0)
        elif re.search('Chan [0-9]{1,2}',curtabname):
            self.ToolsTab.setCurrentIndex(1)

    ########################################################################################################
            
    def closeEvent(self,  *event):
        ''' reimplementation of the closeEvent that closes the h5file before killing the window'''
        if self.H5FileLoaded == True:
            self.h5file.close()
        self.deleteLater()

    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

if __name__ == '__main__':
    if not QtGui.QApplication.instance():
        app = QtGui.QApplication(sys.argv)
    else:
        app = QtGui.QApplication.instance()
    ss = pySpikeSorter()
    ss.show()
    app.exec_()

############################################################################################################
