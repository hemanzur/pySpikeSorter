########## MERGE UNITS WIDGET #######################################################

class MergeUnitsWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Merge Units')
        hlay = QtGui.QHBoxLayout()
        items = ['casa','perro','arbol','tia','ksaka','blablabal']
        self.list1 = QtGui.QListWidget()
        self.list1.addItems(items)
        self.list1.setMaximumWidth(100)
        self.list2 = QtGui.QListWidget()
        self.list2.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-next')
        btnRight = QtGui.QPushButton(icon,'')
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
       
    def Accept_proc(self):
        units = [self.list2.item(k).text() for k in range(self.list2.count())]
        self.close()
        return units

    def toRight(self):
        if self.list1.currentRow() == -1: return
        item = self.list1.currentItem()
        self.list1.takeItem(self.list1.currentRow())
        self.list2.addItem(item)

    def toLeft(self):
        if self.list2.currentRow() == -1: return
        item = self.list2.currentItem()
        self.list2.takeItem(self.list2.currentRow())
        self.list1.addItem(item)

########## MOVE UNITS WIDGET #######################################################

class MoveUnitsWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setWindowTitle('Move Units')
        hlay = QtGui.QHBoxLayout()
        self.list = QtGui.QListWidget()
        items = ['casa','perro','arbol','tia','ksaka','blablabal']
        self.list.addItems(items)
        self.list.setMaximumWidth(100)
        icon = QtGui.QIcon.fromTheme('go-up')
        btnUp = QtGui.QPushButton(icon,'')
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
       
    def Accept_proc(self):
        units = [self.list.item(k).text() for k in range(self.list.count())]
        self.close()
        return units

    def toUp(self):
        if self.list.currentRow() in [-1,0]: return
        place = self.list.currentRow()
        item = self.list.takeItem(place)
        self.list.insertItem(place-1, item)
        self.list.setCurrentRow(place-1)

    def toDown(self):
        if self.list.currentRow() == -1: return
        if self.list.count()-1 == self.list.currentRow(): return
        place = self.list.currentRow()
        item = self.list.takeItem(place)
        self.list.insertItem(place+1, item)
        self.list.setCurrentRow(place+1)
        
###############################################################################
        
class Settings(dt.DataSet):
    WorkingDir   = di.DirectoryItem('Select a Working Dir', default = pth)
    FiguresDir   = di.DirectoryItem('Path to save images', default = pth)
    Figurescolor = di.ColorItem('Fig color', default='black')
    AxesColor    = di.ColorItem('Axes color', default='gray').set_pos(col=1)
    LassoColor   = di.ColorItem('Lasso color', default='gray')
    DensityCM    = di.ChoiceItem('Density Color Map',
                                 tuple(colormaps),
                                 default = colormaps.index('gist_heat')).set_pos(col=1)

###############################################################################
                                 
class NevFilesPth(dt.DataSet):
    filename = di.FileOpenItem ('Select a NEV files', formats=['nev'])
    outdir   = di.DirectoryItem('Select a dir to output the files')

###############################################################################
    
class AutocorrOpts(dt.DataSet):
    Mode        = di.ChoiceItem('Mode', ['ephys','fft','time'])
    BinSize     = di.IntItem('BinSize', default = 1).set_pos(col=1)
    Win0         = di.IntItem('Win 0', default = -150)
    Win1         = di.IntItem('Win 1', default = 150).set_pos(col=1)

###############################################################################
    
class ReadParams(dt.DataSet):
    filename = di.FileOpenItem('NEV file', formats=['nev'])
    outdir   = di.DirectoryItem('Output Dir')