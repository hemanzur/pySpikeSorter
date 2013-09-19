# this is to force traits to use pyqt4 as backend
import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'

import sip
sip.setapi('QString', 2)
sip.setapi('QVariant', 2)

# then import PyQt4

from PyQt4 import QtGui

from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor


################################################################################
#The actual visualization class
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        pass

    # the layout of the dialog created
    view = View(Item('scene',
                      editor = SceneEditor(scene_class=MayaviScene),
                      height = 100, width = 100, show_label = False),
                resizable = True)

################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.
class MayaviQWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        layout = QtGui.QVBoxLayout(self)
        layout.setMargin(0)
        layout.setSpacing(1)
        self.visualization = Visualization()
        
        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.setLayout(layout)
        #self.ui.setParent(self)