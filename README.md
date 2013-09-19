pySpikeSorter
=============

A pure python Spike sorter.
A(nother) spike sorting utility capable of helping in the analisis of extracelullar action potential recordings.

Installation Instructions:
--------------------------

Requirements:
-------------
* Numpy
* Scipy
* Matplotlib
* PyQt4
* Guidata
* traits
* traitsui
* mayavi
* pytables
* optional: Klustakwik

In Debian based linux distributions i.e. **ubuntu** you can do:

    sudo apt-get install python-numpy python-scipy python-tables python-guidata

In **Windows** it is recomemnded to install a python bundle, for example [PythonXY](https://code.google.com/p/pythonxy/).

Under **MacOS** it is recommended to install the [Enthought](https://www.enthought.com/products/canopy/) package, which bundles all the requirements

There is no installation needed. Just download the files and from the command line run:

    cd /folder/that/contains/the/files
    ipython -i pySpikeSorter.py

Notes:
*It is not required to use ipython, but is more convenient if one wants to explore and intereact with the Spyke Sorter object.
*If you are going to use ipython, please make sure that it doesn't load any PyQt4 library in advance, for example, calling:

    ipython --pylab=qt -ipySpikeSorter.py

  will not work because in order to use MayaVi for the 3d visualizatons, there are some tweaks needed to be made before loading PyQt4


Tutorial:
Video Tutorial in progress

H5File Definition:
------------------
pySpikeSorter allows to sort spike events saved in an H5File.
The package [Neo](http://pythonhosted.org/neo/), makes possible to read data in various formats and transform those into an h5file.
The excellent package [Pytables](http://www.pytables.org/moin) provides a convenient way to create and manipulate these files.

To facilitate the process, a sample h5file with the basic structure is provided.
The h5file basic structure is as follows:
* A group or folder called 'Spikes' (capital 'S')
* Inside, one grou/folder named 'Chan_000', where the zeros are replaced by the channel number. Note that the script asumes there are three integers.
* Inside each channel folder a:
    * 'Waveforms' numpy array, with the shape (events, points)
    * 'Unsorted': numpy array containing the indices of the unsorted events.
    * 'isMultiunit': boolean
    * 'TimeStamp': a numpy vector containing the timestamps (in milliseconds) of the neural events.
    * 'isTrash': Boolean.
