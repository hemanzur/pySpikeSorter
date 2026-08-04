# pySpikeSorter
* A(nother) spike sorting utility capable of helping in the analisis of extracelullar action potential recordings.
* It is inspired by some of the features of the offline sorter software from plexon.
* It is intended to be a completely "Manual" solution for spike sorting (old school neurophysiologists don't like automated stuf,
 although there is some work to provide and/or interface with automatic clustering algorithms.
* It is far from being complete, but at this moment it is functional.
* I'm not a programmer, so I think there are lots of things that can be improved, that is why any comments aiming at improving it are welcome.

## Installation Instructions:  

### Requirements:
* Python 3
* NumPy
* SciPy
* Matplotlib
* PyQt5
* PyOpenGL
* pyqtgraph
* PyTables
* scikit-learn
* optional: Klustakwik

To meet all the requirements, create a virtual environment and install the
Python dependencies:

```
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

### Running pySpikeSorter:

There is no installation needed. Just download the files to a known location, and from the command line run:

```
cd /folder/that/contains/the/files
python pySpikeSorter.py
```

#### Notes:
* It is not required to use IPython, but it can still be convenient if one wants to explore and interact with the Spike Sorter object.


## Tutorial:
* A "test.h5" file is provided, which can be loaded in the appplication clicking on the green button.
* It should display an overview of the channels
* By clicking a particular channel it gets selected. Clicking the "Plot Chan" button plots the data in a specific Tab.
* Under that tab there are aeveral analisys and tools that can be used to sort the data.

### Video Tutorial:
in progress ...

## Screenshots:
[![Screenshot 1](screenshots/thumbnails/03.jpg)](screenshots/03.png)
[![Screenshot 2](screenshots/thumbnails/04.jpg)](screenshots/04.png)
[![Screenshot 3](screenshots/thumbnails/06.jpg)](screenshots/06.png)
[![Screenshot 4](screenshots/thumbnails/12.jpg)](screenshots/12.png)
[![Screenshot 5](screenshots/thumbnails/13.jpg)](screenshots/13.png)

##H5File Definition:
pySpikeSorter allows to sort spike events saved in an H5File.
The package [Neo](http://pythonhosted.org/neo/), makes possible to read data in various formats and transform those into an h5file.
The excellent package [Pytables](http://www.pytables.org/moin) provides a convenient way to create and manipulate these files.

To facilitate the process, a sample h5file with the basic structure is provided. You can examine the internal structure of the files
using the [viTables](http://vitables.org/) application.
The h5file basic structure is as follows:
* A group or folder called 'Spikes' (capital 'S')
* Inside, one grou/folder named 'Chan_000', where the zeros are replaced by the channel number. Note that the script asumes there are three integers.
* Inside each channel folder a:
    * 'Waveforms' numpy array, with the shape (events, points)
    * 'Unsorted': numpy array containing the indices of the unsorted events.
    * 'isMultiunit': boolean
    * 'TimeStamp': a numpy vector containing the timestamps (in milliseconds) of the neural events.
    * 'isTrash': Boolean.
    
## To do List:
* Tetrodes implementation.
