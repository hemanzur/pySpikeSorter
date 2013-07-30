from PyQt4 import QtGui
import tables, os, nev, pickle
from glob import glob

##########################################################################################

def bin2h5(pth = None):
    '''Helper function to transform binary files to an h5 file.'''

    # if not path provided open a dialog to select one
    if not pth:
        pth = QtGui.QFileDialog.getExistingDirectory(parent = None)

    if not pth: return

    # read the list of filenames of fragments
    pth = str(pth)
    files = glob(os.path.join(pth,'channel*'))
    files.sort()
    title = os.path.split(pth)[1]

    # create a new h5 file
    filename = os.path.join(pth, title)+'.h5'
    h5file   = tables.openFile(filename, mode = 'w', title = title)

    # load the pickled headers:
    tmp = pickle.load(open(os.path.join(pth,'headers.p'),'rb'))
    bas_header = tmp[0]
    ext_header = tmp[1]

    # read some parameters from the headers
    channel_info_dict         = ext_header['neural event waveform'][1]
    bytes_in_data_packets     = bas_header['bytes in data packets']
    bytes_per_waveform_sample = channel_info_dict['bytes per waveform sample']
    Fs = float(bas_header['time stamp resolution Hz'])
    Ts = 1.0/Fs

    # set parameter to read the fragments
    waveform_format = 'h'
    waveform_size   = (bytes_in_data_packets - 8)/bytes_per_waveform_sample

    # create and display a progression bar
    pd = QtGui.QProgressDialog('Processing Files', 'Cancel', 0, len(files))
    pd.setWindowTitle('Processing Files ...')
    pd.setGeometry(500,500,500,100)
    pd.show()

    # create basic structure inside the h5file
    h5file.createGroup('/','Header')
    h5file.createArray('/Header', 'WaveformSize', waveform_size)
    h5file.createArray('/Header', 'NChans', len(files))
    h5file.createArray('/Header', 'TimeStamp_Res', Fs)
    h5file.createArray('/Header', 'Date', np.array(bas_header['time origin']))
    h5file.createGroup('/','Spikes')

    # read the non neural data
    Timestamps, code = nev.read_frag_nonneural_digital(pth, bas_header)
    if np.any(code):
        binCode = np.int8([ map(int, np.binary_repr(k, width = 16)) for k in code], ndmin=2)
        tmp = np.ones( shape = (1, binCode.shape[1]), dtype = np.int8)
        tmp = np.append(tmp, binCode, axis = 0)
        dx = np.diff(tmp, n=1, axis=0)
        L = dx.shape[1]
        ton  = []
        toff = []
        for k in range(1, L+1):
            ton.append(Timestamps[dx[:,-k]<0])
            toff.append(Timestamps[dx[:,-k]>0])

        # create the group and leaves to hold the non neural data
        h5file.createGroup('/','Non_Neural_Events')
        h5file.createGroup('/Non_Neural_Events', 'ton')
        h5file.createGroup('/Non_Neural_Events', 'toff')
        for j, k in enumerate(ton):
            if len(k)>0:
                h5file.createArray('/Non_Neural_Events/ton', 'ton_%02d' % j, k)

        for j, k in enumerate(toff):
            if len(k)>0:
                h5file.createArray('/Non_Neural_Events/toff', 'toff_%02d' % j, k)

    
    # iterate over the list of fragments
    for n,f in enumerate(files):

        # animate progression bar
        pd.setLabelText('%s' % f)
        pd.setValue(n+1)

        # open binary fragment
        fid      = open(f,'rb')

        waveform = []
        TS       = []
        #waveform = np.empty((0, waveform_size), dtype = np.int16)
        #TS       = np.empty(0, dtype = np.float)

        while 1:
            # read "bytes_in_data_packets" from the current position
            w = fid.read(bytes_in_data_packets)

            # exit if the number of read bytes is less than expected
            if len(w) < bytes_in_data_packets: break

            # transform the waveform and timestamp data into a number
            waveform.append(np.array(struct.unpack(waveform_size*waveform_format, w[8:]),
                                     dtype = np.int16))
            TS.append(struct.unpack('I', w[0:4])[0] * Ts * 1000)

        # close fragment
        fid.close()
        
        # if the channel doesen't contain waveforms or timestamps:
        if not waveform or not TS:
            continue

        # transform waveforms and timestamps into 16 bits integer array
        waveform = np.array(waveform, dtype=np.int16)
        TS = np.array(TS)
        Unsorted = np.arange(len(TS))

        # create the groups and arrays inside the h5file
        # to host the information
        chName = 'Chan_%03d' % (n+1)
        h5file.createGroup('/Spikes', chName)
        h5file.createArray('/Spikes/'+chName,'Waveforms',waveform)
        h5file.createArray('/Spikes/'+chName,'TimeStamp',TS)
        h5file.createArray('/Spikes/'+chName,'Unsorted',Unsorted)
        h5file.createArray('/Spikes/'+chName,'isMultiunit', False)
        h5file.createArray('/Spikes/'+chName,'isTrash', False)

    ###### ADD THE LFP DATA #####################################################

    nsx.addLFP2h5(h5file)
    
    # close the h5file
    h5file.close()

##########################################################################################
    
def ext_fragments(filename=None, outdir=None):

    fid = open(filename, 'rb')
    bas_header = nev.read_basic_header(fid)
    ext_header = nev.read_extended_header(fid, bas_header)

    fname = os.path.split(filename)[1]
    outdir = os.path.join(outdir,fname[0:fname.find('.')])

    # create a directory if it doesn't exists
    if not os.path.isdir(outdir): os.mkdir(outdir)
    headers = open(os.path.join(outdir,'headers.p'),'wb')
    pickle.dump([bas_header, ext_header], headers)
    headers.close()
    nev.fragment(fid,
                 bas_header,
                 ext_header,
                 channel_list = np.arange(1,65),
                 frag_dir = outdir,
                 ignore_spike_sorting = True)
    fid.close()
    return outdir