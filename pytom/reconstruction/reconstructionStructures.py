'''
started on Mar 10, 2011

@author: luiskuhn, foerster
'''

from pytom.basic.structures import PyTomClass
from pytom.agnostic.tools import convert_operation_order_str2list, convert_operation_order_list2str

# TODO what should be the default operation order? => default order seems to be [0, 1, 2]
# TODO this was in previous versions already available as parameter to basic.transformations -> general_transform()


def init_pool_processes(shared_stack_):
    """
    Make a shared array inheritable for all the workers in multiprocessing Pool.
    @param shared_stack_:
    @type shared_stack_:
    @return:
    @rtype:
    """
    global shared_stack
    shared_stack = shared_stack_


def put_in_shared_stack(stack_dtype, stack_shape, image, index):
    """
    For accessing shared arrays in python multiprocessing.

    @param stack_dtype:
    @type stack_dtype:
    @param stack_shape:
    @type stack_shape:
    @param image:
    @type image:
    @param index:
    @type index:
    @return:
    @rtype:
    """
    import numpy as np
    from pytom_numpy import vol2npy
    # create a np view to global shared_stack
    shared_stack_np_ = np.asfortranarray(np.frombuffer(shared_stack.get_obj(),
                                                       dtype=stack_dtype).reshape(stack_shape))
    shared_stack_np_[:, :, index] = vol2npy(image).copy(order='F')


def align_single_projection(index, projection_list, tilt_angles, shared_dtype, shared_stack_shape, weighting, binning,
                            low_pass_freq, apply_circle_filter, scale_factor_particle, angle_specimen, verbose, imdim,
                            imdimX, imdimY, sliceWidth):
    """
    For parallel alignment of tilt series.

    @param index: index of a projection in the list
    @type index: L{int}
    @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
    @type weighting: L{int}
    @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
    pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
    @type binning: L{int}
    @param low_pass_freq: lowpass filter (in Nyquist)
    @type low_pass_freq: L{float}
    @param apply_circle_filter: optional circular filter in fourier space
    @type circleFilter: L{bool}
    @param scale_factor_particle: scaling factor for particles
    @type scale_factor_particle: L{float}
    @param angle_specimen: additional rotation of the specimen
    @type angle_specimen: L{float}
    @param verbose: print output default=False
    @type verbose: L{bool}
    @param imdim: x and y dimension of final image
    @type imdim: L{int}
    @param imdimX: input projection x dimension size
    @type imdimX: L{int}
    @param imdimY: input projection y dimension size
    @type imdimY: L{int}

    @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
    and an offsetStack (user specified x, y offset)
    @rtype: L{list}

    @author: Marten Chaillet
    """
    from pytom_volume import complexRealMult, vol, pasteCenter, mean
    from pytom.basic.functions import taper_edges
    from pytom.basic.transformations import general_transform2d, resize
    from pytom.basic.fourier import ifft, fft
    from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
        fourierFilterShift_ReducedComplex
    from pytom.basic.files import read
    import pytom_freqweight

    projection = projection_list[index]
    tilt_angle = tilt_angles[index]

    # if verbose:
    print(f'start on projection {index}')

    # pre-determine analytical weighting function and lowpass for speedup
    if weighting == -1:
        weightSlice = fourierFilterShift(rampFilter(imdim, imdim))
    elif weighting == 1:
        weightSlice = fourierFilterShift(exactFilter(tilt_angles, tilt_angle,
                                                     imdim, imdim, sliceWidth))

    if apply_circle_filter:
        circleFilterRadius = imdim // 2
        circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
    else:
        circleSlice = vol(imdim, imdim // 2 + 1, 1)
        circleSlice.setAll(1.0)

    # design lowpass filter
    if low_pass_freq > 0:
        lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                      low_pass_freq / 5. * imdim)

    if projection._filename.split('.')[-1] == 'st':
        image = read(file=projection.getFilename(),
                     subregion=[0, 0, index - 1, imdim, imdim, 1],
                     sampling=[0, 0, 0], binning=[0, 0, 0])
        if not (binning == 1):
            image = resize(volume=image, factor=1 / float(binning))[0]
    else:
        image = read(projection.getFilename())
        # image = rotate(image,180.,0.,0.)
        image = resize(volume=image, factor=1 / float(binning))[0]

    # normalize to contrast - subtract mean and norm to mean
    immean = mean(image)
    image = (image - immean) / immean

    # smoothen borders to prevent high contrast oscillations
    image = taper_edges(image, imdim // 30)[0]

    # transform projection according to tilt alignment
    transX = projection.getAlignmentTransX() / binning
    transY = projection.getAlignmentTransY() / binning
    rot = projection.getAlignmentRotation()
    mag = projection.getAlignmentMagnification() * scale_factor_particle
    order = projection.getOperationOrder()

    # 3 -- square if needed
    if imdimY != imdimX:
        newImage = vol(imdim, imdim, 1)
        newImage.setAll(0)
        pasteCenter(image, newImage)
        image = newImage

    # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
    # so IMOD default order is scaling > rotation > translation, that seems a weird order...
    center = (image.sizeX() / 2 - 0.5, image.sizeY() / 2 - 0.5, 0) if order == [1, 2, 0] else None

    if not (rot == 0 and transX == 0 and transY == 0 and mag == 0):
        image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order,
                                    crop=True, center=center)

    if low_pass_freq > 0:
        filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
        image = filtered[0]

    # smoothen once more to avoid edges
    image = taper_edges(image, imdim // 30)[0]

    # apply weighting
    if weighting in [-1, 1]:
        image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

    if verbose:
        print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

    # assign in self._shared_stack_np
    put_in_shared_stack(shared_dtype, shared_stack_shape, image, index)

    print(f'process on projection {index} finished')

    return index, (projection.getAlignmentAxisAngle(), projection.getTiltAngle() - angle_specimen,
                   (projection.getOffsetX(), projection.getOffsetY()))


class Projection(PyTomClass):
    """
    Projection: 
    """
    def __init__(self, filename=None, tiltAngle=None, offsetX=0, offsetY=0, alignmentTransX=0.,
                 alignmentTransY=0., alignmentRotation=0., alignmentMagnification=1.,
                 alignmentAxisAngle=0., operationOrder=[0, 1, 2], index=0):
        """
        Initialize projection class

        @param filename: The filename
        @type filename: string
        @param tiltAngle: Tilt angle of a projection, rotation around microscope y-axis
        @type tiltAngle: float
        @param offsetX: offset X is used for subtomogram reconstruction from cut out projections
        @type offsetX: float
        @param offsetY: offset Y is used for subtomogram reconstruction from cut out projections
        @type offsetY: float
        @param alignmentTransX: X Shift determined from alignment method
        @type alignmentTransX: float
        @param alignmentTransY: Y Shift determined from alignment method
        @type alignmentTransY: float
        @param alignmentRotation: in plane rotation from alignment, around microscope z-axis
        @type alignmentRotation: float
        @param alignmentMagnification: Magnification determined from fine alignment method
        @type alignmentMagnification: float
        @param alignmentAxisAngle: rotation angle around microscope system x-axis
        @type alignmentAxisAngle: float
        @param operationOrder: tuple of alignment operation order (r, t, s) => (2, 0, 1) means translation first,
        then scaling, then rotation
        @type operationOrder: list
        @param index: Index of projection in filename (default: 1)
        @type index: int
        """

        self._filename = filename
        self._checkValidFile()

        self._xSize = None
        self._ySize = None

        if filename and tiltAngle is None:
            self._loadTiltAngle()
        else:
            self._tiltAngle = tiltAngle

        self.setIndex(index)
        self._offsetX = offsetX
        self._offsetY = offsetY
        self._alignmentTransX = alignmentTransX
        self._alignmentTransY = alignmentTransY
        self._alignmentRotation = alignmentRotation
        self._alignmentMagnification = alignmentMagnification
        self._alignmentAxisAngle = alignmentAxisAngle
        if isinstance(operationOrder, str):
            self._operationOrder = convert_operation_order_str2list(operationOrder)
        else:
            self._operationOrder = operationOrder  # default = (2, 0, 1) translation first, then magnification,
                                                   # then rotation

    def info(self):
        tline = (" %3d: " %self.getIndex())
        tline = tline + ("Filename: %20s" %self._filename)
        tline = tline + (", TiltAngle= %6.2f" %self._tiltAngle)
        return tline
        
    def _loadTiltAngle(self):
        """
        _loadTiltAngle: Loads projection EM file and parses header for tilt angle. Assigns it to the object. Works for EM files only!
        """
        import struct
        from pytom.gui.mrcOperations import read_mrc_header

        try:
            if self._filename.split('.')[-1] == 'em':
                fh = open(self._filename,'rb')
                fh.seek(168)
                self._tiltAngle = float(struct.unpack('i',fh.read(4))[0])/1000
                fh.close()
            else:
                fh = read_mrc_header(self._filename)
                self._tiltAngle = fh['user'][0]/1000.  # TODO find correct location of the tiltangle and autowrite it.
        except:
            print(f'could not read tilt angle for projection {self._filename}')
            self._tiltAngle = None

    def _checkValidFile(self):
        """
        Class private method that is called whenever a Projection() is initialized or the filename is set,
        it will warn the user if the current filename is valid
        """
        from pytom.tools.files import checkFileExists
        if not checkFileExists(self._filename):
            print('Warning: Projection() initialized without valid filename')

    def getDimensions(self):
        """
        get dimensions of projection

        @rtype: array (3-dim)
        """
        if not self._filename:
            "Projection Filename not set :("
            return None
        from pytom.basic.files import readProxy as read
        
        proj = read(self._filename)
        return proj.sizeX(), proj.sizeY(), proj.sizeZ()

    def getFilename(self):
        return self._filename
        
    def getTiltAngle(self):
        return self._tiltAngle

    def getOffsetX(self):
        return self._offsetX
    
    def getOffsetY(self):
        return self._offsetY
    
    def getAlignmentTransX(self):
        return self._alignmentTransX
    
    def getAlignmentTransY(self):
        return self._alignmentTransY
    
    def getAlignmentRotation(self):
        return self._alignmentRotation
    
    def getAlignmentMagnification(self):
        return self._alignmentMagnification

    def getAlignedFilename(self):
        """
        get filename of aligned projection

        @rtype: string
        """
        return self._alignedFilename

    def setIndex(self, index):
        """
        set index of projection

        @param index: index of projection (as indicated in filename)
        @type index: int
        """
        self._index = index

    def getIndex(self):
        """
        get index of projection

        @rtype: int
        """
        return self._index

    def getAlignmentAxisAngle(self):
        """
        tilt around microscope X axis

        @rtype: float
        """
        return self._alignmentAxisAngle

    def getOperationOrder(self, as_string=False):
        """
        get order of applying operations

        @rtype: tuple -> (r, t, s)
        """
        if as_string:
            return convert_operation_order_list2str(self._operationOrder)
        else:
            return self._operationOrder
    
    def returnProjectionImage(self,rescale=1.):
        """
        @param rescale: rescale image by factor
	    @type rescale: int
        """
        from pytom.basic.files import readProxy as read
        return read(self.getFilename(),0,0,0,0,0,0,0,0,0,rescale,rescale,1)
    
    def setFilename(self, filename):
        """
        set filename of projection

        @param filename: Filename
        @type filename: float
        """
        self._filename = filename
        self._checkValidFile()
    
    def setTiltAngle(self, tiltAngle):
        """
        set tilt angle in projection

        @param tiltAngle: tilt angle (deg)
        @type tiltAngle: float
        """
        self._tiltAngle = tiltAngle
    
    def setAlignmentTransX(self,value):
        """
        set x-translation in projection

        @param value: translation X
        @type value: float
        """
        self._alignmentTransX = float(value)
    
    def setAlignmentTransY(self,value):
        """
        set y-translation in projection

        @param value: translation Y
        @type value: float
        """
        self._alignmentTransY = float(value)
    
    def setAlignmentRotation(self,value):
        """
        set Rotation value in projection

        @param value: rotation angle (deg)
        @type value: float
        """
        self._alignmentRotation = float(value)
    
    def setAlignmentMagnification(self,value):
        """
        set Magnification value in projection

        @param value: magnification
        @type value: float
        """
        self._alignmentMagnification = float(value)

    def setAlignedFilename(self, alignedFilename):
        """
        set filename of aligned projection

        @type alignedFilename: string
        """
        self._alignedFilename = alignedFilename

    def setAlignmentAxisAngle(self, alignmentAxisAngle):
        """
        set tilt around microscope X axis

        @type axisAngle: float
        """
        self._alignmentAxisAngle = alignmentAxisAngle

    def setOperationOrder(self, operationOrder):
        """
        set operation order for aligning tilts as a list

        @type operationOrder: list -> [r, t, s]
        """
        if isinstance(operationOrder, str):
            self._operationOrder = convert_operation_order_str2list(operationOrder)
        else:
            self._operationOrder = operationOrder

    def _setSize(self):
        """
        set X- and Y-size of proj according to file
        """
        from pytom.basic.files import readProxy as read
        p = read(self._filename)
        self._xSize = p.sizeX()
        self._ySize = p.sizeY()
        
    def getXSize(self):
        """
        get X-dimension of proj

        @return: x-dim
        @rtype: float
        """
        if not self._xSize:
            self._setSize()
            
        return self._xSize
    
    def getYSize(self):
        """
        get Y-dimension of proj

        @return: y-dim
        @rtype: float
        """
        if not self._ySize:
            self._setSize()
            
        return self._ySize
    
    def toXML(self):
        """
        print projection metadata into xml
        """
        from lxml import etree
        projection_element = etree.Element('Projection', Filename=str(self._filename),
                                           TiltAngle=str(float(self._tiltAngle)),OffsetX=str(float(self._offsetX)),
                                           OffsetY=str(float(self._offsetY)),
                                           AlignmentTransX=str(float(self._alignmentTransX)),
                                           AlignmentTransY=str(float(self._alignmentTransY)),
                                           AlignmentAngle=str(float(self._alignmentRotation)),
                                           AlignmentMagnification=str(float(self._alignmentMagnification)))
        return projection_element        
        
    def fromXML(self,xmlObj):
        """
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'Projection':
            projection_element = xmlObj
        else:
            Exception('Is not a Projection! You must provide a valid Projection object.')
        
        self._filename = projection_element.get('Filename')
        self._checkValidFile()
        self._tiltAngle = float(projection_element.get('TiltAngle'))
        self._offsetX = float(projection_element.get('OffsetX'))
        self._offsetY = float(projection_element.get('OffsetY'))
        self._alignmentMagnification = float(projection_element.get('AlignmentMagnification'))
        
        
        if projection_element.get('AlignmentTransX') == [] or projection_element.get('AlignmentTransX') == None:
            self._alignmentTransX = float(projection_element.get('AlignmentOffsetX'))
        else: 
            self._alignmentTransX = float(projection_element.get('AlignmentTransX'))
            
        if projection_element.get('AlignmentTransY') == [] or projection_element.get('AlignmentTransY') == None:
            self._alignmentTransY = float(projection_element.get('AlignmentOffsetY'))
        else: 
            self._alignmentTransY = float(projection_element.get('AlignmentTransY'))
            
        self._alignmentRotation = float(projection_element.get('AlignmentAngle'))
        
        
class ProjectionList(PyTomClass):
    """
    ProjectionList: List of projections
    """
    def __init__(self, projection_list=None):
        self._list = projection_list._list if projection_list is not None else []
        self._tilt_angles = projection_list._tilt_angles if projection_list is not None else []

    def _init_tilt_angles(self):
        self._tilt_angles = [p.getTiltAngle for p in self._list]

    def info(self):
        tline = ""
        for projection in self._list:
            tline = tline + str(projection) + "\n"
        return tline

    def loadDirectory(self, directory, metafile='', prefix=''):
        """
        loadDirectory: Will take all projections in a directory, determine their respective tilt angle from the
        header and populate this list object.

        @param directory: directory where files are located
        @type directory: L{str}
        @param metafile: metafile with information about tilt series
        @type metafile: L{str}
        @param prefix: prefix of projection files in dir, will be used to filter files
        @type prefix: L{str}
        """
        from pytom.tools.files import checkDirExists
        from pytom.basic.files import loadstar
        from pytom.basic.datatypes import DATATYPE_METAFILE
        import os

        if metafile:
            metadata = loadstar(metafile, dtype=DATATYPE_METAFILE)
            tiltAngles = metadata['TiltAngle']

        if not checkDirExists(directory):
            raise RuntimeError('Directory ' + directory + ' does not exist!')

        files = [line for line in sorted(os.listdir(directory)) if prefix in line]

        # if metafile is provided we could just append directory with filename in metafile ...
        # metafile could also contain the full path

        for file in files:
            if os.path.splitext(file)[1] in ['.em', '.mrc']:
                i = int(file.split('_')[-1].split('.')[0])  # TODO this is hackish => but leave it for now
                if metafile:
                    projection = Projection(filename=os.path.join(directory, file), tiltAngle=tiltAngles[i], index=i)
                else:
                    projection = Projection(filename=os.path.join(directory, file), index=i)
                self.append(projection)

        if metafile or not any([p.getTiltAngle() is None for p in self._list]):  # with metafile we have the
            # tilts, otherwise tilts might have been readable from the em/mrc files but we will need to them all
            self.sort()
        else:
            self.sort(key=lambda p: p.getIndex())

    def load_alignment(self, alignment_file):
        """
        Load the alignment parameters for the projections from an alignment results file.
        If the Projections have not yet been initialized the filename from the alignment results file will be used
        to initialize the projections.

        @param alignment_file: alignment file txt
        @type alignment_file: L{str}
        """
        from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO, DATATYPE_ALIGNMENT_RESULTS
        from pytom.basic.files import loadstar

        try:
            alignment_results = loadstar(alignment_file, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
        except:
            alignment_results = loadstar(alignment_file, dtype=DATATYPE_ALIGNMENT_RESULTS)

        if len(self) == 0:
            for i in range(len(alignment_results)):
                # set offsetx and offsety ??
                projection = Projection(filename=alignment_results[i]['FileName'],
                                        tiltAngle=alignment_results[i]['TiltAngle'],
                                        alignmentTransX=alignment_results[i]['AlignmentTransX'],
                                        alignmentTransY=alignment_results[i]['AlignmentTransY'],
                                        alignmentRotation=alignment_results[i]['InPlaneRotation'],
                                        alignmentMagnification=alignment_results[i]['Magnification'],
                                        index=i)
                if 'OperationOrder' in alignment_results.dtype.names:
                    projection.setOperationOrder(alignment_results[i]['OperationOrder'])
                if 'AxisAngle' in alignment_results.dtype.names:
                    projection.setAlignmentAxisAngle(alignment_results[i]['AxisAngle'])
                self.append(projection)
        else:
            assert len(self) == len(alignment_results), print('Current number of projections in ProjectionList does '
                                                              'not correspond with the amount of projections in '
                                                              'the alignment parameters file. Exiting.')
            for i in range(len(self)):
                self[i].setTiltAngle(alignment_results[i]['TiltAngle'])
                self[i].setAlignmentTransX(alignment_results[i]['AlignmentTransX'])
                self[i].setAlignmentTransY(alignment_results[i]['AlignmentTransY'])
                self[i].setAlignmentRotation(alignment_results[i]['InPlaneRotation'])
                self[i].setAlignmentMagnification(alignment_results[i]['Magnification'])
                if 'OperationOrder' in alignment_results.dtype.names:
                    self[i].setOperationOrder(alignment_results[i]['OperationOrder'])
                if 'AxisAngle' in alignment_results.dtype.names:
                    self[i].setAlignmentAxisAngle(alignment_results[i]['AxisAngle'])

        # sort by tilt angle
        self.sort()

    def toXML(self):
        """
        """
        from lxml import etree

        projectionList_element = etree.Element('ProjectionList')

        ##self._list = sorted(self._list, key=lambda Projection: Projection._tiltAngle)

        for projection in self._list:
            projectionList_element.append(projection.toXML())

        return projectionList_element

    def fromXML(self, xmlObj):
        """
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'ProjectionList':
            projectionList_element = xmlObj
        else:
            Exception('Is not a ProjectionList! You must provide a valid ProjectionList object.')

        for projection in projectionList_element:
            newProjection = Projection()
            newProjection.fromXML(projection)
            self.append(newProjection)

        # sort the projections by tilt angle
        if not any([p.getTiltAngle() is None for p in self._list]):
            self.sort()
            self._init_tilt_angles()
        else:
            self.sort(key=lambda p: p.getIndex())

    def sort(self, key=lambda p: p.getTiltAngle()):
        """
        Sort the ProjectionList based on key.

        @param key: lambda expression to select which item of the Projection to sort on
        @type key: L{function}
        """
        self._list.sort(key=key)
        if not any([p.getTiltAngle() is None for p in self._list]):
            self._init_tilt_angles()
        else:
            print('Warning: ProjectionList could not initialize any or all of the tilt angles')

    def append(self, projection):
        """
        append projection to projection list

        @param projection: a projection
        @type projection: L{pytom.reconstructionStructures.Projection}
        """
        self._list.append(projection)
        
    def __getitem__(self, key):
        """
        __getitem__: Retrieve projection at position definded by key
        @param key:
        @type key: int
        @rtype: L{pytom.reconstructionStructures.Projection}. 
        """
        if not isinstance(key, int):
            raise TypeError('Provide an int')
        
        if len(self) <= key:
            raise IndexError('Index out of range.')
        
        return self._list[key]
        
    def __len__(self):
        return len(self._list)
    
    def generateVolumes(self, particles, cubeSize, binning=1, 
            applyWeighting=False, showProgressBar=False, verbose=False,
            preScale=1, postScale=1):
        """
        generateVolumes (old function name, now replaced by reconstructVolumes)
        @deprecated: 
        """
        print('Generate Volumes is deprecated, use reconstructVolumes instead!')
        self.reconstructVolumes(particles, cubeSize, binning, applyWeighting, showProgressBar,verbose,preScale,
                                postScale)

    def reconstructVolume(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0], read_only=False,
                          binning=1, applyWeighting=False, alignResultFile='', gpu=-1, specimen_angle=0, order=[2,1,0]):
        """ reconstruct a a single 3D volume from weighted and aligned projections on either a gpu or a cpu

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param applyWeighting: apply weighting
        @type applyWeighting: bool
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float

        @return: reconstructed volume
        @rtype: L{pytom_volume.vol} or L{numpy.ndarray}
        """

        from pytom.gpu.initialize import device

        if 'gpu' in device:
            vol = self.reconstructVolumeGPU(dims=dims, reconstructionPosition=reconstructionPosition, binning=binning, read_only=read_only,
                                      applyWeighting=applyWeighting, alignResultFile=alignResultFile, gpu=gpu, specimen_angle=specimen_angle, order=order)

        else:
            vol = self.reconstructVolumeCPU(dims=dims, reconstructionPosition=reconstructionPosition, binning=binning, read_only=read_only,
                                      applyWeighting=applyWeighting, alignResultFile=alignResultFile, gpu=gpu, specimen_angle=specimen_angle, order=order)
        return vol

    def reconstructVolumeGPU(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0], interpolation='filt_bspline',
                          binning=1, applyWeighting=False, alignResultFile='', gpu=-1, specimen_angle=0, read_only=False):
        """
        reconstruct a single 3D volume from weighted and aligned projections on a gpu

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param interpolation: interpolation method used (default='filt_bspline')
        @type reconstructionPosition: str
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param applyWeighting: apply weighting
        @type applyWeighting: bool
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool

        @return: reconstructed volume
        @rtype: L{numpy.ndarray}

        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom_volume import vol, backProject
        from pytom.gpu.initialize import device, xp
        from pytom.reconstruction.reconstructionFunctions import alignImagesUsingAlignmentResultFile

        from pytom.agnostic.reconstruction_functions import backProjectGPU as backProject
        from pytom.agnostic.io import write

        import time
        s = time.time()

        reconstruction = xp.zeros((dims[0], dims[1], dims[2]), dtype=xp.float32)

        # stacks for images, projections angles etc.
        if alignResultFile == '':  # TODO this part does not even work
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.toProjectionStack(
                binning=binning, applyWeighting=applyWeighting, showProgressBar=False,
                verbose=False)
            #arrays = [vol2npy(vol_img).copy(), vol2npy(vol_phi).copy(), vol2npy(vol_the).copy(),vol2npy(vol_offsetProjections).copy() ]
        else:
            arrays = alignImagesUsingAlignmentResultFile(alignResultFile, binning=binning, weighting=applyWeighting,
                                                         circleFilter=True, angle_specimen=specimen_angle)
        [projections, vol_phi, proj_angles, vol_offsetProjections] = arrays



        # volume storing reconstruction offset from center (x,y,z)
        recPosVol = xp.zeros((projections.shape[2], 3), dtype=xp.int32)

        for ii in range(0,3):
            recPosVol[:, ii] = float(reconstructionPosition[ii] // binning)

        print(time.time()-s)
        s = time.time()
        # finally backproject into volume
        vol_bp = backProject(projections, reconstruction, vol_phi, proj_angles, recPosVol, vol_offsetProjections,
                             interpolation).get()

        print(f'backproject time: {time.time()-s}')

        return vol_bp


    def reconstructVolumeCPU(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0], read_only=False,
                          binning=1, applyWeighting=False, alignResultFile='', gpu=-1, specimen_angle=0, order=[2,1,0]):
        """
        reconstruct a single 3D volume from weighted and aligned projections on a CPU

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param applyWeighting: apply weighting
        @type applyWeighting: bool
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float

        @return: reconstructed volume
        @rtype: L{pytom_volume.vol}

        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom_volume import vol, backProject
        import time

        vol_bp = vol(dims[0], dims[1], dims[2])
        vol_bp.setAll(0.0)

        self.angle_specimen = specimen_angle

        # stacks for images, projections angles etc.
        if not alignResultFile or read_only:
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.toProjectionStack(
                binning=binning, applyWeighting=applyWeighting, showProgressBar=False,
                verbose=False, order=order)
        else:
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.toProjectionStackFromAlignmentResultsFile(
                alignResultFile, binning=binning, weighting=applyWeighting, circleFilter=True, order=order)

        # volume storing reconstruction offset from center (x,y,z)
        recPosVol = vol(3, vol_img.sizeZ(), 1)
        recPosVol.setAll(0.0)
        for iproj in range(0, vol_img.sizeZ()):
            for ii in range(0, 3):
                recPosVol.setV(float(reconstructionPosition[ii] / binning), ii, iproj, 0)

        # finally backproject into volume
        print('cpu reconstruction')
        s = time.time()
        backProject(vol_img, vol_bp, vol_phi, vol_the, recPosVol, vol_offsetProjections)
        print(f'backproject time: {time.time()-s}')

        return vol_bp

    def reconstructVolumes(self, particles, cubeSize, binning=1, applyWeighting=False,
                           showProgressBar=False, verbose=False, preScale=1, postScale=1, num_procs=5, ctfcenter=None,
                           alignResultFile='', polishResultFile='', scaleFactorParticle=1, gpuIDs=None, specimen_angle=0):
        """
        reconstructVolumes: reconstruct a subtomogram given a particle object.

        @param particles: A list of particles
        @type particles: L{pytom.basic.structures.ParticleList}
        @param cubeSize: 3D Size of reconstructed particle (here, x == y == z)
        @param binning: binning of projections
        @param applyWeighting: Will apply weighting for weighted backprojection to projection. Default is off (False)!
        @param showProgressBar: False by default
        @param verbose: False by default
        @param preScale: Scale (bin) the projections BEFORE reconstruction
        @param postScale: Scale the tomogram AFTER reconstruction
        @param num_procs: Number of parallel processes
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param polishResultFile: result file generate by the particle polishing (see pytom.basic.datatypes).
        @type polishResultFile: txt
        @param scaleFactorParticle: scaling factor for the reconstruction default=1
        @type scaleFactorParticle: float
        @param gpuIDs: List of gpu-ids on which reconstruction will take place.
        @type gpuIDs: list
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        """

        from pytom_volume import vol, paste, backProject, complexRealMult, rescaleSpline
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift
        import time, os
        from pytom_numpy import vol2npy
        from pytom.basic.files import write_em
        from multiprocessing import Process, set_start_method
        import time

        # important to set to allow child processes on multiple GPU's (otherwise cupy init is not working)
        set_start_method('spawn')

        t = time.time()
        if len(particles) == 0:
            raise RuntimeError('ParticleList is empty!')

        if cubeSize.__class__ != int:
            raise TypeError('reconstructVolumes: Parameter cubeSize must be of type int!')

        if showProgressBar:
            progressBar = FixedProgBar(0, len(particles), 'Particle volumes generated ')
            progressBar.update(0)
            numberParticleVolumes = 0

        # stacks for images, projections angles etc.
        if gpuIDs is None:
            # Running on the CPU

            # Are the aligned images already written out?
            if not alignResultFile:
                # YES ALIGNED IMAGES ARE WRITTEN OUT
                resultProjstack = self.toProjectionStack(binning=binning, applyWeighting=applyWeighting,
                                                         showProgressBar=False,
                                                         verbose=False, num_procs=num_procs, scaleFactorParticle=scaleFactorParticle)
            else:
                # NO ALIGNED IMAGES ARE NOT WRITTEN OUT, THUS HAVE TO BE COMPUTED IN MEMORY FROM ALIGNMENT RESULTS FILE
                # Do you want to align (in memory) in parallel?
                if num_procs > 1:
                    # YES
                    resultProjstack = self.toProjectionStackParallel(alignResultFile, weighting=applyWeighting,
                                                                     num_procs=num_procs, binning=binning,
                                                                     circleFilter=True,
                                                                     scaleFactorParticle=scaleFactorParticle)
                else:
                    # NO
                    resultProjstack = self.toProjectionStackFromAlignmentResultsFile(alignResultFile, weighting=applyWeighting,
                                                                                 num_procs=num_procs, binning=binning,
                                                                                 circleFilter=True, scaleFactorParticle=scaleFactorParticle)
        else:
            # Running on the GPU
            from pytom.reconstruction.reconstructionFunctions import alignImagesUsingAlignmentResultFile as align
            from pytom.agnostic.io import write as write_em

            resultProjstack = align( alignResultFile, binning=binning, weighting=applyWeighting, circleFilter=True, angle_specimen=specimen_angle)
            resultProjstack = [x.get() for x in resultProjstack]

        # if verbose: return

        procs = []
        print(time.time()-t)

        if num_procs:

            for n, result in enumerate(resultProjstack):
                if not os.path.exists(os.path.dirname(particles[0].getFilename())): os.mkdir(
                    os.path.dirname(particles[0].getFilename()))
                if gpuIDs is None: write_em(f'{os.path.dirname(particles[0].getFilename())}/.temp_{n}.mrc', result)
                results = resultProjstack
            if gpuIDs is None:
                results = 1

            extract = self.extract_particles_on_gpu if not gpuIDs is None else self.extract_single_particle

            if gpuIDs is None: gpuIDs=[None,]*num_procs

            for procIndex in range(num_procs):
                ps = particles[procIndex::num_procs]
                if (len(ps)) < 1: continue
                #extract(ps, procIndex, verbose, binning, postScale, cubeSize, polishResultFile, gpuIDs[procIndex], results)
                proc = Process(target=extract, args=(ps, procIndex, verbose, binning, postScale, cubeSize, polishResultFile, gpuIDs[procIndex], results))
                procs.append(proc)
                proc.start()

            # print('Subtomogram reconstruction for particle {} has started (batch mode).'.format(particleIndex))
        else:
            vol_bp = vol(cubeSize, cubeSize, cubeSize)
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = resultProjstack

            for particleIndex in range(len(particles)):
                t = time.time()
                p = particles[particleIndex]
                # self.extract_single_particle(p, particleIndex, verbose, binning, postScale, cubeSize)
                # print('Subtomogram reconstruction for particle {} has started.'.format(particleIndex))
                if showProgressBar:
                    progressBar.update(particleIndex)

                if verbose:
                    print(p)

                vol_bp.setAll(0.0)

                reconstructionPosition = vol(3, vol_img.sizeZ(), 1)
                reconstructionPosition.setAll(0.0)

                # adjust coordinates of subvolumes to binned reconstruction
                for i in range(vol_img.sizeZ()):
                    reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)

                if 1:
                    print((p.getPickPosition().getX() / binning, p.getPickPosition().getY() / binning,
                           p.getPickPosition().getZ() / binning))

                backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

                if postScale > 1:
                    volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                    rescaleSpline(vol_bp, volumeRescaled)
                    volumeRescaled.write(p.getFilename())
                else:
                    vol_bp.write(p.getFilename())
                print(time.time()-t)

        while len(procs):
            time.sleep(0.1)
            procs = [proc for proc in procs if proc.is_alive()]
        # for n, result in enumerate(resultProjstack):

        if gpuIDs is None: os.system('rm -f {}'.format(os.path.dirname(particles[0].getFilename()) + '/.temp_*.em'))
        #progressBar.update(len(particles))
        print('\n Subtomogram reconstructions have finished.\n\n')

    def extract_particles_on_gpu(self, particles, procID, verbose, binning, postScale, cubeSize, filename_ppr='', gpuID=None, results=0):
        import os
        from pytom.agnostic.io import read, write
        from pytom.gpu.initialize import xp, device
        import numpy as np
        import time
        from pytom.agnostic.reconstruction_functions import backProjectGPU as backProject

        xp.cuda.Device(gpuID).use()

        ts = time.time()

        try:
            folder = os.path.dirname(particles[0].getFilename())

            if not results:
                results = []

                for index in range(4):
                    results.append(read('{}/.temp_{}.mrc'.format(folder, index)))

            [projections, vol_phi, vol_the, vol_offsetProjections] = [xp.asarray(x) for x in results]
            num_projections = projections.shape[2]

            vol_bp = xp.zeros((cubeSize, cubeSize, cubeSize), dtype=xp.float32)
            reconstructionPosition = xp.zeros((num_projections,3), dtype=xp.float32)

            interpolation = 'filt_bspline'

            proj_angles = vol_the.squeeze()

            for p in particles:
                reconstructionPosition *= 0
                vol_bp *=0

                # adjust coordinates of subvolumes to binned reconstruction

                for i in range(num_projections):
                    reconstructionPosition[i,0] = float(p.getPickPosition().getX() / binning)
                    reconstructionPosition[i,1] = float(p.getPickPosition().getY() / binning)
                    reconstructionPosition[i,2] = float(p.getPickPosition().getZ() / binning)

                vol_bp = backProject(projections, vol_bp, vol_phi.squeeze(), proj_angles, reconstructionPosition,
                                     vol_offsetProjections, interpolation).get()
                write(p.getFilename(), vol_bp)


            #for a in [vol_img, vol_phi, vol_the, vol_offsetProjections]: del a
            #del results


        except Exception as e:
            print('Caught exception in worker thread (x = %d):' % procID)
            print()
            raise e

        print(f'recon time: {time.time()-ts:.3f} sec')

    def extract_single_particle(self, particles, pid, verbose, binning, postScale, cubeSize, filename_ppr='', gpuID=None, rr=0):
        import os
        from pytom.basic.files import read
        from pytom_volume import vol, backProject, rescaleSpline
        import numpy as np
        import time

        ts = time.time()

        if 1:
            folder = os.path.dirname(particles[0].getFilename())

            results = []
            for index in range(4):
                results.append(read('{}/.temp_{}.mrc'.format(folder, index)))

            [vol_img, vol_phi, vol_the, vol_offsetProjections] = results
            num_projections = vol_img.sizeZ()

            vol_bp = vol(cubeSize, cubeSize, cubeSize)
            vol_bp.setAll(0.0)

            reconstructionPosition = vol(3, num_projections, 1)
            reconstructionPosition.setAll(0.0)
            start_index = pid * len(self)

            for p in particles:
                reconstructionPosition.setAll(0.0)
                vol_bp.setAll(0.0)
                # adjust coordinates of subvolumes to binned reconstruction
                if not filename_ppr:
                    for i in range(num_projections):
                        reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)
                else:

                    x, y, z = float(p.getPickPosition().getX() / binning), float(
                        p.getPickPosition().getY() / binning), float(p.getPickPosition().getZ() / binning)

                    for i in range(len(self)):
                        from pytom.gui.guiFunctions import LOCAL_ALIGNMENT_RESULTS, loadstar
                        import pytom.basic.combine_transformations as ct

                        particle_polish_file = loadstar(filename_ppr, dtype=LOCAL_ALIGNMENT_RESULTS)

                        offsetX = float(particle_polish_file['AlignmentTransX'][start_index + i])
                        offsetY = float(particle_polish_file['AlignmentTransY'][start_index + i])

                        pick_position = np.matrix([x, y, z]).T

                        rotation = ct.matrix_rotate_3d_y(self._list[i].getTiltAngle())

                        shift = np.matrix([offsetX, offsetY, 0]).T

                        new_position = rotation * pick_position

                        new_position_shifted = new_position + shift
                        reposition = np.linalg.inv(rotation) * new_position_shifted
                        reposition = np.array(reposition.T)[0]

                        reconstructionPosition(float(reposition[0] / binning), 0, i, 0)
                        reconstructionPosition(float(reposition[1] / binning), 1, i, 0)
                        reconstructionPosition(float(reposition[2] / binning), 2, i, 0)


                backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

                if postScale > 1:
                    volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                    rescaleSpline(vol_bp, volumeRescaled)
                    volumeRescaled.write(p.getFilename())
                else:
                    vol_bp.write(p.getFilename())

            for a in [vol_img, vol_phi, vol_the, vol_offsetProjections]: del a
            del results


        # except Exception as e:
        #     print('Caught exception in worker thread (x = %d):' % pid)
        #     print()
        #     raise e

        print(f'recon time: {time.time()-ts:.3f} sec')

    def saveParticleProjections(self, particles, projectionSize, binning=1,
                                applyWeighting=False, showProgressBar=False, verbose=False, outputScale=1):
        """
        saveParticleProjections: Saves all small projections of particle.

        @param particles:
        @type particles: L{pytom.basic.structures.ParticleList}
        @param projectionSize: x,y size of the resulting small projections
        @param binning: Necessary for positions in projections. Were L{pytom.basic.structures.PickLocation} determined in a binned volume and do you want to reconstruct a unbinned subtomogram? Use libtomc binning notation!
        @param applyWeighting: Will apply weighting for weighted backprojection to projection. Default is off!
        @param showProgressBar:
        @param verbose:
        @param outputScale: Scale the output by a factor. 1 will leave as it is, 2 will scale by 2, 3 by 3, 4 by 4...
        """
        from pytom_volume import vol, complexRealMult, rescaleSpline
        from pytom.basic.files import readProxy as read
        from pytom.reconstruction.reconstructionFunctions import positionInProjection
        from pytom.tools.files import writeSpider, checkDirExists
        from pytom.tools.ProgressBar import FixedProgBar
        from os import mkdir
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift
        from pytom.tools.maths import ZRotationMatrix
        from pytom.basic.normalise import mean0std1

        if showProgressBar:
            progressBar = FixedProgBar(0, len(particles) * len(self._list), 'Particle projections generated ')
            progressBar.update(0)
            numberParticleProjections = 0

        tmpImage = read(self._list[0].getFilename())
        imgSizeX = tmpImage.sizeX()
        imgSizeY = tmpImage.sizeY()

        if applyWeighting:
            weightSlice = fourierFilterShift(rampFilter(projectionSize, projectionSize))
            circleFilterRadius = tmpImage.sizeX() / 2
            circleSlice = fourierFilterShift(circleFilter(projectionSize, projectionSize, circleFilterRadius))

        for p in particles:

            particleName = p.getFilename()
            directoryName = particleName[0:len(particleName) - 3] + '/'

            if not checkDirExists(directoryName):
                mkdir(directoryName)

            particleProjectionList = ProjectionList([])
            particleProjectionListSpider = ProjectionList([])

            if verbose:
                print(p)

            vector = [p.getPickPosition().getX() * binning, p.getPickPosition().getY() * binning,
                      p.getPickPosition().getZ() * binning]

            if verbose:
                print('Original position: ', vector)

            for i in range(len(self._list)):

                if verbose:
                    print(self._list[i])

                # position in aligned projection
                # the angle MUST be negative because we want the inverse rotation! The matrix within the function is the forward rotation!
                position2D = positionInProjection(vector, -self._list[i].getTiltAngle(), 'Y')
                if verbose:
                    print('Position in projection: ', position2D)

                # position in unaligned projection

                vector2 = [position2D[0] * 1 / self._list[i].getAlignmentMagnification(),
                           position2D[1] * 1 / self._list[i].getAlignmentMagnification(), 0.0]
                if verbose:
                    print('Position after magnification: ', vector2)

                # position2D was extended to [x,y,0] -> 3D, can now be multiplied to a 3D rotation matrix. only the 2D part is used
                rotationMatrix = ZRotationMatrix(-self._list[i].getAlignmentRotation())

                newPosition = rotationMatrix * vector2

                if verbose:
                    print('Position after inplane rotation: ', newPosition)

                position2D[0] = newPosition[0] - self._list[i].getAlignmentTransX()
                position2D[1] = newPosition[1] - self._list[i].getAlignmentTransY()

                if verbose:
                    print(position2D)

                # catch io error (RuntimeError) exception when cut out piece is outside of big file.
                # in case of error, omit this projection, print warning
                try:
                    # changing origin, after that in upper left corner
                    # x and y are integers!
                    x = (int(position2D[0]) - (projectionSize / 2) + imgSizeX / 2)
                    y = (int(position2D[1]) - (projectionSize / 2) + imgSizeY / 2)

                    if verbose:
                        print('Final position in projection: ', x, y)

                    # determining error from pick center -> offset
                    offsetX = position2D[0] - float(int(position2D[0]))
                    offsetY = position2D[1] - float(int(position2D[1]))

                    if verbose:
                        print('Position offset: ', offsetX, offsetY)

                    projection = read(self._list[i].getFilename(), x, y, 0, projectionSize, projectionSize, 1, 0, 0, 0,
                                      0, 0, 0)

                    if applyWeighting:
                        projection = ifft(complexRealMult(complexRealMult(fft(projection), weightSlice), circleSlice))

                    try:
                        mean0std1(projection)
                    except ValueError as e:
                        print(str(e.message))
                        continue

                    if outputScale > 1:
                        projectionRescaled = vol(projectionSize / outputScale, projectionSize / outputScale, 1)
                        rescaleSpline(projection, projectionRescaled)
                        projection = projectionRescaled

                    projection.write(directoryName + 'projection_' + str(i) + '.em')
                    writeSpider(projection, directoryName + 'projection_' + str(i) + '.spi')

                    particleProjection = Projection(directoryName + 'projection_' + str(i) + '.em',
                                                    self._list[i].getTiltAngle(), offsetX, offsetY,
                                                    self._list[i].getAlignmentTransX(),
                                                    self._list[i].getAlignmentTransY(),
                                                    self._list[i].getAlignmentRotation(),
                                                    self._list[i].getAlignmentMagnification())
                    particleProjectionList.append(particleProjection)

                    particleProjectionSpider = Projection(directoryName + 'projection_' + str(i) + '.spi',
                                                          self._list[i].getTiltAngle(), offsetX, offsetY,
                                                          self._list[i].getAlignmentTransX(),
                                                          self._list[i].getAlignmentTransY(),
                                                          self._list[i].getAlignmentRotation(),
                                                          self._list[i].getAlignmentMagnification())
                    particleProjectionListSpider.append(particleProjectionSpider)

                except RuntimeError:
                    print('Warning : Particle out of bounds (' + str(x) + ',' + str(y) + ') for projection!')
                    print(p)
                    print(self._list[i])
                # assert False

                if showProgressBar:
                    numberParticleProjections = numberParticleProjections + 1
                    progressBar.update(numberParticleProjections)

            particleProjectionList.toXMLFile(directoryName + 'projectionList.xml')
            particleProjectionListSpider.toXMLFile(directoryName + 'projectionListSpider.xml')

    # TODO new general to_projection_stack functions

    def to_projection_stack(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                            scale_factor_particle=1, angle_specimen=0, show_progress_bar=False, verbose=False):
        """
        Create a projection stack from current Projections in the ProjectionList by weighting the projections.
        Alignment is applied if aligned parameters have been initialized on the Projections.

        @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
        @type weighting: L{int}
        @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
        pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
        @type binning: L{int}
        @param low_pass_freq: lowpass filter (in Nyquist)
        @type low_pass_freq: L{float}
        @param apply_circle_filter: optional circular filter in fourier space
        @type circleFilter: L{bool}
        @param scale_factor_particle: scaling factor for particles
        @type scale_factor_particle: L{float}
        @param angle_specimen: additional rotation of the specimen
        @type angle_specimen: L{float}
        @param show_progress_bar: write progress to terminal
        @type show_progress_bar: L{bool}
        @param verbose: print output default=False
        @type verbose: L{bool}

        @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
        and an offsetStack (user specified x, y offset)
        @rtype: L{list}

        @author: Marten Chaillet
        """
        from pytom_volume import complexRealMult, vol, paste, pasteCenter, mean
        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d, resize
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.basic.files import read
        from pytom.tools.ProgressBar import FixedProgBar
        import pytom_freqweight

        if low_pass_freq > 1.:
            low_pass_freq = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        imdimX = self[0].getXSize()
        imdimY = self[0].getYSize()
        imdim = max(imdimX, imdimY)

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)
            imdim = int(float(imdim) / float(binning) + .5)

        # prepare stacks
        stack = vol(imdim, imdim, len(self))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(self))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self))
        offsetStack.setAll(0.0)

        # design filters
        sliceWidth = imdim  # this should be relative to the particle size as that determines the crowther freq

        # pre-determine analytical weighting function and lowpass for speedup
        if weighting == -1:
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if apply_circle_filter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if low_pass_freq > 0:
            lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          low_pass_freq / 5. * imdim)

        if show_progress_bar:
            progressBar = FixedProgBar(0, len(self), 'Creating aligned and weighted projections')
            progressBar.update(0)

        for i, (tilt_angle, projection) in enumerate(zip(self._tilt_angles, self._list)):

            if verbose:
                print(projection)

            if projection._filename.split('.')[-1] == 'st':
                image = read(file=projection.getFilename(),
                             subregion=[0, 0, i - 1, imdim, imdim, 1],
                             sampling=[0, 0, 0], binning=[0, 0, 0])
                if not (binning == 1):
                    image = resize(volume=image, factor=1 / float(binning))[0]
            else:
                image = read(projection.getFilename())
                # image = rotate(image,180.,0.,0.)
                image = resize(volume=image, factor=1 / float(binning))[0]

            # normalize to contrast - subtract mean and norm to mean
            immean = mean(image)
            image = (image - immean) / immean

            # smoothen borders to prevent high contrast oscillations
            image = taper_edges(image, imdim // 30)[0]

            # transform projection according to tilt alignment
            transX = projection.getAlignmentTransX() / binning
            transY = projection.getAlignmentTransY() / binning
            rot = projection.getAlignmentRotation()
            mag = projection.getAlignmentMagnification() * scale_factor_particle
            order = projection.getOperationOrder()

            # 3 -- square if needed
            if imdimY != imdimX:
                newImage = vol(imdim, imdim, 1)
                newImage.setAll(0)
                pasteCenter(image, newImage)
                image = newImage

            # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
            # so IMOD default order is scaling > rotation > translation, that seems a weird order...
            center = (image.sizeX() / 2 - 0.5, image.sizeY() / 2 - 0.5, 0) if order == [1, 2, 0] else None

            if not (rot == 0 and transX == 0 and transY == 0 and mag == 0):
                image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order,
                                            crop=True, center=center)

            if low_pass_freq > 0:
                filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            # smoothen once more to avoid edges
            image = taper_edges(image, imdim // 30)[0]

            # weighting
            if weighting == -1:  # ramp filter
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)
            elif weighting == 1:  # exact filter
                weightSlice = fourierFilterShift(exactFilter(self._tilt_angles, tilt_angle, imdim, imdim, sliceWidth))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            phiStack(projection.getAlignmentAxisAngle(), 0, 0, i)
            thetaStack(projection.getTiltAngle() - angle_specimen, 0, 0, i)  # TODO why round + int?
            offsetStack(int(round(projection.getOffsetX())), 0, 0, i)
            offsetStack(int(round(projection.getOffsetY())), 0, 1, i)
            paste(image, stack, 0, 0, i)

            if show_progress_bar:
                progressBar.update(i)

            if verbose:
                print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

        return [stack, phiStack, thetaStack, offsetStack]

    def _create_shared_stack(self, imdim):
        import numpy as np
        import ctypes
        import multiprocessing as mp
        # create shared arrays
        self._shared_dtype = np.float32
        self._shared_stack_shape = (imdim, imdim, len(self))
        self._shared_stack = mp.Array(ctypes.c_float, imdim * imdim * len(self))
        self._shared_stack_np = np.frombuffer(self._shared_stack.get_obj(),
                                                dtype=self._shared_dtype).reshape(self._shared_stack_shape)
        self._shared_stack_np[:] = np.zeros(self._shared_stack_shape, dtype=self._shared_dtype)

    def _put_in_shared_stack(self, image, index):
        import numpy as np
        from pytom_numpy import vol2npy
        # create a np view to global shared_stack
        shared_stack_np_ = np.asfortranarray(np.frombuffer(shared_stack.get_obj(),
                                                                dtype=self._shared_dtype).reshape(
            self._shared_stack_shape))
        shared_stack_np_[:, :, index] = vol2npy(image).copy(order='F')

    def _retrieve_shared_stack(self):
        from pytom_numpy import npy2vol
        copy = self._shared_stack_np.copy(order='F')
        del self._shared_stack_np, self._shared_stack
        return npy2vol(copy, len(self._shared_stack_shape))

    def to_projection_stack_parallel(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                                     scale_factor_particle=1, angle_specimen=0, num_procs=1,
                                     show_progress_bar=False, verbose=False):
        """
        Parallel version.
        Create a projection stack from current Projections in the ProjectionList by weighting the projections.
        Alignment is applied if aligned parameters have been initialized on the Projections.

        @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
        @type weighting: L{int}
        @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
        pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
        @type binning: L{int}
        @param low_pass_freq: lowpass filter (in Nyquist)
        @type low_pass_freq: L{float}
        @param apply_circle_filter: optional circular filter in fourier space
        @type circleFilter: L{bool}
        @param scale_factor_particle: scaling factor for particles
        @type scale_factor_particle: L{float}
        @param angle_specimen: additional rotation of the specimen
        @type angle_specimen: L{float}
        @param show_progress_bar: write progress to terminal
        @type show_progress_bar: L{bool}
        @param verbose: print output default=False
        @type verbose: L{bool}
        @param num_procs: number of available processes for parallelization
        @type num_procs: L{int}

        @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
        and an offsetStack (user specified x, y offset)
        @rtype: L{list}

        @author: Marten Chaillet
        """
        from pytom_volume import vol
        from functools import partial
        import multiprocessing as mp
        from contextlib import closing

        if low_pass_freq > 1.:
            low_pass_freq = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        imdimX = self[0].getXSize()
        imdimY = self[0].getYSize()
        imdim = max(imdimX, imdimY)

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)
            imdim = int(float(imdim) / float(binning) + .5)

        # design filters
        sliceWidth = imdim  # this should be relative to the particle size as that determines the crowther freq

        # TODO better use Pool and map the function on args?
        # TODO filters could also be shared arrrays
        # procs = []
        # manager = multiprocessing.Manager()
        # return_dict = manager.dict()
        # iter_args = []
        #
        # for i, _ in enumerate(zip(self._tilt_angles, self._list)):
        #
        #     # while len(procs) >= num_procs:
        #     #     time.sleep(1)
        #     #     procs = [proc for proc in procs if proc.is_alive()]
        #
        #     args = (i, weighting, binning, low_pass_freq, apply_circle_filter, scale_factor_particle,
        #             angle_specimen, verbose, imdim, imdimX, imdimY, sliceWidth)
        #
        #     # ssing.Process(target=self.align_single_projection, args=args)
        #     # procs.append(proc)
        #     # proc.start()
        #
        #     iter_args.append(args)
        # manager = multiprocessing.Manager()
        # results_dict = manager.dict()

        self._create_shared_stack(imdim)

        print([i for i, _ in enumerate(zip(self._tilt_angles, self._list))])

        with closing(mp.Pool(num_procs, initializer=init_pool_processes, initargs=(self._shared_stack, ))) as p:
            # pool = mp.Pool(processes=num_procs)
            results = p.map(partial(align_single_projection, projection_list=self._list,
                                    tilt_angles=self._tilt_angles, shared_dtype=self._shared_dtype,
                                    shared_stack_shape=self._shared_stack_shape, weighting=weighting, binning=binning,
                                    low_pass_freq=low_pass_freq, apply_circle_filter=apply_circle_filter,
                                    scale_factor_particle=scale_factor_particle, angle_specimen=angle_specimen,
                                    verbose=verbose, imdim=imdim, imdimX=imdimX, imdimY=imdimY,
                                    sliceWidth=sliceWidth),
                               [i for i, _ in enumerate(zip(self._tilt_angles, self._list))])
        p.join()
        # map could also iter over range(len(self._list)), however i explicitly zip _tilt_angles and the _list to
        # make sure that they are both present

        if show_progress_bar:
            print('========== All processes finished alignment ============')

        # while len(procs) > 0:
        #     time.sleep(1)
        #     procs = [proc for proc in procs if proc.is_alive()]

        # prepare stacks
        stack = self._retrieve_shared_stack()

        phiStack = vol(1, 1, len(self))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self))
        offsetStack.setAll(0.0)

        # fill the stacks
        for i, result in results:  # why round(int())?
            phiStack(result[0], 0, 0, i)
            thetaStack(result[1], 0, 0, i)
            offsetStack(int(round(result[2][0])), 0, 0, i)
            offsetStack(int(round(result[2][1])), 0, 1, i)

        return [stack, phiStack, thetaStack, offsetStack]

    def align_single_projection(self, index, weighting, binning, low_pass_freq, apply_circle_filter,
                                scale_factor_particle, angle_specimen, verbose, imdim, imdimX, imdimY, sliceWidth,
                                results):
        """
        Align Projection at index of the ProjectionList.
        Create a projection stack from current Projections in the ProjectionList by weighting the projections.
        Alignment is applied if aligned parameters have been initialized on the Projections.

        @param index: index of a projection in the list
        @type index: L{int}
        @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
        @type weighting: L{int}
        @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
        pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
        @type binning: L{int}
        @param low_pass_freq: lowpass filter (in Nyquist)
        @type low_pass_freq: L{float}
        @param apply_circle_filter: optional circular filter in fourier space
        @type circleFilter: L{bool}
        @param scale_factor_particle: scaling factor for particles
        @type scale_factor_particle: L{float}
        @param angle_specimen: additional rotation of the specimen
        @type angle_specimen: L{float}
        @param verbose: print output default=False
        @type verbose: L{bool}
        @param imdim: x and y dimension of final image
        @type imdim: L{int}
        @param imdimX: input projection x dimension size
        @type imdimX: L{int}
        @param imdimY: input projection y dimension size
        @type imdimY: L{int}

        @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
        and an offsetStack (user specified x, y offset)
        @rtype: L{list}

        @author: Marten Chaillet
        """
        from pytom_volume import complexRealMult, vol, pasteCenter, mean
        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d, resize
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.basic.files import read
        import pytom_freqweight

        projection = self._list[index]
        tilt_angle = self._tilt_angles[index]

        if verbose:
            print(projection)

        # pre-determine analytical weighting function and lowpass for speedup
        if weighting == -1:
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))
        elif weighting == 1:
            weightSlice = fourierFilterShift(exactFilter(self._tilt_angles, tilt_angle, imdim, imdim, sliceWidth))

        if apply_circle_filter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if low_pass_freq > 0:
            lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          low_pass_freq / 5. * imdim)

        if projection._filename.split('.')[-1] == 'st':
            image = read(file=projection.getFilename(),
                         subregion=[0, 0, index - 1, imdim, imdim, 1],
                         sampling=[0, 0, 0], binning=[0, 0, 0])
            if not (binning == 1):
                image = resize(volume=image, factor=1 / float(binning))[0]
        else:
            image = read(projection.getFilename())
            # image = rotate(image,180.,0.,0.)
            image = resize(volume=image, factor=1 / float(binning))[0]

        # normalize to contrast - subtract mean and norm to mean
        immean = mean(image)
        image = (image - immean) / immean

        # smoothen borders to prevent high contrast oscillations
        image = taper_edges(image, imdim // 30)[0]

        # transform projection according to tilt alignment
        transX = projection.getAlignmentTransX() / binning
        transY = projection.getAlignmentTransY() / binning
        rot = projection.getAlignmentRotation()
        mag = projection.getAlignmentMagnification() * scale_factor_particle
        order = projection.getOperationOrder()

        # 3 -- square if needed
        if imdimY != imdimX:
            newImage = vol(imdim, imdim, 1)
            newImage.setAll(0)
            pasteCenter(image, newImage)
            image = newImage

        # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
        # so IMOD default order is scaling > rotation > translation, that seems a weird order...
        center = (image.sizeX() / 2 - 0.5, image.sizeY() / 2 - 0.5, 0) if order == [1, 2, 0] else None

        if not (rot == 0 and transX == 0 and transY == 0 and mag == 0):
            image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order,
                                        crop=True, center=center)

        if low_pass_freq > 0:
            filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
            image = filtered[0]

        # smoothen once more to avoid edges
        image = taper_edges(image, imdim // 30)[0]

        # apply weighting
        if weighting in [-1, 1]:
            image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

        if verbose:
            print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

        # assign in self._shared_stack_np
        self._put_in_shared_stack(image, index)

        return index, (projection.getAlignmentAxisAngle(), projection.getTiltAngle() - angle_specimen,
                       (projection.getOffsetX(), projection.getOffsetY()))

    def to_projection_stack_gpu(self):
        return None

    # TODO is the scaleFactorParticle really needed here? has not fuctionality in the GUI yet
    def toProjectionStack(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                          scale_factor_particle=1, angle_specimen=0, showProgressBar=False, verbose=False):
        """
        toProjectionStack:

        @param binning: binning factor
        @type binning: int
        @param applyWeighting: applyWeighting
        @type applyWeighting: bool
        @param tiltAngle: optional tiltAngle of image
        @type tiltAngle: float
        @param showProgressBar: show pretty bar
        @type showProgressBar: bool
        @param verbose: talkative
        @type verbose: bool
        @param num_procs: number of parallel processes (for the future)
        @type num_procs: int
        @param scaleFactorParticle: scaling factor for the reconstruction default=1 (not used)
        @type scaleFactorParticle: float

        @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]
        """
        from pytom_volume import vol, paste, complexRealMult, pasteCenter, mean
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.agnostic.io import read_size
        import pytom_freqweight

        # determine image dimensions according to first image in projection list
        imdimX = read_size(self._list[0].getFilename(), 'x')
        imdimY = read_size(self._list[0].getFilename(), 'y')
        imdim = max(imdimX, imdimY)

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)
            imdim = int(float(imdim) / float(binning) + .5)

        # prepare stacks
        stack = vol(imdim, imdim, len(self._list))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(self._list))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self._list))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self._list))
        offsetStack.setAll(0.0)

        # design filters
        sliceWidth = imdim  # this should be relative to the particle size as that determines the crowther freq

        if weighting == -1:
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if apply_circle_filter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if low_pass_freq > 0:
            print(f'Creating low-pass filter for projections. Fraction of Nyquist {lowpassFilter}')
            if low_pass_freq > 1.:
                low_pass_freq = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
            lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          low_pass_freq / 5. * imdim)

        if showProgressBar:
            progressBar = FixedProgBar(0, len(self._particleList), 'Stack creation')
            progressBar.update(0)

        self.tilt_angles = []  # TODO is this also not already initiated at the class init?

        for ii, projection in enumerate(self._list):
            self.tilt_angles.append(projection._tiltAngle)  # TODO are these correctly assigned????

        self.tilt_angles = sorted(self.tilt_angles)  # they need to be sorted, and the specific angle is got later

        for ii, projection in enumerate(self._list):

            if verbose:
                print(projection)

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

            # normalize to contrast - subtract mean and norm to mean
            immean = mean(image)
            image = (image - immean) / immean

            if imdimY != imdimX:
                newImage = vol(imdim, imdim, 1)
                newImage.setAll(0)
                pasteCenter(image, newImage)
                image = newImage

            if low_pass_freq > 0:
                filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            # weighting
            if weighting == -1:  # ramp filter
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)
            elif weighting == 1:  # exact filter
                print('weighting == 1, exact filter, needs tilt angle infomration, currently not provided')
                weightSlice = fourierFilterShift(exactFilter(self.tilt_angles, projection._tiltAngle,
                                                             imdim, imdim, sliceWidth))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            thetaStack(int(round(projection.getTiltAngle() - angle_specimen)), 0, 0, ii)
            offsetStack(projection.getOffsetX(), 0, 0, ii)
            offsetStack(projection.getOffsetY(), 0, 1, ii)
            paste(image, stack, 0, 0, ii)

            if showProgressBar:
                progressBar.update(ii)

            if verbose:
                print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

        return [stack, phiStack, thetaStack, offsetStack]

    def toProjectionStackFromAlignmentResultsFile(self, alignmentResultsFile, weighting=0,
                                                  binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                                                  scale_factor_particle=1, order=[2,1,0],
                                                  angle_specimen=0, showProgressBar=False, verbose=False):
        # TODO angle_specimen is not used
        """read image and create aligned projection stack, based on the results described in the alignmentResultFile.

           @param alignmentResultsFile: filename of alignment result file (see pytom.basic.datatypes)
           @type alignmentResultsFile: txt
           @param weighting: weighting (<0: analytical weighting, >1: exact weighting, 0/None: no weighting )
           @type weighting: float
           @param lowpassFilter: lowpass filter (in Nyquist)
           @type lowpassFilter: float
           @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.
           @type binnning: int
           @param circleFilter: optional circular filter in fourier space
           @type circleFilter: bool
           @param num_procs: number of parallel processes default=1
           @type num_procs: int
           @param verbose: print output default=False
           @type verbose: bool
           @param scaleFactorParticle: scaling factor for the reconstruction default=1
           @type scaleFactorParticle: float

           @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]

           @author: GS
        """
        from pytom_volume import complexRealMult, vol, paste, pasteCenter, mean
        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d, resize
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO  # TODO why two types of alignment result file
        from pytom.gui.guiFunctions import datatypeAR, loadstar
        from pytom.agnostic.io import read_size  # TODO is this needed?
        from pytom.tools.ProgressBar import FixedProgBar
        import pytom_freqweight

        print(f"Create aligned images from {alignmentResultsFile}")

        try:
            alignmentResults = loadstar(alignmentResultsFile, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        except:
            alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)

        imageList = alignmentResults['FileName']  # TODO do not read these, use self._list instead
        tilt_angles = alignmentResults['TiltAngle']

        imdimX = read_size(imageList[0], 'x')
        imdimY = read_size(imageList[0], 'y')
        imdim = max(imdimX, imdimY)

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)
            imdim = int(float(imdim) / float(binning) + .5)

        # prepare stacks
        stack = vol(imdim, imdim, len(imageList))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(imageList))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(imageList))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(imageList))
        offsetStack.setAll(0.0)

        # design filters
        sliceWidth = imdim  # this should be relative to the particle size as that determines the crowther freq

        # pre-determine analytical weighting function and lowpass for speedup
        if weighting == -1:
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if apply_circle_filter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if low_pass_freq > 0:
            print(f'Creating low-pass filter for projections. Fraction of Nyquist {lowpassFilter}')
            if low_pass_freq > 1.:
                low_pass_freq = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
            lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          low_pass_freq / 5. * imdim)

        projectionList = ProjectionList()
        for n, image in enumerate(imageList):
            atx = alignmentResults['AlignmentTransX'][n]
            aty = alignmentResults['AlignmentTransY'][n]
            rot = alignmentResults['InPlaneRotation'][n]
            mag = alignmentResults['Magnification'][n]
            projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                    alignmentRotation=rot, alignmentMagnification=mag)
            projectionList.append(projection)

        if showProgressBar:
            progressBar = FixedProgBar(0, len(self._particleList), 'Stack creation')
            progressBar.update(0)

        for ii, projection in enumerate(projectionList):

            if verbose:
                print(projection)

            if projection._filename.split('.')[-1] == 'st':
                from pytom.basic.files import EMHeader, read
                idx = projection._index
                image = read(file=projection._filename,
                             subregion=[0, 0, idx - 1, imdim, imdim, 1],
                             sampling=[0, 0, 0], binning=[0, 0, 0])
                if not (binning == 1) or (binning == None):
                    image = resize(volume=image, factor=1 / float(binning))[0]
            else:
                # read projection files
                from pytom.basic.files import EMHeader, read, read_em_header
                image = read(str(projection._filename))
                # image = rotate(image,180.,0.,0.)
                image = resize(volume=image, factor=1 / float(binning))[0]

            tiltAngle = projection._tiltAngle

            # normalize to contrast - subtract mean and norm to mean
            immean = mean(image)
            image = (image - immean) / immean

            # smoothen borders to prevent high contrast oscillations
            image = taper_edges(image, imdim // 30)[0]

            # transform projection according to tilt alignment
            transX = projection._alignmentTransX / binning
            transY = projection._alignmentTransY / binning
            rot = float(projection._alignmentRotation)
            mag = float(projection._alignmentMagnification) * scale_factor_particle

            # 3 -- square if needed
            if imdimY != imdimX:
                newImage = vol(imdim, imdim, 1)
                newImage.setAll(0)
                pasteCenter(image, newImage)
                image = newImage

            if 'OperationOrder' in alignmentResults.dtype.names:
                from pytom.agnostic.tools import convert_operation_order_str2list as str2list
                operation_string = alignmentResults['OperationOrder'][ii]
                order = str2list(operation_string)

            # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
            center = (image.sizeX()/2-0.5, image.sizeY()/2-0.5, 0) if order == [1,2,0] else None

            image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order, crop=True,
                                        center=center)

            if low_pass_freq > 0:
                filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            # smoothen once more to avoid edges
            image = taper_edges(image, imdim // 30)[0]

            # weighting
            if weighting == -1:  # ramp filter
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)
            elif weighting == 1:  # exact filter
                weightSlice = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            thetaStack(int(round(projection.getTiltAngle()))-self.angle_specimen, 0, 0, ii)
            offsetStack(int(round(projection.getOffsetX())), 0, 0, ii)
            offsetStack(int(round(projection.getOffsetY())), 0, 1, ii)
            paste(image, stack, 0, 0, ii)

            if showProgressBar:
                progressBar.update(ii)

            if verbose:
                print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

        return [stack, phiStack, thetaStack, offsetStack]

    def toProjectionStackParallel(self, alignmentResultsFile, weighting=None, lowpassFilter=0.9, binning=1,
                                  circleFilter=False, scaleFactorParticle=1, order=[2,1,0], angle_specimen=0,
                                  num_procs=1, verbose=False):
        """read image and create aligned projection stack in parallel, based on the results described in the alignmentResultFile.

        @param alignmentResultsFile: result file generate by the alignment script.
        @type alignmentResultsFile: str
        @param weighting: weighting (<0: analytical weighting, >1: exact weighting, 0/None: no weighting )
        @type weighting: float
        @param lowpassFilter: lowpass filter (in Nyquist)
        @type lowpassFilter: float
        @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.
        @param circleFilter: optional circular filter in fourier space
        @type circleFilter: bool
        @param num_procs: number of parallel processes default=1
        @type num_procs: int
        @param verbose: print output default=False
        @type verbose: bool
        @param scaleFactorParticle: scaling factor for the reconstruction default=1
        @type scaleFactorParticle: float

        @return Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]

        @author: GS
        """
        import numpy
        from pytom_numpy import vol2npy
        from pytom.basic.files import read
        from pytom.basic.filter import circleFilter as cf, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom_volume import complexRealMult, vol, paste, pasteCenter
        import pytom_freqweight
        from pytom.basic.transformations import resize, rotate
        from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
        from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO
        from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
        from pytom_numpy import vol2npy
        import time, os
        from multiprocessing import Process
        from pytom.agnostic.io import read_size

        if 1: print(f"Create aligned images from {alignmentResultsFile}")
        try:
            alignmentResults = loadstar(alignmentResultsFile, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        except:
            alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
        imageList = alignmentResults['FileName']
        tilt_angles = alignmentResults['TiltAngle']

        imdimX = read_size(imageList[0], 'x')
        imdimY = read_size(imageList[0], 'y')
        imdim = max(imdimX, imdimY)

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)
            imdim = int(float(imdim) / float(binning) + .5)
        else:
            imdim = imdim

        sliceWidth = imdim

        projectionList = ProjectionList()
        for n, image in enumerate(imageList):
            atx = alignmentResults['AlignmentTransX'][n]
            aty = alignmentResults['AlignmentTransY'][n]
            rot = alignmentResults['InPlaneRotation'][n]
            mag = alignmentResults['Magnification'][n]
            projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                    alignmentRotation=rot, alignmentMagnification=mag)
            projectionList.append(projection)

        stack = vol(imdim, imdim, len(imageList))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(imageList))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(imageList))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(imageList))
        offsetStack.setAll(0.0)

        procs = []

        rnum = numpy.random.randint(0, 1000000)

        for (ii, projection) in enumerate(projectionList):
            while len(procs) >= num_procs:
                time.sleep(1)
                procs = [proc for proc in procs if proc.is_alive()]

            if 'OperationOrder' in alignmentResults.dtype.names:
                from pytom.agnostic.tools import convert_operation_order_str2list as str2list
                operation_string = alignmentResults['OperationOrder'][ii]
                order = str2list(operation_string)
                if verbose: print(order)

            else:
                order =[2,1,0]

            args = (projection._filename, projection._index, projection._tiltAngle, tilt_angles, sliceWidth,
                    projection._alignmentTransX, projection._alignmentTransY, projection._alignmentRotation,
                    projection._alignmentMagnification, order,
                    scaleFactorParticle, weighting, imdim, imdimX, imdimY, binning, lowpassFilter, circleFilter,
                    f'TEMP{rnum}_{ii}.mrc')

            proc = Process(target=self.align_single_image, args=args)
            procs.append(proc)
            proc.start()

        while len(procs) > 0:
            time.sleep(1)
            procs = [proc for proc in procs if proc.is_alive()]

        for (ii, projection) in enumerate(projectionList):
            image = read(f'TEMP{rnum}_{ii}.mrc')
            thetaStack(int(round(projection.getTiltAngle())) - self.angle_specimen, 0, 0, ii)
            offsetStack(int(round(projection.getOffsetX())), 0, 0, ii)
            offsetStack(int(round(projection.getOffsetY())), 0, 1, ii)
            paste(image, stack, 0, 0, ii)
            os.system(f'rm TEMP{rnum}_{ii}.mrc')

        return [stack, phiStack, thetaStack, offsetStack]

    def align_single_image(self, filename, index, tiltAngle, tilt_angles, sliceWidth, alignmentTransX, alignmentTransY,
                           alignmentRotation, alignmentMagnification, order,
                           scaleFactorParticle, weighting, imdim, imdimX, imdimY, binning, lowpassFilter, circleFilter,
                           fname):
        """read image and create alignment images using the given params.

           @param fname: filename of tilt image
           @type fname: str
           @param index: index of image
           @type index: int
           @param tiltAngle: tilt angle of tilt image
           @type tiltAngle: float
           @param tilt_angles: all angles of tiltseries
           @type tilt_angles: numpy.ndarray or list
           @param sliceWidth: width of a slice in fourier space
           @type sliceWidth: int
           @param alignmentTransX: translation of tilt image in x-direction
           @type alignmentTransX: float
           @param alignmentTransY: translation of tilt image in y-direction
           @type alignmentTransY: float
           @param alignmentRotation: in-plane rotation of tilt image
           @type alignmentRotation: float
           @param alignmentMagnification: magnificaiton of tilt image
           @type alignmentMagnification: float
           @param scaleFactorParticle: scaling factor for the reconstruction default=1
           @type scaleFactorParticle: float
           @param weighting: weighting (<0: analytical weighting, >1: exact weighting, 0/None: no weighting )
           @type weighting: float
           @param imdim: max(imdimX, imdimY)
           @type imdim: int
           @param imdimX: num pixels in x-direction
           @type imdimX: int
           @param imdimY: num pix in y-dimension
           @type imdimY: int
           @param lowpassFilter: lowpass filter (in Nyquist)
           @type lowpassFilter: float
           @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.
           @param circleFilter: optional circular filter in fourier space
           @type circleFilter: bool
           @param fname: output filename, used for storage of alignment image. If np2vol would work, writing would not be needed.
           @type fname: str

           @author: GS
        """


        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter as filterFunction, bandpassFilter
        from pytom.basic.filter import circleFilter as cf, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom_volume import complexRealMult, vol, paste, pasteCenter
        import pytom_freqweight
        from pytom.basic.transformations import resize
        from pytom_numpy import vol2npy
        import time

        s = time.time()

        if (weighting != None) and (float(weighting) < -0.001):
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if circleFilter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(cf(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if lowpassFilter:
            print(f'Creating low-pass filter for projections. Fraction of Nyquist {lowpassFilter}')
            if lowpassFilter > 1.:
                lowpassFilter = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
            # weighting filter: arguments: (angle, cutoff radius, dimx, dimy,
            lpf = pytom_freqweight.weight(0.0, lowpassFilter * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          lowpassFilter / 5. * imdim)

        if filename.split('.')[-1] == 'st':
            from pytom.basic.files import EMHeader, read
            idx = index
            image = read(file=filename,
                         subregion=[0, 0, idx - 1, imdim, imdim, 1],
                         sampling=[0, 0, 0], binning=[0, 0, 0])
            if not (binning == 1) or (binning == None):
                image = resize(volume=image, factor=1 / float(binning))[0]
        else:
            # read projection files
            from pytom.basic.files import EMHeader, read, read_em_header
            image = read(str(filename))
            # image = rotate(image,180.,0.,0.)
            image = resize(volume=image, factor=1 / float(binning))[0]

        tiltAngle = tiltAngle

        # normalize to contrast - subtract mean and norm to mean
        immean = vol2npy(image).mean()
        image = (image - immean) / immean

        # smoothen borders to prevent high contrast oscillations
        image = taper_edges(image, imdim // 30)[0]

        # transform projection according to tilt alignment
        transX = alignmentTransX / binning
        transY = alignmentTransY / binning
        rot = float(alignmentRotation)
        mag = float(alignmentMagnification) * scaleFactorParticle

        # 3 -- square if needed
        if imdimY != imdimX:
            newImage = vol(imdim, imdim, 1)
            newImage.setAll(0)
            pasteCenter(image, newImage)
            image = newImage

        # IMOD Rotates around size/2 -0.5, wheras pytom rotates around size/2

        center = (image.sizeX()/2-1, image.sizeY()/2-1, 0) if order == [1,2,0] else None

        image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order, crop=True, center=center)

        if lowpassFilter:
            filtered = filterFunction(volume=image, filterObject=lpf, fourierOnly=False)
            image = filtered[0]

        # smoothen once more to avoid edges
        image = taper_edges(image, imdim // 30)[0]

        # analytical weighting
        if (weighting != None) and (weighting < 0):
            image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

        elif (weighting != None) and (weighting > 0):
            weightSlice = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
            image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

        image.write(fname)

        print(f'time to calculate {fname}: {time.time()-s:.2f}')

    # TODO add toProjectionStackGPU and toProjectionStackFr...FileGPU
    # TODO tell user that order is not used if not yet implemented
    def toProjectionStackGPU(self, binning=1, applyWeighting=False, scaleFactorParticle=1, showProgressBar=False,
                             verbose=False,):
        """
        toProjectionStack for the GPU

        @param binning: binning factor
        @type binning: int
        @param applyWeighting: applyWeighting
        @type applyWeighting: bool
        @param tiltAngle: optional tiltAngle of image
        @type tiltAngle: float
        @param showProgressBar: show pretty bar
        @type showProgressBar: bool
        @param verbose: talkative
        @type verbose: bool
        @param num_procs: number of parallel processes (for the future)
        @type num_procs: int
        @param scaleFactorParticle: scaling factor for the reconstruction default=1 (not used)
        @type scaleFactorParticle: float

        @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]
        """
        from pytom_numpy import vol2npy
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, exactFilter, rotateFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.agnostic.io import read_size

        # determine image dimensions according to first image in projection list
        imdimX = read_size(self._list[0].getFilename(), 'x')
        imdimY = read_size(self._list[0].getFilename(), 'y')
        imdim = max(imdimX, imdimY)

        stack = vol(imdim, imdim, len(self._list))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(self._list))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self._list))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self._list))
        offsetStack.setAll(0.0)

        if int(applyWeighting):
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))

        if showProgressBar:
            progressBar = FixedProgBar(0, len(self._particleList), 'Particle volumes generated ')
            progressBar.update(0)

        self.tilt_angles = []

        for (i, projection) in enumerate(self._list):
            self.tilt_angles.append(projection._tiltAngle)

        self.tilt_angles = sorted(self.tilt_angles)

        for (i, projection) in enumerate(self._list):

            if verbose:
                print(projection)

            if int((applyWeighting)) >= 1:

                if int(applyWeighting) != 2:
                    weightSlice = fourierFilterShift(exactFilter(self.tilt_angles, projection._tiltAngle,
                                                                 imdim, imdim, imdim))

                else:
                    weightSlice = fourierFilterShift(rotateFilter(self.tilt_angles, projection._tiltAngle,
                                                                  imdim, imdim, imdim))

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

            if int(applyWeighting):
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            if verbose:
                print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())

            thetaStack(int(round(projection.getTiltAngle() - self.angle_specimen)), 0, 0, i)
            offsetStack(projection.getOffsetX(), 0, 0, i)
            offsetStack(projection.getOffsetY(), 0, 1, i)
            paste(image, stack, 0, 0, i)

            if showProgressBar:
                progressBar.update(i)

        return [stack, phiStack, thetaStack, offsetStack]

    def toProjectionStackFromAlignmentResultsFileGPU(self, alignmentResultsFile, weighting=None, lowpassFilter=0.9,
                                                     binning=1, circleFilter=True, scaleFactorParicle=1, order=[2,1,0],
                                                     angle_specimen=0, verbose=False):
        from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
        from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
        from pytom.agnostic.io import read, write, read_size
        from pytom.agnostic.tools import taper_edges, create_circle, paste_in_center as pasteCenter
        from pytom.agnostic.filter import circle_filter, ramp_filter, exact_filter, ellipse_filter
        import pytom.voltools as vt
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.transform import resize

        print("Create aligned images from alignResults.txt")

        alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
        imageList = alignmentResults['FileName']
        tilt_angles = alignmentResults['TiltAngle']

        imdimX = read_size(imageList[0], 'x')
        imdimY = read_size(imageList[0], 'y')

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)

        imdim = max(imdimY, imdimX)
        sliceWidth = imdim

        if (weighting != None) and (float(weighting) < -0.001):
            cfreq = abs(tilt_angles[1:] - tilt_angles[:-1])
            cfreq = float(1 / xp.sin(cfreq.min() * xp.pi / 180)) // 1
            weightSlice = xp.fft.fftshift(ramp_filter(imdim, imdim, cfreq, len(tilt_angles)))

        if circleFilter:
            circleFilterRadiusX = imdim // 2
            circleFilterRadiusY = imdim // 2
            circleSlice = xp.fft.fftshift(ellipse_filter(imdim, imdim, circleFilterRadiusX, circleFilterRadiusY))
        else:
            circleSlice = xp.ones((imdim, imdim))

        # design lowpass filter
        if lowpassFilter:
            if lowpassFilter > 1.:
                lowpassFilter = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        projectionList = ProjectionList()
        for n, image in enumerate(imageList):
            atx = alignmentResults['AlignmentTransX'][n]
            aty = alignmentResults['AlignmentTransY'][n]
            rot = alignmentResults['InPlaneRotation'][n]
            mag = 1 / (alignmentResults['Magnification'][n])
            projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                    alignmentRotation=rot, alignmentMagnification=mag)
            projectionList.append(projection)

        stack = xp.zeros((imdim, imdim, len(imageList)), dtype=xp.float32)
        phiStack = xp.zeros((len(imageList)), dtype=xp.float32)
        thetaStack = xp.zeros((len(imageList)), dtype=xp.float32)
        offsetStack = xp.zeros((len(imageList), 2), dtype=xp.float32)

        for (ii, projection) in enumerate(projectionList):
            # print(f'Align {projection._filename}')
            image = read(str(projection._filename)).squeeze()

            if binning > 1:
                image = resize(image, 1 / binning)

            tiltAngle = projection._tiltAngle

            # 1 -- normalize to contrast - subtract mean and norm to mean
            immean = image.mean()
            image = (image - immean) / immean

            # 2 -- smoothen borders to prevent high contrast oscillations
            if ii == 0:
                image, taper_mask = taper_edges(image, imdim // 30)

            else:
                image *= taper_mask  # taper_edges(image,imdim//30,taper_mask)[0]

            # 3 -- square if needed
            if imdimY != imdimX:
                newImage = xp.zeros((imdim, imdim), dtype=xp.float32)
                image = pasteCenter(image, newImage)

            # 4 -- transform projection according to tilt alignment
            transX = projection._alignmentTransX / binning
            transY = projection._alignmentTransY / binning
            rot = float(projection._alignmentRotation)
            mag = float(projection._alignmentMagnification)

            inputImage = xp.expand_dims(image, 2)

            if ii == 0:
                outputImage = xp.zeros_like(inputImage, dtype=xp.float32)
            else:
                outputImage *= 0

            # TODO provide operation order parameter to voltools
            vt.transform(inputImage.astype(xp.float32), rotation=[0, 0, rot], rotation_order='rxyz', output=outputImage,
                         device=device, translation=[transX, transY, 0], scale=[mag, mag, 1],
                         interpolation='filt_bspline')
            del inputImage
            image = outputImage.squeeze()

            # 5 -- lowpass filter (optional)
            if lowpassFilter:
                from pytom.agnostic.filter import bandpass_circle
                image = bandpass_circle(image, high=lowpassFilter * imdim // 2, sigma=lowpassFilter / 5. * imdim)

            # 6 -- smoothen once more to avoid edges
            if imdimY == imdimX:
                image *= taper_mask
            else:
                image = taper_edges(image, imdim // 30)[0]

            # 7 -- weighting
            if (weighting != None) and (weighting < 0):
                # image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
                image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real
            elif (weighting != None) and (weighting > 0):
                weightSlice = xp.fft.fftshift(exact_filter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
                image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real

            thetaStack[ii] = float(round(projection.getTiltAngle() - angle_specimen))
            offsetStack[ii, :] = xp.array([int(round(projection.getOffsetX())), int(round(projection.getOffsetY()))])
            stack[:, :, ii] = image

        return [stack, phiStack, thetaStack, offsetStack]

    # TODO Bottom three functions are currently useless
    def read_projections(self, pid, num_procs, imgDim, applyWeighting, verbose, binning):
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex

        reduced_list = self._list[pid::num_procs]

        for (i, projection) in enumerate(reduced_list):
            print(i * num_procs + pid)

            if int(applyWeighting) >= 1:
                print(self.tilt_angles, projection._tiltAngle, imgDim)

                self.temp_weightSlice = fourierFilterShift(exactFilter(self.tilt_angles, projection._tiltAngle,
                                                                       imgDim, imgDim, imgDim))

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

            if int(applyWeighting):
                image = ifft(complexRealMult(complexRealMult(fft(image), self.temp_weightSlice), self.temp_circleSlice))

            self.temp_thetaStack(int(round(projection.getTiltAngle())), 0, 0, pid + num_procs * i)
            self.temp_offsetStack(projection.getOffsetX(), 0, 0, i * num_procs + pid)
            self.temp_offsetStack(projection.getOffsetY(), 0, 1, i * num_procs + pid)
            paste(image, self.temp_stack, 0, 0, i * num_procs + pid)

    def saveAsProjectionStack(self,filename,scale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        saveAsProjectionStack: Saves this projection list as a stack of projections
        """
        stack = self.toProjectionStack(scale,applyWeighting,showProgressBar,verbose)
        stack.write(filename)

    def setAlignmentParameters(self,alignmentXMLFile,verbose=False):
        """
        setAlignmentParameters:

        @param alignmentXMLFile:
        @param verbose: verbose level
	    @type verbose: bool
        """
        from lxml import etree
        from pytom.tools.files import readStringFile
        
        string = readStringFile(alignmentXMLFile)
        alignXML = etree.fromstring(string)
        
        if verbose:
            print(etree.tostring(alignXML,pretty_print = True))
        
        for projection in self:
            
            projectionName = projection.getFilename()
            
            names = alignXML.xpath('/root/projection/name') 
            
            for name in names:
                id = -1
                
                if str(name.text).find(projectionName) > -1 or projectionName.find(str(name.text)) > -1:
                    id = int(name.get('idx'))
                else:
                    continue
                
                if verbose:
                    print(etree.tostring(name,pretty_print = True))
                
                query = '/root/projection/alpha[@idx=' +str(id)+ ']'
                if verbose:
                    print(query)
                angle = alignXML.xpath(query)
                projection.setAlignmentRotation(float(angle[0].text))
                
                query = '/root/projection/tx[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                x = alignXML.xpath(query)
                projection.setAlignmentTransX(float(x[0].text))
                
                query = '/root/projection/ty[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                y = alignXML.xpath(query)
                projection.setAlignmentTransY(float(y[0].text))
                
                query = '/root/projection/s[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                magnification = alignXML.xpath(query)
                projection.setAlignmentMagnification(float(magnification[0].text))
