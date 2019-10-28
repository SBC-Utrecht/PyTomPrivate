'''
started on Mar 10, 2011

@author: luiskuhn, foerster
'''

from pytom.basic.structures import PyTomClass


class Projection(PyTomClass):
    """
    Projection: 
    """
    def __init__(self, filename=None, tiltAngle=None, offsetX=None,offsetY=None,
            alignmentTransX=None,alignmentTransY=None,
            alignmentRotation=None,alignmentMagnification=1.,
            alignedFilename=None, index=1):
        """
        Initialize projection class

        @param filename: The filename
        @type filename: string
        @param tiltAngle: Tilt angle of a projection
        @type tiltAngle: float
        @param offsetX: Offset X is used for subtomogram reconstruction from cut out projections. Takes care of offset due to positions determined in binned tomograms.
        @type offsetX: float
        @param offsetY:
        @type offsetY: float
        @param alignmentTransX: X Shift determined from alignment method
        @type alignmentTransX: float
        @param alignmentTransY: Y Shift determined from alignment method
        @type alignmentTransY: float
        @param alignmentRotation: Rotation determined from fine alignment method
        @type alignmentRotation: float
        @param alignmentMagnification: Magnification determined from fine alignment method      
        @type alignmentMagnification: float
        @param alignedFilename: filename of aligned projection
        @type alignedFilename: string
        @param index: Index of projection in filename (default: 1)
        @type index: int
    
        """
        self._filename = filename
        self._xSize = None
        self._ySize = None
        
        if filename and tiltAngle == None:
            self._loadTiltAngle()
#            print "TiltAngle = "+str(self._tiltAngle)+" deg"
        else:
            self._tiltAngle = tiltAngle

        self.setIndex( index)
        
        if offsetX == None:    
            self._offsetX = 0
        else:
            self._offsetX = offsetX    
        
        if offsetX == None:    
            self._offsetY = 0
        else:
            self._offsetY = offsetY
        
        if alignmentTransX:
            self._alignmentTransX = alignmentTransX 
        else:
            self._alignmentTransX = 0.
            
        if alignmentTransY:
            self._lignmentTransY = alignmentTransY
        else:
            self._alignmentTransY = 0.
        
        if alignmentRotation:
            self._alignmentRotation = alignmentRotation
        else:
            self._alignmentRotation = 0.
        
        self._alignmentMagnification = alignmentMagnification

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

        fh = open(self._filename,'rb')
        
        fh.seek(168)
        
        self._tiltAngle = float(struct.unpack('i',fh.read(4))[0])/1000
        
        fh.close()


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
    #use ParticleList as template class
    
    def __init__(self, list = None):
        self._list = list or []

    def info(self):
        tline = ""
        for projection in self._list:
            tline = tline + str(projection) + "\n"
        return tline

        
    def loadDirectory(self, directory):
        """
        loadDirectory: Will take all projections in a directory, determine their respective tilt angle from the header and populate this list object.
        @param directory: directory where files are located
        @type directory: L{str}
        """    
        from pytom.tools.files import checkDirExists
        
        if not checkDirExists(directory):
            raise RuntimeError('Directory ' + directory +' does not exist!')
        
        import os
        
        if directory[-1] != os.sep:
            #append / to dir
            directory += os.sep
        
        
        files = os.listdir(directory)
        
        for file in files:
            if file[len(file)-3:len(file)] == '.em':
                projection = Projection(directory + file)
            
                self.append(projection)
    
    def sort(self):
        self._list.sort(key=lambda p: p.getTiltAngle())
        
    def append(self, projection):
        """
        append projection to projection list

        @param projection: a projection
        @type projection: L{pytom.reconstructionStructures.Projection}
        """
        self._list.append(projection)
        
    def __getitem__(self,key):
        """
        __getitem__: Retrieve projection at position definded by key
        @param key:
        @type key: int
        @rtype: L{pytom.reconstructionStructures.Projection}. 
        """
        if not isinstance(key, (int, long)):
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
        print 'Generate Volumes is deprecated, use reconstructVolumes instead!'
        self.reconstructVolumes(particles, cubeSize, binning, applyWeighting,
	        showProgressBar ,verbose,preScale,postScale)

    def reconstructVolume(self, dims=[512,512,128], reconstructionPosition=[0,0,0], 
            binning=1, applyWeighting=False):
        """
        reconstruct a single 3D volume from weighted and aligned projections

        @param dims: 3D Size of reconstructed volume
        @type dims: 3d list
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: 3-dim list
        @param applyWeighting: apply weighting
        @type applyWeighting: bool

        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom_volume import vol,backProject

        vol_bp = vol(dims[0], dims[1], dims[2])
        vol_bp.setAll(0.0)

        # stacks for images, projections angles etc.
        [vol_img, vol_phi, vol_the, vol_offsetProjections] =  self.toProjectionStack(
                binning=binning, applyWeighting=applyWeighting, showProgressBar=False,
                verbose=False)

        # volume storing reconstruction offset from center (x,y,z)   
        recPosVol = vol(3,vol_img.sizeZ(),1)
        recPosVol.setAll(0.0)
        for iproj in range(0,vol_img.sizeZ()):
            for ii in range(0,3):
                recPosVol.setV(float(reconstructionPosition[ii]/binning),ii,iproj,0)
        
        # finally backproject into volume
        backProject(vol_img, vol_bp, vol_phi, vol_the, recPosVol, vol_offsetProjections)

        return vol_bp


    def reconstructVolumes(self, particles, cubeSize, binning=1, applyWeighting=False,
                           showProgressBar=False, verbose=False, preScale=1, postScale=1,
                           alignResultFile='', num_procs=10, num_procs_read=10, particle_polish_file=None, dimz=None, notpolished=False):
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

        @changes: Made it use multiple threads [Douwe, 2019]
        """
        from pytom_volume import vol, paste, backProject, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift

        if len(self) == 0:
            print RuntimeWarning('This ProjectionList contains no projections!Abort!')
            return

        if len(particles) == 0:
            raise RuntimeError('ParticleList is empty!')

        if cubeSize.__class__ != int:
            raise TypeError('reconstructVolumes: Parameter cubeSize must be of type int!')

        if showProgressBar:
            progressBar = FixedProgBar(0, len(particles), 'Particle volumes generated ')
            progressBar.update(0)

        if not alignResultFile and particle_polish_file is None:
            self.projectionstacks = self.toProjectionStack(
                binning=binning, applyWeighting=applyWeighting, showProgressBar=False,
                verbose=False)
        elif particle_polish_file is None:
            from pytom.reconstruction.generateAlignedTiltImagesInMemory import \
                toProjectionStackFromAlignmentResultsFile
            self.projectionstacks = toProjectionStackFromAlignmentResultsFile(
                alignResultFile, binning=binning, weighting=applyWeighting, verbose=False)

        self.reconstruct_volume_particles = particles
        self.particle_polish_file = particle_polish_file

        print('start subtomogram reconstructions!' )
        from pytom_numpy import vol2npy
        #print(vol2npy(self.projectionstacks[0]))
        from multiprocessing import Process
        import time
        procs = []

        for pid in range(num_procs):
            args = (pid, num_procs, verbose, binning, postScale, cubeSize, dimz, notpolished)
            if particle_polish_file is None:
                proc = Process(target=self._worker_reconstruct_volumes, args=args)
            else:
                proc = Process(target=self._worker_reconstruct_volumes_polished, args=args)
            procs.append(proc)
            proc.start()

        while len(procs):
            time.sleep(0.4)
            procs = [proc for proc in procs if proc.is_alive()]

    def _worker_reconstruct_volumes(self, pid, num_procs, verbose, binning, postScale, cubeSize, dimz, notpolished):
        from pytom_volume import vol, backProject, rescaleSpline
        from pytom.basic.files import read

        [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.projectionstacks
        particles = self.reconstruct_volume_particles

        index = 0
        assert num_procs >= 1
        while True:
            pos = index*num_procs+pid
            if pos >= len(particles): break
            index += 1
            p = particles[pos]
            x = p.getPickPosition().getX()
            y = p.getPickPosition().getY()
            z = p.getPickPosition().getZ()
            filename = p.getFilename()

            vol_bp = vol(cubeSize, cubeSize, cubeSize)
            vol_bp.setAll(0.0)

            reconstructionPosition = vol(3, len(self), 1)
            reconstructionPosition.setAll(0.0)

            # adjust coordinates of subvolumes to binned reconstruction
            for i in range(len(self)):
                reconstructionPosition(float(x / binning), 0, i, 0)
                reconstructionPosition(float(y / binning), 1, i, 0)
                reconstructionPosition(float(z / binning), 2, i, 0)

            if verbose:
                print(x / binning, y / binning, z / binning)

            backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

            if postScale > 1:
                volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                rescaleSpline(vol_bp, volumeRescaled)
                volumeRescaled.write(filename)
            else:
                vol_bp.write(filename)

            del vol_bp

    def _worker_reconstruct_volumes_polished(self, pid, num_procs, verbose, binning, postScale, cubeSize, dimz, notpolished):
        from pytom_volume import vol, backProject, rescaleSpline, pasteCenter, paste, complexRealMult
        from pytom.basic.files import read, read_size
        from pytom.basic.transformations import general_transform2d
        from math import cos, sin, pi, floor
        import os
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift_ReducedComplex

        # For proper results from WBP to avoid interpolation effects at the edges, the final result will be the original cubeSize
        cubeSize *= 2

        particles = self.reconstruct_volume_particles
        print("polishfile:", self.particle_polish_file)

        dimx, dimy, _ = read_size(self._list[0].getFilename())
        dimz = dimx if dimz is None else dimz
        print(dimx, dimy)
        print(self._list)
        print(len(self._list))

        weightSlice = fourierFilterShift_ReducedComplex(rampFilter(cubeSize, cubeSize))
        circleFilterRadius = cubeSize / 2
        circleSlice = fourierFilterShift_ReducedComplex(circleFilter(cubeSize, cubeSize, circleFilterRadius))

        index = 0
        assert num_procs >= 1
        while True:
            pos = index * num_procs + pid
            if pos >= len(particles): break
            index += 1
            p = particles[pos]
            x = p.getPickPosition().getX()
            y = p.getPickPosition().getY()
            z = p.getPickPosition().getZ()

            stack = vol(cubeSize, cubeSize, len(self._list))
            stack.setAll(0.0)

            phiStack = vol(1, 1, len(self._list))
            phiStack.setAll(0.0)

            thetaStack = vol(1, 1, len(self._list))
            thetaStack.setAll(0.0)

            offsetStack = vol(1, 2, len(self._list))
            offsetStack.setAll(0.0)

            reconstructionPosition = vol(3, len(self), 1)
            reconstructionPosition.setAll(0.0)

            start_index = pos * len(self._list)
            raw_size = int(cubeSize*1.5)
            print(raw_size)

            for n, proj in enumerate(self._list):
                tiltangle = int(self.particle_polish_file['TiltAngle'][start_index + n]) # Int because that seems to give better results, it is still unclear why (Douwe and Gijs)

                # print(tiltangle, start_index, n)

                # Retrieve the right file based on the tiltangle, assumes the same ordering in the particlepolishfile
                # and numbering of the projections, and assumes a naming scheme for the projections
                folder = os.path.dirname(proj.getFilename())
                filename = '{}/sorted_aligned_{}.em'.format(folder, n + 1)

                offsetX = 0 if notpolished else (self.particle_polish_file['AlignmentTransX'][start_index + n])
                offsetY = 0 if notpolished else float(self.particle_polish_file['AlignmentTransY'][start_index + n])

                # Set the coordinates to absolute (top left) and to 2D on this projection
                yy = y #+ dimy / 2
                xx = cos(tiltangle * pi / 180) * (x - dimx / 2) - sin(tiltangle * pi / 180) * (z - dimz / 2) + dimx / 2

                x_top = int(round(xx) - raw_size / 2 + round(offsetX))
                y_top = int(round(yy) - raw_size / 2 + round(offsetY))
                x_bot = dimx - x_top - raw_size
                y_bot = dimy - y_top - raw_size

                # The difference between raw_size and the resulting size of the image because of edges
                x_diff = min(0, x_top) + min(0, x_bot)
                y_diff = min(0, y_top) + min(0, y_bot)

                try:
                    img = read(filename, [max(0, x_top), max(0, y_top), 0, max(0, raw_size + x_diff),
                                         max(0, raw_size + y_diff), 1], [0, 0, 0], [binning, binning, 1])
                except Exception as e:
                    raise Exception("Exception in read {:s}\ntop {:d},{:d} bot {:d},{:d} size {:d},{:d} 3d {:f},{:f},{:f} 2d {:f},{:f} offset {:f},{:f}".format(e.message, x_top, y_top, x_bot, y_bot, max(0, raw_size + x_diff), max(0, raw_size + y_diff), x, y, z, xx, yy, offsetX, offsetY))

                x_shift = xx - round(xx) + offsetX - round(offsetX) + (-min(0, x_top) + min(0, x_bot)) / 2.
                y_shift = yy - round(yy) + offsetY - round(offsetY) + (-min(0, y_top) + min(0, y_bot)) / 2.

                print(x_shift, y_shift, offsetX - round(offsetX), offsetY - round(offsetY))

                interpolated_temp = general_transform2d(v=img, shift=[x_shift, y_shift], rot=0, scale=1, order=[2,1,0], crop=True)

                # interpolated_temp = general_transform(img, shift=[
                #     xx - floor(xx) - offsetX + floor(offsetX) + (-min(0, x_top) + min(0, x_bot)) / 2.,
                #     yy - floor(yy) - offsetY + floor(offsetY) + (-min(0, y_top) + min(0, y_bot)) / 2., 0])

                interpolated = vol(cubeSize, cubeSize, 1)
                interpolated.setAll(0.0)
                pasteCenter(interpolated_temp, interpolated)

                # Apply weighting
                interpolated = ifft(complexRealMult(complexRealMult(fft(interpolated), weightSlice), circleSlice)) / interpolated.numelem()

                thetaStack(tiltangle, 0, 0, n)
                paste(interpolated, stack, 0, 0, n)
                #reconstructionPosition(cubeSize / 2, 0, n, 0)
                #reconstructionPosition(cubeSize / 2, 1, n, 0)
                #reconstructionPosition(cubeSize / 2, 2, n, 0)

                from pytom_numpy import vol2npy
                print(vol2npy(img).mean(), vol2npy(interpolated).mean(), tiltangle)

                # import pylab as pl
                # from pytom_numpy import vol2npy
                # f, ax = pl.subplots(1, 2)
                # im0 = ax[0].imshow(vol2npy(img))
                # f.colorbar(im0)
                # im1 = ax[1].imshow(vol2npy(interpolated))
                # f.colorbar(im1)
                # pl.show()

            vol_bp = vol(cubeSize/2, cubeSize/2, cubeSize/2)
            vol_bp.setAll(0.0)
            filename = str(self.particle_polish_file['FileName'][start_index])

            # print(vol2npy(thetaStack).mean(), vol2npy(stack).mean(), 0xff)
            # import pylab as pl

            # f, ax = pl.subplots(5, 5)

            # for i in range(len(self._list)):
            #     ax[i//5][i%5].imshow(vol2npy(stack)[:, :, i])
            # pl.show()

            if verbose:
                print(x / binning, y / binning, z / binning)

            # Create custom projectionstack per particle to have all the aligned images as given by the particle polish
            # file, and the position is the center and the offsets are 0
            # use general_transform (basic.transformations) to interpolate the remainder of the shift after getting it
            # to integer position with the help of read (basic.files)

            from pytom.basic.files import write_em
            write_em("stack.em", stack)

            print(os.getcwd())

            backProject(stack, vol_bp, phiStack, thetaStack, reconstructionPosition, offsetStack)

            if postScale > 1:
                volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                rescaleSpline(vol_bp, volumeRescaled)
                volumeRescaled.write(filename)
            else:
                vol_bp.write(filename)

            del vol_bp

    def saveParticleProjections(self, particles, projectionSize,binning=1,
            applyWeighting = False,showProgressBar = False,verbose=False,outputScale=1):
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
        from pytom_volume import vol,complexRealMult,rescaleSpline
        from pytom.basic.files import readProxy as read
        from pytom.reconstruction.reconstructionFunctions import positionInProjection
        from pytom.tools.files import writeSpider,checkDirExists
        from pytom.tools.ProgressBar import FixedProgBar
        from os import mkdir
        from pytom.basic.fourier import fft,ifft
        from pytom.basic.filter import circleFilter,rampFilter,fourierFilterShift
        from pytom.tools.maths import ZRotationMatrix
        from pytom.basic.normalise import mean0std1
        
        if showProgressBar:
            progressBar = FixedProgBar(0,len(particles)*len(self._list),'Particle projections generated ')
            progressBar.update(0)
            numberParticleProjections = 0
            
        tmpImage = read(self._list[0].getFilename())
        imgSizeX = tmpImage.sizeX()
        imgSizeY = tmpImage.sizeY()
        
        if applyWeighting:
            weightSlice = fourierFilterShift(rampFilter(projectionSize,projectionSize))
            circleFilterRadius = tmpImage.sizeX()/2
            circleSlice = fourierFilterShift(circleFilter(projectionSize,projectionSize,circleFilterRadius))
         
        for p in particles:
            
            particleName = p.getFilename()
            directoryName = particleName[0:len(particleName)-3]+ '/'
            
            if not checkDirExists(directoryName):
                mkdir(directoryName)
            
            particleProjectionList = ProjectionList([])
            particleProjectionListSpider = ProjectionList([])
            
            if verbose:
                print p    
                        
            vector = [p.getPickPosition().getX()*binning, p.getPickPosition().getY()*binning, p.getPickPosition().getZ()*binning]
            
            if verbose:
                print 'Original position: ' , vector
            
            for i in range(len(self._list)):
                               
                if verbose:
                    print self._list[i]
                
                #position in aligned projection
                #the angle MUST be negative because we want the inverse rotation! The matrix within the function is the forward rotation!
                position2D = positionInProjection(vector,-self._list[i].getTiltAngle(),'Y')
                if verbose:
                    print 'Position in projection: ' ,position2D
                
                #position in unaligned projection
                
                vector2 = [position2D[0] * 1/self._list[i].getAlignmentMagnification() ,position2D[1] * 1/self._list[i].getAlignmentMagnification(),0.0]
                if verbose:
                    print 'Position after magnification: ', vector2
                
                #position2D was extended to [x,y,0] -> 3D, can now be multiplied to a 3D rotation matrix. only the 2D part is used     
                rotationMatrix = ZRotationMatrix(-self._list[i].getAlignmentRotation())
                
                newPosition = rotationMatrix * vector2 
                
                if verbose:
                    print 'Position after inplane rotation: ', newPosition
                    
                position2D[0] = newPosition[0] - self._list[i].getAlignmentTransX()
                position2D[1] = newPosition[1] - self._list[i].getAlignmentTransY()
                
                if verbose:
                    print position2D
                    
                #catch io error (RuntimeError) exception when cut out piece is outside of big file. 
                #in case of error, omit this projection, print warning
                try:                    
                    #changing origin, after that in upper left corner
                    #x and y are integers!
                    x = (int(position2D[0]) -(projectionSize/2) + imgSizeX/2)
                    y = (int(position2D[1]) -(projectionSize/2) + imgSizeY/2)                    
                    
                    if verbose:
                        print 'Final position in projection: ',x,y
                
                    #determining error from pick center -> offset
                    offsetX = position2D[0] - float(int(position2D[0]))
                    offsetY = position2D[1] - float(int(position2D[1])) 
                    
                    if verbose:
                        print 'Position offset: ',offsetX,offsetY
                        
                    projection = read(self._list[i].getFilename(), x, y, 0, projectionSize, projectionSize, 1, 0, 0, 0, 0, 0, 0)
                    
                    
                    if applyWeighting:
                        projection = ifft( complexRealMult( complexRealMult( fft(projection), weightSlice), circleSlice) )
                    
                    try:
                        mean0std1(projection)
                    except ValueError,e:
                        print str(e.message)
                        continue
                     
                    if outputScale > 1:
                        projectionRescaled = vol(projectionSize/outputScale,projectionSize/outputScale,1)
                        rescaleSpline(projection,projectionRescaled)
                        projection = projectionRescaled
                       
                    projection.write(directoryName + 'projection_' + str(i) + '.em')
                    writeSpider(projection,directoryName + 'projection_' + str(i) + '.spi')
                    
                    particleProjection = Projection(directoryName + 'projection_' + str(i) + '.em', self._list[i].getTiltAngle(),offsetX,offsetY,self._list[i].getAlignmentTransX(),self._list[i].getAlignmentTransY(),self._list[i].getAlignmentRotation(),self._list[i].getAlignmentMagnification())
                    particleProjectionList.append(particleProjection)                    
                    
                    particleProjectionSpider = Projection(directoryName + 'projection_' + str(i) + '.spi', self._list[i].getTiltAngle(),offsetX,offsetY,self._list[i].getAlignmentTransX(),self._list[i].getAlignmentTransY(),self._list[i].getAlignmentRotation(),self._list[i].getAlignmentMagnification())
                    particleProjectionListSpider.append(particleProjectionSpider)
                    
                except RuntimeError:
                    print 'Warning : Particle out of bounds ('+ str(x) + ',' + str(y) + ') for projection!'
                    print p
                    print self._list[i]
                #assert False
                
                if showProgressBar:
                    numberParticleProjections = numberParticleProjections + 1
                    progressBar.update(numberParticleProjections)
                           
            particleProjectionList.toXMLFile(directoryName + 'projectionList.xml')
            particleProjectionListSpider.toXMLFile(directoryName + 'projectionListSpider.xml')
    
    def toXML(self):
        """
        """
        from lxml import etree
                
        projectionList_element = etree.Element('ProjectionList')
        
        ##self._list = sorted(self._list, key=lambda Projection: Projection._tiltAngle)
        
        for projection in self._list:
            projectionList_element.append(projection.toXML())
        
        return projectionList_element   
        
    def fromXML(self,xmlObj):
        """
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
      
        if xmlObj.tag == 'ProjectionList':
            projectionList_element = xmlObj
        else:
            Exception('Is not a ProjectionList! You must provide a valid ProjectionList object.')
            
        for projection in projectionList_element:
            
            newProjection = Projection()
            newProjection.fromXML(projection)
            self.append(newProjection)
        
        self._list = sorted(self._list, key=lambda Projection: Projection._tiltAngle)

    def toProjectionStack(self, binning=1, applyWeighting=False, showProgressBar=False, verbose=False):
        """
        toProjectionStack:
        @param binning: binning factor
        @type binning: int
        @param applyWeighting: applyWeighting
        @type applyWeighting: bool
        @param showProgressBar: show pretty bar
        @type showProgressBar: bool
        @param verbose: talkative
        @type verbose: bool
        @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]
        """
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift_ReducedComplex

        # determine image dimensions according to first image in projection list
        imgDim = read(self._list[0].getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1).sizeX()

        stack = vol(imgDim, imgDim, len(self._list))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(self._list))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self._list))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self._list))
        offsetStack.setAll(0.0)

        if applyWeighting:
            weightSlice = fourierFilterShift_ReducedComplex(rampFilter(imgDim, imgDim))
            circleFilterRadius = imgDim / 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imgDim, imgDim, circleFilterRadius))

        if showProgressBar:
            progressBar = FixedProgBar(0, len(self._particleList), 'Particle volumes generated ')
            progressBar.update(0)

        for (i, projection) in enumerate(self._list):
            if verbose:
                print projection

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)
            # Shift the projection based on the particle polish file

            if applyWeighting:
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice)) / image.numelem()

            thetaStack(projection.getTiltAngle(), 0, 0, i)
            offsetStack(projection.getOffsetX(), 0, 0, i)
            offsetStack(projection.getOffsetY(), 0, 1, i)
            paste(image, stack, 0, 0, i)

            if showProgressBar:
                progressBar.update(i)

        return [stack, phiStack, thetaStack, offsetStack]
               
    def toProjectionStackMT(self,binning=1,applyWeighting=False,showProgressBar=False,verbose=False, num_procs=1):
        """
        toProjectionStack:

        @param binning: binning factor
        @type binning: int
        @param applyWeighting: applyWeighting
        @type applyWeighting: bool
        @param showProgressBar: show pretty bar
        @type showProgressBar: bool
        @param verbose: talkative
        @type verbose: bool

        @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]
        """
        from pytom_volume import vol,paste,complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft,ifft
        from pytom.basic.filter import circleFilter,rampFilter,fourierFilterShift
        
        # determine image dimensions according to first image in projection list
        imgDim = read(self._list[0].getFilename(),0,0,0,0,0,0,0,0,0,binning,binning,1).sizeX()

        self.projectionstack_stack = vol(imgDim, imgDim, len(self._list))
        self.projectionstack_stack.setAll(0.0)
        
        phiStack = vol(1, 1, len(self._list))
        phiStack.setAll(0.0)
        if num_procs > 1: self.projectionstack_phiStack = phiStack
        
        thetaStack = vol(1, 1, len(self._list))
        thetaStack.setAll(0.0)
        if num_procs > 1: self.projectionstack_thetaStack = thetaStack
        
        offsetStack = vol(1, 2, len(self._list))
        offsetStack.setAll(0.0)
        if num_procs > 1: self.projectionstack_offsetStack = offsetStack
        
        if applyWeighting:
            self.projectionstack_weightSlice = fourierFilterShift(rampFilter(imgDim,imgDim))
            circleFilterRadius = imgDim/2
            self.projectionstack_circleSlice = fourierFilterShift(circleFilter(imgDim,imgDim, circleFilterRadius))
        
        if showProgressBar:
            progressBar = FixedProgBar(0,len(self._particleList),'Particle volumes generated ')
            progressBar.update(0)

        if num_procs > 1:
            from multiprocessing import Process
            import time
            procs = []

            for pid in range(num_procs):
                args = (pid, num_procs, binning, applyWeighting, verbose)
                proc = Process(target=self._worker_to_projection_stack, args=args)
                procs.append(proc)
                proc.start()

            while len(procs):
                time.sleep(0.4)
                procs = [proc for proc in procs if proc.is_alive()]

        else:
            for (i, projection) in enumerate(self._list):
                if verbose:
                    print projection

                image = read(projection.getFilename(),0,0,0,0,0,0,0,0,0,binning,binning,1)

                if applyWeighting:
                    image = ifft( complexRealMult( complexRealMult( fft(image), self.projectionstack_weightSlice), self.projectionstack_circleSlice) )

                thetaStack(projection.getTiltAngle(), 0, 0, i)
                offsetStack(projection.getOffsetX(),0,0,i)
                offsetStack(projection.getOffsetY(),0,1,i)
                paste(image, self.projectionstack_stack, 0, 0, i)

                if showProgressBar:
                    progressBar.update(i)

        print("Finished reconstruction")
        return [self.projectionstack_stack,self.projectionstack_phiStack,self.projectionstack_thetaStack,self.projectionstack_offsetStack]

    def _worker_to_projection_stack(self, pid, num_procs, binning, applyWeighting, verbose):
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.basic.fourier import fft, ifft

        index = 0
        assert num_procs >= 1
        while True:
            i = index * num_procs + pid # Each worker has its own group of tasks
            if i >= len(self._list): break
            index += 1

            projection = self._list[i]

            if verbose:
                print projection

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

            if applyWeighting:
                image = ifft(complexRealMult(complexRealMult(fft(image), self.projectionstack_weightSlice),
                                             self.projectionstack_circleSlice))

            from pytom_numpy import vol2npy
            print("Image")
            print(vol2npy(image))

            self.projectionstack_thetaStack(projection.getTiltAngle(), 0, 0, i)
            self.projectionstack_offsetStack(projection.getOffsetX(), 0, 0, i)
            self.projectionstack_offsetStack(projection.getOffsetY(), 0, 1, i)
            paste(image, self.projectionstack_stack, 0, 0, i)

            print("Stack")
            print(vol2npy(self.projectionstack_stack))

            #del image

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
            print etree.tostring(alignXML,pretty_print = True)
        
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
                    print etree.tostring(name,pretty_print = True)
                
                query = '/root/projection/alpha[@idx=' +str(id)+ ']'
                if verbose:
                    print query
                angle = alignXML.xpath(query)
                projection.setAlignmentRotation(float(angle[0].text))
                
                query = '/root/projection/tx[@idx="'+str(id)+'"]'
                if verbose:
                    print query
                x = alignXML.xpath(query)
                projection.setAlignmentTransX(float(x[0].text))
                
                query = '/root/projection/ty[@idx="'+str(id)+'"]'
                if verbose:
                    print query
                y = alignXML.xpath(query)
                projection.setAlignmentTransY(float(y[0].text))
                
                query = '/root/projection/s[@idx="'+str(id)+'"]'
                if verbose:
                    print query
                magnification = alignXML.xpath(query)
                projection.setAlignmentMagnification(float(magnification[0].text))

            
class Reconstruction(PyTomClass):            

                
    def __init__(self,projectionList,particleList):
        """
        __init__
        """

        self._projectionList = None
        self._particleList = None

        if projectionList:
            self._projectionList = projectionList
                
        if particleList:
            self._particleList = particleList    
        
    
    def apply(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1):
        """
        apply: Apply reconstruction method
        @param cubeSize: Size of resulting cube
        @param preScale: Downscale projections BEFORE reconstruction. Default = 1. Can be 1x, 2x
        @param postScale:   Downscale projections AFTER reconstruction. Default = 1. Can be 1x, 2x
        @param target: Target volume if no particle list is given.
        @param coordinatesScale: Do the determined coordinates stem from a "binned" volume?  
        """
        
        raise RuntimeError('You must run this function on a child of Reconstruction')


class WeightedBackprojection(PyTomClass):
    """
    WeightedBackprojection: A class that can do everything related to WB
    """          
                
    def apply(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        @param applyWeighting: Are the projections weighted. Apply weighting if not. Is False by default 
        """

        from pytom.tools.ProgressBar import FixedProgBar
        from pytom_volume import vol, backProject
        
        if not self._projectionList or len(self._projectionList) == 0:
            raise RuntimeError('This reconstruction object does not have a valid ProjectionList')
        
        if not self._particleList:
            from pytom.basic.structures import ParticleList,Particle
            
            p = Particle(target)
            pl = ParticleList('.')
            
            self._particleList = pl
            
        targetVolume = vol(cubeSize, cubeSize, cubeSize)
        targetVolume.setAll(0.0)
        
        [projectionStack,phiStack,thetaStack,projectionOffsetStack] = self._projectionStack.toProjectionStack() 

        if showProgressBar:
            progressBar = FixedProgBar(0,len(self._particleList),'Particle volumes generated ')
            progressBar.update(0)
            numberParticleVolumes = 0
        
        
        particleOffsetStack = vol(3, len(self) ,1)
        particleOffsetStack.setAll(0.0)
        
        
        for particleIndex in xrange(len(self._particleList)):    
            p = self._particleList[particleIndex]
            
            if verbose:
                print p
            
            for i in xrange(len(self)):
                particleOffsetStack.setV(p.getPickPosition().getX()*coordinatesScale, 0, i, 0) 
                particleOffsetStack.setV(p.getPickPosition().getY()*coordinatesScale, 1, i, 0)
                particleOffsetStack.setV(p.getPickPosition().getZ()*coordinatesScale, 2, i, 0)   
            
            backProject(projectionStack, targetVolume, phiStack, thetaStack, projectionOffsetStack,projectionOffsetStack)
            targetVolume.write(p.getFilename())

            if showProgressBar:
                numberParticleVolumes = numberParticleVolumes + 1
                progressBar.update(numberParticleVolumes)
            
        
    def reconstructTomogram(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        reconstructTomogram: Will reconstruct a WB tomogram
        """
        raise Exception("Not implemented: reconstructTomogram")
