def xcc(volume,template,mask=None, volumeIsNormalized=False):
    """
    xcc: Calculates the cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: only used for compatibility with nxcc - not used
    @type volumeIsNormalized: L{bool}
    @return: A unscaled value
    @raise exception: Raises a runtime error if volume and template have a different size.  
    @author: Thomas Hrabe 
    """
    from pytom_volume import sum
    from pytom.tools.macros import volumesSameSize
    
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')
   
    if mask: # mask is given
        result = mask * volume * template
    else:
        result = volume * template
    
    cc = sum(result)
    
    cc = cc / float(volume.numelem())
    
    return cc 
    
def nxcc(volume, template, mask=None, volumeIsNormalized=False):
    """
    nxcc: Calculates the normalized cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: speed up if volume is already normalized
    @type volumeIsNormalized: L{bool}
    @return: A value between -1 and 1
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe 
    @change: flag for pre-normalized volume, FF
    """

    from pytom_volume import vol,sum,limit
    from pytom.tools.macros import volumesSameSize
    
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')
    
    if not mask:
        from pytom.basic.normalise import mean0std1
        if not volumeIsNormalized:
           v = mean0std1(volume, True)
        t = mean0std1(template, True)
        p = volume.numelem()
        result = v*t
    else:
        from pytom_numpy import vol2npy
        from pytom.basic.normalise import normaliseUnderMask
        if not volumeIsNormalized:
            (v,p) = normaliseUnderMask(volume, mask)
            (t,p) = normaliseUnderMask(template, mask, p)
            t = t * mask # multiply with the mask
            result = v * t
        else:
            (t,p) = normaliseUnderMask(template,mask)
            t = t * mask # multiply with the mask
            result = volume * t
    
    ncc = sum(result)
    ncc = ncc / float(p)

    return ncc 

    
def xcf(volume, template, mask=None, stdV=None):
    """
    XCF: returns the non-normalised cross correlation function. The xcf 
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: changed: will be used if specified
    @type mask: L{pytom_volume.vol}
    @param stdV: Will be unused, only for compatibility reasons with FLCF 
    @return: XCF volume
    @rtype: L{pytom_volume.vol}
    @author: Thomas Hrabe
    """
    import pytom_volume
    from pytom.basic import fourier

    if mask != None:
        volume = volume * mask
        template = template * mask

    #determine fourier transforms of volumes
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
    
    if template.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        ftemplate = fft(template)
    else:
        ftemplate = template
    
    #perform element wise - conjugate multiplication    
    pytom_volume.conjugate(ftemplate)
    fresult = fvolume * ftemplate

    #transform back to real space
    result = fourier.ifft(fresult)
    
    fourier.iftshift(result)

    n = result.numelem()
    if mask:
        n1 = pytom_volume.sum(mask)
        # real -> FFT requires n1, FFT -> real n
        result.shiftscale(0,1/float(n1*n))
    else:
        result.shiftscale(0,1/float(n*n))
    
    return result

def nXcf(volume,template,mask=None, stdV=None):
    """
    nXCF: returns the normalised cross correlation function. Autocorrelation 
    of two equal objects would yield a max nxcf peak of 1.

    @param volume: The search volume
    @param template: The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the same size with the given template.
    @type mask: L{pytom_volume.vol}
    @param stdV: Will be unused, only for compatibility reasons with FLCF
    @return: the calculated nXcf volume
    @rtype: L{pytom_volume.vol}
    @author: Thomas Hrabe
    @change: masking of template implemented
    """
    from pytom.basic.normalise import mean0std1

    if mask == None:
        result = xcf(mean0std1(volume,True),mean0std1(template,True), mask=None, stdV=None)
    else:
        from pytom.basic.normalise import normaliseUnderMask
        result = xcf(normaliseUnderMask(volume=volume, mask=mask, p=None)[0],
                     normaliseUnderMask(volume=template, mask=mask, p=None)[0],
                     mask=mask, stdV=None)
    #n = result.numelem()
    #result.shiftscale(0,1/float(n*n))

    return result


def meanValueUnderMask(volume, mask=None, p=None):
    """
    meanValueUnderMask: Determines the mean value under a mask
    @param volume: The volume
    @type volume:  L{pytom_volume.vol}
    @param mask:  The mask
    @type mask:  L{pytom_volume.vol}
    @param p: precomputed number of voxels in mask
    @type p: float
    @return: A value (scalar)
    @rtype: single
    @change: support None as mask, FF 08.07.2014
    """
    from pytom_volume import vol, sum
  
    if not volume.__class__ == vol:
        raise TypeError("meanValueUnderMask: Volume parameter must be of type vol!")
    
    if mask:
        if not p:
            p = sum(mask)
        resV = volume*mask
        res = sum(resV)/p
    else:
        res = sum(volume)/(volume.sizeX()*volume.sizeY()*volume.sizeZ())
    
    return res

def stdValueUnderMask(volume, mask, meanValue, p=None):
    """
    stdValueUnderMask: Determines the std value under a mask

    @param volume: input volume
    @type volume:  L{pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{float} or L{int}
    @return: A value
    @rtype: L{float}
    @change: support None as mask, FF 08.07.2014
    """
    from pytom_volume import sum
    from pytom_volume import vol, power, variance
    
    assert volume.__class__ == vol
    if mask:
        if not p:
            p = sum(mask)
    else:
        p = volume.sizeX()*volume.sizeY()*volume.sizeZ()
    
    squareM = meanValue**2
    squareV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
    squareV.copyVolume(volume)
    power(squareV, 2)
    res = meanValueUnderMask(squareV, mask, p)
    res = res - squareM
    if res >= 0:
        res = res**0.5
    else:
        print("Res = %.6f < 0 => standard deviation determination fails :(")
        print("   something went terribly wrong and program has to stop")
        raise ValueError('Program stopped in stdValueUnderMask')

    return res

def meanUnderMask(volume, mask, p):
    """
    meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
    @param volume: input volume
    @type volume:  L{pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{int} or L{float}
    @return: the calculated mean volume under mask
    @rtype:  L{pytom_volume.vol}
    @author: Yuxiang Chen
    """
    size = volume.numelem()
    
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom_volume import conjugate
    # for some reason, this has to be conjugated. (Otherwise the asym mask won't work)
    fMask = fft(mask)
    conjugate(fMask)
    
    result = iftshift(ifft(fMask*fft(volume)))/(size*p)

    return result

def stdUnderMask(volume, mask, p, meanV):
    """
    stdUnderMask: calculate the std volume under the given mask
    @param volume: input volume
    @type volume:  L{pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{int}
    @param meanV: mean volume under mask, which should already been caculated
    @type meanV:  L{pytom_volume.vol}
    @return: the calculated std volume under mask
    @rtype:  L{pytom_volume.vol}
    @author: Yuxiang Chen
    """
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom_volume import vol, power, limit
    
    copyV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
    copyV.copyVolume(volume)
    power(copyV, 2) #calculate the square of the volume
    
    copyMean = vol(meanV.sizeX(), meanV.sizeY(), meanV.sizeZ())
    copyMean.copyVolume(meanV)
    power(copyMean, 2)

    result = meanUnderMask(copyV, mask, p) - copyMean

#    from pytom_volume import abs
#    abs(result)
    limit(result, 1e-09, 1, 0, 0, True, False) # this step is needed to set all those value (close to 0) to 1
    power(result, 0.5)

    return result


def FLCF(volume, template, mask=None, stdV=None, wedge=1):
    '''
    Created on Apr 13, 2010
    FLCF: compute the fast local correlation function
    This functions works only when the mask is in the middle.
    
    @param volume: target volume
    @type volume: L{pytom_volume.vol}
    @param template: template to be searched. It can have smaller size then target volume.
    @type template: L{pytom_volume.vol}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the same size with the given template.
    @type mask: L{pytom_volume.vol}
    @param stdV: standard deviation of the target volume under mask, which do not need to be calculated again when the mask is identical.
    @type stdV: L{pytom_volume.vol}
    @return: the local correlation function
    @rtype: L{pytom_volume.vol}
    
    @author: Yuxiang Chen
    '''
    from pytom_volume import vol, pasteCenter
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom_volume import conjugate
    from pytom.basic.structures import Mask
    from pytom_volume import sum
    from pytom.basic.files import write_em

    if volume.__class__ != vol and template.__class__ != vol:
        raise RuntimeError('Wrong input type!')
    
    if volume.sizeX()<template.sizeX() or volume.sizeY()<template.sizeY() or volume.sizeZ()<template.sizeZ():
        raise RuntimeError('Template size is bigger than the target volume')

    # generate the mask 
    if mask.__class__ != vol:
        from pytom_volume import initSphere
        mask = vol(template.sizeX(), template.sizeY(), template.sizeZ())
        mask.setAll(0)
        initSphere(mask, template.sizeX()/2,0,0,template.sizeX()/2, template.sizeY()/2, template.sizeZ()/2)
    else:
        if template.sizeX()!=mask.sizeX() and template.sizeY()!=mask.sizeY() and template.sizeZ()!=mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')
    
    # calculate the non-zeros
    p = sum(mask)

    # normalize the template under mask
    meanT = meanValueUnderMask(template, mask, p)
    stdT = stdValueUnderMask(template, mask, meanT, p)

    temp = (template - meanT)/stdT
    temp = temp * mask

    # construct both the template and the mask which has the same size as target volume
    tempV = temp
    if volume.sizeX() != temp.sizeX() or volume.sizeY() != temp.sizeY() or volume.sizeZ() != temp.sizeZ():
        tempV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        tempV.setAll(0)
        pasteCenter(temp, tempV)
    
    maskV = mask
    if volume.sizeX() != mask.sizeX() or volume.sizeY() != mask.sizeY() or volume.sizeZ() != mask.sizeZ():
        maskV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        maskV.setAll(0)
        pasteCenter(mask, maskV)
    
    # calculate the mean and std of volume
    if stdV.__class__ != vol:
        meanV = meanUnderMask(volume, maskV, p)
        stdV = stdUnderMask(volume, maskV, p, meanV)

    size = volume.numelem()
    fT = fft(tempV)
    conjugate(fT)
    result = iftshift(ifft(fT*fft(volume)))/stdV

    result.shiftscale(0,1/(size*p))
    
    return result


def bandCC(volume,reference,band,verbose = False):
    """
    bandCC: Determines the normalised correlation coefficient within a band
    @param volume: The volume
    @type volume: L{pytom_volume.vol}
    @param reference: The reference
    @type reference: L{pytom_volume.vol}
    @param band: [a,b] - specify the lower and upper end of band.
    @return: First parameter - The correlation of the two volumes in the specified band. 
             Second parameter - The bandpass filter used.
    @rtype: List - [float,L{pytom_freqweight.weight}]
    @author: Thomas Hrabe    
    """
    import pytom_volume
    from pytom.basic.filter import bandpassFilter
    from pytom.basic.correlation import xcf
    from math import sqrt
    
    if verbose:
        print('lowest freq : ', band[0],' highest freq' , band[1])
        
    vf = bandpassFilter(volume,band[0],band[1],fourierOnly=True)
    rf = bandpassFilter(reference,band[0],band[1],vf[1],fourierOnly=True)
    
    ccVolume = pytom_volume.vol_comp(rf[0].sizeX(),rf[0].sizeY(),rf[0].sizeZ())
    ccVolume.copyVolume(rf[0])
    
    pytom_volume.conj_mult(ccVolume,vf[0])
    
    cc = pytom_volume.sum(ccVolume)

    cc = cc.real
    
    v = vf[0]
    r = rf[0]
    
    absV = pytom_volume.abs(v)
    absR = pytom_volume.abs(r)
    
    pytom_volume.power(absV,2)
    pytom_volume.power(absR,2)
    
    sumV = pytom_volume.sum(absV)
    sumR = pytom_volume.sum(absR)
    
    sumV = abs(sumV)
    sumR = abs(sumR)
    
    if sumV == 0:
        sumV =1
        
    if sumR == 0:
        sumR =1
        
    cc = cc / (sqrt(sumV*sumR)) 
    
    #numerical errors will be punished with nan
    if abs(cc) > 1.1 :
        cc = float('nan')
    
    return [cc,vf[1]];
    
def weightedXCC(volume,reference,numberOfBands,wedgeAngle=-1):
        """
        weightedXCC: Determines the band weighted correlation coefficient for a volume and reference. Notation according Steward/Grigorieff paper
        @param volume: A volume
        @type volume: L{pytom_volume.vol}
        @param reference: A reference of same size as volume
        @type reference: L{pytom_volume.vol}
        @param numberOfBands: Number of bands
        @param wedgeAngle: A optional wedge angle
        @return: The weighted correlation coefficient
        @rtype: float  
        @author: Thomas Hrabe   
        """    
        import pytom_volume
        from pytom.basic.fourier import fft
        from math import sqrt
        import pytom_freqweight
        result = 0
        numberVoxels = 0
        
        #volume.write('vol.em');
        #reference.write('ref.em');
        fvolume = fft(volume)
        freference = fft(reference)
        numelem = volume.numelem()
        
        fvolume.shiftscale(0,1/float(numelem))
        freference.shiftscale(0,1/float(numelem))
        from pytom.basic.structures import WedgeInfo
        wedge = WedgeInfo(wedgeAngle)
        wedgeVolume = wedge.returnWedgeVolume(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        
        increment = int(volume.sizeX()/2 * 1/numberOfBands)
        band = [0,100]
        for i in range(0,volume.sizeX()/2, increment):
        
            band[0] = i
            band[1] = i + increment
    
            r = bandCC(volume,reference,band)
            cc = r[0]
            
            #print cc;
            filter = r[1]
            #get bandVolume
            bandVolume = filter.getWeightVolume(True)
            
            filterVolumeReduced = bandVolume * wedgeVolume
            filterVolume = pytom_volume.reducedToFull(filterVolumeReduced)
            
            #determine number of voxels != 0    
            N = pytom_volume.numberSetVoxels(filterVolume)
            
            w = sqrt(1/float(N))
            
            #add to number of total voxels
            numberVoxels=numberVoxels + N
            #print 'w',w;
            #print 'cc',cc;
            #print 'N',N;
            
            cc2 = cc*cc
            #print 'cc2',cc2;
            if cc <= 0.0:
                cc = cc2
            else:
                cc = cc2/(cc+w)
            
            #print 'cc',cc;
            cc = cc *cc *cc; #no abs
            #print 'cc',cc;
            
            #add up result
            result = result + cc*N
        
        return result*(1/float(numberVoxels))
    
    
    
def FSCSum(volume,reference,numberOfBands,wedgeAngle=-1):
    """
    FSCSum: Determines the sum of the Fourier Shell Correlation coefficient for a volume and reference. 
    @param volume: A volume
    @type volume: L{pytom_volume.vol}
    @param reference: A reference of same size as volume
    @type reference: L{pytom_volume.vol}
    @param numberOfBands: Number of bands
    @param wedgeAngle: A optional wedge angle
    @return: The sum FSC coefficient
    @rtype: float  
    @author: Thomas Hrabe   
    """    
    
    from pytom.basic.correlation import bandCC
    import pytom_volume
    from pytom.basic.fourier import fft
    from math import sqrt
    import pytom_freqweight
    result = 0
    numberVoxels = 0
    
    #volume.write('vol.em');
    #reference.write('ref.em');
    fvolume = fft(volume)
    freference = fft(reference)
    numelem = volume.numelem()
    
    fvolume.shiftscale(0,1/float(numelem))
    freference.shiftscale(0,1/float(numelem))

    #print '-----'
    for i in range(numberOfBands):
        #process bandCorrelation
        band = []
        band[0] = i*volume.sizeX()/numberOfBands 
        band[1] = (i+1)*volume.sizeX()/numberOfBands
        
        r = bandCC(fvolume,freference,band)
        cc = r[0]
        #print cc
        result = result + cc
    #print '-----'
    
    return result*(1/float(numberOfBands))

def bandCF(volume,reference,band=[0,100]):
    """
    bandCF:
    @param volume: The volume
    @param reference: The reference
    @param band: [a,b] - specify the lower and upper end of band. [0,1] if not set.
    @return: First parameter - The correlation of the two volumes in the specified ring. 
             Second parameter - The bandpass filter used.
    @rtype: List - [L{pytom_volume.vol},L{pytom_freqweight.weight}]
    @author: Thomas Hrabe   
    @todo: does not work yet -> test is disabled
    """
    import pytom_volume
    from math import sqrt
    from pytom.basic import fourier
    from pytom.basic.filter import bandpassFilter
    from pytom.basic.correlation import nXcf
    
    vf = bandpassFilter(volume,band[0],band[1],fourierOnly=True)
    rf = bandpassFilter(reference,band[0],band[1],vf[1],fourierOnly=True)
    
    v = pytom_volume.reducedToFull(vf[0])
    r = pytom_volume.reducedToFull(rf[0])
    
    absV = pytom_volume.abs(v)
    absR = pytom_volume.abs(r)
    
    pytom_volume.power(absV,2)
    pytom_volume.power(absR,2)
    
    sumV = abs(pytom_volume.sum(absV))
    sumR = abs(pytom_volume.sum(absR))
    
    if sumV == 0:
        sumV = 1
        
    if sumR == 0:
        sumR = 1
    
    pytom_volume.conjugate(rf[0])
    
    fresult = vf[0] * rf[0]

    #transform back to real space
    result = fourier.ifft(fresult)
    
    fourier.iftshift(result)
    
    result.shiftscale(0,1/float(sqrt(sumV*sumR)))
    
    return [result,vf[1]]

def weightedXCF(volume,reference,numberOfBands,wedgeAngle=-1):
    """
    weightedXCF: Determines the weighted correlation function for volume and reference
    @param volume: A volume 
    @param reference: A reference 
    @param numberOfBands:Number of bands
    @param wedgeAngle: A optional wedge angle
    @return: The weighted correlation function 
    @rtype: L{pytom_volume.vol} 
    @author: Thomas Hrabe 
    @todo: does not work yet -> test is disabled
    """
    from pytom.basic.correlation import bandCF
    import pytom_volume
    from math import sqrt
    import pytom_freqweight
    
    result = pytom_volume.vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
    result.setAll(0)
    cc2 = pytom_volume.vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
    cc2.setAll(0)
    q = 0
    
    if wedgeAngle >=0:
        wedgeFilter = pytom_freqweight.weight(wedgeAngle,0,volume.sizeX(),volume.sizeY(),volume.sizeZ())
        wedgeVolume = wedgeFilter.getWeightVolume(True)
    else:
        wedgeVolume = pytom_volume.vol(volume.sizeX(), volume.sizeY(), int(volume.sizeZ()/2+1))
        wedgeVolume.setAll(1.0)
        
    w = sqrt(1/float(volume.sizeX()*volume.sizeY()*volume.sizeZ()))
    
    numberVoxels = 0
    
    for i in range(numberOfBands):
        """
        notation according Steward/Grigorieff paper
        """
        band = [0,0]
        band[0] = i*volume.sizeX()/numberOfBands 
        band[1] = (i+1)*volume.sizeX()/numberOfBands
        
        r = bandCF(volume,reference,band)
        
        cc = r[0]
                
        filter = r[1]
        #get bandVolume
        bandVolume = filter.getWeightVolume(True)
            
        filterVolumeReduced = bandVolume * wedgeVolume
        filterVolume = pytom_volume.reducedToFull(filterVolumeReduced)
        #determine number of voxels != 0    
        N = pytom_volume.numberSetVoxels(filterVolume)
            
        #add to number of total voxels
        numberVoxels = numberVoxels + N
                 
        cc2.copyVolume(r[0])
                
        pytom_volume.power(cc2,2)
        
        cc.shiftscale(w,1)
        ccdiv = cc2/(cc)
                
        pytom_volume.power(ccdiv,3)
        
        #abs(ccdiv); as suggested by grigorief
        ccdiv.shiftscale(0,N)
        
        result = result + ccdiv
    
    result.shiftscale(0,1/float(numberVoxels))
    
    return result

def FSC(volume1,volume2,numberBands,mask=None,verbose=False, filename=None):
    """
    FSC - Calculates the Fourier Shell Correlation for two volumes
    @param volume1: volume one
    @type volume1: L{pytom_volume.vol}
    @param volume2: volume two
    @type volume2: L{pytom_volume.vol}
    @param numberBands: number of shells for FSC
    @type numberBands: int
    @param filename: write FSC to ascii file if specified
    @type filename: string

    @return: Returns a list of cc values 
    @author: Thomas Hrabe  
    @rtype: list[floats]
    """
    
    from pytom.basic.correlation import bandCC
    from pytom.tools.macros import volumesSameSize
    from pytom.basic.structures import Mask
    import pytom_volume
    
    if not volumesSameSize(volume1, volume2):
        raise RuntimeError('Volumes must have the same size!')
    
    if mask:
        if mask.__class__ == pytom_volume.vol:
            volume1 = volume1 * mask
            volume2 = volume2 * mask
          
        elif mask.__class__ == Mask:
            mask = mask.getVolume()
            volume1 = volume1 * mask
            volume2 = volume2 * mask
        elif mask.__class__ == str:
            mask = pytom_volume.read(mask)
            volume1 = volume1 * mask
            volume2 = volume2 * mask 
        else:
            raise RuntimeError('FSC: Mask must be a volume OR a Mask object OR a string path to a mask!')  
        
    fscResult = []
    band = [-1,-1]
    if type(numberBands) == tuple:
        numberBands = tuple[0]
    
    increment = int(volume1.sizeX()/2 * 1/numberBands)
    
    for i in range(0,volume1.sizeX()//2, increment):
        
        band[0] = i
        band[1] = i + increment
        
        if verbose:
            print('Band : ' ,band)
            
        res = bandCC(volume1,volume2,band,verbose)
        
        if i == 0 and increment == 1:
            #force a 1 for correlation of the zero frequency 
            res[0] = 1
  
        if verbose:
            print('Correlation ' ,res[0])

        fscResult.append(res[0])

    if filename:
        f = open(filename,'w')
        for item in fscResult:
            f.write("%s\n" % item)
        f.close()

    return fscResult

def determineResolution(fsc,resolutionCriterion,verbose=False):
    """
    determineResolution: Determines frequency and band where correlation drops below the resolutionCriterion. Uses linear interpolation between two positions
    @param fsc: The fsc list determined by L{pytom.basic.correlation.FSC}
    @param resolutionCriterion: A value between 0 and 1
    @return: [resolution,interpolatedBand,numberBands] 
    @author: Thomas Hrabe 
    @todo: Add test! 
    """
    numberBands = len(fsc)
    
    band = numberBands

    for i in range(numberBands):
        if fsc[i] < resolutionCriterion:     
            band = i-1  #select the band that is still larger than criterion
            break
    
    if verbose:
        print('Band detected at ', band)
    
    if band == -1:
        raise RuntimeError("Please check your resolution criterion or you FSC!")
            
    elif band < numberBands:
        fsc1 = fsc[band]
        fsc2 = fsc[band+1]
        
        try:
            interpolatedBand = (resolutionCriterion-fsc1)/(fsc2-fsc1)+band            
        
        except ZeroDivisionError:
            interpolatedBand = band
        
    else:
        interpolatedBand = band
        
    if verbose:
        print('Band interpolated to ', interpolatedBand)
        
    resolution = (interpolatedBand+1) / float(numberBands)
    
    if resolution < 0 :
        resolution = 1
        interpolatedBand = numberBands
        print('Warning: PyTom determined a resolution < 0 for your data. Please check "mass" in data is positive or negative for all cubes.')
        print('Warning: Setting resolution to 1 and ',interpolatedBand)
        print('')
        
    return [resolution,interpolatedBand,numberBands]


def soc(volume,reference,mask=None, stdV=None):
    """
    soc : Second Order Correlation. Correlation of correlation peaks.
    @param volume: The volume
    @type volume:  L{pytom_volume.vol}
    @param reference: The reference / template
    @type reference:  L{pytom_volume.vol}
    @author: Thomas Hrabe   
    """
    referencePeak = FLCF(reference,reference,mask)
    peaks = FLCF(volume,reference,mask)
    
    return FLCF(peaks,referencePeak,mask)


def subPixelPeakParabolic(scoreVolume, coordinates, verbose=False):
    """
    quadratic interpolation of three adjacent samples
    @param scoreVolume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param verbose: be talkative
    @type verbose: bool
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy
    """
    from pytom_volume import vol
    assert type(scoreVolume) == vol, "scoreVolume must be vol!"
    if (coordinates[0] == 0 or coordinates[0] == scoreVolume.sizeX()-1 or
                coordinates[1] == 0 or coordinates[1] == scoreVolume.sizeY()-1 or
                coordinates[2] == 0 or coordinates[2] == scoreVolume.sizeZ()-1 ):
        if verbose:
            print("subPixelPeakParabolic: peak near borders - no interpolation done")
        return [scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]), coordinates]
    peakCoordinates = coordinates
    (x, p1, a1) = qint(ym1=scoreVolume.getV(coordinates[0]-1, coordinates[1], coordinates[2]), 
                     y0=scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                     yp1=scoreVolume.getV(coordinates[0]+1, coordinates[1], coordinates[2]))
    (y, p2, a2) = qint(ym1=scoreVolume.getV(coordinates[0], coordinates[1]-1, coordinates[2]), 
                     y0=scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                     yp1=scoreVolume.getV(coordinates[0], coordinates[1]+1, coordinates[2]))
    if scoreVolume.sizeZ() != 1:
        (z, p3, a3) = qint(ym1=scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]-1), 
                         y0=scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                         yp1=scoreVolume.getV(coordinates[0], coordinates[1], coordinates[2]+1))
        peakCoordinates[0] += x
        peakCoordinates[1] += y
        peakCoordinates[2] += z
        peakValue = p1 + a2*y**2 + a3*z**2
    else:
        peakCoordinates[0] += x
        peakCoordinates[1] += y
        peakValue = p1 + a2*y**2
    return [peakValue,peakCoordinates]


def qint(ym1, y0, yp1):
    """
    quadratic interpolation of three adjacent samples
     
    [p,y,a] = qint(ym1,y0,yp1)
     
    returns the extremum location p, height y, and half-curvature a
    of a parabolic fit through three points.
    Parabola is given by y(x) = a*(x-p)^2+b,
    where y(-1)=ym1, y(0)=y0, y(1)=yp1.
    @param ym1: y(-1)
    @type ym1: float
    @param y0: y(0)
    @type y0: float
    @param yp1: y(+1)
    @type yp1: float
    @return: peak-ordinate, peak-value
    """
    p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1))
    y = y0 - 0.25*(ym1-yp1)*p
    a = 0.5*(ym1 - 2*y0 + yp1)
    #b = y0 - a*(p**2)
    return p, y, a
            
def subPixelPeak(scoreVolume, coordinates, cubeLength=8, interpolation='Spline', verbose=False):
    """
    subPixelPeak: Will determine the sub pixel area of peak. Utilizes spline, fourier or parabolic interpolation.

    @param verbose: be talkative
    @type verbose: L{str}
    @param scoreVolume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param cubeLength: length of cube - only used for Spline and Fourier interpolation
    @type cubeLength: int (even)
    @param interpolation: interpolation type: 'Spline', 'Quadratic', or 'Fourier'
    @type interpolation: str
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy

    last change: 02/07/2013 FF: 2D functionality added
    """
    assert type(interpolation) == str, 'subPixelPeak: interpolation must be str'
    if (interpolation.lower() == 'quadratic') or (interpolation.lower() == 'parabolic'):
        (peakValue,peakCoordinates) = subPixelPeakParabolic(scoreVolume=scoreVolume, coordinates=coordinates, verbose=verbose)
        return [peakValue,peakCoordinates]
    from pytom_volume import vol,subvolume,rescaleSpline,peak
    from pytom.basic.transformations import resize
  
    #extend function for 2D
    twoD = False
    if scoreVolume.sizeZ() == 1:
        twoD = True

    cubeStart = cubeLength//2
    sizeX = scoreVolume.sizeX()
    sizeY = scoreVolume.sizeY()
    sizeZ = scoreVolume.sizeZ()
    
    if twoD:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1) or\
            (coordinates[0]-cubeStart + cubeLength >= sizeX or coordinates[1]-cubeStart + cubeLength >= sizeY):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(scoreVolume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,0,cubeLength,cubeLength,1)
    else:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1 or coordinates[2]-cubeStart < 1) or \
                (coordinates[0]-cubeStart + cubeLength >= sizeX or coordinates[1]-cubeStart + cubeLength >= sizeY or \
                 coordinates[2]-cubeStart + cubeLength >= sizeZ):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(scoreVolume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,coordinates[2]-cubeStart,
                              cubeLength,cubeLength,cubeLength)
    
    #size of interpolated volume
    scaleSize = 10*cubeLength

    #ratio between interpolation area and large volume
    scaleRatio = 1.0 * cubeLength / scaleSize
    
    #resize into bigger volume
    if interpolation=='Spline':
        if twoD:
            subVolumeScaled = vol(scaleSize,scaleSize,1)
        else:
            subVolumeScaled = vol(scaleSize,scaleSize,scaleSize)
        rescaleSpline(subVolume,subVolumeScaled)
    else:
        subVolumeScaled = resize(volume=subVolume, factor=10)[0]
    
    peakCoordinates = peak(subVolumeScaled)
    
    peakValue = subVolumeScaled(peakCoordinates[0],peakCoordinates[1],peakCoordinates[2])
    
    #calculate sub pixel coordinates of interpolated peak
    peakCoordinates[0] = peakCoordinates[0]*scaleRatio - cubeStart + coordinates[0]
    peakCoordinates[1] = peakCoordinates[1]*scaleRatio - cubeStart + coordinates[1]
    if twoD:
        peakCoordinates[2] = 0
    else:
        peakCoordinates[2] = peakCoordinates[2]*scaleRatio - cubeStart + coordinates[2]
    if ( peakCoordinates[0] > scoreVolume.sizeX() or peakCoordinates[1] > scoreVolume.sizeY() or
            peakCoordinates[2] > scoreVolume.sizeZ() ):
        if verbose:
            print ("SubPixelPeak: peak position too large :( return input value")
        #something went awfully wrong here. return regular value 
        return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]
    
    return [peakValue,peakCoordinates]


def dev(volume, template, mask=None, volumeIsNormalized=False):
    """
    dev: Calculates the squared deviation of volume and template in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: speed up if volume is already normalized
    @type volumeIsNormalized: L{bool}
    @return: deviation
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: FF
    """
    from pytom_volume import vol,sum
    from pytom.tools.macros import volumesSameSize

    assert type(volume) == vol, "dev: volume has to be of type vol!"
    assert type(template) == vol, "dev: template has to be of type vol!"
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')

    if not mask:
        p = volume.numelem()
        result = volume-template
    else:
        assert type(mask) == vol, "dev: mask has to be of type vol!"
        p = sum(mask)
        result = (volume-template)*mask

    deviat = sum(result**2)
    deviat = deviat / float(p)

    return deviat

def POF(volume, template, mask=None, stdV=None):
    """
    POF: Phase only filter correlation function
    @param volume: The search volume
    @type volume: L{pytom_volume.vol}
    @param template: the template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: Only used if specified
    @type mask: L{pytom_volume.vol}
    @param stdV: Will be unused, only for compatibility with FLCF
    @return: POF volume
    @type return: L{pytom_volume.vol}
    @author: Maria Cristina Trueba & Marten Chaillet
    """
    import pytom_volume
    from pytom_volume import vol, vol_comp, pasteCenter, conjugate, sum, real, abs, complexDiv, power, limit
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom.tools.macros import volumesSameSize

    if volume.__class__ != vol and template.__class__ != vol:
        raise RuntimeError('Wrong input type!')

    if volume.sizeX() < template.sizeX() or volume.sizeY() < template.sizeY() or volume.sizeZ() < template.sizeZ():
        raise RuntimeError('Template size is bigger than the target volume')

    # generate the mask
    if mask.__class__ != vol:
        from pytom_volume import initSphere
        mask = vol(template.sizeX(), template.sizeY(), template.sizeZ())
        mask.setAll(0)
        initSphere(mask, template.sizeX() / 2, 0, 0, template.sizeX() / 2, template.sizeY() / 2, template.sizeZ() / 2)
    else:
        if template.sizeX() != mask.sizeX() and template.sizeY() != mask.sizeY() and template.sizeZ() != mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')

    # Normalize template with ampli =1
    ftemplate = fft(template)
    template = ifft(ftemplate / abs(ftemplate))

    # calculate the non-zeros
    p = sum(mask)

    # normalize the template under mask
    meanT = meanValueUnderMask(template, mask, p)
    stdT = stdValueUnderMask(template, mask, meanT, p)

    temp = (template - meanT) / stdT
    temp = temp * mask

    # construct both the template and the mask which has the same size as target volume
    tempV = temp
    if volume.sizeX() != temp.sizeX() or volume.sizeY() != temp.sizeY() or volume.sizeZ() != temp.sizeZ():
        tempV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        tempV.setAll(0)
        pasteCenter(temp, tempV)

    maskV = mask
    if volume.sizeX() != mask.sizeX() or volume.sizeY() != mask.sizeY() or volume.sizeZ() != mask.sizeZ():
        maskV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        maskV.setAll(0)
        pasteCenter(mask, maskV)

    size = volume.numelem()

    # Normalize volume with ampli = 1
    fvolume = fft(volume)
    fvolume = (fvolume / abs(fvolume))

    # calculate the mean and std of volume
    if stdV.__class__ != vol:
        fMask = fft(maskV)
        conjugate(fMask)
        meanV = iftshift(ifft(fMask*fvolume))/(size*p)

        volume = ifft(fvolume)
        stdV = stdUnderMask(volume, maskV, p, meanV)


    fT = fft(tempV)
    conjugate(fT)
    result = iftshift(ifft(fT * fft(volume))) / stdV

    result.shiftscale(0, 1 / (size * p))

    return result

def FPOF(volume, template, mask=None, stdV=None):
    import pytom_volume
    from pytom_volume import vol, pasteCenter, conjugate, sum, real, abs, complexDiv
    from pytom.basic.fourier import fft, ifft, iftshift

    fvolume = fft(volume)
    ftemplate = fft(template)

    amp1_volume = ifft(fvolume/abs(fvolume))
    amp1_template = ifft(ftemplate/abs(ftemplate))

    return FLCF(amp1_volume, amp1_template)

def MCF(volume, template, mask=None, stdV=None):
    """
    POF: Mutual correlation function
    @param volume: The search volume
    @type volume: L{pytom_volume.vol}
    @param template: the template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: Only used if specified
    @type mask: L{pytom_volume.vol}
    @param stdV: Will be unused, only for compatibility with FLCF
    @return: POF volume
    @type return: L{pytom_volume.vol}
    @author: Maria Cristina Trueba
    """

    import pytom_volume
    from pytom_volume import vol, vol_comp, pasteCenter, conjugate, sum, real, abs, complexDiv, power, limit
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom.tools.macros import volumesSameSize

    if volume.__class__ != vol and template.__class__ != vol:
        raise RuntimeError('Wrong input type!')

    if volume.sizeX() < template.sizeX() or volume.sizeY() < template.sizeY() or volume.sizeZ() < template.sizeZ():
        raise RuntimeError('Template size is bigger than the target volume')

    # generate the mask
    if mask.__class__ != vol:
        from pytom_volume import initSphere
        mask = vol(template.sizeX(), template.sizeY(), template.sizeZ())
        mask.setAll(0)
        initSphere(mask, template.sizeX() / 2, 0, 0, template.sizeX() / 2, template.sizeY() / 2, template.sizeZ() / 2)
    else:
        if template.sizeX() != mask.sizeX() and template.sizeY() != mask.sizeY() and template.sizeZ() != mask.sizeZ():
            raise RuntimeError('Template and mask size are not consistent!')

    # Normalize template with ampli =1
    ftemplate = fft(template)
    amp_template = abs(ftemplate)
    power(amp_template, 0.5)
    template = iftshift(ifft(ftemplate/amp_template))

    # calculate the non-zeros
    p = sum(mask)

    # normalize the template under mask
    meanT = meanValueUnderMask(template, mask, p)
    stdT = stdValueUnderMask(template, mask, meanT, p)

    temp = (template - meanT) / stdT
    temp = temp * mask

    # construct both the template and the mask which has the same size as target volume
    tempV = temp
    if volume.sizeX() != temp.sizeX() or volume.sizeY() != temp.sizeY() or volume.sizeZ() != temp.sizeZ():
        tempV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        tempV.setAll(0)
        pasteCenter(temp, tempV)

    maskV = mask
    if volume.sizeX() != mask.sizeX() or volume.sizeY() != mask.sizeY() or volume.sizeZ() != mask.sizeZ():
        maskV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        maskV.setAll(0)
        pasteCenter(mask, maskV)

    size = volume.numelem()

    # Normalize volume with ampli = 1
    fvolume = fft(volume)
    amp_volume = abs(fvolume)
    power(amp_volume, 0.5)
    fvolume = (fvolume/amp_volume)

    # calculate the mean and std of volume
    if stdV.__class__ != vol:
        fMask = fft(maskV)
        conjugate(fMask)
        meanV = iftshift(ifft(fMask * fvolume)) / (size * p)

        volume = iftshift(ifft(fvolume))
        stdV = stdUnderMask(volume, maskV, p, meanV)

    fT = fft(tempV)
    conjugate(fT)
    result = iftshift(ifft(fT * fft(volume))) / stdV

    result.shiftscale(0, 1 / (size * p))

    return result