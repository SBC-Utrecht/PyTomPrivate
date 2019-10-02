#!/usr/bin/env pytom
"""
Created on Jul 20, 2013

@author: FF
"""

if __name__ == '__main__':
    import sys
    #from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct

    options=[ScriptOption(['--tiltSeriesName'], 'Name tilt series - either prefix of sequential tilt series files \
             expected as "tiltSeriesName_index.em/mrc" or full name of stack "tiltSeriesName.st"',
                          'has arguments', 'required'),
             ScriptOption(['--tiltSeriesFormat'], 'Format of tilt series (series of "em" or "mrc" images or "st" stack).',
                          'has arguments', 'optional'),
             ScriptOption(['--firstIndex'], 'Index of first projection.', 'has arguments', 'optional'),
             ScriptOption(['--lastIndex'], 'Index of last projection.', 'has arguments', 'required'),
             ScriptOption(['--tltFile'], 'tltFile containing tilt angles.', 'has arguments', 'optional'),
             ScriptOption(['--prexgFile'], 'prexgFile containing pre-shifts from IMOD.', 'has arguments', 'optional'),
             ScriptOption(['--preBin'], 'pre-Binning in IMOD prior to marker determination.', 'has arguments', 'optional'),
             ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', 'has arguments',
                          'required'),
             ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.',
                          'has arguments', 'required'),
             ScriptOption(['--referenceMarkerIndex'], 'Index of reference marker to set up coordinate system.',
                          'has arguments', 'required'),
             ScriptOption(['--handflip'], 'Is your tilt series outside of 0-180deg (Specify if yes).', 'no arguments',
                          'optional'),
             ScriptOption(['--projectionTargets'],
                          'Relative or absolute path to the aligned projections that will be generated + file prefix.\
                          default: "align/myTilt"', 'has arguments', 'optional'),
             ScriptOption(['--fineAlignFile'],
                          'Relative or absolute path to the file with fineAlign parameters (type should be *.dat).',
                          'has arguments', 'optional'),
             ScriptOption(['--projectionBinning'], 'Binning of projections during read - default: 1.', 'has arguments',
                          'optional'),
             ScriptOption(['--lowpassFilter'], 'Lowpass filter in Nyquist after binning.', 'has arguments', 'required'),
             ScriptOption(['--tomogramFile'],
                          'Relative or absolute path to final tomogram (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--fileType'], 'File type (can be em or mrc - no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--tomogramSizeX'], 'Size of tomogram in x (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--tomogramSizeY'], 'Size of tomogram in y (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--tomogramSizeZ'], 'Size of tomogram in z (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--reconstructionCenterX'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--reconstructionCenterY'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--reconstructionCenterZ'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          'has arguments', 'optional'),
             ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', 'has arguments',
                          'optional'),
             ScriptOption(['--verbose'], 'Enable verbose mode', 'no arguments', 'optional')]
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Align and weight projections, save them and reconstruct tomogram (optional). \n\
                                      See http://pytom.org/doc/pytom/reconstructTomograms.html for documentation.',
                          authors='Friedrich Foerster',
                          options = options)

    try:
        tiltSeriesName, tiltSeriesFormat, firstProj, lastProj, \
        tltFile, prexgFile, preBin, referenceIndex, markerFileName, referenceMarkerIndex, handflip, \
        projectionTargets, fineAlignFile, projBinning, lowpassFilter, \
        volumeName, filetype, \
        tomogramSizeX, tomogramSizeY, tomogramSizeZ, \
        reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ, \
        weightingType, verbose = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print sys.version_info
        print e
        sys.exit()

    # input parameters
    #tiltSeriesName = tiltSeriesPath + tiltSeriesPrefix  # "../projections/tomo01_sorted" # ending is supposed to be tiltSeriesName_index.em (or mrc)
    if not tiltSeriesFormat:
        tiltSeriesFormat = 'em'
    if firstProj:
        firstProj = int(firstProj)  # 1 # index of first projection
    else:
        firstProj = 1
    lastProj = int(lastProj)  # 41 # index of last projection
    ireftilt = int(referenceIndex)  # 21 # reference projection (used for alignment)
    if referenceMarkerIndex:
        irefmark = int(referenceMarkerIndex)  # reference marker (defines 3D coordinate system)
    else:
        irefmark = 1

    handflip = handflip is not None  # False # is your tilt axis outside 0-180 deg?
    # output parameters
    if projectionTargets:
        # weighted and aligned projections are stored as alignedTiltSeriesName_index.em
        alignedTiltSeriesName = projectionTargets
    else:
        alignedTiltSeriesName = 'align/myTilt'
    if projBinning:
        projBinning = int(projBinning)  # binning factor
    else:
        projBinning = 1
    if lowpassFilter:
        lowpassFilter = float(lowpassFilter)  # lowpass filter in Nyquist (post-binning)
    else:
        lowpassFilter = 1.

    if weightingType is None:
        weightingType = -1
    else:
        weightingType = int(weightingType)
    # only write projections and do NOT reconstruct tomogram (following parameters would be obsolete)
    onlyWeightedProjections = False
    if not volumeName:
        onlyWeightedProjections = True
         
    if filetype is None:
        if volumeName:
            filetype = volumeName.split('.')[-1]
    else:
        filetype = 'em'

    if tomogramSizeX is not None or tomogramSizeY is not None or tomogramSizeZ is not None:
        # dimensions of reconstructed tomogram
        voldims = [int(tomogramSizeX), int(tomogramSizeY), int(tomogramSizeZ)]
    else:
        onlyWeightedProjections = True
        voldims = [0, 0, 0]       # offset from center of volume - for example choose z!=0 to shift in z (post-binning coordinates)
    if reconstructionCenterX:
        reconstructionCenterX = int(reconstructionCenterX)
    else:
        reconstructionCenterX = 0
    if reconstructionCenterY:
        reconstructionCenterY = int(reconstructionCenterY)
    else:
        reconstructionCenterY = 0
    if reconstructionCenterZ:
        reconstructionCenterZ = int(reconstructionCenterZ)
    else:
        reconstructionCenterZ = 0
    reconstructionPosition = [reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ]
    if preBin:
        preBin=int(preBin)

    outMarkerFileName = 'MyMarkerFile.em'
    if verbose:
        print "Tilt Series: "+str(tiltSeriesName)+", "+str(firstProj)+"-"+str(lastProj)
        print "Index of Reference Projection: "+str(referenceIndex)
        print "Marker Filename: "+str(markerFileName)
        print "TltFile: "+str(tltFile)
        print "prexgFile: "+str(prexgFile)
        print "Index of Reference Marker: "+str(referenceMarkerIndex)
        print "Handflip: "+str(handflip)
        print "Projection Targets: "+str(projectionTargets)
        print "FineAlignmentFile: "+str(fineAlignFile)
        print "Binning Factor of Projections: "+str(projBinning)+", lowpass filter (in Ny): "+str(lowpassFilter)
        print "Name of Reconstruction Volume: "+str(volumeName)+" of Filetype: "+str(filetype)
        print "Reconstruction size: "+str(voldims)
        print "Reconstruction center: "+str(reconstructionPosition)
        print "write only aligned projections out: "+str(onlyWeightedProjections)

    alignWeightReconstruct(tiltSeriesName=tiltSeriesName, markerFileName=markerFileName, lastProj=lastProj,
                           tltfile=tltFile, prexgfile=prexgFile, preBin=preBin,
                           volumeName=volumeName, volumeFileType=filetype,
                           voldims=voldims, recCent=reconstructionPosition,
                           tiltSeriesFormat=tiltSeriesFormat, firstProj=firstProj, irefmark=irefmark, ireftilt=ireftilt,
                           handflip=handflip,
                           alignedTiltSeriesName=alignedTiltSeriesName,
                           weightingType=weightingType,
                           lowpassFilter=lowpassFilter, projBinning=projBinning,
                           outMarkerFileName=outMarkerFileName, verbose=True)


