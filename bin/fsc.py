#!/usr/bin/env pytom

'''
Created on Jul 21, 2011

@author: hrabe
'''


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.correlation import FSC,determineResolution
    from pytom_volume import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Determine resolution by FSC.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption('--v1', 'First volume path.', 'string', 'optional'),
                                   ScriptOption('--v2', 'Second volume path.', 'string', 'optional'),
                                   ScriptOption('--pl', 'A particleList if v1 and v2 are not available.', 'string', 'optional'),
                                   ScriptOption('--fsc', 'The FSC criterion. Value between 0.0 and 1.0. Standard values are 0.5 or 0.3', 'float', 'optional', 0.5),
                                   ScriptOption('--numberBands', 'Number of bands (optional). If not set, numberBands = cubesize/4.', 'int', 'optional'),
                                   ScriptOption(['-m','--mask'], 'Mask (optional, but recomended).', 'string', 'optional'),
                                   ScriptOption('--pixelsize', 'Pixelsize in Angstrom (optional). Will return resolution in Angstrom.', 'float', 'required'),
                                   ScriptOption('--xml', 'Output in XML. (optional) ', 'no arguments', 'optional'),
                                   ScriptOption(['-v','--verbose'], 'Verbose data. (optional) ', 'no arguments', 'optional'),
                                   ScriptOption(['--noShift'], 'Removes the shift from the particlelist.', 'no arguments', 'optional')])


    v1Filename, v2Filename, particleList, fscCriterion, numberBands, mask, pixelSize, xml, verbose, noshift = parse_script_options(sys.argv[1:], helper)
    
    try:
        mask= read(mask) 
    except:
        mask = None

    
        
    if v1Filename and v2Filename:    
        v1 = read(v1Filename)
        v2 = read(v2Filename)
        
        if not numberBands:
            numberBands = int(v1.sizeX()/2)
        
        f = FSC(v1,v2,numberBands,mask,verbose)
        print 'FSC:\n', f
        
        r = determineResolution(f,fscCriterion,verbose)
    elif particleList:
        from pytom.basic.structures import ParticleList
        
        pl = ParticleList('.')
        pl.fromXMLFile(particleList)
        # Remove all shifts, useful for particle polishing
        if noshift:
            print('set shift to zero')
            for particle in pl:
                particle.getShift().scale(0)


        if len(pl) <= 1:
            raise RuntimeError('There is only 1 or less particles in the particle list. Need at least two! Abort!')
        
        if not numberBands:
            p = pl[0]
            pv = p.getVolume()
            numberBands = int(pv.sizeX()/4)
            
        r = pl.determineResolution(fscCriterion,numberBands,mask,verbose=verbose,plot='',keepHalfsetAverages = True,halfsetPrefix='plFSC')
        print 'Even and odd halfsets were written into current directory and stored as plFSCeven / odd .em!'
        
    else:
        raise RuntimeError('You must provide either two files or a particlelist as parameters for determining resolution!')
        
    if not xml:
        if v1Filename and v2Filename:
            print 'Resolution determined for ', v1Filename, ' and ', v2Filename
        elif particleList:
            print 'Resolution determined for ', particleList, ' '
            
        print ''
        print 'FSC Criterion:   ', fscCriterion
        print 'Number of Bands: ', numberBands
        print ''
        print 'Nyquist: ', r[0]
        print 'Band:    ', r[1]
        
        
        if pixelSize:
            from pytom.basic.resolution import bandToAngstrom
            
            print 'Resolution determined for pixelsize : ', pixelSize , ' at ', bandToAngstrom(r[1],pixelSize,numberBands), ' Angstrom'
    
    else:
        print 'XML'
