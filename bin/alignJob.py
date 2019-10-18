#!/usr/bin/env pytom
'''
Created on Jan 30, 2013

@author: Thomas Hrabe
'''

if __name__ == '__main__':
    # parse command line arguments

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.alignment.ExMaxAlignment import ExMaxJob
    from pytom.basic.structures import ParticleList, Reference, Mask, SampleInformation, PointSymmetry
    from pytom.score.score import FLCFScore
    from pytom.frontend.serverpages.createAlignmentJob import createRunscripts
    from pytom.angles.localSampling import LocalSampling
    from pytom.alignment.preprocessing import Preprocessing
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], 
                          description='Create an alignment job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/alignment.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-p', '--particleList'], 'Particle list : xml file storing information to all subvolumes', 'string', 'required'),
                                   ScriptOption(['-r', '--reference'], 'Reference : the alignment reference', 'string', 'required'),
                                   ScriptOption(['-m', '--mask'], 'Mask : a mask ', 'string', 'required'),
                                   ScriptOption(['--angleShells'], '# angle shells used for angular refinement.', 'has arguments', 'required'),
                                   ScriptOption(['--angleIncrement'], 'Angular increment for refinement.', 'has arguments', 'required'),
                                   ScriptOption(['--symmetry'], 'PointSymmetry : specify n-fold symmetry (n)', 'uint', 'optional'),
                                   ScriptOption(['--symmetryAngleZ'], 'PointSymmetry axis tilt around Z axis', 'float', 'optional'),
                                   ScriptOption(['--symmetryAngleX'], 'PointSymmetry axis tilt around X axis', 'float', 'optional'),
                                   ScriptOption(['-l', '--lowestFrequency'], 'Highest frequency band used : in pixels', 'float', 'required'),
                                   ScriptOption(['-h', '--highestFrequency'], 'Lowest frequency band used : in pixels', 'float', 'required'),
                                   ScriptOption(['-d', '--destination'], 'Destination : destination directory', 'string', 'required'),
                                   ScriptOption(['-n', '--numberIterations'], 'Number of iterations', 'has arguments', 'required'),
                                   ScriptOption(['-b', '--binning'], 'Perform binning (downscale) of subvolumes: 1: no binning, 2: 2 pixels -> 1, 3: 3 -> 1 ...', 'uint', 'required'),
                                   ScriptOption(['--pixelSize'], 'Pixelsize in Angstrom', 'float', 'required'),
                                   ScriptOption(['--particleDiameter'], 'Particle diameter in Angstrom', 'float', 'required'),
                                   ScriptOption(['-j', '--jobName'], 'Specify job.xml filename', 'string', 'required')])

    particleList, reference, mask, angShells,angleInc,symmetryN,symmetryAxisZ,symmetryAxisX,\
    lowestFrequency,highestFrequency,destination,numberIterations,binning,\
    pixelSize,diameter,jobName = parse_script_options(sys.argv[1:], helper)
        
    if not checkFileExists(particleList):
        raise RuntimeError('ParticleList file ' + volume + ' does not exist!')
    
    if not checkFileExists(reference):
        raise RuntimeError('Reference file ' + reference + ' does not exist!')
    
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')
    
    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')    

    p      = ParticleList()
    p.fromXMLFile(particleList)
    r      = Reference(reference)
    m      = Mask(mask)
    a      = LocalSampling(angShells,angleInc)
    pre    = Preprocessing(lowestFrequency = lowestFrequency,highestFrequency = highestFrequency)
    sample = SampleInformation(pixelSize = pixelSize, particleDiameter = diameter)
    
    if symmetryN is None or symmetryAxisZ is None or symmetryAxisX is None:
        sym = None
    else:
        sym = PointSymmetry(nfold=symmetryN,z2=symmetryAxisZ,x=symmetryAxisX)
    
    job = ExMaxJob(particleList=p,destination=destination,reference=r,score=FLCFScore(),
                   rotations=a,mask=m,symmetry=sym,
                   numberRefinementRounds=None,numberIterations=numberIterations,preprocessing=pre,
                   excludeThreshold=-1,binning=binning,sampleInformation=sample,fscCriterion=0.5,
                   adaptiveResolution=True,adaptiveOffset=0.1,angleFactor=0.5)

    
    job.toXMLFile(jobName)
    
    createRunscripts(jobName[:-3] + 'sh',jobName)
    
    
    
    
    
