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
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], 
                          description='Create an EXMX alignment job. Documentation is available at\n \
                          http://www.pytom.org/doc/pytom/classification.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-p','--particleList'], 'Particle list : xml file storing information to all subvolumes', 'has arguments', 'required'),
                                   ScriptOption(['-m','--mask'], 'Mask : a mask ', 'has arguments', 'required'),
                                   ScriptOption(['-c','--classes'], 'Number of classes', 'has arguments', 'required'),
                                   ScriptOption(['-t','--endThreshold'], 'End threshold', 'has arguments', 'required'),
                                   ScriptOption(['--wedge1'], 'Wedge : first tilt angle. Must be 90-tilt!', 'has arguments', 'required'),
                                   ScriptOption(['--wedge2'], 'Wedge : second tilt angle.  Must be 90-tilt!', 'has arguments', 'required'),
                                   ScriptOption(['--symmetry'], 'PointSymmetry : specify n-fold symmetry (n)', 'has arguments', 'optional'),
                                   ScriptOption(['--symmetryAngleZ'], 'PointSymmetry axis tilt around Z axis', 'has arguments', 'optional'),
                                   ScriptOption(['--symmetryAngleX'], 'PointSymmetry axis tilt around X axis', 'has arguments', 'optional'),
                                   ScriptOption(['-l','--lowestFrequency'], 'Highest frequency band used : in pixels', 'has arguments', 'required'),
                                   ScriptOption(['-h','--highestFrequency'], 'Lowest frequency band used : in pixels', 'has arguments', 'required'),
                                   ScriptOption(['-d','--destination'], 'Destination : destination directory', 'has arguments', 'required'),
                                   ScriptOption(['--startTemperature'], 'Start temperature : ', 'has arguments', 'optional'),
                                   ScriptOption(['--annealingStep'], 'Annealing step : temperature drop ', 'has arguments', 'optional'),
                                   ScriptOption(['-n','--numberRefinementIterations'], 'Number of local refinement iterations', 'has arguments', 'required'),
                                   ScriptOption(['-b','--binning'], 'Perform binning (downscale) of subvolumes: 1: no binning, 2: 2 pixels -> 1, 3: 3 -> 1 ...', 'has arguments', 'required'),
                                   ScriptOption(['--pixelSize'], 'Pixelsize in Angstrom', 'has arguments', 'required'),
                                   ScriptOption(['--particleDiameter'], 'Particle diameter in Angstrom', 'has arguments', 'required'),
                                   ScriptOption(['-j','--jobName'], 'Specify job.xml filename', 'has arguments', 'required')])

    try:
        particleList, mask, numberClasses, endThreshold,wedge1,wedge2,symmetryN,symmetryAxisZ,symmetryAxisX,\
        lowestFrequency,highestFrequency,destination,\
        startTemperature,annealingStep,numberRefinementIterations,binning,\
        pixelSize,diameter,jobName = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
   
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')
    
    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')    

    from pytom.cluster.mcoACStructures import MCOACJob,MetropolisCriterion,SigmaTemperature 
    from pytom.basic.structures import ParticleList,Reference,Mask,Wedge,SampleInformation,PointSymmetry
    from pytom.score.score import FLCFScore
    from pytom.frontend.serverpages.createMCOACJob import createRunscripts
    from pytom.alignment.preprocessing import Preprocessing
     
    p       = ParticleList()
    p.fromXMLFile(particleList)
    
    m       = Mask(mask)
    w       = Wedge([float(wedge1),float(wedge2)])
    pre     = Preprocessing(lowestFrequency = float(lowestFrequency),highestFrequency = float(highestFrequency))
    sample  = SampleInformation(pixelSize = float(pixelSize), particleDiameter = float(diameter))
    temperature = SigmaTemperature(startTemperature,annealingStep)
    
    if symmetryN is None or symmetryAxisZ is None or symmetryAxisX is None:
        sym = None
    else:
        sym = PointSymmetry(nfold=int(symmetryN),z2=float(symmetryAxisZ),x=float(symmetryAxisX))

    job = MCOACJob(particleList=p,destinationDirectory=destination,mask=mask,\
                   score=FLCFScore(),preprocessing=pre,wedgeInfo=w,binning=int(binning),\
                   sampleInformation=sample,numberClasses=int(numberClasses),\
                   temperature=temperature,criterion=MetropolisCriterion(),localSearchIncrement=numberRefinementIterations,\
                   endThreshold=float(endThreshold),symmetry = sym)
    
    job.toXMLFile(jobName)
    
    createRunscripts(jobName[:-3] + 'sh',jobName)
    
    
    
    
    
