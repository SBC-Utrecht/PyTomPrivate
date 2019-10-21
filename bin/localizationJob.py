#!/usr/bin/env pytom
'''
Created on Jan 28, 2013

@author: thrabe
'''


if __name__ == '__main__':
    # parse command line arguments

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists,checkDirExists
    from pytom.basic.structures import Mask, Reference, Wedge, BandPassFilter
    from pytom.localization.structures import Volume
    from pytom.angles.globalSampling import GlobalSampling
    from pytom.score.score import FLCFScore
    from pytom.localization.peak_job import PeakJob
    from pytom.frontend.serverpages.createLocalizationJob import createRunscripts

    helper = ScriptHelper(
        sys.argv[0].split('/')[-1],
        description='Create a localization job. Documentation is available at\n\
        http://www.pytom.org/doc/pytom/localization.html',
        authors='Thomas Hrabe',
        options=[ScriptOption(['-v','--volume'], 'Volume : the big volume', 'string', 'required'),
                 ScriptOption(['-r','--reference'], 'Reference : the molecule searched', 'string', 'required'),
                 ScriptOption(['-m','--mask'], 'Mask : a mask ', 'string', 'required'),
                 ScriptOption(['--wedge1'], 'Wedge : first tilt angle. Must be 90-tilt!', 'float', 'required'),
                 ScriptOption(['--wedge2'], 'Wedge : second tilt angle.  Must be 90-tilt!', 'float', 'required'),
                 ScriptOption(['-a','--angles'], '''Angles : name of angle list. Either : 
                  angles_50_100.em
                  angles_38.53_256.em
                  angles_35.76_320.em
                  angles_25.25_980.em
                  angles_19.95_1944.em
                  angles_18_3040.em    
                  angles_12.85_7112.em    
                  angles_11_15192.em    
                  angles_07_45123.em
                  angles_3_553680.em
                  ''', 'string', 'required'),
                 ScriptOption(['-d','--destination'], 'Destination : destination directory', 'string', 'required'),
                 ScriptOption(['-b','--band'], 'Lowpass filter : band - in pixels', 'float', 'required'),
                 ScriptOption(['--splitX'], 'Into how many parts do you want to split volume (X dimension)', 'int', 'required'),
                 ScriptOption(['--splitY'], 'Into how many parts do you want to split volume (Y dimension)', 'int', 'required'),
                 ScriptOption(['--splitZ'], 'Into how many parts do you want to split volume (Z dimension)', 'int', 'required'),
                 ScriptOption(['-j','--jobName'], 'Specify job.xml filename', 'string', 'required')])

    volume, reference, mask, wedge1,wedge2,angles,destination,band,sx,sy,sz,jobName = parse_script_options(sys.argv[1:], helper)
    
    if not checkFileExists(volume):
        raise RuntimeError('Volume file ' + volume + ' does not exist!')
    
    if not checkFileExists(reference):
        raise RuntimeError('Reference file ' + reference + ' does not exist!')
    
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')
    
    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')
    
    v = Volume(volume)
    r = Reference(reference)
    m = Mask(mask)
    w = Wedge([wedge1, wedge2])
    a = GlobalSampling(angles)
    
    job = PeakJob(volume=v, reference=r, mask=m, wedge=w, rotations=a, score=FLCFScore(), jobID=0, members=1, dstDir=destination, bandpass=BandPassFilter(0,band,0))
    
    job.toXMLFile(jobName)
    
    createRunscripts(jobName[:-3] + 'sh', jobName, sx, sy, sz)
