#!/usr/bin/env pytom

def startLocalizationJob(filename, splitX=0, splitY=0, splitZ=0, doSplitAngles=False):
    """
    @author: chen
    """
    verbose=True
    from pytom.localization.peak_job import PeakJob
    job = PeakJob()
    job.fromXMLFile(filename)
    job.check()

    if doSplitAngles:
        print 'Ignore split volume parameters ...'
        from pytom.localization.parallel_extract_peaks import PeakManager
        manager = PeakManager()
        manager.parallelStart_splitAng(job, verbose)
    else:
        from pytom.localization.parallel_extract_peaks import PeakLeader
        leader = PeakLeader()
        leader.parallelRun(job, splitX, splitY, splitZ, verbose)

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run a localization job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/localization.html',
                          authors='Yuxiang Chen, Thomas Hrabe',
                          options=[ScriptOption(['-j', '--jobName'], 'Specify job.xml filename', 'string', 'required'),
                                   ScriptOption(['-x', '--splitX'], 'Parts you want to split the volume in X dimension', 'int', 'optional', 0),
                                   ScriptOption(['-y', '--splitY'], 'Parts you want to split the volume in Y dimension', 'int', 'optional', 0),
                                   ScriptOption(['-z', '--splitZ'], 'Parts you want to split the volume in Z dimension', 'int', 'optional', 0)])
    
    try:
        jobName, splitX, splitY, splitZ = parse_script_options(sys.argv[1:], helper)
        
        if jobName is None:
            raise RuntimeError()
        
    except:  # backward compatibility, does not work anymore because of the new scripthelper (2019)
        
        jobName = sys.argv[1]
        
        if len(sys.argv) == 5:
            splitX = int(sys.argv[2])
            splitY = int(sys.argv[3])
            splitZ = int(sys.argv[4])

    
    from pytom.tools.timing import Timing
    t = Timing(); t.start()
    
    startLocalizationJob(jobName, splitX, splitY, splitZ, doSplitAngles=False)
    
    time = t.end(); print 'The overall execution time: %f' % time
