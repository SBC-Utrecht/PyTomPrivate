#!/usr/bin/env pytom

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.cluster.mcoACStructures import MCOACJob 
    from pytom.cluster.mcoAC import mcoAC 
    import pytom_mpi
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run a mcoAC job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/classification.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-j','--job'], 'Job', 'has arguments', 'optional'),
                                   ScriptOption(['-v','--verbose'], 'Verbose', 'no arguments', 'required'),
                                   ScriptOption(['-h', '--help'], 'Help.', 'no arguments', 'optional')])

    verbose = False

    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        jobFile, verbose ,helpme = parse_script_options(sys.argv[1:], helper)
    except Exception:
        #print e
        sys.exit()
    if helpme is True:
        print helper
        sys.exit()

    job = MCOACJob(0,0,0,0,0,0,0,0,0,0,0,0)
    job.fromXMLFile(jobFile)
    
    verbose = verbose is True
    try:
        mcoAC(job,True,verbose)
    except:
        if pytom_mpi.rank() == 0:
            from pytom.alignment.ExMaxAlignment import ExMaxManager
            manager = ExMaxManager(annealingJob.getExMaxJob())
            manager.parallelEnd()
            
        pytom_mpi.finalise()