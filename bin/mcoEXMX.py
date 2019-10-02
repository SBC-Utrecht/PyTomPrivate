#!/usr/bin/env pytom
if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.cluster.mcoEXMXStructures import MCOEXMXJob 
    from pytom.cluster.mcoEXMX import mcoEXMX 
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run a mcoEXMX job. Documentation is available at\n \
                          http://www.pytom.org/doc/pytom/classification.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-j','--job'], 'Job', 'has arguments', 'required'),
                                   ScriptOption(['-v','--verbose'], 'Verbose', 'has arguments', 'required')])


    try:
        jobFile, verbose = parse_script_options(sys.argv[1:], helper)
    except Exception:
#        print e
        sys.exit()

    verbose = verbose == 'True'
    
    job = MCOEXMXJob(0,0,0,0,0,0,0,0,0,0,0,0)
    job.fromXMLFile(jobFile)

    mcoEXMX(job,True,verbose)