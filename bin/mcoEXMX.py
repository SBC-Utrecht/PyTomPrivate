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
                          options=[ScriptOption(['-j', '--job'], 'Jobfile', 'string', 'required'),
                                   ScriptOption(['-v', '--verbose'], 'Verbose', 'no arguments', 'optional', False)])

    jobFile, verbose = parse_script_options(sys.argv[1:], helper)
    
    job = MCOEXMXJob(0,0,0,0,0,0,0,0,0,0,0,0)
    job.fromXMLFile(jobFile)

    mcoEXMX(job, True, verbose)