#!/usr/bin/env pytom
'''
Created on Nov 3, 2011

@author: hrabe
'''

if __name__ == '__main__':
    # parse command line arguments
    
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options 
    from pytom.alignment.ExMaxAlignment import parallelStart,ExMaxJob 
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run an alignment job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/alignment.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-j','--job'], 'Job', 'has arguments', 'required'),
                                   ScriptOption(['-v','--verbose'], 'Verbose', 'no arguments','optional', False)])

    verbose = False

    jobFile, verbose = parse_script_options(sys.argv[1:], helper)

    #from pytom.tools.timing import Timing
    #t = Timing()
    #t.start()

    exj = ExMaxJob()
    exj.fromXMLFile(jobFile)

    parallelStart(exj,verbose = verbose)
    
    #print 'Overall execution time: %f s.' % t.end()