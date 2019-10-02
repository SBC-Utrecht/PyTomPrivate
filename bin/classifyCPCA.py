#!/usr/bin/env pytom
'''
Created on Jan 21, 2013

@author: FF
'''

from pytom.classification.CPCAfunctions import *

def doClassification(pl, cpl, ccc, neig, nclass, cName):
    """
    """
    subTomoClust(particleListFilename=pl, classifiedParticleListFilename=cpl,
            cccName=ccc, neig=neig, nclass=nclass)
    averageClasses(particleListFilename=cpl, avName=cName)



if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Classification of particle list using CPCA. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/classification.html',
                          authors='Friedrich Foerster',
                          options=[ScriptOption(['-p','--particleList'], 
			               'ParticleList', 'has arguments', 'optional'),
                                   ScriptOption(['-o','--outputParticleList'], 
				       'classified Particle List.', 'has arguments', 'optional'),
                                   ScriptOption(['-c','--ccc'], 
				       'constrained correlation matrix.', 'has arguments', 'optional'),
                                   ScriptOption(['-e','--neig'], 
				       'number of eigenvectors.', 'has arguments', 'optional'),
                                   ScriptOption(['-n','--nclass'], 
				       'number of classes.', 'has arguments', 'optional'),
                                   ScriptOption(['-a','--average'], 
				       'name for class averages.', 'has arguments', 'optional')])

    try:
        pl, cpl, ccc, neig, nclass, cName = parse_script_options(sys.argv[1:], helper)
    except Exception:
        sys.exit()

    neig = int(neig)
    nclass = int(nclass)
    if not cName:
        cName = 'class'
    doClassification(pl=pl, cpl=cpl, ccc=ccc, neig=neig, nclass=nclass, cName=cName)
 
