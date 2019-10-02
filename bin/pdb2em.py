#!/usr/bin/env pytom
'''
Created on Nov 8, 2013

@author: thrabe
'''

if __name__ == '__main__':
 # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options 
    from pytom.basic.files import pdb2em
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Compile a electron density from PDB file\n\
                          http://pytom.org/doc/pytom/files.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-p','--pdbFile'], 'A PDB file', 'has arguments', 'required'),
                                   ScriptOption(['-c','--chain'], 'A Chain', 'has arguments', 'optional'),
                                   ScriptOption(['-s','--pixelSize'], 'Pixel size of output volume (in Angstrom)', 'has arguments', 'optional'),
                                   ScriptOption(['-v','--volumeSize'], 'Volume length (size) in all dimensions', 'has arguments', 'optional'),
                                   ScriptOption(['-o','--outputVolumePath'], 'Path to output volume ', 'has arguments', 'required'),
                                   ScriptOption(['-i','--invertDensity'],'Set if density should be negative', 'no arguments', 'required')])

    try:
        pdbFile, chain, pixelSize, cubeSize, volumePath ,densityNegative = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    volume = pdb2em(pdbFile, float(pixelSize), int(cubeSize), chain = chain,densityNegative = densityNegative)
    
    volume.write(volumePath)