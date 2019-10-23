#!/usr/bin/env pytom

"""
Created on Aug 24, 2011

@author: hrabe
"""



if __name__ == '__main__':
# parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList, Rotation, PointSymmetry
    from pytom_volume import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Generate a symmetrized density.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-p', '--particleList'], 'Aligned XML particle list.', 'string', 'optional'),
                                   ScriptOption(['-v', '--volume'], 'Volume to be symmetrized.', 'string', 'optional'),
                                   ScriptOption(['-r', '--result'], 'Resulting symmetry filename.', 'string', 'required'),
                                   ScriptOption(['-s', '--symmetries'], 'How many symmetries', 'int', 'required'),
                                   ScriptOption(['--z1Rotation'], 'First rotation around z axis', 'float', 'optional', 0),
                                   ScriptOption(['--xRotation'], 'Rotation around x axis', 'float', 'optional', 0),
                                   ScriptOption(['--z2Rotation'], 'Second rotation around z axis', 'float', 'optional', 0)])

    particleListName, volumeName, result, numberSymmetries, z1, x, z2 = parse_script_options(sys.argv[1:], helper)
        
    if particleListName:
        pl = ParticleList('.')
        pl.fromXMLFile(particleListName)
    elif volumeName:
        volume = read(volumeName)
    else:        
        raise RuntimeError('You must specify either a particle list or a volume file for symmetrization.')
        
    symmetry = PointSymmetry(numberSymmetries, z2, x)
    
    if particleListName:
        newParticleList = symmetry.apply(pl)
        newParticleList.average(result)
    elif volumeName:
        newVolume = symmetry.applyToParticle(volume, Rotation(0, z2, x))  # the 0 could maybe be meant to be z1, remark in passing the code, I could not know, Douwe Schulte October 2019
        newVolume.write(result)
