#!/usr/bin/env pytom
"""
    little script to create new particle list of particles that belong to selected classes
    GvdS May 2019
"""
from pytom.basic.structures import ParticleList, Rotation

from copy import deepcopy
import random
import numpy
#from pytom_volume import read
#from pytom_numpy import vol2npy
import os
import lxml.etree as et


def remove_element(el):
    parent = el.getparent()
    if el.tail.strip():
        prev = el.getprevious()
        if prev:
            prev.tail = (prev.tail or '') + el.tail
        else:
            parent.text = (parent.text or '') + el.tail
    parent.remove(el)

def extractParticleListsByTomoNameFromXML(xmlfile, directory='./', query='all', prefix='particleList_'):

    if not query:
        query_list = []
    else:
        query_list = query.split(',')

    excludeList = []
    outfiles = []
    pL = [1]
    while len(pL):
        tree = et.parse(xmlfile)
        tomogram = ''
        for n , particle in enumerate( tree.xpath("Particle") ):
            remove = True
            origin = particle.xpath('PickPosition')[0].get('Origin')
            if not origin:
                origin = particle.xpath('InfoTomogram')[0].get('TomoName')

            if not tomogram and not origin in excludeList: 
                tomogram = origin
                
            if not tomogram == origin or origin in excludeList:
                remove_element(particle)

        excludeList.append(tomogram)
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            tomoname = os.path.basename(tomogram).replace('.mrc','').replace('.em','')
            outfile = "{}{}_{}.xml".format(prefix, os.path.basename(xmlfile)[:-4], tomoname)
            outfile = os.path.join(directory, outfile)
            if False == (tomoname in query_list or 'all' in query_list or not query_list):
                tomoname=''
            if tomoname:
                tree.write(outfile, pretty_print=True)
                outfiles.append(outfile)
                print(outfile)
        except Exception as e:
            print(e)
            print('writing {} failed.'.format(xmlfile))
            print('No file written.')

        pL = tree.xpath('Particle')
        #print(os.path.basename(tomogram), len(pL))

    return outfiles

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-t', '--directory'], 'Target Directory', True, True),
               ScriptOption(['--prefix'], 'Prefix that is prepended to your output file name (in the target dir)', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert particle list to n particle list based on tomoname.',
                          authors='Gijs van der Schot',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, target, prefix, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    if not os.path.exists(plName):
        sys.exit()

    target = './' if target is None else target
    prefix = '' if prefix is None else prefix
    if not os.path.exists(target) or not os.path.isdir(target):
        sys.exit()

    
    extractParticleListsByTomoNameFromXML(plName, target, prefix=prefix)
