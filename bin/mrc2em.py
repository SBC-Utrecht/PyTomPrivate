#!/usr/bin/env pytom

from pytom.basic.files import mrc2em


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert mrc file to em.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-f','--file'], 'Filename', 'has arguments', 'optional'),
                                   ScriptOption(['-d','--directory'], 'A directory of files.', 'has arguments', 'optional'),
                                   ScriptOption(['-t','--targetPath'], 'Path to new file.', 'has arguments', 'optional')])

    try:
        filename, directory, target = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    
    if filename:
        #convert only one file
        mrc2em(filename,target)
    elif directory:
        import os
        
        fileList = os.listdir(directory)
        for file in fileList:
            if file[-4:] == '.mrc':
                print directory + os.sep + file , target
                mrc2em(directory + os.sep + file,target)
