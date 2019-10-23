#!/usr/bin/env pytom

"""
Created on October, 2019

@author: Douwe Schulte

To join the conversion scripts scattered all around the bin folder into one single script.
"""

import pytom.basic.files as f
import os

# The lookup tables for the conversion functions
extensions = ["em", "mrc", "ccp4"]
lookuptable = [
    ["Same",        f.em2mrc,  f.em2ccp4],
    [f.mrc2em,      "Same",    f.mrc2ccp4],
    [f.ccp42em,     f.ccp42em, "Same"],
]
special_extensions = ["pdb", "mmCIF"]
special_lookuptable = [f.pdb2em, f.mmCIF2em]


def special(file, format, target, chaindata):
    """Converts a special fileformat to the given format"""
    in_ex = file.split('.')[-1]

    func = special_lookuptable[special_extensions.index(in_ex)]

    volume = func(file, chaindata[0], chaindata[1], invertDensity=chaindata[2], chain=chaindata[3])

    # Write it as an em file
    em_name = f.name_to_format(file, target, "em")
    f.write_em(em_name, volume)

    if format != "em":
        # Convert em file to target
        convertfile(em_name, target)
        os.remove(em_name)


def convertfile(file, format, target, chaindata):
    """Convert a file to the given format, returns (statuscode, errormessage)"""
    in_ex = file.split('.')[-1]

    if in_ex in special_extensions:
        if chaindata is None:
            return 1, "For conversion from special format ({:s} and {:s}) chaindata is needed".format(", ".join([a for a in special_extensions[:-1]]), special_extensions[-1])
        special(file, target, chaindata)
        return 0, ""

    try:
        index_from = extensions.index(in_ex)
    except:
        return -1, "The given input file format is not supported: {:s}".format(in_ex)

    index_to = extensions.index(format)

    func = lookuptable[index_from][index_to]

    if func == "Same":
        return -1, "Converting from {:s} to the same format seems not needed".format(in_ex)
    else:
        func(file, target)
        return 0, ""


def print_errors(message):
    """To print to stderr and use fancy styling, for multiple errors"""
    import sys
    color1 = ""
    color2 = ""

    if sys.stdout.isatty():
        color1 = "\033[91m"
        color2 = "\033[0m"

    sys.stderr.write(color1 + message + color2 + "\nSee --help for help on how to use this function\n")
    sys.exit()


def print_error(message, exit=True, warning=False):
    """To print to stderr and use fancy styling, on a single line"""
    import sys
    color = "93" if warning else "91"

    color1 = ""
    color2 = ""

    if sys.stdout.isatty():
        color1 = "\033[" + color + "m"
        color2 = "\033[0m"

    sys.stderr.write("{:s}{:s}{:s}\n".format(color1, message, color2))
    if exit:
        sys.exit()


def test_validity(filename, directory, target, format, chaindata):
    """Test the validity of the input arguments"""
    exceptions = ""

    # Determine the validity of the call
    if filename and not os.path.isfile(filename):
        exceptions += "The filename specified is not an existing file.\n"

    if directory and not os.path.isdir(directory):
        exceptions += "The directory specified is not an existing directory.\n"

    if not os.path.isdir(target):
        exceptions += "The target specified is not an existing directory.\n"

    if chaindata and chaindata[2] > 1:
        exceptions += "The only valid values for invertDensity are 0 and 1\n"

    if filename and directory:
        exceptions += "Only use one input at a time, only --file or --directory\n"

    if not format in extensions:
        exceptions += "The output format is not valid, choose between: {:s} and {:s}\n".format(
            ", ".join([a for a in extensions[:-1]]), extensions[-1])

    if not filename and not directory:
        exceptions += "There needs to be at least one input, so provide --file or --directory\n"

    if exceptions != "":
        print_errors(exceptions)


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    helper = ScriptHelper(
        sys.argv[0].split('/')[-1],  # script name
        description='Converts file(s) to the given output file, converts freely between em, mrc and ccp4 and converts '
                    'from pdb and mmCIF to the other three formats.',
        authors='Douwe Schulte',
        options=[ScriptOption(
                     ['-f', '--file'], 'Filename will convert this single file to the target format, not to be used in '
                                       'combination with --directory', 'string', 'optional'),
                 ScriptOption(
                     ['-d', '--directory'], 'A directory of files, will convert all possible files to the outputFormat '
                                            'specified, not te be used in combination with --file', 'string', 'optional'),
                 ScriptOption(
                     ['-t', '--targetPath'], 'Path to new file, as a relative or absolute path to a directory, the file'
                                             'name will be the same as the original file', 'string', 'optional', '.'),
                 ScriptOption(
                     ['-o', '--outputFormat'], 'The format of the output file, will keep the same name but only change '
                                               'the extension.', 'string', 'required'),
                 ScriptOption(
                     ['-c', '--chainData'], 'Data needed for conversion from the special formats ({:s} and {:s}), in the format \
                     pixelsSize,cubeSize,invertDensity,chain (invertDensity is 1 or 0)'.format(", ".join([a for a in special_extensions[:-1]]), special_extensions[-1]),
                     'float,uint,uint,string', 'optional')])

    filename, directory, target, format, chaindata = parse_script_options(sys.argv[1:], helper)

    # Test for validity of the arguments passed, will stop execution if an error is found
    test_validity(filename, directory, target, format, chaindata)

    if filename:
        num, ex = convertfile(filename, format, target, chaindata)

        # Print possible errors
        if num != 0:
            print_errors(ex)

    elif directory:
        import os

        fileList = os.listdir(directory)
        for filename in fileList:
            # Try for every file if the conversion is possible
            newname = f.name_to_format(filename, target, format)
            if os.path.exists(newname):
                print_error("File: {:s} already exists will be overwritten".format(newname), exit=False, warning=True)

            num = 1
            try:
                num, ex = convertfile(filename, format, target, chaindata)
            except Exception as e:
                ex = "CONVERSION EXCEPTION " + e.message

            # Print the result, ignores status code -1
            if num == 0:
                print("Converted {:s} to {:s}".format(filename, newname))
            elif num == 1:
                print_error("File: {:s} gave error: {:s}".format(filename, ex), exit=False)