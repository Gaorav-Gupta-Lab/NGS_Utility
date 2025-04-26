
import collections
import os
import csv
import re
from collections import namedtuple
from contextlib import suppress


def options_file(options_parser):
    """
    this function parses the options file and adds the data to the argparse object.
    :param: options_parser
    :return:
    """
    count = 0
    config_file = options_parser.parse_args().options_file

    if not os.path.isfile(config_file):
        print("\033[1;31mWARNING:\n\tOptions_File {} Not Found.  Check File Name and Path.".format(config_file))
        raise SystemExit(1)

    options = csv.reader(open(config_file), delimiter='\t')

    for line in options:
        count += 1
        if len(line) < 2 or "#" in str(line[0]):  # Skip lines that are comments or blank.
            continue

        try:
            value = re.sub(r"[\s]", "", line[1].split("#")[0])  # Strip out any end of line comments and whitespace.

        except IndexError:
            raise SystemExit("There is a syntax error in the options file on line ".format(count))

        key = line[0].strip('--')
        options_parser.add_argument(line[0], dest=key, default=value)
    args = options_parser.parse_args()
    return args


def options_file1(options_file):
    """
    This function parses the file and returns an object.
    :return:
    """
    count = 0
    options_dictionary = collections.defaultdict(str)

    if not os.path.isfile(options_file):
        print("\033[1;31mWARNING:\n\tOptions_File {} Not Found.  Check File Name and Path.".format(options_file))
        raise SystemExit(1)

    options = csv.reader(open(options_file), delimiter='\t')

    for line in options:
        count += 1
        if len(line) < 2 or "#" in str(line[0]):  # Skip lines that are comments or blank.
            continue

        try:
            value = re.sub(r"[\s]", "", line[1].split("#")[0])  # Strip out any end of line comments and whitespace.

        except IndexError:
            raise SystemExit("There is a syntax error in the options file on line " .format(count))

        key = line[0].strip('--')
        options_dictionary[key] = value

    return namedtuple('options_file', options_dictionary.keys())(**options_dictionary)


def sample_indices(input_file):
    """
    Parse the index file or target file and return a list of values.
    :return:
    """
    # log.info("Parsing {}".format(input_file))
    if not os.path.isfile(input_file):
        # log.error("{} Not Found.  Check File Name and Path.".format(input_file))
        raise SystemExit("{} Not Found.  Check File Name and Path.".format(input_file))

    index_list = []
    line_num = 0
    index_file = list(csv.reader(open(input_file), delimiter='\t'))
    for line in index_file:
        line_num += 1
        col_count = len(line)

        if col_count > 0 and len(line[0].split("#")[0]) > 0:  # Skip any lines that are blank or comments.
            tmp_line = []

            for i in range(col_count):
                try:
                    line[i] = line[i].split("#")[0]  # Strip out end of line comments and white space.
                except IndexError:
                    raise SystemExit("There is a syntax error in file {0} on line {1}, column {2} "
                                    .format(input_file, str(line_num), str(i)))

                line[i] = re.sub(",", '', line[i])  # Strip out any commas.

                tmp_line.append(line[i])
            index_list.append(tmp_line)

    # log.debug("Parsing Complete for  {}".format(input_file))

    return index_list

def dual_indices(index_file):
    master_index_dict = {}
    with open(index_file) as f:
        for l in f:
            if "#" in l or not l:
                continue
            l_list = [x for x in l.strip("\n").split("\t")]
            master_index_dict[l_list[0]] = [l_list[1], l_list[2]]
    return master_index_dict

def compress_files(file):
    """
    This function will compress our files.  Takes a very long time.
    :param file:
    """
    # ToDo: Allow user access to gzip compression level.
    if os.path.isfile(file):
        exists_list = [file + ".gz"]
        delete(exists_list)  # if the compressed file already exists we need to delete it first.
        cmd = "gzip -8 " + file
        os.system(cmd)
        print("\t-->{0} Compressed".format(file))
    else:
        print("{0} not found.".format(file))


def delete(file_list):
    """
    Delete one or more files.
    :param file_list:
    :return:
    """
    for file in file_list:
        if os.path.isfile(file):
            with suppress(FileNotFoundError):
                os.remove(str(file))
