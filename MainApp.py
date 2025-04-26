import argparse
import re
import csv
import datetime
import collections
import shutil
import Valkyries.FASTQ_Tools
import Valkyries.InputFileParser
from Valkyries import SequenceIndexMatching
from Valkyries import Alignment_Launcher
import os

__author__ = 'Dennis A. Simpson'
__version__ = '0.0.5'
__package__ = 'NGS_Utility'


def main():
    date_format = "%a %b %d %Y %H:%M:%S"
    run_start = datetime.datetime.today().strftime(date_format)

    parser = argparse.ArgumentParser(description="A package to process FASTQ Files.\n {0} v{1}"
                                     .format(__package__, __version__),
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    # options_parser = options_file(parser)
    args = Valkyries.InputFileParser.options_file(parser)

    dual_index_dict = Valkyries.InputFileParser.dual_indices(args.Master_Index_File)
    sample_index = Valkyries.InputFileParser.sample_indices(args.SampleManifest)
    read1_blocks = Valkyries.FASTQ_Tools.FASTQ_Reader(args.FASTQ1)
    read2_blocks = Valkyries.FASTQ_Tools.FASTQ_Reader(args.FASTQ2)
    data_dict = collections.defaultdict(list)

    fq_outfile_dir = "{}{}fq_output_files".format(args.Working_Folder, os.sep)
    if os.path.exists(fq_outfile_dir):
        shutil.rmtree(fq_outfile_dir, ignore_errors=True)
    os.makedirs(fq_outfile_dir, exist_ok=True)

    read_counter = 0
    indexed_counter = 0

    file_end = False
    not1_found = False
    single_found = False
    dual_found = False

    while not file_end:
        for r1_line, r2_line in zip(read1_blocks.seq_read(), read2_blocks.seq_read()):
            if r1_line.name == "EOF":
                file_end = True
                break

            data_dict, index_key = SequenceIndexMatching.index_matching(args, fastq1_name=r1_line.name, fastq2_name=None,
                                                                        index_dict=dual_index_dict,
                                                                        read_count_dict=data_dict)

            r1 = open("{}{}{}_R1.fq".format(fq_outfile_dir, os.sep, index_key), "a")
            r2 = open("{}{}{}_R2.fq".format(fq_outfile_dir, os.sep, index_key), "a")

            r1.write("@{}\n{}\n+\n{}\n".format(r1_line.name, r1_line.seq, r1_line.qual))
            r2.write("@{}\n{}\n+\n{}\n".format(r2_line.name, r2_line.seq, r2_line.qual))
            r1.close()
            r2.close()

            if index_key != "unidentified":
                indexed_counter += 1

            if args.Target_Sequence in r1_line.seq:
                # SNV Count
                data_dict[index_key][1] += 1
                single_found = True

            if args.Target_Sequence2 in r1_line.seq:
                # PAM Count
                data_dict[index_key][2] += 1
                if single_found:
                    data_dict[index_key][5] += 1
                    dual_found = True

            if args.Target_Sequence4 in r1_line.seq:
                # Not1 Count
                data_dict[index_key][4] += 1
                not1_found = True

            read_counter += 1

            if args.Target_Sequence3 in r1_line.seq:
                data_dict[index_key][3] += 1
                if dual_found:
                    data_dict[index_key][6] += 1

            if not1_found and dual_found:
                data_dict[index_key][7] += 1

            if single_found and dual_found and not1_found:
                data_dict[index_key][8] += 1

            single_found = False
            dual_found = False
            not1_found = False

        if read_counter >= 10000:
            file_end = True
            print("Limiting Reads")

        if read_counter % int(args.Read_Limit) == 0:
            print("Processed {} Reads".format(read_counter))
        if file_end:
            break

    outstring_header = \
        (("Package:\t{}\tVersion:\t{}\nAnalysis Run:\t{}\nTotal Read Count:\t{}\nReads With Proper Indices:\t{}\n"
        "SNV Pattern:\tGAACATGTC\nPAM Pattern:\tTTGTGGTTAG\nMarker Pattern:\tTTTATCTAGGT\n"
        "NotI Pattern:\tAGAATGCGGCCGCAATCT\n\n"
        "Index\tRead Count\tSNV\tPAM\tMarker\tNotI\tSNV+PAM\tSNV+PAM+Marker\tSNV+PAM+NotI\tSNV+PAM+Marker+NotI\n")
        .format(__package__, __version__, run_start, read_counter, indexed_counter))

    outstring_data = ""
    for index_key in data_dict:
        read_count = data_dict[index_key][0]
        snv = data_dict[index_key][1]
        pam = data_dict[index_key][2]
        marker = data_dict[index_key][3]
        not1 = data_dict[index_key][4]
        snv_pam = data_dict[index_key][5]
        snv_pam_marker = data_dict[index_key][6]
        snv_pam_not1 = data_dict[index_key][7]
        snv_pam_marker_not1 = data_dict[index_key][8]

        outstring_data += ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                           .format(index_key, read_count, snv, pam, marker, not1, snv_pam, snv_pam_marker,
                                   snv_pam_not1, snv_pam_marker_not1))
    outfile = open("{}{}{}.csv".format(args.Working_Folder, os.sep, args.Job_Name), 'w')
    outfile.write(outstring_header+outstring_data)
    outfile.close()

    print("Begin Compressing FASTQ Files.")
    outfile_list_dict = collections.defaultdict(list)
    for index in sample_index:
        Valkyries.InputFileParser.compress_files("{}{}{}_R1.fq".format(fq_outfile_dir, os.sep, index[0]))
        Valkyries.InputFileParser.compress_files("{}{}{}_R2.fq".format(fq_outfile_dir, os.sep, index[0]))

        outfile_list_dict[index[0]] = ["{}{}{}_R1.fq.gz".format(fq_outfile_dir, os.sep, index[0]),
                                            "{}{}{}_R2.fq.gz".format(fq_outfile_dir, os.sep, index[0])]

    print("All FASTQ Files Compressed.\n")

    if args.AlignDemultiplexedFASTQ == "True":
       print("Begin Aligning FASTQ Files.")
       paired_end = True
       aligner = Alignment_Launcher.AlignmentLauncher(args, paired_end)
       aligner.run_aligner(outfile_list_dict)

    print("All FASTQ Files Aligned.  Exiting Program.  Thank you for using NGS_Utility.  Goodbye.\n\n")
    exit(0)

def ptions_file(options_parser):
        """
        This function parses the options file and adds the data to the argparse object.
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
                raise SystemExit("There is a syntax error in the options file on line " .format(count))

            key = line[0].strip('--')
            options_parser.add_argument(line[0], dest=key, default=value)

        return options_parser


if __name__ == "__main__":
    main()