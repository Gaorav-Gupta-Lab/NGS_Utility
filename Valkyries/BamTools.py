
import pysam
import os

def bamfile_sort(args, unsorted_bamfile_list):
    """
    Sort and index the BAM file using pysam samtools.
    :param args:
    :param unsorted_bamfile_list:
    :return:
    """
    sorted_bamfile_dict = {}
    for unsorted_bamfile in unsorted_bamfile_list:
        name, extension = os.path.splitext(unsorted_bamfile)
        sorted_bamfile = unsorted_bamfile.replace(extension, "_sorted.bam")
        sorted_bamfile_dict[name] = sorted_bamfile

        #  "-l", compression_level,

        pysam.samtools.sort("-@", str(args.Spawn), '-o', sorted_bamfile, unsorted_bamfile, catch_stdout=False)
        pysam.samtools.index(sorted_bamfile, catch_stdout=False)

    return sorted_bamfile_dict


def total_align_count(input_bam, chromosome=None):
    """Returns count of all mapped alignments in input BAM (based on index)
    :param input_bam:
    :param chromosome:
    :return:
    """
    count = 0

    for line in pysam.samtools.idxstats(input_bam).split('\n'):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')

            if chrom != '*' and not chromosome:
                count += int(mapped) + int(unmapped)
            elif chrom == chromosome:
                count += int(mapped) + int(unmapped)
    return count


def coverage(logger, region):
    """
    Determines sequence coverage for sub-regions.  Should be cleaned up some
    :param logger:
    :param region:
    :return:
    """
    logger.info("\t-->Determining read coverage and depth for \033[1;35m{0}\033[m.".format(region[0]))

    # data_file = self.data_file
    pysam_depth = pysam.depth("-r{0}:1-{1}".format(region[0], region[1][0]), data_file, split_lines=True)
    depth_list = []
    depth_counts = 0
    breadth_counts = 0
    for line in pysam_depth:
        depth_list.append(line.split("\t")[2])
        depth_counts += int(line.split("\t")[2])
        breadth_counts += 1

    logger.info("\t-->Read coverage and depth analysis complete for \033[1;35m{0}\033[m.".format(region[0]))

    return depth_counts / int(region[1][0]), breadth_counts / int(region[1][0])


