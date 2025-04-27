"""
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27517
@copyright: 2025
"""
import os
import subprocess

__author__ = 'Dennis A. Simpson'
__version__ = "0.5.0"


class AlignmentLauncher:
    """
    Builds the commands and launches the aligner.
    """

    def __init__(self, args, logger, paired_end):

        self.args = args
        self.logger = logger
        self.paired_end = paired_end

    def run_bwa_aligner(self, fq1_name, fq2_name, sam_file, bam_file):
        threads = int(self.args.Spawn)
        aligner_options = getattr(self.args, "Aligner_Options", "")

        if self.args.BWA_Method == "mem":
            self.logger.info("\t-->Begin alignment of {0} and {1} with BWA mem.".format(fq1_name, fq2_name))
            cmd = "bwa mem -t {0} {1} {2} {3} {4} > {5}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq1_name, fq2_name, sam_file)
            self.logger.debug(cmd)
            subprocess.run([cmd], shell=True)

            self.logger.info("\t-->Alignment complete, begin SAM to BAM conversion.")
            cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
            subprocess.run([cmd], shell=True)
            self.logger.debug(cmd)

        elif self.args.BWA_Method == "aln":
            print("Begin alignment of {0} and {1} with BWA aln.".format(fq1_name, fq2_name))
            sai_file1 = fq1_name.replace(".fastq", ".sai")
            sai_file2 = fq2_name.replace(".fastq", ".sai")

            cmd1 = "bwa aln {0} {1} {2} {3} > {4}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq1_name, sai_file1)
            cmd2 = "bwa aln {0} {1} {2} {3} > {4}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq2_name, sai_file2)

            subprocess.run([cmd1], shell=True)
            subprocess.run([cmd2], shell=True)
            print("Alignment complete, begin SAI to SAM conversion.")

            cmd3 = "bwa sampe {0} {1} {2} {3} {4} > {5}" \
                    .format(self.args.Aligner_RefSeq, sai_file1, sai_file2, fq1_name, fq2_name, sam_file)

            if not self.paired_end:
                cmd3 = "bwa {0} samse {1} {2} > {3}".format(self.args.Aligner_RefSeq, sai_file1, fq1_name,
                                                            sam_file)

            subprocess.run([cmd3], shell=True)

            print("SAI to SAM conversion complete. Begin SAM to BAM conversion.")
            cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
            subprocess.run([cmd], shell=True)

        else:
            self.logger.error("No valid BWA alignment method provided")
            raise SystemExit(1)

    def run_bowtie(self, fq1_name, fq2_name, sam_file, bam_file):
        if self.paired_end:
            sub_cmd = " -1 {0} -2 {1}".format(fq1_name, fq2_name)
            message = "file {0} and {1}".format(fq1_name, fq2_name)
        else:
            sub_cmd = " -U {0}".format(fq1_name)
            message = "file {0}".format(fq1_name)

        if self.args.local:
            bowtie2_local = " --local --ma {0}".format(self.args.ma)
        else:
            bowtie2_local = ''

            self.logger.info("Begin alignment of {0} with Bowtie2".format(message))

        # Construct command for aligner and call subprocess to run.
        cmd = "bowtie2 {0} --trim5 {1} --trim3 {2} -x {3} {4} -S {5}" \
            .format(bowtie2_local, self.args.trim5, self.args.trim3, self.args.Aligner_RefSeq, sub_cmd,
                    sam_file)
        subprocess.run([cmd], shell=True)

        self.logger.info("Alignment complete. Converting SAM format to BAM format")

        # Construct command for samtools and call subprocess to run.
        cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
        subprocess.run([cmd], shell=True)

    def run_aligner(self, outfile_list_dict):
        """
        This will run our aligner of choice.  The final output is a BAM file.
        :param outfile_list_dict:
        :return:
        """
        bamfile_list = []
        samfile_list = []
        for sample_index in outfile_list_dict:
            bam_file = "{}{}{}.bam".format(self.args.WorkingFolder, os.sep, sample_index)
            sam_file = "{}{}{}.sam".format(self.args.WorkingFolder, os.sep, sample_index)
            fq1_name = outfile_list_dict[sample_index][0]
            fq2_name = None
            if self.paired_end:
                fq2_name = outfile_list_dict[sample_index][1]

            if self.args.Aligner == "BWA":
                self.run_bwa_aligner(fq1_name, fq2_name, sam_file, bam_file)

            elif self.args.Aligner == "Bowtie2":
                self.run_bowtie(fq1_name, fq2_name, sam_file, bam_file)

            else:
                self.logger.error("\033[1;31m***Warning:\033[m {0} not allowed.  Currently only BWA or Bowtie2 supported."
                               .format(self.args.Aligner))
                raise SystemExit(1)

            bamfile_list.append(bam_file)
            samfile_list.append(sam_file)

        return bamfile_list, samfile_list
