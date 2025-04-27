import gzip
import ntpath
import magic
import collections

from Valkyries import Sequence_Magic


class FASTQ_Reader:
    """
    Main class that creates FASTQ reads using a generator
    """
    __slots__ = ['input_file', 'log', 'name', 'seq', 'index', 'qual', 'read_block', 'file_name', 'fq_file']

    def __init__(self, input_file, log=None):
        """
        Splits the FASTQ read list from the FASTQ Iterator into the lines to be manipulated.  Also does a check to make
        sure the sequence length = quality string length.

        :param input_file:
        :return:
        """

        self.name = None
        self.seq = None
        self.index = None
        self.qual = None
        self.input_file = input_file
        self.log = log
        self.read_block = []
        self.file_name = ntpath.basename(input_file)
        self.fq_file = self.__fastq_file()

    def __fastq_file(self):
        """
        Create a FASTQ file object.
        :return:
        """

        """
        if len(self.input_file) < 3:
            self.log.warning("FASTQ file parameter missing from options file. Correct error and try again.")
            raise SystemExit(1)

        if not pathlib.Path(self.input_file).is_file():
            self.log.warning("FASTQ file {} not found.  Correct error and run again.".format(self.input_file))
            raise SystemExit(1)
        """
        try:
            mime_type = magic.from_file(self.input_file, mime=True).decode()
        except AttributeError:
            mime_type = magic.from_file(self.input_file, mime=True)

        if "text" in mime_type:
            fq_file = open(self.input_file, 'rU')
        elif "gzip" in mime_type:
            fq_file = gzip.open(self.input_file, 'rt', encoding='utf-8')
        else:
            self.log.warning("Unsupported file-type for {}.  Only TEXT or GZIP Allowed.".format(self.input_file))
            raise SystemExit(1)
        return fq_file

    def line_reader(self):
        """
        Part of the generator for FASTQ reads
        """
        for line in self.fq_file:
            while True:
                yield line

    def seq_read(self):
        """
        Generator reads FASTQ file creating read block.
        """
        read_block = []
        count = 0
        eof = False
        try:
            while count < 4:
                read_block.append(next(FASTQ_Reader.line_reader(self)))
                count += 1
        except StopIteration:
            eof = True
            self.name = "EOF"
            yield self

        if len(read_block) == 4 and not eof:

            self.name = read_block[0].strip("\n").strip("@")
            self.seq = read_block[1].strip("\n").strip()
            self.index = read_block[2].strip("\n").strip()
            self.qual = read_block[3].strip("\n").strip()

            if len(self.seq) != len(self.qual):
                raise ValueError("Sequence and quality scores of different lengths! \n{0:s}\n{1:s}\n{2:s}\n{3:s}"
                                 .format(self.name, self.seq, self.index, self.qual))
            yield self



def find_targets(argv, fq_group):
    """
    Does sgRNA search.  Called by multiprocessor.  Given one block of reads, returns count data, sgRNA position,
    anchor data and no anchor data.
    :param argv:
    :param fq_group:
    :return:
    """
    def __anchor_search(fastq_read):
        """
        Looks for the anchor sequence and returns the start position of the sgRNA.
        :param fastq_read:
        :return:
        """
        anchored = False
        start_pos = args.AnchorStart
        sgrna_start = ""
        while not anchored:
            sgrna_start = start_pos + len(args.AnchorSeq)
            mismatch_index = Sequence_Magic.match_maker(args.AnchorSeq, fastq_read[start_pos:][:len(args.AnchorSeq)])
            # self.args.AnchorSeq, fastq_read.seq[start_pos:][:len(self.args.AnchorSeq)])

            # If we do not find the anchor sequence exit the loop and go to the next read.
            if start_pos > args.AnchorStop:
                sgrna_start = args.Expected_Position - args.Target_Padding
                data_dict["NoAnchorCount"][0] += 1
                break

            elif mismatch_index <= args.AnchorMismatch:
                anchored = True
                data_dict["AnchorCount"][0] += 1

            start_pos += 1
        return anchored, sgrna_start

    def __target_match(fastq_read, unknown_seq_start):
        target_mismatch = args.Target_Mismatch
        targets_found_dict = collections.defaultdict(list)

        for line in target_file:
            sgrna_name = line[0]
            target_seq = line[1]
            target_length = len(line[1])
            '''
            if args.RevComp:
                target_seq = Sequence_Magic.rcomp(target_seq)
            '''
            # Go through targets.
            unknown_seq = fastq_read[unknown_seq_start:][:target_length]
            mismatch_index = Sequence_Magic.match_maker(target_seq, unknown_seq)

            if mismatch_index <= 1:
                # found_name = sgrna_name
                return sgrna_name, mismatch_index

            elif mismatch_index <= target_mismatch:
                targets_found_dict[mismatch_index].append((sgrna_name, target_seq))

        if targets_found_dict:
            for i in range(2, target_mismatch + 1):
                if i in targets_found_dict:
                    found_name = targets_found_dict[i][0][0]

                    return found_name, i

        return False, False

    data_dict = collections.defaultdict(list)
    args = argv[0]
    log = argv[1]
    target_file = FileParser.indices(log, args.Target_File)

    data_dict["FoundPosition"] = []
    data_dict["AnchorCount"] = [0]*(args.Target_Mismatch+1)
    data_dict["NoAnchorCount"] = [0]*(args.Target_Mismatch+1)

    for targets in target_file:
        target_name = targets[0]
        data_dict[target_name] = []

    for read in fq_group[0]:
        fq_read = read[1]

        # Find the first position of the sgRNA
        anchor_found, sgrna_seq_start = __anchor_search(fq_read)

        # Compare the sgRNA sequence to the targets.
        if anchor_found:
            target_name, mismatch_found = __target_match(fq_read, sgrna_seq_start)

            # Count our targets.
            if target_name:
                data_dict[target_name].append(mismatch_found)
                data_dict["FoundPosition"].append(sgrna_seq_start)

    return data_dict

