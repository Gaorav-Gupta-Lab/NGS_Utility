
from Valkyries import Sequence_Magic

def index_matching(args, read1, read2=None, index_dict=None, read_count_dict=None):
    """
    This matches an index sequence with the index found in the sequence reads.
    @param read1:
    @param read2:
    @return:
    """

    match_found = False
    index_key = 'unidentified'
    left_match = 5
    right_match = 5

    # Set stringency of index match.
    if args.Platform == "Illumina":
        mismatch = 1
    elif args.Platform == "TruSeq":
        mismatch = 0
    elif args.Platform == "Ramsden":
        mismatch = 3

    for index_key in index_dict:

        left_index = index_dict[index_key][0]
        right_index = index_dict[index_key][1]

        if args.Platform == "Illumina":
            # The indices are after the last ":" in the header.
            right_match = Sequence_Magic.match_maker(right_index, read1.name.split(":")[-1].split("+")[1])
            left_match = Sequence_Magic.match_maker(left_index, read1.name.split(":")[-1].split("+")[0])

        elif args.Platform == "TruSeq":
            # The indices are the first 6 and last 6 nucleotides of the consensus read.
            left_match = Sequence_Magic.match_maker(right_index, read1.seq[:6])
            right_match = Sequence_Magic.match_maker(left_index, read1.seq[-6:])

        elif args.Platform == "Ramsden":
            if args.PEAR:
                left_match = \
                    Sequence_Magic.match_maker(left_index, read1.seq[-len(left_index):])
            else:
                left_match = \
                    Sequence_Magic.match_maker(Sequence_Magic.rcomp(left_index), read2.seq[:len(left_index)])
            right_match = \
                Sequence_Magic.match_maker(right_index, read1.seq[:len(right_index)])

        if index_key not in read_count_dict:
            read_count_dict[index_key] = [0]*9

        if left_match <= mismatch and right_match <= mismatch:
            read_count_dict[index_key][0] += 1

            if args.Platform == "TruSeq":
                pass
                # right_seq = fastq1_read.seq[6:-6]

            match_found = True
            if read2 is None:
                break
        '''
        if match_found and read2:
            # iSeq runs generally have low quality reads on the 3' ends.  This does a blanket trim to remove them.
            left_seq = read2.seq[:-5]
            right_seq = read1.seq[:-5]
            break
        '''

    if not match_found:
        index_key = 'unidentified'

        if 'unidentified' not in read_count_dict:
            read_count_dict['unidentified'] = [0]*9
        read_count_dict['unidentified'][0] += 1

    # return match_found, left_seq, right_seq, index_key, fastq1_read, fastq2_read
    return read_count_dict, index_key
