#!/usr/bin/python3
__author__ = 'Przemek Decewicz'
__date__ = '2023.05.5'
__version__= '0.1.1'

import binascii
import gzip
import os
import pandas as pd
import re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sys import argv
from typing import Iterable, List, Union


def is_gzip_file(
        f : str,
    ) -> bool:
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(f, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'

def get_span(
        matches : Union [List, Iterable]
    ) -> List:
    """
    Returns the span of the match
    :param matches: list of matches
    :return: list of tuples with the start and end of the match
    """

    return [x.span() for x in matches]

def correct_for_bs_position(
        apos : int, 
        bs_len : int, 
        r_len: int, 
        forward : bool
    ) -> int:
    """
    Corrects the insert position
    :param apos: border index of the adapter
    :param bs_len: maximum length of the binding sequence
    :param r_len: length of the read to avoiid going beyond the index
    :param forward: if the read was forward or not
    :return new_index: new index of insertion
    """
    if forward:
        new_index = apos + bs_len
        if new_index > r_len:
            new_index = r_len
    else:
        new_index = apos - bs_len
        if new_index < 0:
            new_index = 0
            
    return new_index

def identify_insert(
        seq : SeqRecord, 
        i5 : str, 
        i7 : str,
        bs : str,
    ) -> List[int]:

    """
    The function returns the positions of unwanted sequences of i5 and i7 adaptors
    :param seq: str or Seq.Sequnce 
    :param i5: first adapter sequence
    :param i7: second adapter sequence
    :param bs: binding sequence
    :return: list of positions
    """

    i5 = i5.lower()
    i7 = i7.lower()
    bs = bs.lower()
    i5rc = str(Seq(i5).reverse_complement()).lower()
    i7rc = str(Seq(i7).reverse_complement()).lower()
    bsrc = str(Seq(bs).reverse_complement()).lower() # TODO: use it????
    seq = str(seq.seq).lower()

    # FORWARD
    # i5 adapter
    i5r = re.compile(i5)
    i5r_pos_fw = get_span(i5r.finditer(seq))
    # i7 adapter
    i7r = re.compile(i7)
    i7r_pos_fw = get_span(i7r.finditer(seq))

    # REVERSE-COMPLIMENT
    i5r_rc = re.compile(i5rc)
    i5r_pos_rc = get_span(i5r_rc.finditer(seq))
    # reverse-compliment
    i7r_rc = re.compile(i7rc)
    i7r_pos_rc = get_span(i7r_rc.finditer(seq))

    # if both i5 and i7 are present in the forward direction
    if len(i5r_pos_fw) > 0 and len(i7r_pos_fw) > 0:
        # if i5 is before i7
        if i5r_pos_fw[0][1] < i7r_pos_fw[0][0]:
            insert_start = i5r_pos_fw[0][0]
            insert_end = i7r_pos_fw[0][1]
        # if i5 is after i7
        else:
            insert_start = i7r_pos_fw[0][0]
            insert_end = i5r_pos_fw[0][1]
    # if both i5 and i7 are present in the reverse direction
    elif len(i5r_pos_rc) > 0 and len(i7r_pos_rc) > 0:
        # if i5 is before i7
        if i5r_pos_rc[0][1] < i7r_pos_rc[0][0]:
            insert_start = i5r_pos_rc[0][0]
            insert_end = i7r_pos_rc[0][1]
        # if i5 is after i7
        else:
            insert_start = i7r_pos_rc[0][0]
            insert_end = i5r_pos_rc[0][1]
        
    # if only i5 is present in the forward direction
    elif len(i5r_pos_fw) > 0:
        insert_start = i5r_pos_fw[0][0]
        insert_end = i5r_pos_fw[0][1]
        # fix end
        insert_end = correct_for_bs_position(
            apos = insert_end, 
            bs_len = len(bs), 
            r_len = len(seq),
            forward = True
        )
    
    # if only i5 is present in the reverse direction
    elif len(i5r_pos_rc) > 0:
        insert_start = i5r_pos_rc[0][0]
        insert_end = i5r_pos_rc[0][1]
        # fix end
        insert_end = correct_for_bs_position(
            apos = insert_end, 
            bs_len = len(bs), 
            r_len = len(seq),
            forward = True
        )
    
    # if only i7 is present in the forward direction
    elif len(i7r_pos_fw) > 0:
        insert_start = i7r_pos_fw[0][0]
        insert_end = i7r_pos_fw[0][1]
        # fix start
        insert_start = correct_for_bs_position(
            apos = insert_start, 
            bs_len = len(bs), 
            r_len = len(seq),
            forward = False # because of the i7
        )
    
    # if only i7 is present in the reverse direction
    elif len(i7r_pos_rc) > 0:
        insert_start = i7r_pos_rc[0][0]
        insert_end = i7r_pos_rc[0][1]
        # fix start
        insert_start = correct_for_bs_position(
            apos = insert_start, 
            bs_len = len(bs), 
            r_len = len(seq),
            forward = False # because of the i7
        )
    
    # if neither i5 nor i7 are present
    else:
        return []

    return [insert_start, insert_end]


def up_and_down(r : SeqRecord, indexes : List[int]) -> List[SeqRecord]:

    if indexes:
        upstream_r = r[ : indexes[0]]
        downstream_r = r[indexes[1] : ]
    else:
        upstream_r = r[:0]
        downstream_r = r

    return [upstream_r, downstream_r]

def main() -> None:
    args = ArgumentParser(
        description="This scripts splits reads into upstream and downstream parts of R1/R2 reads based on i5 and i7 indexes and the suspected binding sequence.",
    )

    args.add_argument('-r', '--read_1',
                        type = str,
                        help = 'Path to input directory.',
                        required = True)
    
    args.add_argument('-R', '--read_2',
                        type = str,
                        help = 'Path to input directory.',
                        required = True)
    
    args.add_argument('-i', '--index_i5',
                        type = str,
                        help = 'The i5 index sequence.',
                        required = True)
    
    args.add_argument('-I', '--index_i7',
                        type = str,
                        help = 'The i7 index sequence.',
                        required = True)
    
    args.add_argument('-b', '--binding_sequence',
                        type = str,
                        help = 'Binding sequence - we need its length at the moment only. Default: %(default)s',
                        default= 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                        required = False)
    
    args.add_argument('-o', '--outdir',
                        type = str,
                        help = 'Path to output directory.',
                        required = True)

    args.add_argument('-s', '--record_stats',
                        action = 'store_true',
                        help = 'Record stats about the splitting process. These will be written into a tsv-file in the output directory.')
    
    args.add_argument('-q', '--quiet',
                        action = 'store_true',
                        help = 'Do not print anything to stdout.')
    
    args.add_argument('-v', '--version',
                        action = 'version',
                        version = f"%(prog)s v{__version__} ({__date__})")
    

    

    if len(argv[1:]) == 0:
        args.print_usage()
        args.exit()
    try:
        args = args.parse_args()
    except:
        args.exit()

    r1_path = args.read_1
    r2_path = args.read_2
    i5 = args.index_i5
    i7 = args.index_i7
    bs = args.binding_sequence
    outdir = args.outdir
    quiet = args.quiet
    record_stats = args.record_stats

    if not os.path.exists(outdir): os.makedirs(outdir)
    # check if input exists and are not empty
    for infile in [r1_path, r2_path]:
        if not os.path.exists(infile):
            print(f"Input file {infile} does not exist.")
            exit(1)
        if os.stat(infile).st_size == 0:
            print(f"Input file {infile} is empty.")
            exit(1)
    
    # READ FASTQ
    if not quiet: print('Reading fastq files...')
    # read R1
    handle = gzip.open(r1_path, 'rt') if is_gzip_file(r1_path) else open(r1_path, 'r')
    r1_s = SeqIO.parse(handle, 'fastq')
    # read R2
    handle = gzip.open(r2_path, 'rt') if is_gzip_file(r2_path) else open(r2_path, 'r')
    r2_s = SeqIO.parse(handle, 'fastq')
    
    # PREPARE OUTPUT FILES
    if r1_path.endswith('.gz'):
        r1_path_up = os.path.join(outdir, os.path.basename(r1_path).rsplit('.', 2)[0] + '.up.fastq.gz')
        r1_path_down = os.path.join(outdir, os.path.basename(r1_path).rsplit('.', 2)[0] + '.down.fastq.gz')
    else:
        r1_path_up = os.path.join(outdir, os.path.basename(r1_path).rsplit('.', 1)[0] + '.up.fastq.gz')
        r1_path_down = os.path.join(outdir, os.path.basename(r1_path).rsplit('.', 1)[0] + '.down.fastq.gz')

    if r2_path.endswith('.gz'):
        r2_path_up = os.path.join(outdir, os.path.basename(r2_path).rsplit('.', 2)[0] + '.up.fastq.gz')
        r2_path_down = os.path.join(outdir, os.path.basename(r2_path).rsplit('.', 2)[0] + '.down.fastq.gz')
    else:
        r2_path_up = os.path.join(outdir, os.path.basename(r2_path).rsplit('.', 1)[0] + '.up.fastq.gz')
        r2_path_down = os.path.join(outdir, os.path.basename(r2_path).rsplit('.', 1)[0] + '.down.fastq.gz')
    df_path = os.path.join(outdir, os.path.basename(r1_path).rsplit('.', 2)[0] + '.stats.tsv')

    # OPEN OUTPUT FASTQ FILES
    r1_up_handle = gzip.open(r1_path_up, 'wt')
    r1_down_handle = gzip.open(r1_path_down, 'wt')
    r2_up_handle = gzip.open(r2_path_up, 'wt')
    r2_down_handle = gzip.open(r2_path_down, 'wt')

    # FIND INSERTS
    stats = []
    for i, reads in enumerate(zip(r1_s, r2_s), 1):
        r1, r2 = reads
        r1_insert = identify_insert(r1, i5 = i5, i7 = i7, bs = bs)
        r2_insert = identify_insert(r2, i5 = i5, i7 = i7, bs = bs)
        
        r1_ins_len = 0 if not r1_insert else r1_insert[1] - r1_insert[0]
        r2_ins_len = 0 if not r2_insert else r2_insert[1] - r2_insert[0]

        r1_up, r1_down = up_and_down(r1, r1_insert)
        r2_up, r2_down = up_and_down(r2, r2_insert)

        # record some stats
        if record_stats:
            stats.append({  
                'order'         : i,
                'r_id'          : r1.id,
                'r1_len'        : len(r1.seq),
                'r2_len'        : len(r2.seq),
                'r1_ins_len'    : r1_ins_len,
                'r2_ins_len'    : r2_ins_len,
                'r1_up_len'     : 0 if not r1_up else len(r1_up),
                'r1_down_len'   : 0 if not r1_down else len(r1_down),
                'r2_up_len'     : 0 if not r2_up else len(r2_up),
                'r2_down_len'   : 0 if not r2_down else len(r2_down),
            })

        if r1_up: SeqIO.write(r1_up, r1_up_handle, 'fastq')
        if r1_down: SeqIO.write(r1_down, r1_down_handle, 'fastq')
        if r2_up: SeqIO.write(r2_up, r2_up_handle, 'fastq')
        if r2_down: SeqIO.write(r2_down, r2_down_handle, 'fastq')

        if not quiet:
            print(f' -- processed {i} read pairs', end = '\r', flush = True)
    
    if not quiet:
        print()

    # close output files
    r1_up_handle.close()
    r1_down_handle.close()
    r2_up_handle.close()
    r2_down_handle.close()

    if record_stats:
        df = pd.DataFrame(stats)
        df['ins_both'] = [True if (row['r1_ins_len'] > False and row['r2_ins_len'] > 0) else 0 for i, row in df.iterrows()]
        df.to_csv(df_path, sep = '\t', index = False)

    if not quiet:
        print(f'Done!')


if __name__ == '__main__':
    main()