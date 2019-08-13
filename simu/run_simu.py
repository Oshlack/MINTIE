'''
Module      : run_simu
Description : Create simulated cryptic variants
Copyright   : (c) Marek Cmero, 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''
import os
import numpy as np
import pandas as pd
import pybedtools
import simu
import subprocess
import configparser
import logging
import sys
from argparse import ArgumentParser
from pybedtools import BedTool

PROGRAM_NAME = 'RUN_SIMU'
def exit_with_error(message, exit_status):
    '''
    function from https://github.com/bionitio-team/bionitio-python
    Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)

def init_logging(log_filename):
    '''
    function from https://github.com/bionitio-team/bionitio-python
    If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
    logging.info('command line: %s', ' '.join(sys.argv))

def parse_args():
    '''
    function from https://github.com/bionitio-team/bionitio-python
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Run simulations'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='params',
                        metavar='PARAMS',
                        type=str,
                        help='''Parameters file.''')
    return parser.parse_args()

def preproc(case_fasta, control_fasta, out_prefix):
    '''
    Cleanup and make outdir
    '''
    if os.path.exists(control_fasta):
        os.remove(control_fasta)
    if os.path.exists(case_fasta):
        os.remove(case_fasta)
    outdir = '/'.join(out_prefix.split('/')[:-1])
    subprocess.call(['mkdir', '-p', outdir])

def simulate_fusions(simp, params, available_genes, all_exons, gene_trees, valid_txs, paths):
    fuscols = ['loc1', 'gene1', 'tx1', 'insert', 'loc2', 'gene2', 'tx2', 'fusion_type']
    n_vars = int(simp['n_cfus']) \
             + int(simp['n_ee_fus']) \
             + int(simp['n_ne_fus']) \
             + int(simp['n_ins_fus']) \
             + int(simp['n_ufus'])

    fus_genes, available_genes = simu.pick_genes(n_vars, available_genes)
    partner_genes, available_genes = simu.pick_genes(n_vars - int(simp['n_ufus']), available_genes)
    partner_genes = list(partner_genes) + [''] * int(simp['n_ufus'])
    add = [''] * int(simp['n_cfus']) \
        + ['EE'] * int(simp['n_ee_fus'])\
        + ['NE'] * int(simp['n_ne_fus']) \
        + ['INS'] * int(simp['n_ins_fus']) \
        + [''] * int(simp['n_ufus'])

    fusions = pd.DataFrame({'gene1': fus_genes,
                            'gene2': partner_genes,
                            'add': add})
    fus_info = []
    for idx, row in fusions.iterrows():
        vartype = 'canonical' if row['add'] == '' else row['add']
        vartype = 'unpartnered' if row['gene2'] == '' else vartype

        gene1, gene2 = row['gene1'], row['gene2'] if row['gene2'] != '' else None
        fusname = '%s:%s' % (gene1, gene2) if gene2 else gene1
        logging.info('Generating %s %s fusion' % (fusname, vartype))

        gene_ref = all_exons.filter(lambda x: simu.get_gene_name(x) in [gene1, gene2]).saveas()
        tx1 = simu.get_transcripts(gene1, gene_ref, valid_txs=valid_txs)[0]
        tx2 = simu.get_transcripts(gene2, gene_ref, valid_txs=valid_txs)[0] if gene2 else None

        add = row['add'] if row['add'] != '' else None
        fus_parts = simu.write_fusion((tx1, tx2), (gene1, gene2), gene_ref, paths['genome_fasta'], \
                                      params, gene_trees, all_exons=all_exons, add=add)
        fus_info.append(fus_parts)

    fus_info = pd.DataFrame(fus_info, columns=fuscols)
    fus_info.to_csv('%s_fusions_simulated.tsv' % paths['out_prefix'], index=False, sep='\t')
    return available_genes

def simulate_tsvs_and_splicevars(simp, params, available_genes, all_exons,
                                 gene_trees, junc_ref, valid_txs, paths):
    vartypes = ['DEL'] * int(simp['n_del']) \
             + ['INS'] * int(simp['n_ins']) \
             + ['ITD'] * int(simp['n_itd']) \
             + ['PTD'] * int(simp['n_ptd']) \
             + ['INV'] * int(simp['n_inv']) \
             + ['EE'] * int(simp['n_ee']) \
             + ['NE'] * int(simp['n_ne']) \
             + ['RI'] * int(simp['n_ri']) \
             + ['NEJ'] * int(simp['n_nej']) \
             + ['US'] * int(simp['n_us'])
    var_genes, available_genes = simu.pick_genes(len(vartypes), available_genes)
    varcols = ['loc', 'tx', 'gene', 'size', 'exon', 'vartype']
    indel_range, block_range, exons_range = params['ins_range'], params['block_range'], params['exons_range']
    genome_fasta, out_prefix = paths['genome_fasta'], paths['out_prefix']

    var_df = pd.DataFrame({'gene': var_genes, 'vartype': vartypes})
    var_info = []

    for idx, row in var_df.iterrows():
        vartype, gene = row['vartype'], row['gene']
        logging.info('Generating %s in gene %s' % (vartype, gene))

        gene_ref = all_exons.filter(lambda x: simu.get_gene_name(x) == gene).saveas()
        tx = simu.get_transcripts(gene, gene_ref, valid_txs=valid_txs)[0]

        if vartype in ['DEL', 'INS', 'ITD']:
            varsize, exon = simu.write_indel(tx, gene_ref, genome_fasta,
                                             indel_range, out_prefix, vartype=vartype)
        elif vartype in ['PTD', 'INV']:
            varsize, exon = simu.write_large_tsv(tx, gene_ref, genome_fasta,
                                                 out_prefix, exons_range, vartype=vartype)
        elif vartype == 'NEJ':
            varsize, exon = simu.write_trunc_exons(tx, all_exons, genome_fasta, out_prefix, block_range)
            while varsize == '':
                logging.info('Gene contains no valid transcripts, reselecting...')

                genes, available_genes = simu.pick_genes(1, available_genes)
                gene = genes[0]
                logging.info('Generating %s in gene %s' % (vartype, gene))

                tx = simu.get_transcripts(gene, all_exons, valid_txs=valid_txs)[0]
                varsize, exon = simu.write_trunc_exons(tx, all_exons, genome_fasta, out_prefix, block_range)
        elif vartype == 'US':
            tx, loc, stats = simu.write_unannot_splice(gene, gene_ref, valid_txs, genome_fasta,
                                                      out_prefix, junc_ref, gene_trees)
            while tx == '':
                logging.info('Gene contains no valid transcripts, reselecting...')

                genes, available_genes = simu.pick_genes(1, available_genes)
                gene = genes[0]
                logging.info('Generating %s in gene %s' % (vartype, gene))
                tx, loc, stats = simu.write_unannot_splice(gene, all_exons, valid_txs, genome_fasta,
                                                           out_prefix, junc_ref, gene_trees)
            varsize, exon = stats['varsize'], stats['exon']
        else:
            # NE, EE or RI
            tx, loc, stats = simu.write_novel_exon(gene, valid_txs, all_exons,
                                           genome_fasta, out_prefix,
                                           block_range, gene_trees, vartype=vartype)
            while tx == '':
                logging.info('Gene contains no valid transcripts, reselecting...')

                genes, available_genes = simu.pick_genes(1, available_genes)
                gene = genes[0]
                logging.info('Generating %s in gene %s' % (vartype, gene))
                tx, loc, stats = simu.write_novel_exon(gene, valid_txs, all_exons,
                                                       genome_fasta, out_prefix,
                                                       block_range, gene_trees, vartype=vartype)
            varsize, exon = stats['varsize'], stats['exon']

        chrom = simu.get_tx_chrom(tx, all_exons)
        loc = simu.get_gene_loc(chrom, gene_trees, gene)
        var_info.append([loc, tx, gene, varsize, exon, vartype])

    var_info = pd.DataFrame(var_info, columns=varcols)
    var_info.to_csv('%s_tsvs_splice_simulated.tsv' % paths['out_prefix'], index=False, sep='\t')
    return available_genes


def write_background_genes(bg_genes, simp, paths, valid_txs, all_exons):
    control_fasta = '%s-control.fasta' % paths['out_prefix']
    case_fasta = '%s-case.fasta' % paths['out_prefix']

    for gene in bg_genes:
        logging.info('Writing out %s' % gene)
        tx = simu.get_transcripts(gene, all_exons, valid_txs=valid_txs)[0]
        exons = all_exons.filter(lambda x: x['transcript_id'] == tx).saveas()
        tx_seq, strand = simu.get_seq(exons, paths['genome_fasta'])

        simu.write_wildtype_sequence(tx_seq, strand, control_fasta, tx)
        simu.write_wildtype_sequence(tx_seq, strand, case_fasta, tx)
        simu.write_wildtype_sequence(tx_seq, strand, case_fasta, tx) # write case twice (two alleles)

def simulate_reads(fasta, readp, paths, control=False):
    fold = str(int(readp['fold']) * 2) if control else readp['fold']
    sample = 'control' if control else 'case'
    logging.info('Generating simulated reads for %s' % sample)
    try:
        seed = int(readp['art_seed_%s'] % sample)
        assert seed > 0
    except KeyError:
        seed = np.random.randint(99999)
    logging.info('ART-Illumina random seed is %d' % seed)
    cmd = [paths['art_illumina'],
            '-ss', 'HS25',
            '-i', fasta,
            '-p',
            '-l', readp['read_len'],
            '-f', fold,
            '-m', readp['frag_size'],
            '-s', readp['frag_sd'],
            '-rs', str(seed),
            '-o', '%s-%s_R' % (paths['out_prefix'], sample)]
    logging.info('Running command: ' + ' '.join(cmd))
    subprocess.call(cmd)

def gzip_files(sample, out_prefix):
    for r in range(2):
        outf = open('%s-%s_R%d.fastq.gz' % (out_prefix, sample, (r+1)), 'w')
        cmd = ['gzip', '-c', '%s-%s_R%d.fq' % (out_prefix, sample, (r+1))]
        logging.info('Running command: ' + ' '.join(cmd))
        subprocess.call(cmd, stdout=outf)
        outf.close()

def simulate(args):
    config = configparser.ConfigParser()
    config.read(args.params)

    simp = config['SimParams']
    readp = config['ReadParams']
    paths = config['Paths']

    control_fasta = '%s-control.fasta' % paths['out_prefix']
    case_fasta = '%s-case.fasta' % paths['out_prefix']
    preproc(case_fasta, control_fasta, paths['out_prefix'])

    logging.info('Building references...')
    gr = BedTool(paths['gtf_ref'])
    junc_ref = simu.build_junc_ref(paths['junc_ref'])

    # make gene start/end reference
    gene_trees = simu.get_gene_features(gr)

    # get exons for records that have a gene name
    all_exons = gr.filter(lambda x: x[2] == 'exon').saveas()
    all_genes = np.unique([simu.get_gene_name(ex) for ex in all_exons if simu.get_gene_name(ex)!=''])

    logging.info('Extracting valid transcripts')
    valid_txs, valid_genes = simu.get_valid_txs(all_exons, int(simp['min_exons']))
    valid_txs = np.unique([tx for tx, gn in valid_txs])
    available_genes = [gene for gene in all_genes if gene in valid_genes]

    # set up parameters
    block_range = (int(simp['block_min']),
                   int(simp['block_max']))
    indel_range = (int(simp['indel_min']),
                   int(simp['indel_max']))
    exons_range = (int(simp['exon_min']),
                   int(simp['exon_max']))
    params = {'n_exons': int(simp['n_exons']),
              'ins_range': indel_range,
              'block_range': block_range,
              'exons_range': exons_range,
              'out_prefix': paths['out_prefix']}

    try:
        seed = int(simp['seed_init'])
        if seed >= 0:
            np.random.seed(seed)
            logging.info('Set seed to %d' % seed)
    except KeyError:
        pass

    logging.info('Simulating fusions')
    available_genes = simulate_fusions(simp, params, available_genes, all_exons,
                                       gene_trees, valid_txs, paths)

    logging.info('Simulating TSVs and splice variants')
    available_genes = simulate_tsvs_and_splicevars(simp, params, available_genes, all_exons,
                                                   gene_trees, junc_ref, valid_txs, paths)

    logging.info('Writing background genes')
    bg_genes, available_genes = simu.pick_genes(int(simp['n_background_genes']), available_genes)
    write_background_genes(bg_genes, simp, paths, valid_txs, all_exons)

    simulate_reads(case_fasta, readp, paths)
    simulate_reads(control_fasta, readp, paths, control=True)

    logging.info('Gzipping files')
    gzip_files('case', paths['out_prefix'])
    gzip_files('control', paths['out_prefix'])

def main():
    args = parse_args()
    init_logging(args.log)
    simulate(args)

if __name__ == '__main__':
    main()
