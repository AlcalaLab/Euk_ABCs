import argparse, sys

def capture_args():
    # create the top-level parser
    gen_descript = 'Something goes here...\n\nUsage:\n     UPDATE! phyg.py <module> <options>\n\n' \
        'To see module options:\n     phyg.py <module> -h'

    parser = argparse.ArgumentParser(description = gen_descript,
            usage = argparse.SUPPRESS, add_help = False,
            formatter_class = argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)
    parser.add_argument('--version', action = 'store_true', help = argparse.SUPPRESS)

    contam = parser.add_argument_group('General Contamination Options', description = (
    '''--in (-i)             input FASTA file\n'''
    '''--out-dir (-o)        output directory name\n'''
    '''--taxon-name (-n)     taxon name (genus species)\n'''
    '''--taxon-code (-c)     additional taxonomic code included\n'''
    '''--db (-d)             diagnostic gene-family database\n'''
    '''--kmer                k-mer based clustering\n'''
    '''--phylo               phylogenetic based evaluation (default)\n'''
    '''--threads (-p)        number of CPU threads to use (default = 4)\n'''
    '''--clean               remove intermediate files\n'''
    '''--quiet (-q)          no console output\n'''
    '''--gzip (-gz)          tar and gzip output\n'''
    '''--help (-h)           show this help message and exit\n'''))

    contam_phylo = parser.add_argument_group('Phylogenetic Tree Based Options', description = (
    '''--ref-msa (-m)            diagnostic mutli-sequence alignments\n'''
    '''--ref-trees (-t)         reference diagnostic phylogenies\n'''))

    contam_kmer = parser.add_argument_group('K-mer Based Options', description = (
    '''--cluster-threshold   distance threshold for calling clusters (default = 3)\n\n'''
    '''--eval-threshold      representative proportion of a clade [e.g. Rhizaria'] to
                      assign a clade-name to a cluster (default = 0.6)\n\n'''))

    contam.add_argument('--help', '-h', action = "help", help = argparse.SUPPRESS)

    contam.add_argument('--input', '--in', '-i', action = 'store', metavar = '[FASTA-file or FASTA-DIR]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--out-dir','-o', action = 'store', metavar = '[output-directory]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--taxon-name', '-n', action = 'store', metavar = '[taxon-name]',
        nargs='+', type = str, help = argparse.SUPPRESS)

    contam.add_argument('--taxon-code', action = 'store', metavar = '[taxon-code]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--db', '-d', action = 'store', metavar = '[Diagnostic Database]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--kmer', '-k', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--cluster-threshold', action = 'store', default = 3,
        metavar = '[cluster-thresh]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--eval-threshold', action = 'store', default = 0.6,
        metavar = '[eval-thresh]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--phylo', action = 'store_false', help = argparse.SUPPRESS)

    contam.add_argument('--ref-msa', '-m', action = 'store', metavar = '[diagnostic msas]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--ref-trees', '-t', action = 'store', metavar = '[diagnostic trees]',
        type = str, help = argparse.SUPPRESS)

    contam.add_argument('--blen-mode', '-bl', action = 'store', default = 'median',
        metavar = '[branch-length-eval]', type = str, help = argparse.SUPPRESS)

    contam.add_argument('--threads','-p', action = 'store', default = 4,
        metavar = '[Threads]', type = int, help = argparse.SUPPRESS)

    contam.add_argument('--clean', action = 'store_true', help = argparse.SUPPRESS)

    contam.add_argument('--quiet', '-q', action = 'store_false', help = argparse.SUPPRESS)

    contam.add_argument('--gzip', '-gz', action = 'store_true', help = argparse.SUPPRESS)

    if len(sys.argv) == 1:
#     print(ascii_logo_vsn())
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    if args.kmer:
        args.approach = 'kmer'
        args.phylo = False
    else:
        args.approach = 'phylo'

    return args


def double_check_args(args):
    req_args = []

    if not args.input:
        req_args.append('    --in <input-[FASTA-file or FASTA-DIR]>')

    if not args.db:
        req_args.append('    --db <gene-family-database>')

    if args.phylo:
        if not args.ref_msa:
            req_args.append('    --ref-msa <reference-msa-directory>')
        if not args.ref_trees:
            req_args.append('    --ref-trees <reference-tree-directory>')
        if args.blen_mode.lower() not in ['median','average']:
            req_args.append('    --blen-mode <"median" or "average">')

    if req_args:
        print('\nERROR: Missing the following required arguments:')
        print("\n".join(req_args) + '\n')
        return False
    else:
        return True
