"""truml.truml: defines the entry point for the TRuML translator"""


__version__ = "0.1.0"


import argparse as ap
import logging
import rbexceptions
import readers
import re


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-b', '--bngl_files', nargs='+', metavar='model.bngl', help='BNGL files to convert to Kappa')
    parser.add_argument('-k', '--kappa_files', nargs='+', metavar='model.ka', help='Kappa files to convert to BNGL')
    parser.add_argument('-nf', '--no_print_funcs', action='store_false',
                        help='writes dynamic quantities as variables instead of observables (Kappa)')
    parser.add_argument('-dp', '--dot_and_plus', action='store_true',
                        help='outputs BNGL files with both intra- and intermolecular rules')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('-d', '--debug', action='store_true', help='verbose output useful for debugging')
    parser.add_argument('-l', '--log_file', metavar='convert.log', help='Specify log file')
    args = parser.parse_args()

    log_format = "%(name)s\t%(levelname)s\t%(message)s"

    logging.info("Running the TRuML translator version %s" % __version__)

    # Assign logging level
    ll = logging.WARNING
    if args.debug:
        ll = logging.DEBUG
    elif args.verbose:
        ll = logging.INFO

    if args.log_file:
        logging.basicConfig(format=log_format, filename=args.log_file, level=ll)
    else:
        logging.basicConfig(format=log_format, level=ll)

    if args.bngl_files is None and args.kappa_files is None:
        raise rbexceptions.NoModelsException("\nNo models specified.  See `truml.py -h` for information\n")
    elif args.bngl_files is None:
        for kf in args.kappa_files:
            kr = readers.KappaReader(kf)
            model = kr.parse()
            model.write_as_bngl(re.sub('ka', 'bngl', kf), args.dot_and_plus)
    else:
        for bf in args.bngl_files:
            br = readers.BNGLReader(bf)
            model = br.parse()
            model.write_as_kappa(re.sub('bngl', 'ka', bf), args.no_print_funcs)
