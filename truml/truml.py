"""Defines the entry point for the TRuML translator"""


__version__ = "0.0.0"


import argparse as ap
import logging
import rbexceptions
import readers
import re
import pyparsing as pp


def main():
    """Entry point for translation"""

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

    log_format = "%(levelname)s\t%(message)s"

    ll = logging.WARNING
    if args.debug:
        ll = logging.DEBUG
    elif args.verbose:
        ll = logging.INFO

    if args.log_file:
        logging.basicConfig(format=log_format, filename=args.log_file, level=ll, filemode='w')
    else:
        logging.basicConfig(format=log_format, level=ll, filemode='w')

    logging.info("Running the TRuML translator version %s" % __version__)

    if args.bngl_files is None and args.kappa_files is None:
        raise rbexceptions.NoModelsException("\nNo models specified.  See `truml.py -h` for information\n")
    elif args.bngl_files is None:
        for kf in args.kappa_files:
            kr = readers.KappaReader(kf)
            try:
                model = kr.parse()
            except pp.ParseException as pe:
                logging.error("Could not parse line: %s" % pe.line)
                print "Could not parse line: '%s'.  Exiting" % pe.line
                exit(1)
            model.write_as_bngl(re.sub('ka', 'bngl', kf), args.dot_and_plus)
    else:
        for bf in args.bngl_files:
            br = readers.BNGLReader(bf)
            try:
                model = br.parse()
            except pp.ParseException as pe:
                logging.error("Could not parse line: %s" % pe.line)
                print "Could not parse line: '%s'.  Exiting" % pe.line
                exit(1)
            model.write_as_kappa(re.sub('bngl', 'ka', bf), args.no_print_funcs)
