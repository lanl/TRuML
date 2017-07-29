import argparse as ap
import logging
import rbexceptions
import readers
import re


# List of arguments to be implemented:
# add flag for outputting Kappa file with/without x{y} rate notation
def main():
    parser = ap.ArgumentParser()
    parser.add_argument('-b', '--bngl_files', nargs='+', metavar='model.bngl', help='BNGL files to convert to Kappa')
    parser.add_argument('-k', '--kappa_files', nargs='+', metavar='model.ka', help='Kappa files to convert to BNGL')
    parser.add_argument('-nf', '--no_print_funcs', action='store_false',
                        help='writes dynamic quantities as variables instead of observables (Kappa)')
    parser.add_argument('-v', '--verbose', action='store_true', help='prints information useful for debugging')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    if args.bngl_files is None and args.kappa_files is None:
        raise rbexceptions.NoModelsException("\nNo models specified.  See `convert.py -h` for information\n")
    elif args.bngl_files is None:
        for kf in args.kappa_files:
            kr = readers.KappaReader(kf)
            model = kr.parse()
            model.write_as_bngl(re.sub('ka', 'bngl', kf))
    else:
        for bf in args.bngl_files:
            br = readers.BNGLReader(bf)
            model = br.parse()
            model.write_as_kappa(re.sub('bngl', 'ka', bf), args.no_print_funcs)


if __name__ == "__main__":
    main()
