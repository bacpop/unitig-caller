# Copyright 2019 John Lees

'''Wrapper around mantis to detect presence of sequence elements'''

import os, sys

def get_options():
    import argparse

    description = 'Call unitigs in a population'
    parser = argparse.ArgumentParser(description=description,
                                     prog='unitig-caller')

    parser.add_argument('--mode',
                        choices=['index', 'call'],
                        required=True
                        help='\'index\' sequences, or \'call\' on '
                             'an indexed dataset.')
    parser.add_argument('--strains',
                        help='List of strains to index')
    parser.add_argument('--unitigs',
                        help='List of unitigs to call')



    other = parser.add_argument_group('Other')
    other.add_argument('--cpus',
                        type=int,
                        default=1,
                        help='Number of CPUs to use.')
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    sys.exit(0)

if __name__ == "__main__":
    main()
