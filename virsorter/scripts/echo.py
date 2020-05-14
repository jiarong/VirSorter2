#!/usr/bin/env python

import sys
import os
import logging
import click


@click.command()
@click.option('--msg', default=None, help='message to echo')
@click.option('--outfile', default=None, help='save message to file')
@click.option('--level', default='info', 
        type=click.Choice(['info', 'warning', 'error'], case_sensitive=False),
        help='levelname in logging module')
def main(msg, outfile, level):
    '''Print message in similar way to echo in bash with python logging style
    '''
    #numeric_level = getattr(logging, level.upper())
    logging.basicConfig(
        filename=outfile,
        level=logging.INFO, 
        datefmt='%Y-%m-%d %H:%M',
        format='[%(asctime)s %(levelname)s] %(message)s',
    )
    if msg == None:
        infile = '/dev/stdin'
        msg = open(infile).read().rstrip('\n')
    if level == 'info':
        logging.info(msg)
    elif level == 'warning':
        logging.warning(msg)
    elif level == 'error':
        logging.error(msg)

if __name__ == '__main__':
    main()
