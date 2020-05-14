#!/usr/bin/env python
# split cluters in .mmmsa file into different files

import sys
import os

def main():
    if len(sys.argv) != 3:
        mes = '*** Usage: python {} <file.mmmsa> <outdir>\n'
        sys.stderr.write(mes.format(os.path.basename(sys.argv[0])))
        sys.exit(1)
        
    msa_f = sys.argv[1]
    outdir = sys.argv[2].rstrip('/')

    os.makedirs(outdir)

    with open(msa_f) as fp:
        n = 1
        fw = open('{}/cluster-{}.afa'.format(outdir, n), 'w')
        try:
            for line in fp:
                if line.startswith('\0'):
                    fw.close()
                    line=line[1:]
                    if not line:
                        continue
                    n += 1
                    fw = open('{}/cluster-{}.afa'.format(outdir, n), 'w')

                fw.write(line)
        finally:
            fw.close()

    sys.stderr.write('*** clusters generated: {}\n'.format(n))

if __name__ == '__main__':
    main()
