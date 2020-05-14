#!/usr/bin/env python
# No need to run this when running for the first time; NCBITaxa() will create
#   a new database automatically

from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
