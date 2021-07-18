import sys
import os
from setuptools import setup, find_packages


# read the contents of your README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()


def get_version(relpath):
    '''Read version info from a file without importing it'''
    for line in open(os.path.join(os.path.dirname(__file__), relpath)):
        if '__version__' in line:
            v = line.strip().split('=', 1)[1]
            v = v.strip().strip('"\' ')
            return v
        else:
            sys.stderr.write('*** No version info found..\n')
            return ''

setup(
    name='virsorter',
    version='2.2.3',
    url='https://github.com/jiarong/VirSorter2',
    license='GPL-2',
    author='Jiarong Guo',
    author_email='guojiaro@gmail.com',
    description=('VirSorter2: A multi-classifier, expert-guided approach to '
        'detect diverse DNA and RNA virus genomes'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['virsorter'],
    # use MANIFEST.in instead
    #package_data={
    #    #'': ["virsorter/*", ]  
    #    # include anything under virsorter; 
    #    'virsorter': ['*']
    #},
    # include_package_data flag controls whether non-python files 
    #  distributed in the sdist AND existing in a package 
    #  (directory with __init__.py) install to the site-packages location. 
    #  So in order for non-python files to install with your code 
    #  they need to a) be in the sdist (controlled by MANIFEST.in) and 
    #  b) exist inside an installable python package.
    include_package_data=True,  # include all files in MANIFEST.in
    data_files=[],
    zip_safe=False,
    install_requires= [], # install via conda instead
    entry_points={
          'console_scripts': [
              'virsorter = virsorter.virsorter:cli'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bioinformatics"],
)
