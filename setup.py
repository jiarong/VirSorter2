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
            sys.write('*** No version info found..\n')
            return ''

setup(
    name='virsorter',
    version=get_version("virsorter/__init__.py"),
    url='https://github.com/jiarong/VirSorter2',
    license='GPL-2',
    author='Guo et al.',
    author_email='guojiaro@gmail.com',
    description=('VirSorter2: A multi-classifier, expert-guided approach to '
        'detect diverse DNA and RNA virus genomes'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['virsorter'],
    package_data={
        '': ["virsorter/*", ]  
        # include anything under virsorter; 
        # same as virsorter: ['*']
    },
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires= [
    ],
    # install via conda: click, ruamel.yaml, snakemake
    entry_points={
          'console_scripts': [
              'virsorter = virsorter.virsorter:cli'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bioinformatics"],
)
