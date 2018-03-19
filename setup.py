from __future__ import print_function
from setuptools import setup
from setuptools.command.install import install
import subprocess as subp
import os


config = {
'include_package_data': True,
'description': '3DChromatin_ReplicateQC',
'download_url': 'https://github.com/kundajelab/3DChromatin_ReplicateQC',
'version': '0.0.1',
'packages': ['3DChromatin_ReplicateQC'],
'setup_requires': [],
'install_requires': ['numpy>=1.9', 'matplotlib>=1.5.0','h5py','hifive==1.5.6','genomedisco>=1.0.0'],
'scripts': [],
'entry_points': {'console_scripts': ['3DChromatin_ReplicateQC = 3DChromatin_ReplicateQC.__main__:main']},
'name': '3DChromatin_ReplicateQC',
}
if __name__== '__main__':
    setup(**config)
