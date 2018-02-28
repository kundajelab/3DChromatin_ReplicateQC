from __future__ import print_function
from setuptools import setup
from setuptools.command.install import install
import subprocess as subp
import os

#from https://stackoverflow.com/questions/42384294/how-to-run-multiple-bash-sh-scripts-in-setup-py-file
class MyInstall(install):
    def run(self):
        repo_dir=os.path.dirname(os.path.realpath(__file__))
        software_dir=repo_dir+'/software'
        #make init files across the path to genomedisco
        open(software_dir+'/__init__.py', 'a').close()
        open(software_dir+'/genomedisco/__init__.py', 'a').close()
        install.run(self)


config = {
'include_package_data': True,
'description': '3DChromatin_ReplicateQC',
'download_url': 'https://github.com/kundajelab/3DChromatin_ReplicateQC',
'version': '0.0.1',
'packages': ['3DChromatin_ReplicateQC'],
'setup_requires': [],
'install_requires': ['numpy>=1.9', 'matplotlib>=1.5.0','h5py','hifive==1.5.6','genomedisco>=0.0.0'],
#'dependency_links': ["https://github.com/kundajelab/deeplift/tarball/v0.5.1-theano#egg=deeplift-0.5.1-theano",
# "https://github.com/kundajelab/simdna/tarball/0.3#egg=simdna-0.3"],
'scripts': [],
'entry_points': {'console_scripts': ['3DChromatin_ReplicateQC = 3DChromatin_ReplicateQC.__main__:main']},
'name': '3DChromatin_ReplicateQC',
'cmdclass': {'install': MyInstall},
#'package_data': {'': extra_files},
}
if __name__== '__main__':
    setup(**config)
