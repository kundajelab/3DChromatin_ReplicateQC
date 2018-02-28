from __future__ import print_function
from setuptools import setup
from setuptools.command.install import install
import subprocess as subp
import os

#from https://stackoverflow.com/questions/42384294/how-to-run-multiple-bash-sh-scripts-in-setup-py-file
class MyInstall(install):
    def run(self):
        repo_dir=os.path.dirname(os.path.realpath(__file__))
        #===========================
        # read the paths from a file
        #===========================
        paths={}
        pathsfile=repo_dir+'/install_scripts/paths.txt'
        for line in open(pathsfile,'r'):
            items=line.strip().split('\t')
            k,v=items[0],items[1]
            if k in ['python','R']:
                paths[k]=v
        #==========================
        #get other methods
        #==========================
        software_dir=repo_dir+'/software'
        #get genomedisco
        '''
        subp.check_output(['bash','-c','git clone https://github.com/kundajelab/genomedisco '+software_dir+'/genomedisco'])
        subp.check_output(['bash','-c','pip install --editable '+software_dir+'/genomedisco/'])
        #get hicrep (and r-essentials for running it)
        
        hicrep_name='http://genome.cshlp.org/content/suppl/2017/10/06/gr.220640.117.DC1/Supplemental_hicrep_1.0.1.tar.gz'
        subp.check_output(['bash','-c','wget '+hicrep_name])
        subp.check_output(['bash','-c',paths['R']+'script '+repo_dir+'/install_scripts/install_R_packages.R'])
        subp.check_output(['bash','-c','rm '+os.path.basename(hicrep_path)])
        
        #Supplemental_hicrep_1.0.1.tar.gz",dependencies="logical")
        #get hic-spector
        subp.check_output(['bash','-c','git clone https://github.com/gersteinlab/HiC-spector '+software_dir+'/HiC-spector'])
        '''
        #make init files across the path to genomedisco
        open(software_dir+'/__init__.py', 'a').close()
        open(software_dir+'/genomedisco/__init__.py', 'a').close()
        install.run(self)


config = {
'include_package_data': True,
'description': '3DChromatin_ReplicateQC',
'download_url': 'https://github.com/kundajelab/3DChromatin_ReplicateQC',
'version': '0.0.1',
'packages': ['3DChromatin_ReplicateQC','3DChromatin_ReplicateQC/software/genomedisco/'],
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
