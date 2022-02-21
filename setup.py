"""Install the pynmdc Package"""
import os

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(BASE_DIR, 'README.md'), encoding='utf-8').read()
VERSION = open(os.path.join(BASE_DIR, 'VERSION'), encoding='utf-8').read()

INSTALL_REQUIRES = [
    pkg for pkg in open('requirements.txt').readlines()
]

PYTHON_REQUIRES = '>=3.6.*'


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)


setup(
    name='variant_viz',
    version=VERSION,
    description="SARS-CoV-2 visualization tool",
    long_description=README,
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: POSIX :: BSD',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords='EDGE bioinformatics, covid19, sars-cov-2, bioinformatics',
    author='Po-E Li',
    author_email='po-e@lanl.gov',
    url='https://github.com/poeli/variant_viz',
    license='BSD 3 "Clause"',
    # packages=find_packages('src'),
    packages=['variant_viz', 'variant_viz.scripts', 'variant_viz.surv_viz', 'variant_viz.surv_viz.data'],
    package_dir={'': 'src'},
    package_data={'variant_viz.surv_viz.data': ['*.tsv']},
    include_package_data=False,
    zip_safe=False,
    install_requires=INSTALL_REQUIRES,
    python_requires=PYTHON_REQUIRES,
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand
        },
    entry_points={
        'console_scripts': ['ec19_varviz = variant_viz.scripts.__main__:vizcli']
        }
)
