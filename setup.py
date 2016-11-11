# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os
from os import environ, path
import warnings

from setuptools import setup, find_packages, Extension
from distutils.core import Extension

readme_dir = os.path.dirname(__file__)
readme_filename = os.path.join(readme_dir, 'README.md')

try:
    with open(readme_filename, 'r') as f:
        readme = f.read()
except:
    readme = ""

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except:
    print(
        "Conversion of long_description from MD to reStructuredText failed...")

seq_align_path = environ.get('SEQ_ALIGN_PATH')
extensions = []
if seq_align_path:
    seq_align_dirs = [path.join(seq_align_path, 'src'),
                      path.join(seq_align_path, 'libs'),
                      path.join(seq_align_path, 'libs/string_buffer')]
    pmbec_align = Extension('pmbecalign',
                            include_dirs=seq_align_dirs,
                            library_dirs=seq_align_dirs,
                            libraries=[
                                'align',
                                'strbuf',
                                'pthread',
                                'z'
                            ],
                            sources=['pmbec_align.c'])
    extensions.append(pmbec_align)
else:
    warnings.warn("seq-align is not installed and/or SEQ_ALIGN_PATH is not set, "
                  "so fast Smith-Waterman alignment is currently disabled.")

if __name__ == '__main__':
    setup(
        name='topeology',
        version="0.0.2",
        description="Compare epitope homology",
        author="Tavi Nathanson",
        author_email="tavi {dot} nathanson {at} gmail {dot} com",
        url="https://github.com/hammerlab/topeology",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'topeology = topeology.shell:run'
            ],
        },
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'numpy >=1.10.1, <2.0',
            'pandas >=0.17.0',
            'pepdata >=0.6.7',
            'statsmodels >=0.6.1',
            'scikit-bio >=0.4.0',
            'six >=1.10.0',
            'nose >=1.3.7',
            'pylint >=1.4.4',
        ],
        long_description=readme,
        ext_modules=extensions,
        packages=find_packages(exclude="test")
    )
