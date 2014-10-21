"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from setuptools import setup
from pip.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements("./requirements.txt")

# reqs is a list of requirement
reqs = [str(ir.req) for ir in install_reqs]

if __name__ == '__main__':
    # https://docs.python.org/2/distutils/setupscript.html#additional-meta-data
    setup(
          name='pcml',
          version='0.1',
          description='The parallel cartographic modeling language (PCML) provides spatial operations while hiding the implementation complexities of parallelism.',
          url='https://github.com/HPCGISLab/pcml',
          author='Eric Shook',
          author_email='eshook@kent.edu',
          license='Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved. Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.',
          long_description=
"""
PCML
====

The parallel cartographic modeling language (PCML) is a multi-institutional 
collaborative project aiming to create a computing language for 
cyberGIScientists that is designed for (1) usability, (2) programmability, and 
(3) scalability. PCML provides multi-core parallel processing for spatial 
operations while hiding the implementation complexities of parallelism.

**Author**

Eric Shook <eshook@kent.edu>

**Contributors**

* Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
""",
          #platform=[''],
          install_requires = reqs,
          packages=['pcml', 'pcml.core', 'pcml.lib', 'pcml.util'],
          test_suite = 'tests',
          )
