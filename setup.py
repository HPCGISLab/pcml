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
    setup(
          name='pcml',
          version='0.1',
          license='Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved. Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.',
          install_requires = reqs,
          packages=['pcml', 'pcml.core', 'pcml.lib', 'pcml.util'],
          test_suite = 'tests',
          )
