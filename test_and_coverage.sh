#Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
#Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
#Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)

# Check for coverage.
# coverage: Code coverage measurement for Python.
command -v coverage >/dev/null 2>&1 || {
	echo >&2 '`coverage` is required, but not found.'
	echo >&2 "Install with: "
	echo >&2 "    pip install coverage"
	exit 1
}

# Remove old coverage data.
rm -Rfi .coverage

# Measure, and report results.
coverage run --source pcml -m unittest discover -b -v
coverage report -m

