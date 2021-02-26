# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages
from setuptools.extension import Extension
# from distutils.extension import Extension
from Cython.Build import cythonize
# from codecs import open
from os import path
import glob
# import re
import sys
import numpy as np
from scipy._build_utils import numpy_nodepr_api


here = path.abspath(path.dirname("__file__"))

# with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
#     long_description = description.read()

# version = {}
# with open(path.join(here, "Mikado", "version.py")) as fp:
#     exec(fp.read(), version)
# version = version["__version__"]

version = "0.1alpha"

# if version is None:
#     print("No version found, exiting", file=sys.stderr)
#     sys.exit(1)

if sys.version_info.major != 3:
    raise EnvironmentError("""Mikado is a pipeline specifically programmed for python3,
    and is not compatible with Python2. Please upgrade your python before proceeding!""")

extensions = [Extension("pytritex.graph_utils.k_opt_tsp",
                        sources=[path.join("pytritex", "graph_utils", "k_opt_tsp.pyx")],
                        include_dirs=[np.get_include()],
                        language="c++",
                        **numpy_nodepr_api),
              Extension("pytritex.sequencing_coverage.collapse_bins",
                        sources=[path.join("pytritex", "sequencing_coverage", "collapse_bins.pyx")],
                        include_dirs=[np.get_include()],
                        language="c++",
                        **numpy_nodepr_api),
              Extension("pytritex.graph_utils.insert_nodes",
                        sources=[path.join("pytritex", "graph_utils", "insert_nodes.pyx")],
                        include_dirs=[np.get_include()],
                        language="c++",
                        **numpy_nodepr_api),
              Extension("pytritex.graph_utils.node_relocation",
                        sources=[path.join("pytritex", "graph_utils", "node_relocation.pyx")],
                        include_dirs=[np.get_include()],
                        language="c++",
                        **numpy_nodepr_api),
              ]

setup(
    name="pytritex",
    version=version,
    description="Python3-port of Tritex",
    # long_description=long_description,
    url="https://github.com/lucventurini/pytritex",
    author="Luca Venturini",
    author_email="lucventurini@gmail.com",
    license="LGPL3",
    tests_require=["pytest"],
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX :: Linux",
        "Framework :: Pytest",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        'Programming Language :: Python :: 3.7'
    ],
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    # zip_safe=False,
    keywords="wheat genomics",
    packages=find_packages(),
    # scripts=glob.glob("util/*.py"),
    # entry_points={"console_scripts": ["mikado = Mikado:main",
    #                                   "daijin = Mikado.daijin:main",
    #                                   ]},
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    # extras_require={
    #     "postgresql": ["psycopg2"],
    #     "mysql": ["mysqlclient>=1.3.6"],
    #     "bam": ["pysam>=0.8"]
    # },
    # test_suite="nose2.collector.collector",
    package_data={
        "pytritex.share": glob.glob(path.join("pytritex", "share", "*pkl*"))
        },
    include_package_data=True
)
