# from distutils.core import setup
from Cython.Build import cythonize
from setuptools import setup, Extension
import numpy as np


extensions = [
    Extension(
        name="extend.extender",         # full dotted module path
        sources=["extend/extender.pyx"],# where the .pyx lives
        include_dirs=[np.get_include()],                # add numpy.get_include() later if you cimport numpy
    ),
    Extension(
        name="extend.chainer",         # full dotted module path
        sources=["extend/chainer.pyx"],# where the .pyx lives
        include_dirs=[np.get_include()],                # add numpy.get_include() later if you cimport numpy
    ),
    Extension(
        name="hashing.hash",         # full dotted module path
        sources=["hashing/hash.pyx"],# where the .pyx lives
        include_dirs=[np.get_include()],                # add numpy.get_include() later if you cimport numpy
    ),
    Extension(
        name ="parallelization.batch_reads",
        sources=["parallelization/batch_reads.pyx"],
        include_dirs=[np.get_include()], 
    ),
    Extension(
        name ="mmm_parser.parser",
        sources=["mmm_parser/parser.pyx"],
        include_dirs=[np.get_include()], 
    ),
    Extension(
        name ="mmm_parser.readParser",
        sources=["mmm_parser/readParser.pyx"],
        include_dirs=[np.get_include()], 
    )
]

setup(
    name="mmm",
    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),
    zip_safe=False,
)