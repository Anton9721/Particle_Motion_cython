from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

compilation_args = ['/openmp', '/std:c++20']
link_args = ['/openmp']

ext_modules = [
    Extension(
        "particle_motion",
        ["particle_motion.pyx"],
        extra_compile_args=compilation_args,
        extra_link_args=link_args,
        language='c++',
        extra_objects=["build/Debug/library.lib"],
        include_dirs=[numpy.get_include()],
    )
]

setup(
    ext_modules=cythonize(ext_modules, annotate=True)
)
