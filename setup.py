from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

compilation_args = ['-fPIC', '-std=c++20', '-fopenmp']  
link_args = ['-fopenmp']

ext_modules = [
    Extension(
        "particle_motion",
        ["particle_motion.pyx"],
        extra_compile_args=compilation_args,
        extra_link_args=link_args,
        language='c++',
        extra_objects=["build/liblibrary.a"],
        include_dirs=[numpy.get_include()],
        library_dirs=['build'],
        runtime_library_dirs=['build'],
    )
]

setup(
    ext_modules=cythonize(ext_modules, annotate=True)
)
