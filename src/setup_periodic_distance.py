from setuptools import Extension, setup
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "periodic_distance",
        ["peridodic_distance.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name = 'hello_periodic_distance',
    ext_modules=cythonize("periodic_distance.pyx", annotate = True),

)
