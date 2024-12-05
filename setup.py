import os
from setuptools import Extension, setup
import sys
import numpy as np



with open("README.md", 'r') as readme:
    if sys.version_info[0] < 3:
        long_description = readme.read().decode('utf-8')
    else:
        long_description = readme.read()

current_dir = os.getcwd()
arrow_dir = os.path.join(current_dir, 'stochastic_arrow')

# Compile the Cython code to C for development builds:
#    USE_CYTHON=1 python -m pip install -e .
# and for building source distribution packages:
#    USE_CYTHON=1 python -m build --sdist
# and *not* when installing a distribution package.
# See http://docs.cython.org/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
USE_CYTHON = 'USE_CYTHON' in os.environ

ext = '.pyx' if USE_CYTHON else '.c'

cython_extensions = [
    Extension('stochastic_arrow.arrowhead',
              sources=[
                  'stochastic_arrow/mersenne.c',
                  'stochastic_arrow/obsidian.c',
                  'stochastic_arrow/arrowhead'+ext,],
              include_dirs=[arrow_dir, np.get_include()],
              define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
              )]

if USE_CYTHON:
    from Cython.Build import cythonize
    cython_extensions = cythonize(
        cython_extensions,
        include_path=[arrow_dir],
        annotate=True,  # to get an HTML code listing
    )

setup(
    name='stochastic-arrow',
    version='1.0.0',
    packages=['stochastic_arrow'],
    author='Ryan Spangler, John Mason, Jerry Morrison, Chris Skalnik, Travis Ahn-Horst, Sean Cheah',
    author_email='ryan.spangler@gmail.com',
    url='https://github.com/CovertLab/arrow',
    license='MIT',
    ext_modules=cython_extensions,
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['numpy>=1.14', 'six'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ])
