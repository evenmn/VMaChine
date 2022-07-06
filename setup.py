import os
from glob import glob
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext


# set version number both for pip and python interface
__version__ = "0.0.1"


this_dir = os.path.dirname(os.path.abspath(__file__))

with open("README.md", "r") as fh:
    long_description = fh.read()


#source_files = glob('src/*.cpp')
#source_files.extend(glob('src/**/*.cpp'))
#source_files.remove('src/main.cpp')
source_files = glob('src/RNG/*.cpp')
source_files.append('python/wrap.cpp')
print(source_files)


ext_modules = [
    Extension(
    'vmc',
        source_files,
        include_dirs=[os.path.join(this_dir, 'include'),
                      os.path.join(this_dir, 'src'),
                ],
    language='c++',
    extra_compile_args = ['-std=c++14'],
    define_macros = [('VERSION_INFO', __version__)],
    ),
]


setup(name='vmc',
      version=__version__,
      description='Python library for variational Monte Carlo simulations with support for neural network wave functions',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/evenmn/VMaChine',
      author='Even Marius Nordhagen',
      author_email='evenmn@mn.uio.no',
      license='GPLv3',
      include_package_data=True,
      zip_safe=False,
      packages=find_packages(),
      cmdclass={'build_ext': build_ext},
      classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: C++",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU General License v3 or later (GPLv3+)",
      ],
      python_requires='>=3.5',
      ext_modules=ext_modules,
)
