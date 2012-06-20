from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(name="pyplinkseq",
                         sources=["pyplinkseq.pyx"],
                         libraries=["pyplinkseqint"],
                         language="c++"
                         )]

setup(
  name = 'pyplinkseq',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
