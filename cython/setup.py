from distutils.core import setup
from Cython.Build import cythonize

setup(
#    ext_modules = cythonize("perform_ids.pyx",language_level=3)
#    ext_modules = cythonize("load_spectra.pyx",language_level=3)
#    ext_modules = cythonize("load_kernel.pyx",language_level=3)
    ext_modules = cythonize(["*.pyx"],language_level=3)
)

