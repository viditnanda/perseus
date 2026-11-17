import os
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

HERE = os.path.dirname(os.path.abspath(__file__))
# headers: ../include in the repository, ./include in a source distribution
INC = os.path.join(HERE, "include") if os.path.isdir(os.path.join(HERE, "include", "perseus")) else os.path.join(HERE, "..", "include")
ext = Pybind11Extension(
    "perseus._perseus",
    [os.path.join("perseus_py.cpp")],
    include_dirs=[INC],
    cxx_std=17,
    extra_compile_args=["-O3"],
)
setup(packages=["perseus"], ext_modules=[ext], cmdclass={"build_ext": build_ext})
