# setup.py with cmake extension following
# https://stackoverflow.com/questions/42585210/

import os
import pathlib
import sys
from setuptools import setup, Extension
from setuptools import find_packages
from setuptools.command.build_ext import build_ext as build_ext_orig

# versioneer (+dependencies) does not work in a pep518/7 context w/o modification here
sys.path.append(os.path.dirname(__file__))
import versioneer  # noqa


class CMakeExtension(Extension):
    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        config = "Debug" if self.debug else "Release"
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + str(extdir.parent.absolute()),
            "-DCMAKE_BUILD_TYPE=" + config,
        ]

        # example of build args
        build_args = ["--config", config, "--", "-j4"]

        os.chdir(str(build_temp))
        self.spawn(["cmake", str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(["cmake", "--build", "."] + build_args)
        os.chdir(str(cwd))


cmdclass = versioneer.get_cmdclass()
cmdclass["build_ext"] = build_ext

setup(
    version=versioneer.get_version(),
    name="pymor_dealii",
    author="pyMOR developers",
    author_email="contact@pymor.org",
    python_requires=">=3.7",
    maintainer="Rene Fritze",
    maintainer_email="rene.fritze@wwu.de",
    install_requires=[],
    package_dir={"": "src"},
    packages=find_packages("src"),
    include_package_data=True,
    ext_modules=[CMakeExtension("pymor_dealii")],
    cmdclass=cmdclass,
)
