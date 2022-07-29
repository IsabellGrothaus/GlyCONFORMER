from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
classifiers = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Education',
    'Programming Language :: Python :: 3'
    ]

setup(name='GlyCONFORMER',
      version='0.0.1',
      description="GlyCONFORMERS is a Python package that assigns conformer strings to N-glycan conformers, based on their torsion angle values.",
      url="https://github.com/IsabellGrothaus/GlyCONFORMER",
      long_description = (here / "README.rst").read_text(encoding="utf-8"),
      author='Isabell Grothaus',
      author_email='grothaus@uni-bremen.de',
      classifiers=classifiers,
      keywords="glycan, conformer, classification, sugar, carbohydrates",
      license='GPL-3.0 license',
      packages=find_packages(),
      install_requires=[
                       "numpy","matplotlib","pandas","plumed","scipy"
                       ],
      )
