from setuptools import setup, find_packages

setup(name='GlyCONFORMER',
      version='0.1',
      description="GlyCONFORMERS is a Python package that assigns conformer strings to N-glycan conformers, based on their torsion angle values.",
      long_description="",
      author='Isabell Grothaus',
      author_email='grothaus@uni-bremen.de',
      license='GPL-3.0 license',
      packages=find_packages(),
      zip_safe=False,
      install_requires=[
          "numpy","matplotlib","pandas","plumed","scipy","json","panedr","os","glob","sys"
          # 'Sphinx',
          # ^^^ Not sure if this is needed on readthedocs.org
          # 'something else?',
          ],
      )
