from codecs import open as codecs_open
from setuptools import setup, find_packages
from warnings import warn
import os


# Get the long description from the relevant file
with codecs_open('README.md', encoding='utf-8') as f:
    long_description = f.read()

try:
    import tskit
except ImportError:
    warn("`tskit` not present and must be installed")

try:
    import kastore
except ImportError:
    warn("`kastore` not present and must be installed")

# After exec'ing this file we have tskit_version defined.
tskit_version = None  # Keep PEP8 happy.
version_file = os.path.join("pyslim", "_version.py")
with open(version_file) as f:
    exec(f.read())

setup(name='pyslim',
      version=pyslim_version,
      description=u"Manipulate tree sequences produced by SLiM.",
      long_description=long_description,
      classifiers=[],
      keywords=['tree sequences', 'tskit'],
      author=u"Peter Ralph",
      author_email='petrel.harp@gmail.com',
      url='https://github.com/tskit-dev/pyslim',
      license='MIT',
      packages=['pyslim'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'msprime>=0.7.0',
          'tskit',
          'kastore',
          'attrs',
          'numpy'],
      extras_require={
          'dev': [],
      },

      setup_requires=[],
      project_urls={
          'Bug Reports': 'https://github.com/tskit-dev/pyslim/issues',
          'Source': 'https://github.com/tskit-dev/pyslim',
      },
      )
