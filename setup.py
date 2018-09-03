from codecs import open as codecs_open
from setuptools import setup, find_packages
from warnings import warn


# Get the long description from the relevant file
with codecs_open('README.md', encoding='utf-8') as f:
    long_description = f.read()

HAVE_MSPRIME = False
try:
    import msprime
except ImportError:
    warn("`msprime` not present and must be installed")


setup(name='pyslim',
      version='0.1',
      description=u"Manipulate tree sequences produced by SLiM.",
      long_description=long_description,
      classifiers=[],
      keywords='',
      author=u"Peter Ralph",
      author_email='petrel.harp@gmail.com',
      url='https://github.com/tskit-dev/pyslim',
      license='MIT',
      packages=find_packages(exclude=[]),
      include_package_data=True,
      zip_safe=False,
      install_requires=['msprime', 'attrs'],
      extras_require={
          'dev': [],
      },
      setup_requires=[],
      )
