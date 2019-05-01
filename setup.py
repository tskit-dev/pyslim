from codecs import open as codecs_open
from setuptools import setup, find_packages
from warnings import warn


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
      install_requires=['msprime>=0.7.0', 'tskit', 'kastore', 'attrs'],
      extras_require={
          'dev': [],
      },
      setup_requires=[],
      )
