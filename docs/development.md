---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(sec_development)=

# Development


To install a particular version of ``pyslim`` from source, e.g., to obtain a recent update:

```bash
git clone https://github.com/tskit-dev/pyslim.git
cd pyslim
python setup.py install --user
```


Then, to run the tests to make sure everything is working, do:

```bash
python -m pytest tests
```

*Note:* if you use ``python3`` you may need to replace ``python`` with ``python3`` above.

If you would like to add some features to ``pyslim``, please read the
following. If you think there is anything missing,
please open an `issue <http://github.com/tskit-dev/pyslim/issues>`_ or
`pull request <http://github.com/tskit-dev/pyslim/pulls>`_ on GitHub!

## Quickstart

- Make a fork of the pyslim repo on `GitHub <http://github.com/tskit-dev/pyslim>`_
- Clone your fork into a local directory:

```bash
git clone git@github.com:YOUR_GITHUB/pyslim.git
```

- Install the development requirements using
  ``python3 -m pip install -r requirements/development.txt``.
- Run the tests to ensure everything has worked: ``python3 -m pytest tests``. These should
  all pass.
- Make your changes in a local branch, and open a pull request on GitHub when you
  are ready. Please make sure that (a) the tests pass before you open the pull request; and
  (b) your code passes PEP8 checks before opening the pull request.

