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

(sec_installation)=

# Installation

Pyslim can be installed using the standard [pip](https://pypi.org/project/pyslim/) distribution method:

```bash
python3 -m pip install pyslim
```

To install a [different version](https://pypi.org/project/pyslim/#history), e.g., 
the 1.0 beta release whose version string is "1.0b1", just append this, like:

```bash
python3 -m pip install pyslim==1.0b1
```

Alternatively, you can install from source as described in [](sec_development).


## Troubleshooting

If you find a bug in ``pyslim`` or want to suggest an improvement, please
[open an issue](https://github.com/tskit-dev/pyslim/issues) on our github page.
If you have a question about using tree sequences,
please ask it at [the tskit discussion page](https://github.com/tskit-dev/tskit/discussions).
Finally, questions about SLiM should be directed to
[the SLiM mailing list](https://groups.google.com/forum/#!forum/slim-discuss).
