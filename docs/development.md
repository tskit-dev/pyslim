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

All contributions, bug reports, documentation improvements and ideas are welcome. If you think
there is anything missing, please open an [issue](https://github.com/tskit-dev/pyslim/issues)
or [pull request](https://github.com/tskit-dev/pyslim/pulls) on GitHub.

See the [tskit developer documentation](https://tskit.dev/tskit/docs/stable/development.html)
for the general development workflow (git, prek, testing, documentation).

Install development dependencies with:

```bash
uv sync
```

Run the tests with:

```bash
uv run pytest
```
