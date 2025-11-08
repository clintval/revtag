# revtag

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/revtag/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/revtag/badges/version.svg)](http://bioconda.github.io/recipes/revtag/README.html)
[![Build Status](https://github.com/clintval/revtag/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/clintval/revtag/actions/workflows/rust.yml)
[![Coverage Status](https://coveralls.io/repos/github/clintval/revtag/badge.svg?branch=main)](https://coveralls.io/github/clintval/revtag?branch=main)
[![Language](https://img.shields.io/badge/language-rust-a72144.svg)](https://www.rust-lang.org/)

Reverse (and complement) array-like SAM tags for reverse alignments.

![Waimea - Kauaʻi](.github/img/cover.jpg)

Install with the Conda or Mamba package manager after setting your [Bioconda channels](https://bioconda.github.io/#usage):

```bash
❯ mamba install revtag
```

### Example Usage

```bash
❯ revtag --rev ad ae aq bd be bq cd ce --revcomp ac bc < input.bam > output.bam
```
