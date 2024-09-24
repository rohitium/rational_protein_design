# Rational Protein Design

[![Build Status](https://github.com/rohitium/rational_protein_design/actions/workflows/ci.yml/badge.svg)](https://github.com/rohitium/rational_protein_design/actions)
[![PyPI version](https://badge.fury.io/py/rational_protein_design.svg)](https://badge.fury.io/py/rational_protein_design)
[![License](https://img.shields.io/github/license/rohitium/rational_protein_design)](LICENSE)

A Python package for rational protein binder design using Biopython.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)

## Installation

You can install the package from PyPI:

```bash
pip install rational_protein_design
```

## Usage

```python
from rational_protein_design import BinderDesigner

designer = BinderDesigner(
    pdb_file="path/to/pdb/file.pdb",
    chain_id="A",
    target_residues_range=(100, 120),
    neighbor_chain_id="B",
    N=50,
    num_seq=10
)
designer.design_binder("output_fasta_file.fasta")
```

## Features

- Parse PDB files and extract chains
- Identify neighboring residues based on distance thresholds
- Extract protein sequences and generate rational binders
- Support for randomization and multiple design attempts

## Contributing

We welcome contributions! Here's how you can help:

### Reporting Bugs

- Search existing issues to avoid duplicates.
- Open a new issue and provide detailed information.

### Feature Requests

- Open an issue with the tag `enhancement`.

### Code Contributions

- Fork the repository.
- Create a new branch for your feature or bugfix.
- Write clear, concise commit messages.
- Submit a pull request.

### Coding Standards

- Follow PEP 8 style guidelines.
- Write docstrings for all public methods and classes.
- Include unit tests for new features or bug fixes.

### Testing

- Run existing tests with `python -m unittest discover tests`.
- Add tests for your contributions.

Thank you for contributing!

## License

This project is licensed under the MIT License.

