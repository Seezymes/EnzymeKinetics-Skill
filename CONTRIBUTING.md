# Contributing to EnzymeKinetics-Skill

Thank you for your interest in contributing to EnzymeKinetics-Skill! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with the following information:

- **Title**: Clear and descriptive title
- **Description**: Detailed description of the bug
- **Steps to Reproduce**: Step-by-step instructions to reproduce the bug
- **Expected Behavior**: What you expected to happen
- **Actual Behavior**: What actually happened
- **Environment**: Python version, operating system, package versions
- **Code Sample**: Minimal code example that reproduces the issue

### Suggesting Enhancements

Enhancement suggestions are welcome! Please open an issue with:

- **Title**: Clear and descriptive title
- **Description**: Detailed description of the proposed enhancement
- **Use Case**: Explain why this enhancement would be useful
- **Proposed Implementation**: If you have ideas on how to implement it

### Pull Requests

1. **Fork the repository** and create your branch from `main`.
2. **Make your changes** following the coding standards below.
3. **Add tests** for any new functionality.
4. **Update documentation** as needed.
5. **Ensure all tests pass** by running `pytest`.
6. **Submit a pull request** with a clear description of the changes.

## Coding Standards

### Python Style Guide

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide
- Use 4 spaces for indentation
- Maximum line length: 100 characters
- Use descriptive variable names

### Documentation

- All functions must have docstrings following [Google Style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
- Include type hints where appropriate
- Update README.md if adding new features

### Testing

- Write unit tests for all new functions
- Maintain test coverage above 80%
- Run tests before submitting PR: `pytest tests/ -v --cov=src`

### Commit Messages

- Use clear and meaningful commit messages
- Start with a verb in present tense (e.g., "Add", "Fix", "Update")
- Keep the first line under 72 characters
- Reference issue numbers when applicable

Example:
```
Add support for substrate inhibition model

- Implement substrate_inhibition() function
- Add unit tests for new model
- Update documentation

Fixes #123
```

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/seezymes/EnzymeKinetics-Skill.git
cd EnzymeKinetics-Skill
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install development dependencies:
```bash
pip install -e ".[dev]"
```

4. Run tests:
```bash
pytest tests/ -v --cov=src
```

## Code Review Process

All submissions require review before being merged:

1. At least one maintainer must approve the PR
2. All CI checks must pass
3. Code coverage must not decrease
4. Documentation must be updated if needed

## Areas for Contribution

We particularly welcome contributions in the following areas:

- **New kinetic models**: Allosteric enzymes, inhibition models, etc.
- **Visualization**: Additional plot types, interactive plots
- **Documentation**: Tutorials, examples, API documentation
- **Performance**: Optimization for large datasets
- **Web interface**: Browser-based user interface
- **Database integration**: Connection to enzyme databases (BRENDA, etc.)

## Questions?

If you have questions about contributing, please:

- Open a [GitHub Discussion](https://github.com/seezymes/EnzymeKinetics-Skill/discussions)
- Email: joan.gao@seezymes.com

## Code of Conduct

This project adheres to a code of conduct. By participating, you are expected to:

- Be respectful and inclusive
- Accept constructive criticism gracefully
- Focus on what is best for the community
- Show empathy towards others

Thank you for contributing to EnzymeKinetics-Skill!
