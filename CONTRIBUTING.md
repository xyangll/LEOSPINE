# Contributing to SPINE

Thank you for your interest in contributing to SPINE! This document provides guidelines and information for contributors.

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [How to Contribute](#how-to-contribute)
3. [Development Workflow](#development-workflow)
4. [Code Standards](#code-standards)
5. [Testing](#testing)
6. [Documentation](#documentation)
7. [Pull Request Process](#pull-request-process)

---

## Getting Started

### Development Environment

1. **Fork and Clone**
   ```bash
   git clone https://github.com/xyangll/LEOSPINE.git
   cd LEOSPINE
   ```

2. **Create Virtual Environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/macOS
   # or: venv\Scripts\activate  # Windows
   ```

3. **Install Dependencies**
   ```bash
   # Option 1: Install from requirements files
   pip install -r requirements.txt
   pip install -r requirements-dev.txt
   
   # Option 2: Install as a package (recommended for development)
   pip install -e .
   # Or with dev dependencies:
   pip install -e ".[dev]"
   ```

4. **Verify Setup**
   ```bash
   pytest tests/ -v
   python main.py
   ```

---

## How to Contribute

### Reporting Bugs

Before creating a bug report, please check existing issues to avoid duplicates.

**Include in your report:**
- SPINE version and Python version
- Operating system
- Steps to reproduce
- Expected vs actual behavior
- Screenshots (if applicable)
- Sample data files (if applicable)

### Suggesting Features

Feature requests are welcome! Please describe:
- The problem you're trying to solve
- Your proposed solution
- Alternative solutions considered
- Any additional context

### Code Contributions

We welcome:
- Bug fixes
- New features
- Documentation improvements
- Test coverage improvements
- Performance optimizations

---

## Development Workflow

### Branch Naming

Use descriptive branch names:

| Type | Format | Example |
|------|--------|---------|
| Feature | `feature/description` | `feature/add-glonass-support` |
| Bug fix | `fix/description` | `fix/positioning-convergence` |
| Documentation | `docs/description` | `docs/update-api-reference` |
| Refactor | `refactor/description` | `refactor/sim-data-cleanup` |

### Workflow Steps

1. **Create Branch**
   ```bash
   git checkout -b feature/your-feature
   ```

2. **Make Changes**
   - Write code
   - Add tests
   - Update documentation

3. **Test Locally**
   ```bash
   pytest tests/ -v
   python main.py  # Manual testing
   ```

4. **Commit**
   ```bash
   git add .
   git commit -m "feat: add your feature description"
   ```

5. **Push and Create PR**
   ```bash
   git push origin feature/your-feature
   ```

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(constellation): add RAAN spread configuration
fix(positioning): correct relativistic effect formula
docs(api): update sim_data documentation
```

---

## Code Standards

### Python Style Guide

- Follow [PEP 8](https://pep8.org/)
- Use [Black](https://black.readthedocs.io/) for formatting
- Maximum line length: 120 characters (configured in `pyproject.toml`)
- Use type hints where practical

### Formatting

```bash
# Format code
black .

# Check style
flake8 core/ app/ tools/
```

### Code Organization

```python
# Standard library imports
import os
import sys

# Third-party imports
import numpy as np
from PySide6.QtWidgets import QWidget

# Local imports
from core.utilities import geodetic_to_ecef
```

### Documentation

- Use docstrings for all public functions/classes
- Follow NumPy docstring format
- Include type hints

```python
def compute_position(observations: np.ndarray, 
                     initial_guess: np.ndarray) -> tuple[np.ndarray, float]:
    """
    Compute receiver position using least squares.
    
    Parameters
    ----------
    observations : np.ndarray
        Observation array with shape (n_obs, 4)
    initial_guess : np.ndarray
        Initial position estimate [x, y, z]
    
    Returns
    -------
    position : np.ndarray
        Estimated position [x, y, z]
    residual : float
        RMS residual
    
    Raises
    ------
    ValueError
        If fewer than 4 observations provided
    
    Examples
    --------
    >>> obs = np.array([[...], [...], [...], [...]])
    >>> pos, res = compute_position(obs, np.array([0, 0, 0]))
    """
    pass
```

---

## Testing

### Running Tests

```bash
# All tests
pytest tests/ -v

# Specific test file
pytest tests/test_imports.py -v

# With coverage
pytest tests/ --cov=core --cov-report=html
```

### Writing Tests

- Place tests in `tests/` directory
- Name test files `test_*.py`
- Name test functions `test_*`
- Use descriptive test names

```python
# tests/test_positioning.py
import pytest
import numpy as np
from core.positioning import run_lsq_positioning

class TestPositioning:
    def test_basic_positioning(self):
        """Test basic positioning with valid observations."""
        pos_obj, result = run_lsq_positioning(
            sim_data_file='tests/data/sample_obs.csv',
            tle_file='InputData/test_con.tle',
            use_doppler=True
        )
        assert result is not None
        assert 'positions' in result
        assert len(result['positions']) > 0
    
    def test_insufficient_satellites(self):
        """Test error handling with too few satellites."""
        pos_obj, result = run_lsq_positioning(
            sim_data_file='tests/data/few_sats.csv',
            use_pseudorange=True,
            use_doppler=False
        )
        # Should handle gracefully or raise appropriate error
        assert result is None or len(result.get('positions', [])) == 0
```

### Test Data

- Store test data in `tests/data/`
- Use minimal data files
- Never commit private or proprietary data

---

## Documentation

### Updating Docs

- API changes → Update `docs/api.md`
- New features → Update `docs/user_guide.md`
- Architecture changes → Update `docs/architecture.md`
- Bug fixes → Update `CHANGELOG.md`

### Building Docs

Documentation is in Markdown format and can be viewed directly on GitHub or with any Markdown viewer.

---

## Pull Request Process

### Before Submitting

- [ ] Code follows style guidelines
- [ ] Tests pass locally
- [ ] Documentation updated (if needed)
- [ ] CHANGELOG.md updated (for features/fixes)
- [ ] No merge conflicts

### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Refactoring

## Testing
Describe testing performed

## Screenshots (if applicable)

## Checklist
- [ ] Code follows project style
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] CHANGELOG updated
```

### Review Process

1. PR created → Automated checks run
2. Maintainer reviews code
3. Address feedback (if any)
4. Approval and merge

### After Merge

- Delete your branch
- Pull latest changes to main
- Celebrate! 🎉

---

## Questions?

- Open a [Discussion](https://github.com/xyangll/LEOSPINE/discussions)
- Create an [Issue](https://github.com/xyangll/LEOSPINE/issues)

---

Thank you for contributing to SPINE! 🛰️
