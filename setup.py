"""Setup script for lotus package."""

from pathlib import Path

from setuptools import find_packages, setup

# Read the README file if it exists
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

# Read requirements from requirements.txt
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, encoding="utf-8") as f:
        install_requires = [
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        ]
else:
    # Fallback to minimal requirements if requirements.txt doesn't exist
    install_requires = [
        "numpy",
        "pandas",
        "anndata",
        "scanpy",
    ]

setup(
    name="lotus",
    version="0.1.0",
    description="Lotus: A single-cell RNA sequencing analysis package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Lotus Team",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=install_requires,
)

