from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzymekinetics-skill",
    version="2.0.0",
    author="Qi Gao",
    author_email="joan.gao@seezymes.com",
    description="An Intelligent Tool for Automated Enzyme Kinetic Parameter Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seezymes/EnzymeKinetics-Skill",
    project_urls={
        "Bug Tracker": "https://github.com/seezymes/EnzymeKinetics-Skill/issues",
        "Documentation": "https://enzymekinetics-skill.readthedocs.io",
        "Source Code": "https://github.com/seezymes/EnzymeKinetics-Skill",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "pandas>=1.3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "sphinx>=4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "enzymekinetics=enzymekinetics.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
