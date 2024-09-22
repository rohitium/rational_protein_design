from setuptools import setup, find_packages

setup(
    name="rational_protein_design",
    version="0.1.0",
    description="A package for rational protein binder design using Biopython",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Rohit Satija",
    author_email="rohitsatija0092@gmail.com",
    url="https://github.com/rohitium/rational_protein_design",
    packages=find_packages(exclude=("tests", "docs")),
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.21.0",
        "setuptools>=49.6.0"
    ],
    python_requires='>=3.7',
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        "console_scripts": [
            "design-binder=rational_protein_design.binder_designer:main",  # Optional command-line script
        ],
    },
)
