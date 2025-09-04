from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

setup(
    name="pairwiseANIviz",
    version="1.3",
    author="Runjia Ji",
    author_email="jirunjia@gmail.com",
    description="Pairwise ANI (Average Nucleotide Identity) visualization tool",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/RunJiaJi/pairwiseANIviz",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
    ],
    entry_points={
        "console_scripts": [
            "pairwiseANIviz=pairwiseANIviz:main",
        ],
    },
    include_package_data=True,
    package_data={
        "pairwiseANIviz": ["static/*"],
    },
    keywords="bioinformatics, ANI, visualization, clustering, heatmap, microbiology",
    project_urls={
        "Bug Reports": "https://github.com/RunJiaJi/pairwiseANIviz/issues",
        "Source": "https://github.com/RunJiaJi/pairwiseANIviz",
        "Documentation": "https://github.com/RunJiaJi/pairwiseANIviz#readme",
    },
)
