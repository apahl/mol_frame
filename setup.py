from setuptools import setup

setup(
    name="mol_frame",
    version="0.1.0",
    description="Chemical structure handling for Pandas dataframes",
    url="https://github.com/apahl/mol_frame",
    author="Axel Pahl",
    author_email="",
    license="MIT",
    packages=["mol_frame"],
    install_requires=[
        "rdkit",
        "pandas",
        "matplotlib",
        "cairocffi",
        "holoviews",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
