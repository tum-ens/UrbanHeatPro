from setuptools import setup, find_packages

setup(
    name='UrbanHeatPro',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/tum-ens/UrbanHeatPro',
    license=' GPL-3.0',
    author='TUM ENS',
    author_email='',
    description='A bottom-up model for the simulation of heat demand profiles of urban areas ',
    classifiers=[
            "Programming Language :: Python :: 3.9",
        ],
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "geopandas"
    ],
)
