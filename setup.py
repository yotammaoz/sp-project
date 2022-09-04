from setuptools import setup, Extension, find_packages

setup(
    name='spkmeansc',
    version='1.0',
    description='kmeans algo in c',
    packages=find_packages(),
    ext_modules=[
        Extension('spkmeansc', ['spkmeansmodule.c'])
    ]
)