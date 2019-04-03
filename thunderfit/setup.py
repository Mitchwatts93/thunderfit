import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='thunderfit',
    python_requires='>3.6',
    version='1.0.4',
    description='Thunderfit fitting code',
    long_description='Routines to allow robust fitting to data. Mainly built for Raman analysis but flexible enough for '
                     'most data types',
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    keywords='',
    url='https://github.com/Mitchwatts93/thunderfit',
    author='https://github.com/Mitchwatts93',
    author_email='',
    license='MIT',
    packages=setuptools.find_packages(),
    install_reqs = ['jsonschema=2.6.0',
        'dill=0.2.9',
        'scipy=1.1.0',
        'numpy=1.15.1',
        'numpy-base=1.15.1',
        'numpydoc=0.8.0',
        'matplotlib=2.2.3',
        'pandas=0.23.4',
        'pandoc=1.19.2.1',
        'pandocfilters=1.4.2',
        'lmfit=0.9.11'],
    include_package_data=True)