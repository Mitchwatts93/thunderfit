import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='thunderfit',
    python_requires='>3.6',
    version='1.0.7.4',
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
    install_requires = ['jsonschema>=2.6.0',
        'dill>=0.2.9',
        'scipy>=1.2.1',
        'numpy>=1.16.2',
        'matplotlib>=2.2.3',
        'pandas>=0.23.4',
        'lmfit>=0.9.11'],
    include_package_data=True)