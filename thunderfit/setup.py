import setuptools

def readme():
    with open("README.rst") as f:
        return f.read()


setuptools.setup(
    name='thunderfit',
    python_requires='>3.6',
    version='1.0.7.14',
    description='Thunderfit fitting code',
    long_description=readme(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    keywords='raman analysis fitting',
    url='https://github.com/Mitchwatts93/thunderfit',
    author='https://github.com/Mitchwatts93',
    author_email='',
    license='MIT',
    packages=['thunderfit'],
    install_requires = ['jsonschema>=2.6.0',
        'dill>=0.2.9',
        'scipy>=1.2.1',
        'numpy>=1.16.2',
        'matplotlib>=2.2.3',
        'pandas>=0.23.4',
        'lmfit>=0.9.11'],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose'],
    entry_points={
        "console_scripts": ['ramananalyse = thunderfit.thundobj:main']},
    )

