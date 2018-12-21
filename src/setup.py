import setuptools


setuptools.setup(
    name='src',
    python_requires='>3.5.2',
    version='1.0',
    description='',
    long_description='',
    classifiers=[
        'Development Status :: 1 - Alpha',
        'License :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Physics :: Plotting package',
    ],
    keywords='',
    url='',
    author='',
    author_email='mitchcwatts@gmail.com',
    license='MIT',
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy>=1.15.1',
        'pandas>=0.23.4',
        'tqdm>=4.26.0',
        'matplotlib>=2.2.3'
    ],
    include_package_data=True,
    zip_safe=False)