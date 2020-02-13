from setuptools import setup

setup(
    name='FASMA',
    version='2.0',
    description='FASMA-synthesis',
    url='http://github.com/MariaTsantaki/FASMA-synthesis',
    author='Maria Tsantaki',
    license='MIT',
    packages=['FASMA'],
    scripts=['bin/fasma'],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=['numpy>=1.7.0', 'scipy>=0.16', 'matplotlib',
    'PyYAML', 'statsmodels', 'patsy', 'pandas>=0.17.0', 'astropy',
    'PyAstronomy', 'periodictable==1.5.2', 'setuptools', 'nose', 'coverage'],
    zip_safe=False,
)
