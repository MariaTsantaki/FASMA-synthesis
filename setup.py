from setuptools import setup

setup(
    name='FASMA',
    version='2.0',
    description='FASMA-synthesis',
    url='http://github.com/MariaTsantaki/FASMA-synthesis',
    author='Maria Tsantaki',
    license='MIT',
    packages=['FASMA'],
    scripts=['bin/fasma', 'bin/fasma_gui'],
    include_package_data=True,
    zip_safe=False,
)
