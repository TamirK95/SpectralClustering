from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    version='0.1.0',
    author="Tamir Kashi",
    author_email="tamirkashi@gmail.com",
    description="spkmeans algorithm implementation",
    install_requires=['invoke'],
    packages=find_packages(),

    license='GPL-2',
    classifiers=[
        'Development status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'spkmeansmodule',
            ['spkmeansmodule.c', 'spkmeans.c'],
        ),
    ]
)
