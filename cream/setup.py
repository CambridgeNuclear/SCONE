from setuptools import setup, find_packages

kwargs = {
    # Basic Info
    'name': 'cream',
    'version': '0.1.0',
    'description': 'Command Line Utility for SCONE MC Code',
    'licence': 'MIT',

    # Packages & Scripts
    'packages': find_packages(exclude=['tests*']),
    'entry_points': {
        'console_scripts': [
            'cream = cream.__main__:main'
        ]
    },

    # Clasifiers
    'classifiers': [
        'Development Status :: 1 -Planning',
        'Programming Language :: Fortran',
    ],


    # Depedencies
    'python_requires': '>=3.6',
    # 'instal_requires': [],
    'extras_require': {
            'test': ['pytest']
    }
}

setup(**kwargs)
