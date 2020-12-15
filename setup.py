
import setuptools

setuptools.setup(
    name='mot',
    version='1.0',
    packages=['mot', 'mot.demo'],
    package_dir={
        'mot':'src/',
        'mot.demo':'src/demo'
    },
    install_requires=[
        'matplotlib',
        'numpy',
        'seaborn',
        'scikit-learn',
        'tqdm',
    ]
)
