from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["numpy>=1.18.0", "scipy>=1.3.3"]

setup(
    name="pyGSW",
    version="0.0.2",
    author="Mironenko Artem",
    author_email="artem.mironenko1998@mail.ru",
    description="FHE GSW",
    long_description=readme,
    long_description_content_type='text/markdown',

    url="https://github.com/Homomorphic-encryption-GSW/pyGSW/",

    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
    ],
)