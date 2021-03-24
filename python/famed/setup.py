import setuptools

exec(open('famed/version.py').read)

setuptools.setup(
    name="famed", # Replace with your own username
    version=__version__,
    author="FAMED Authors",
    author_email="author@example.com",
    description="A python implementation of FAMED",
    long_description=open('README.md').read(),
    #long_description_content_type="text/markdown",
    url="https://github.com/jeanm12/py_famed",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science Research",
        "Topic :: Scientific/Engineering :: Astronomy"
    ],
)
