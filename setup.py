import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bp_tools_dzorine",  # Replace with your own username
    version="0.0.1",
    author="Dmitri Zorine",
    author_email="dzorine@gmail.com",
    description="Tools for fragment assembly runs in protein design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=["bp_tools"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
