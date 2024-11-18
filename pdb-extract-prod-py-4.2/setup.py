from setuptools import find_packages
from setuptools import setup
# from setuptools.command.egg_info import egg_info


# class egg_info_ex(egg_info):
#     """Includes license file into `.egg-info` folder."""

#     def run(self):
#         # don't duplicate license into `.egg-info` when building a distribution
#         if not self.distribution.have_run.get('install', True):
#             # `install` command is in progress, copy license
#             self.mkpath(self.egg_info)
#             self.copy_file('LICENSE.txt', self.egg_info)

#         egg_info.run(self)


setup(
    name="pdb-extract-prod-py",
    version="4.2",
    author="Chenghua Shao",
    author_email="chenghua.shao@rcsb.org",
    description="pdb_extract is a pre-deposition service for assembling structure files for wwPDB OneDep deposition",
    long_description="See README.txt",
    # url="https://github.com/rcsb/test.git",
    license="Apache 2.0",
    license_files = ("LICENSE.txt",),
    # cmdclass = {'egg_info': egg_info_ex},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS",
        "Operating System :: Linux"
    ],
    packages=find_packages(exclude=["tests"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg", "*.sh", "*.csh, *.dic, *.cif, *.pdb"]
    },
    include_package_data=True,
    test_suite="tests",
    python_requires=">=3.6",
)
