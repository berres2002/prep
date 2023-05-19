from setuptools import setup, find_packages

setup(
    packages=find_packages(
      where='.',
      include=['prep*']),
   use_scm_version={
      "write_to": "prep/version.py",
      "write_to_template": "__version__ = '{version}'",
    }
)