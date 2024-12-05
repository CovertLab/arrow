#!/bin/bash
# Tag a new release and publish it on PyPI.

set -eu

version="$(python setup.py --version)"

# Create and push a git tag.
git tag -m "Version v$version" "v$version"
git push --tags

twine upload dist/*

echo "Version v$version has been published on PyPI and has a git tag."
echo "Please make a GitHub Release from that tag."
