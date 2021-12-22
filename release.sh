#!/bin/bash
# Build a new release including a git tag and and publish it to PyPI.
# Do some basic checks to avoid mistakes.
#
# Usage: ./release.sh 0.1.0

set -eu

version=$1

# Check the version number.
setup_py_version="$(python setup.py --version)"
if [ "$setup_py_version" != "$version" ]; then
    echo "setup.py has version `$setup_py_version`, not `$version`."
    echo "Aborting."
    exit 1
fi

# Check that the working directory is clean.
if [ ! -z "$(git status --untracked-files=no --porcelain)" ]; then
    echo "You have uncommitted git changes."
    echo "Aborting."
    exit 1
fi

# Check that we are on master.
branch="$(git rev-parse --abbrev-ref HEAD)"
if [ "$branch" != "master" ]; then
    echo "You are on git branch `$branch` but should release from `master`."
    echo "Aborting."
    exit 1
fi

# Create and push a git tag.
git tag -m "Version v$version" "v$version"
git push --tags

# Compile and test, allowing test failures.
make clean compile
set +e
pytest
set -e

# Create and publish the package.
rm -rf dist build *.egg-info
python setup.py sdist
twine upload dist/*

echo "Version v$version has been published on PyPI and has a git tag."
echo "Please make a GitHub Release from that tag."
