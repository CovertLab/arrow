#!/bin/bash
# Build a new release.
# Do some basic checks to avoid mistakes.
#
# Usage: ./release.sh 0.1.0

set -eu

version=$1

# Compile and test, allowing test failures.
make clean compile
set +e
pytest
set -e

# Create the package.
# TODO: Why can't it find numpy without `--no-isolation`?
rm -rf dist build *.egg-info
python -m build --no-isolation


# Cross-check the version number.
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


echo "Do additional testing, then run push.sh"
