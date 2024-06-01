#!/bin/sh

brew install gsl boost automake autoconf libtool
ln -s $(which glibtoolize) /usr/local/bin/libtoolize
