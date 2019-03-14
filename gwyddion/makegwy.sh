#!/bin/bash

set -e

GWY_VER="2.48" # keep track of current release!
SVN_CMD="svn checkout"
GWY_URL="svn://svn.code.sf.net/p/gwyddion/code/tags/gwyddion-${GWY_VER}"
GWY_DIR="gwyddion"

${SVN_CMD} ${GWY_URL} ${GWY_DIR}

pushd "$GWY_DIR"

set +e
patch --forward < ../Makefile.am.patch # do not complain when the patch has already been applied
set -e

export CFLAGs=-Ofast
./autogen.sh --prefix="`pwd`/build" \
             --enable-maintainer-mode --disable-gtk-doc --without-x \
             --disable-dependency-tracking \
             --without-libiconv-prefix --without-libintl-prefix \
             --disable-schemas-install --disable-desktop-file-update \
             --without-pascal --without-perl --without-python --without-ruby \
             --without-kde4-thumbnailer --disable-pygwy \
             --enable-module-bundling

make
# make -C po update-gmo
make install
popd
