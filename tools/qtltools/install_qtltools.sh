#!/usr/bin/env bash
# Build QTLtools v1.3.1 from source on Linux (x86_64).
#
#   bash tools/qtltools/install_qtltools.sh [install_dir]     (default: tools)
#
# Produces <install_dir>/bin/QTLtools. Then set in your config:
#   qtl:
#     qtltools_binary: "tools/bin/QTLtools"
#
# QTLtools links statically against HTSlib and R's standalone math library
# (libRmath), so both are built here from source; Boost, GSL and compression
# libraries come from the system. On Ubuntu/Debian install them with:
#
#   sudo apt-get install -y build-essential gfortran curl bzip2 \
#     libboost-iostreams-dev libboost-program-options-dev libgsl-dev \
#     zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
#     libpcre2-dev
#
# On CentOS/RHEL: gcc-c++ gcc-gfortran boost-devel boost-static gsl-devel
#   zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel pcre2-devel
#   (then set BOOST_LIB=/usr/lib64 when running this script).
#
# Versions and checksums are pinned below.
set -euo pipefail

INSTALL_DIR="${1:-tools}"
QTLTOOLS_URL="https://qtltools.github.io/qtltools/binaries/QTLtools_1.3.1_source.tar.gz"
QTLTOOLS_SHA256="8d812145a210fe96ae31590423a0eab6bf9e3f46b33a62c8b04b0d428bcfed8d"
HTSLIB_VERSION="1.21"
R_VERSION="4.4.1"
BOOST_INC="${BOOST_INC:-/usr/include}"
BOOST_LIB="${BOOST_LIB:-/usr/lib/x86_64-linux-gnu}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 2)}"

mkdir -p "$INSTALL_DIR/bin"
INSTALL_DIR="$(cd "$INSTALL_DIR" && pwd)"
BUILD="$INSTALL_DIR/build"
mkdir -p "$BUILD"
cd "$BUILD"

echo "== HTSlib $HTSLIB_VERSION"
if [[ ! -f "htslib-$HTSLIB_VERSION/libhts.a" ]]; then
  curl -fsSL -o htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2"
  tar -xjf htslib.tar.bz2
  (cd "htslib-$HTSLIB_VERSION" && ./configure --disable-libcurl --disable-gcs --disable-s3 --without-libdeflate && make -j "$JOBS" lib-static)
fi

echo "== R $R_VERSION standalone math library"
if [[ ! -f "R-$R_VERSION/src/nmath/standalone/libRmath.a" ]]; then
  curl -fsSL -o R.tar.gz "https://cran.r-project.org/src/base/R-4/R-$R_VERSION.tar.gz"
  tar -xzf R.tar.gz
  (cd "R-$R_VERSION" && ./configure --with-readline=no --with-x=no --without-recommended-packages \
      --disable-java --enable-R-shlib=no > configure.log && cd src/nmath/standalone && make -j "$JOBS" static)
fi

echo "== QTLtools 1.3.1"
curl -fsSL -o qtltools.tar.gz "$QTLTOOLS_URL"
echo "$QTLTOOLS_SHA256  qtltools.tar.gz" | sha256sum -c -
rm -rf qtltools && tar -xzf qtltools.tar.gz
cd qtltools
make -j "$JOBS" \
  BOOST_INC="$BOOST_INC" BOOST_LIB="$BOOST_LIB" \
  RMATH_INC="$BUILD/R-$R_VERSION/src/include" RMATH_LIB="$BUILD/R-$R_VERSION/src/nmath/standalone" \
  HTSLD_INC="$BUILD/htslib-$HTSLIB_VERSION" HTSLD_LIB="$BUILD/htslib-$HTSLIB_VERSION"
cp bin/QTLtools "$INSTALL_DIR/bin/QTLtools"

"$INSTALL_DIR/bin/QTLtools" --help 2>&1 | head -5 || true
echo "Installed: $INSTALL_DIR/bin/QTLtools"
