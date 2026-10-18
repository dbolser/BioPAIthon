# Build image for Anthropic's OSS Scanner (https://github.com/anthropics/oss-scanner).
# The scanner builds this with the repository as the context, then works in the
# image with no network, so everything the tests need is installed here.
# threat_model.md, beside this file, says how to use the image.

# Pinned by digest, as .github/workflows pins its actions. Bump it deliberately.
FROM python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1

# The extensions get the flags Python was built with (-O3 -DNDEBUG ...), as in a
# release, plus those in CC. Flags go in CC, as ci.yml does it: setting CFLAGS
# would replace Python's, not add to them.
ENV CC="gcc -g -fno-omit-frame-pointer" \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_ROOT_USER_ACTION=ignore

# gcc compiles the C extensions and brings the ASan and UBSan runtimes.
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# scanner contract: the checkout lives inside the image, at /src
COPY . /src
WORKDIR /src

# Every build below uses --no-build-isolation, so the same commands rebuild Bio
# offline, against this setuptools.
RUN pip install "setuptools>=77"

# A second copy of Bio in /opt/asan, its C extensions built with ASan and UBSan
# as the test_linux job in .github/workflows/ci.yml builds them. Run it with
# python-asan, below. Built first, while /src holds no compiled extensions.
RUN CC="$CC -fsanitize=address -fsanitize=undefined" \
      pip install --no-deps --no-build-isolation --target /opt/asan . \
    && rm -rf build

# Bio itself, an editable install as AGENTS.md sets it up: the extensions are
# built in the source tree, and an edit under /src takes effect at once. The
# "test" extra in pyproject.toml is everything the offline test suite uses.
RUN pip install --no-build-isolation -e ".[test]"

# python-asan: python, with the ASan+UBSan copy of Bio first on the path. The
# settings are ci.yml's, except two. PYTHONMALLOC=malloc sends small
# allocations to malloc too, so ASan sees their bounds. UBSan stops at its
# first report, as ASan does, so a report cannot hide in pytest's captured
# output. -P stops /src/Bio in the current directory from shadowing /opt/asan.
# Written with printf, not a heredoc, which Docker 20.10 (Debian 12) cannot parse.
RUN printf '%s\n' '#!/bin/sh' \
      'export PYTHONPATH="/opt/asan${PYTHONPATH:+:$PYTHONPATH}"' \
      'export PYTHONMALLOC="${PYTHONMALLOC:-malloc}"' \
      'export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=0}"' \
      'export UBSAN_OPTIONS="${UBSAN_OPTIONS:-print_stacktrace=1:halt_on_error=1}"' \
      'export LD_PRELOAD="$(gcc -print-file-name=libasan.so) $(gcc -print-file-name=libstdc++.so.6)"' \
      'exec python -P "$@"' \
      > /usr/local/bin/python-asan \
    && chmod 755 /usr/local/bin/python-asan

# Check that every extension imports in both builds, then run the offline test
# suite from Tests/, as AGENTS.md does. A failing test does not abort the
# build, and a hang cannot use up its time.
RUN python .github/check_wheel.py \
    && python-asan -c "import Bio; assert Bio.__file__.startswith('/opt/asan/'), Bio.__file__" \
    && python-asan .github/check_wheel.py
RUN cd Tests && (timeout 900 python -m pytest --offline -n auto -q -p no:cacheprovider \
    || echo "WARNING: the offline test suite reported failures")
