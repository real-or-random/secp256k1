FROM debian:stable-slim

SHELL ["/bin/bash", "-c"]

RUN dpkg --add-architecture i386 && \
    dpkg --add-architecture s390x && \
    dpkg --add-architecture armhf && \
    dpkg --add-architecture arm64 && \
    dpkg --add-architecture ppc64el

# dkpg-dev: to make pkg-config work in cross-builds
# llvm: for llvm-symbolizer, which is used by clang's UBSan for symbolized stack traces
RUN apt-get update && apt-get install --no-install-recommends -y \
        git ca-certificates wget \
        make automake libtool pkg-config dpkg-dev valgrind qemu-user \
        gcc clang llvm libclang-rt-dev libc6-dbg \
        g++ \
        gcc-i686-linux-gnu libc6-dev-i386-cross libc6-dbg:i386 libubsan1:i386 libasan8:i386 \
        gcc-s390x-linux-gnu libc6-dev-s390x-cross libc6-dbg:s390x \
        gcc-arm-linux-gnueabihf libc6-dev-armhf-cross libc6-dbg:armhf \
        gcc-aarch64-linux-gnu libc6-dev-arm64-cross libc6-dbg:arm64 \
        gcc-powerpc64le-linux-gnu libc6-dev-ppc64el-cross libc6-dbg:ppc64el

WORKDIR /root

# Build and install gcc snapshot
ARG GCC_SNAPSHOT_MAJOR=14
RUN mkdir gcc && cd gcc && \
    wget --progress=dot:giga --https-only --recursive --accept '*.tar.xz' --level 1 --no-directories "https://gcc.gnu.org/pub/gcc/snapshots/LATEST-${GCC_SNAPSHOT_MAJOR}" && \
    wget "https://gcc.gnu.org/pub/gcc/snapshots/LATEST-${GCC_SNAPSHOT_MAJOR}/sha512.sum" && \
    sha512sum --check --ignore-missing sha512.sum && \
    # We should have downloaded exactly one tar.xz file
    ls && \
    [[ $(ls *.tar.xz | wc -l) -eq "1" ]] && \
    tar xf *.tar.xz && \
    mkdir gcc-build && cd gcc-build && \
    apt-get update && apt-get install --no-install-recommends -y libgmp-dev libmpfr-dev libmpc-dev flex && \
    ../*/configure --prefix=/opt/gcc-snapshot --enable-languages=c --disable-bootstrap --disable-multilib --without-isl && \
    make -j $(nproc) && \
    make install && \
    apt-get autoremove -y libgmp-dev libmpfr-dev libmpc-dev flex && \
    apt-get clean && \
    cd ../.. && rm -rf gcc && \
    ln -s /opt/gcc-snapshot/bin/gcc /usr/bin/gcc-snapshot

# Install clang snapshot
RUN wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc && \
    # Add repository for this Debian release
    . /etc/os-release && echo "deb http://apt.llvm.org/${VERSION_CODENAME} llvm-toolchain-${VERSION_CODENAME} main" >> /etc/apt/sources.list && \
    # Install clang snapshot
    apt-get update && apt-get install --no-install-recommends -y clang && \
    # Remove just the "clang" symlink again
    apt-get remove -y clang && \
    # We should have exactly two clang versions now
    ls /usr/bin/clang* && \
    [[    $(ls /usr/bin/clang-?? | sort | wc -l) -eq "2" ]] && \
    # Create symlinks for them
    ln -s $(ls /usr/bin/clang-?? | sort | tail -1) /usr/bin/clang-snapshot && \
    ln -s $(ls /usr/bin/clang-?? | sort | head -1) /usr/bin/clang
