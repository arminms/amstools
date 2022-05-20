# amstools
[![Build and Test](https://github.com/arminms/amstools/actions/workflows/build-n-test.yml/badge.svg)](https://github.com/arminms/amstools/actions/workflows/build-n-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Collection of fast and efficient biological sequence tools developed in mixed
`C` (file i/o) and `C++` (parallel algorithms) with minimal dependencies (only
[zlib](http://www.zlib.net/)).

## List of Tools
* `acgt` – Print residue statistics and optionally GC and AT contents.
* `ngx`  – Print the contiguity statistics (_e.g._ `N50`, `L50`).
* `sc`   – Print sequence and residue counts.

## Building from Source
`amstools` is cross-platform (`Linux`/`macOS`/`Windows`). Using
[CMake](https://cmake.org/) you can build and install all the tools
by the same commands.

### Prerequisites
* `CMake` version 3.16 or higher
* `C++` compiler supporting the `C++11` standard (_e.g._ `gcc 4.8`)
* [zlib](http://www.zlib.net/) – to install `zlib` under `Windows` you can
use [vcpkg](https://vcpkg.io/)

### Build and Install
```
$ git clone https://github.com/arminms/amstools.git
$ cd amstools
$ cmake -S . -B build
$ cmake --build build -j 3 # using 3 concurrent threads
$ cmake --install build
```