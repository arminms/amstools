# amstools
[![Build and Test](https://github.com/arminms/amstools/actions/workflows/build-n-test.yml/badge.svg)](https://github.com/arminms/amstools/actions/workflows/build-n-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Collection of high-performance biological sequence tools developed in mixed
`C` (file i/o) and `C++` (parallel algorithms) with minimal dependencies (only
[zlib](http://www.zlib.net/)).

## List of Tools
* `acgt` – Print residue statistics and optionally GC and AT contents.
* `ngx`  – Print the contiguity statistics (_e.g._ _N50_, _L50_).
* `sc`   – Print sequence and residue counts.

## Building from Source
`amstools` is cross-platform (_Linux_/_macOS_/_Windows_). Using `CMake`
 you can build and install all the tools by the same commands.

### Prerequisites
* [CMake](https://cmake.org/) version 3.16 or higher
* `C++` compiler supporting the `C++11` standard (_e.g._ `gcc 4.8`)
* [zlib](http://www.zlib.net/) – to install `zlib` under _Windows_ you can
use [vcpkg](https://vcpkg.io/)

### Build and Install
```
$ git clone https://github.com/arminms/amstools.git
$ cd amstools
$ cmake -S . -B build
$ cmake --build build -j 3 # using 3 concurrent threads
$ cmake --install build
```
The last command on _Linux_ and _macOS_ needs to be preceded by `sudo`.