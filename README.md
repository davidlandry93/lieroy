# lieroy
Lie algebra library for both Python and C++

## Getting started

### Dependencies

Lieroy has Eigen3 as a dependency.

```
sudo apt install libeigen3-dev
```

### Installation

Build and install the library using cmake.

```
git clone git@github.com:davidlandry93/lieroy
cd lieroy
mkdir build && cd build
cmake ..
make
sudo make install
```

## Python bindings

Lieroy ships with python bindings.
To build them you need Boost >= 1.63 and (consequently) CMake >= 3.7.2.
Boost 1.63 was the first version to include the Boost Python Numpy library which is used to generate the bindings.

You can build the python package with:

```
cmake .. -DBUILD_PYTHON_BINDINGS=1
make
```

After that the python package will be available in `build/python`.
You can install it with:

```
pip install build/python
```
