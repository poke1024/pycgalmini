# Installation

In order to compile this, you need to install:

* gmp
* mpfr
* cmake
* CGAL 5.0
* pybind11

Here's an example how this is done on macOS:

```
brew install boost gmp mpfr cmake
```

Install CGAL 5.0 (see https://doc.cgal.org/latest/Manual/usage.html):
```
cd /path/to/CGAL
cmake .
make install
```

Install pybind11 (see https://pybind11.readthedocs.io/en/stable/compiling.html):
```
git clone https://github.com/pybind/pybind11
cd pybind11
mkdir build
cd build
cmake ..
make install
```

Now you should be able to do:

```
cd /path/to/pycgalmini
cmake .
make
```

At this point, you should be able to run `demo.py`.
