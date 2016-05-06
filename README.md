- [Contact Mode Score (CMS)](#sec-1)
  - [Installation](#sec-1-1)
  - [Usage](#sec-1-2)
  - [Examples](#sec-1-3)
  - [Dependencies](#sec-1-4)
- [eXtended Contact Mode Score (XCMS)](#sec-2)
  - [Installation Dependencies](#sec-2-1)
    - [Python 2.7 and pip](#sec-2-1-1)
    - [OpenBabel and Pybel](#sec-2-1-2)
    - [Pkcombu](#sec-2-1-3)
    - [APoc](#sec-2-1-4)
  - [Installation](#sec-2-2)
  - [Usage](#sec-2-3)

# Contact Mode Score (CMS)<a id="orgheadline5"></a>

Measure of protein-ligand complex (dis)similarity based on intermolecular
contacts

## Installation<a id="orgheadline1"></a>

Download the package, unpack the tarball, go to `contact-mode-score/src`, and
run `make`:

```sh
$ cd src
$ make cms -f makefile
```

If successful, the binary `cms` will show up in `contact-mode-score/bin`:

```sh
$ ls ../bin

cms
```

Add the path of `cms` to the `PATH` environment variable, like:

## Usage<a id="orgheadline2"></a>

```sh
$ cms -h
```

## Examples<a id="orgheadline3"></a>

```sh
$ cd data
$ sh cms.sh
```

## Dependencies<a id="orgheadline4"></a>

The codes have been tested on g++ 4.2.1, but any modern compilers supporting c++
standard libraries should work.

# eXtended Contact Mode Score (XCMS)<a id="orgheadline13"></a>

Compare two protein-ligand binding conformations even when they are
non-identical systems using **eXtended Contact Mode Score (XCMS)**

## Installation Dependencies<a id="orgheadline10"></a>

### Python 2.7 and pip<a id="orgheadline6"></a>

It is recommended to install [Anaconda](https://www.continuum.io/downloads) on your system to set up the python
working environment, which provides most of the python dependencies to run XCMS,
except for Pybel, whose installation shall be explained below.

### OpenBabel and Pybel<a id="orgheadline7"></a>

We use [OpenBabel](http://openbabel.org/wiki/Main_Page) and its python binding [pybel](https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html) to handle various input formats.

See [here](https://openbabel.org/docs/dev/Installation/install.html#compiling-open-babel) for the instructions how you can compile OpenBabel together with its
python binding.

Be aware, you must use the old [Eigen2 ](http://eigen.tuxfamily.org/index.php?title=Eigen2)instead of Eigen3.

### Pkcombu<a id="orgheadline8"></a>

We use [pkcombu](http://strcomp.protein.osaka-u.ac.jp/kcombu/doc/README_pkcombu.html) in the KCOMBU package to compare two ligands.

You may download the KCOMBU program [here](http://strcomp.protein.osaka-u.ac.jp/kcombu/download_src.html).

Please add the `pkcombu` into your system's `PATH` variable.

### APoc<a id="orgheadline9"></a>

We use Apoc to compare two binding pockets.

You may donwload Apoc program [here](http://cssb.biology.gatech.edu/APoc).

Please add the `apoc` into your system's `PATH` variable.

## Installation<a id="orgheadline11"></a>

At the command line:

```sh
$ pip install xcms
```

To install locally only for yourself:

```sh
$ pip install xcms --user
```

## Usage<a id="orgheadline12"></a>

To use it:

```sh
$ x-cms --help
Usage: x-cms [OPTIONS]

calculated extended contact mode score provided the query and template
protein-ligand structures

Options:
  --template-protein PATH  template protein file path, in PDB format
  --query-protein PATH     query protein file path, in PDB format
  --query-ligand PATH      query ligand file path, in SDF format
  --template-ligand PATH   template ligand file path, in SDF format
  --help                   Show this message and exit.
```
