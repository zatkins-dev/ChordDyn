## ChordDyn

ChordDyn is a software library and collection of scripts written in Python 3 and Julia to generate novel chord progressions using chaotic trajectories.

### Installation

#### Requirements
- Julia >= 1.8.0
- Python >= 3.8

#### Installing
First, clone this repository:
```sh
git clone https://github.com/zatkins-dev/ChordDyn.git --recurse-submodules
cd ChordDyn
```

Then, install as a development Julia package:
```sh
julia -e 'import Pkg; Pkg.activate("."); Pkg.instantiate()'
```

Install python deps:
```sh
pip install --user -r requirements.txt
pip install -e ./submodules/TIVlib
```
