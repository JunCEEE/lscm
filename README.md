# LSCM
Large Single Crystal Maker
v1.2

## Installation

Compile

```bash
$ make
```

## Usage

```
Usage: ./LSCM  [-p c]  [-f] [-b] [-d] [-n] [-a lattice_constant] [-r file.data] -x nx -y ny -z nz [--cfg] [--custom] -o output
Example: ./LSCM -f -a 3.615 -x 50 -y 50 -z 50  -o singleCu.custom
Example: ./LSCM -p 5.21033 -a 3.20927 -x 30 -y 50 -z 30 --cfg -o singleMg.cfg
    -s, Spherical shape 
    -f, FCC 
    -b, BCC 
    -d, Diamond 
    -n, Nacl 
    -r, Read and replicate 
    -o, outputfile
    -h, print help
	--cfg, output in .cfg format
	--custom, output in LAMMPS  .custom format
```

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

Copyright (C) 2017 [Juncheng E](mailto:jce@pims.ac.cn)
