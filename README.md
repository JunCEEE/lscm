# LSCM
Large Single Crystal Maker

## Installation

Compile

```bash
$ gcc LSCM.c -std=c99 -o LSCM
```

## Usage

```
Usage: ./LSCM nx ny nz a [structure] -o output
Example: ./LSCM 50 50 50 3.615 -f -o singleCu.custom
    -s, Spherical shape 
    -f, FCC 
    -b, BCC 
    -n, Nacl 
    -o, outputfile
    -h, print help
```

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

Copyright (C) 2017 [Juncheng E](mailto:jce@pims.ac.cn)
