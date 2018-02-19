# Fast Self-forced Inspirals

Software to rapidly compute inspirals trajectories and their associated waveforms for eccentric small mass-ratio inspirals into a Schwarzschild black hole. 
The computed inspirals and waveform include local self-force effects.

Details of the near-identity transformation method used can be found in: https://arxiv.org/abs/1802.05281

### Dependencies

The NIT inspiral code depends upon:

 - The GNU Scientific Library (https://www.gnu.org/software/gsl/)
 - FFTW (http://www.fftw.org/)
 - libconfig (https://hyperrealm.github.io/libconfig/)
 - Scons (http://scons.org/)

### Compile the code

Type `scons` in the  main directory to compile.

### Usage

Running `./NIT_inspiral` without any arguments will give a list of the possible arguments. The various options and how to use the code are given below.

Before computing a NIT inspiral you must first compute the the averaged forcing functions. This is performed in two steps:

1. ./NIT_inspiral -d (to decompose the self-force into Fourier modes)
2. ./NIT_inspiral -c (to compute the averaged forcing functions from the Fourier coefficients)

A NIT inspiral can now be computed with 

```
./NIT_inspiral -n p0 e0 q
```

where (p0, e0) are the initial semi-latus rectum and orbital eccentricity and q is the (small) mass ratio. The inspiral is computed from the initial parameters
until the onset of plunge near the separatrix. The NIT inspiral will be computed in milliseconds.

A Full self-forced inspiral can be computed with

```
./NIT_inspiral -f p0 e0 q
```

This inspiral will take seconds to hours to compute depending on the value of q.

To compute a waveform one must first compute the inspiral using the commands above. The parameters for the waveform (sampling rate etc) are specified in the
configuration file found in config/parameters.cfg. To compute the waveform associated with a NIT inspiral use

```
./NIT -w p0 e0 q -n
```

Similarly to compute the waveform associated with a full self-forced inspiral use:

```
./NIT -w p0 e0 q -f
```

### Licence

The code is licensed under the GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)