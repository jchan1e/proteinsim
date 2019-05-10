## Coarse Grained Protein Simulation

A simple simulation and visualization of basic coarse grained protein folding using C++ and OpenGL.  

Compile all components
```
make
make all
```

compile individual components
```
make sim
make simg
make vis
```

Generate random amino acid sequences
```
python generate.py num_aminos [num_repeats] > input.dat
```

Run simulations
```
./sim input.dat output.bin [num_frames]
./simg input.dat output.bin [num_frames]
```

Run visualizations
```
./vis output.bin
```

Visualization camera/playback controls
```
space - play/pause
,/. - slower/faster
[/] - reverse/forward
z/x - zoom in/out
arrows - rotate camera
0 - reset camera
```
