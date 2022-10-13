### Sensor Placement for 3D Coverage
---
This work provides solutions for sensor locations to have 3D Coverage (e.g. location for UV lights for sanitization in a bus, train or ICU). 

Three Models were studied,
1. Bare Visibility
2. Maximum Quality
3. Cumulative Quality

To run the program, Gurobi and Point Cloud Library are needed. First replace the Gurobi related file path in CMakeList as the actual installation path. Then,
```
$ mkdir build && cd build
$ cmake ../
$ make
$ mv ../icu_model/uniform_out.pcd . # or a pcd file and rename it into uniform_out.pcd
$ ./3d_coverage [parameter1 parameter2 ...] 
```

