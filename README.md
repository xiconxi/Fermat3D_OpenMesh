# Fermat3D
An implementation for Coverage Path Planning(CPP) problem. It can fill any watertight manifold tri mesh with a connected fermat spirals. Simpily, one inward and one outward only.

Fill a mesh with one line only(网格一笔画).
## Demo Time
![](./data/fermatcpp.gif)

## Dependencies

[libigl](https://libigl.github.io/): Heat Geodesic Distance Scalar Feild Generator 

[OpenMesh](https://www.openmesh.org/): iso-countour cut the watertight manifold mesh by the scalar field.

## Algorithm
```
1. Calc heat geodesic distance field.
2. Gen Iso-countour from the mesh and corresponding scalar field.
3. Construct the spaning tree
4. Convert iso-countour into connected fermat spirals.
5. Connect all child connnected fermat spirals to the parent spiral.
```
## References
1. [Connected Fermat Spirals](https://homes.cs.washington.edu/~haisen/CFS/index.html)
2. [Coverage Path Planning for General Terrain Surfaces](http://www.mae.cuhk.edu.hk/~cwang/pubs/ICRA19FermatCPP.pdf)

## TODO
1. curves optimizations.
