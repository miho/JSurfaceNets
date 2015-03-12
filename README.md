# JSurfaceNets

Isosurface extraction & CSG based on SurfaceNets

Based on the slightly scary-looking, compact and efficient JavaScript implementation by Mikola Lysenko which is based on: S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report. The only difference to the method described in the paper is that this algorithms chooses the average point of the edge-crossing points on the surface as cell coordinate for the quad-mesh generation.

This works surprisingly well! For many cases this algorithm is a good substitute for the marching cubes algorithm. This implementation can also be used to create more intelligent dual methods by replacing the "average point" code.

In addition this implementation comes with simple CSG support (union, difference and intersection). Meshes can be saved as *.obj-File.


 



