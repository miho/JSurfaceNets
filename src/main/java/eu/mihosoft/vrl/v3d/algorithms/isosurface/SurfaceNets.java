/*
 * Copyright 2015 Michael Hoffer <info@michaelhoffer.de>. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Michael Hoffer <info@michaelhoffer.de> "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Michael Hoffer <info@michaelhoffer.de> OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those of the
 * authors and should not be interpreted as representing official policies, either expressed
 * or implied, of Michael Hoffer <info@michaelhoffer.de>.
 */
package eu.mihosoft.vrl.v3d.algorithms.isosurface;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * SurfaceNets based isosurface extraction.
 *
 * Based on the slightly scary-looking, compact and efficient JavaScript
 * implementation by Mikola Lysenko. 
 *
 * In addition this implementation comes with simple CSG support (union,
 * difference and intersection). Meshes can be saved as *.obj-File.
 *
 * Based on: S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech
 * Report.
 *
 * The only difference to the method described in the paper is that this
 * algorithms chooses the average point of the edge-crossing points on the
 * surface as cell coordinate for the quad-mesh generation.
 * 
 * This works surprisingly well! For many cases this algorithm is a good
 * substitute for the marching cubes algorithm. This implementation can also
 * be used to create more intelligent dual methods by replacing the
 * "average point" code.
 *
 * @author Michael Hoffer &lt;info@michaelhoffer.de&gt;
 */
public class SurfaceNets {

    private static final int[] cube_edges = new int[24];
    private static final int[] edge_table = new int[256];
    private static int[] vertexBuffer = new int[4096];

    /**
     * Initializes the cube edge indices. This follows the idea of Paul Bourke
     * link: http://paulbourke.net/geometry/polygonise/
     */
    static void initCubeEdges() {
        int k = 0;
        for (int i = 0; i < 8; ++i) {
            for (int j = 1; j <= 4; j <<= 1) {
                int p = i ^ j;
                if (i <= p) {
                    cube_edges[k++] = i;
                    cube_edges[k++] = p;
                }
            }
        }
    }

    /**
     * Initializes the cube edge table. This follows the idea of Paul Bourke
     * link: http://paulbourke.net/geometry/polygonise/
     */
    static void initEdgeTable() {
        for (int i = 0; i < 256; ++i) {
            int em = 0;
            for (int j = 0; j < 24; j += 2) {
                boolean a = bool(i & (1 << cube_edges[j]));
                boolean b = bool(i & (1 << cube_edges[j + 1]));
                em |= a != b ? (1 << (j >> 1)) : 0;
            }
            edge_table[i] = em;
        }
    }

    static {
        initCubeEdges();
        initEdgeTable();
    }

    /**
     * Minimal mesh class.
     */
    public static final class Mesh {

        public final List<double[]> vertices;
        public final List<int[]> faces;

        public Mesh(List<double[]> vertices, List<int[]> faces) {
            this.vertices = vertices;
            this.faces = faces;
        }
    }

    /**
     * Runs the surface nets algorithm.
     *
     * @param data iso data
     * @param dims dimensions
     * @return mesh
     */
    public static Mesh run(double[] data, int[] dims) {

        // location (location[0]=x, location[1]=y, location[2]=z)
        final int[] location = new int[3];
        // layout for one-dimensional data array
        // we use this to reference vertex buffer
        final int[] R = {
            // x
            1,
            // y * width
            dims[0] + 1,
            // z * width * height
            (dims[0] + 1) * (dims[1] + 1)
        };
        // grid cell
        final double grid[] = new double[8];
        
        // TODO: is is the only mystery that is left
        int buf_no = 1;
        

        // Resize buffer if necessary 
        if (R[2] * 2 > vertexBuffer.length) {
            vertexBuffer = new int[R[2] * 2];
        }
        
        // we make some assumptions about the number of vertices and faces
        // to reduce GC overhead
        final List<double[]> vertices = new ArrayList<>(vertexBuffer.length/2);
        final List<int[]> faces = new ArrayList<>(vertices.size()/4);
        
        int n = 0;

        // March over the voxel grid
        for (location[2] = 0; location[2] < dims[2] - 1; ++location[2], n += dims[0], buf_no ^= 1 /*even or odd*/, R[2] = -R[2]) {

            // m is the pointer into the buffer we are going to use.  
            // This is slightly obtuse because javascript does not
            // have good support for packed data structures, so we must
            // use typed arrays :(
            // The contents of the buffer will be the indices of the
            // vertices on the previous x/y slice of the volume
            int m = 1 + (dims[0] + 1) * (1 + buf_no * (dims[1] + 1));

            for (location[1] = 0; location[1] < dims[1] - 1; ++location[1], ++n, m += 2) {
                for (location[0] = 0; location[0] < dims[0] - 1; ++location[0], ++n, ++m) {

                    // Read in 8 field values around this vertex
                    // and store them in an array
                    // Also calculate 8-bit mask, like in marching cubes,
                    // so we can speed up sign checks later
                    int mask = 0, g = 0, idx = n;
                    for (int k = 0; k < 2; ++k, idx += dims[0] * (dims[1] - 2)) {
                        for (int j = 0; j < 2; ++j, idx += dims[0] - 2) {
                            for (int i = 0; i < 2; ++i, ++g, ++idx) {
                                final double p = data[idx];
                                grid[g] = p;
                                mask |= (p < 0) ? (1 << g) : 0;
                            }
                        }
                    }

                    // Check for early termination
                    // if cell does not intersect boundary
                    if (mask == 0 || mask == 0xff) {
                        continue;
                    }

                    // Sum up edge intersections
                    int edge_mask = edge_table[mask];
                    double[] v = {0.0, 0.0, 0.0};
                    int e_count = 0;

                    // For every edge of the cube...
                    for (int i = 0; i < 12; ++i) {

                        // Use edge mask to check if it is crossed
                        if (!bool((edge_mask & (1 << i)))) {
                            continue;
                        }

                        // If it did, increment number of edge crossings
                        ++e_count;

                        // Now find the point of intersection
                        int firstEdgeIndex = i << 1;
                        int secondEdgeIndex = firstEdgeIndex + 1;
                        // Unpack vertices
                        int e0 = cube_edges[firstEdgeIndex];
                        int e1 = cube_edges[secondEdgeIndex];
                        // Unpack grid values
                        double g0 = grid[e0];
                        double g1 = grid[e1];

                        // Compute point of intersection (linear interpolation)
                        double t = g0 - g1;
                        if (Math.abs(t) > 1e-6) {
                            t = g0 / t;
                        } else {
                            continue;
                        }

                        // Interpolate vertices and add up intersections
                        // (this can be done without multiplying)
                        for (int j = 0; j < 3; j++) {
                            int k = 1 << j; // (1,2,4)
                            int a = e0 & k;
                            int b = e1 & k;
                            if (a != b) {
                                v[j] += bool(a) ? 1.0 - t : t;
                            } else {
                                v[j] += bool(a) ? 1.0 : 0;
                            }
                        }
                    }

                    // Now we just average the edge intersections
                    // and add them to coordinate
                    double s = 1.0 / e_count;
                    for (int i = 0; i < 3; ++i) {
                        v[i] = location[i] + s * v[i];
                    }

                    // Add vertex to buffer, store pointer to
                    // vertex index in buffer
                    vertexBuffer[m] = vertices.size();
                    vertices.add(v);

                    // Now we need to add faces together, to do this we just
                    // loop over 3 basis components
                    for (int i = 0; i < 3; ++i) {
                        
                        // The first three entries of the edge_mask
                        // count the crossings along the edge
                        if (!bool(edge_mask & (1 << i))) {
                            continue;
                        }

                        // i = axes we are point along.
                        // iu, iv = orthogonal axes
                        int iu = (i + 1) % 3;
                        int iv = (i + 2) % 3;

                        // If we are on a boundary, skip it
                        if (location[iu] == 0 || location[iv] == 0) {
                            continue;
                        }

                        // Otherwise, look up adjacent edges in buffer
                        int du = R[iu];
                        int dv = R[iv];

                        // finally, the indices for the 4 vertices
                        // that define the face
                        final int indexM = vertexBuffer[m];
                        final int indexMMinusDU = vertexBuffer[m - du];
                        final int indexMMinusDV = vertexBuffer[m - dv];
                        final int indexMMinusDUMinusDV = vertexBuffer[m - du - dv];

                        // Remember to flip orientation depending on the sign
                        // of the corner.
                        if (bool(mask & 1)) {
                            faces.add(new int[]{
                                indexM,
                                indexMMinusDU,
                                indexMMinusDUMinusDV,
                                indexMMinusDV
                            });
                        } else {
                            faces.add(new int[]{
                                indexM,
                                indexMMinusDV,
                                indexMMinusDUMinusDV,
                                indexMMinusDU
                            });
                        }
                    }
                } // end x
            } // end y
        } // end z

        //All done!  Return the result
        return new Mesh(vertices, faces);
    }

    /**
     * Saves the specified mesh as obj file.
     *
     * @param p destination path
     * @param mesh mesh that shall be saved
     */
    public static void saveAsObj(Path p, Mesh mesh) {

        List<double[]> vertices = mesh.vertices;
        List<int[]> faces = mesh.faces;

        StringBuilder txtFile = new StringBuilder();

        txtFile.append("# JIsosurface OBJ File:\n");
        txtFile.append("o Geometry\n\n");

        for (int i = 0; i < vertices.size(); i++) {
            double[] vertex = vertices.get(i);
            txtFile.
                    append("v ").append(vertex[0]).
                    append(" ").append(vertex[1]).append(" ").
                    append(vertex[2]).append("\n");
        }

        for (int i = 0; i < faces.size(); i++) {
            int[] vertIdx = faces.get(i);
            txtFile.append("f ").append(vertIdx[0] + 1).
                    append(" ").append(vertIdx[1] + 1).append(" ").
                    append(vertIdx[2] + 1).append(" ").
                    append(vertIdx[3] + 1).append("\n");
        }

        try {
            Files.write(p, txtFile.toString().getBytes("UTF-8"));
        } catch (IOException ex) {
            Logger.getLogger(SurfaceNets.class.getName()).
                    log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Converts int to bool.
     *
     * @param i integer to convert
     * @return {@code true} if i > 0; {@code false} otherwise
     */
    private static boolean bool(int i) {
        return i > 0;
    }

    /**
     * Test method (tests union, difference and intersection).
     *
     * @param args (are ignored)
     */
    public static void main(String[] args) {

        // stack size
        int sizeX = 256;
        int sizeY = 256;
        int sizeZ = 256;

        // cube CSG
        CSG cube = new CubeCSG(
                new Vector3_double(
                        sizeX / 4 + 1,
                        sizeY / 4 + 1,
                        sizeZ / 4 + 1),
                new Vector3_double(
                        sizeX / 8,
                        sizeX / 8,
                        sizeX / 8));

        // sphere CSG
        CSG sphere = new SphereCSG(
                new Vector3_double(
                        sizeX / 2.5f,
                        sizeX / 2.5f,
                        sizeX / 2.5f),
                sizeX / 8);

        // cube minus sphere
        SampleFunctionData cubeMinusSphere = new SampleFunctionData(
                cube.difference(sphere),
                sizeX, sizeY, sizeZ);

        SurfaceNets.Mesh mesh1
                = SurfaceNets.run(cubeMinusSphere.getRawData(),
                        new int[]{sizeX, sizeY, sizeZ});

        SurfaceNets.saveAsObj(Paths.get("cube-minus-sphere.obj"), mesh1);

        // cube plus sphere
        SampleFunctionData cubePlusSphere = new SampleFunctionData(
                cube.union(sphere),
                sizeX, sizeY, sizeZ);

        SurfaceNets.Mesh mesh2
                = SurfaceNets.run(cubePlusSphere.getRawData(),
                        new int[]{sizeX, sizeY, sizeZ});

        SurfaceNets.saveAsObj(Paths.get("cube-plus-sphere.obj"), mesh2);

        // intersection of cube and sphere
        SampleFunctionData intersectionBetweenCubeAndSphere = new SampleFunctionData(
                cube.intersection(sphere),
                sizeX, sizeY, sizeZ);

        SurfaceNets.Mesh mesh3
                = SurfaceNets.run(intersectionBetweenCubeAndSphere.getRawData(),
                        new int[]{sizeX, sizeY, sizeZ});

        SurfaceNets.saveAsObj(
                Paths.get("instersection-of-cube-and-sphere.obj"), mesh3);
    } // end main()

    static interface CSG {

        default CSG union(CSG other) {
            CSG result = (x, y, z)
                    -> Math.min(CSG.this.get(x, y, z), other.get(x, y, z));
            return result;
        }

        default CSG difference(CSG other) {
            CSG result = (x, y, z) -> {
                double sv1 = CSG.this.get(x, y, z);
                double sv2 = other.get(x, y, z);

                return Math.max(-sv2, sv1);
            };
            return result;
        }

        default CSG intersection(CSG other) {
            CSG result = (x, y, z)
                    -> Math.max(CSG.this.get(x, y, z), other.get(x, y, z));
            return result;
        }

        public double get(double x, double y, double z);
    }

    /**
     * Sphere CSG.
     */
    static class SphereCSG implements CSG {

        private float radius;
        private Vector3_double center;

        /**
         * Constructor.
         *
         * @param center sphere center
         * @param radius sphere radius
         */
        public SphereCSG(Vector3_double center, float radius) {
            this.radius = radius;
            this.center = center;
        }

        @Override
        public double get(double x, double y, double z) {
            double diffX = center.x - x;
            double diffY = center.y - y;
            double diffZ = center.z - z;

            diffX /= radius;
            diffY /= radius;
            diffZ /= radius;

            double sv1 = radius * (Math.sqrt(
                    diffX * diffX
                    + diffY * diffY
                    + diffZ * diffZ) - 1);

            return sv1;
        }
    }

    /**
     * Cube CSG.
     */
    static class CubeCSG implements CSG {

        private Vector3_double size;
        private Vector3_double center;

        /**
         * Constructor.
         *
         * @param center cube center
         * @param size cube dimensions (w, h, d)
         */
        public CubeCSG(Vector3_double center, Vector3_double size) {
            this.size = size;
            this.center = center;
        }

        @Override
        public double get(double x, double y, double z) {

            Vector3_double dist = center.minus(
                    new Vector3_double((float) x, (float) y, (float) z));
            dist = new Vector3_double(
                    Math.abs(dist.x),
                    Math.abs(dist.y),
                    Math.abs(dist.z)).minus(size);

            return Math.max(dist.x, Math.max(dist.y, dist.z));
        }

    }

    /**
     * Data interface used for isosurfe extraction.
     */
    interface IData3_double {

        /**
         * Returns the width of the data stack.
         *
         * @return width of the data stack
         */
        public int getWidth();

        /**
         * Returns the height of the data stack.
         *
         * @return height of the data stack
         */
        public int getHeight();

        /**
         * Returns the depth of the data stack.
         *
         * @return depth of the data stack
         */
        public int getDepth();

        /**
         * Returns the iso value for the specified location (x,y,z).
         *
         * @param x x coordinate
         * @param y y coordinate
         * @param z z coordinate
         * @return iso value for the specified location
         */
        public double get(double x, double y, double z);

        /**
         * Returns the iso value for the specified location (x,y,z).
         *
         * @param loc location
         * @return iso value for the specified location
         */
        public double get(Vector3_double loc);
    }

    /**
     * Samples the specified csg function.
     */
    static class SampleFunctionData implements IData3_double {

        private final CSG csg;
        private int w, h, d;

        /**
         * Constructor.
         *
         * @param csg csg function that shall be sampled
         */
        public SampleFunctionData(CSG csg) {
            this.csg = csg;
        }

        /**
         * Constructor.
         *
         * @param csg csg function that shall be sampled
         * @param w stack width
         * @param h stack height
         * @param d stack depth
         */
        SampleFunctionData(CSG csg, int w, int h, int d) {
            this.csg = csg;
            this.w = w;
            this.h = h;
            this.d = d;
        }

        @Override
        public int getWidth() {
            return w;
        }

        @Override
        public int getHeight() {
            return h;
        }

        @Override
        public int getDepth() {
            return d;
        }

        @Override
        public double get(double x, double y, double z) {
            if (x < 1 || y < 1 || z < 1
                    || x >= getWidth() - 1
                    || y >= getHeight() - 1
                    || z >= getDepth() - 1) {
                return Double.MAX_VALUE;
            }

            return csg.get(x, y, z);
        }

        @Override
        public double get(Vector3_double pos) {
            return get(pos.x, pos.y, pos.z);
        }

        /**
         * Returns the stack data as raw data (double[])
         *
         * @return the stack data as raw data (double[])
         */
        public double[] getRawData() {
            double[] raw = new double[w * h * d];

            for (int z = 0; z < d; z++) {
                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < w; x++) {
                        raw[z * h * w + y * w + x] = get(x, y, z);
                    }
                }
            }

            return raw;
        }
    };

    /**
     * Minimal Vector3 implementation (not optimized for speed).
     */
    static class Vector3_double {

        public double x;
        public double y;
        public double z;

        /**
         * Constructor. Creates the specified vector (x,y,z).
         *
         * @param x x coordinate
         * @param y y coordinate
         * @param z z coordinate
         */
        public Vector3_double(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        /**
         * Returns v1+v2.
         *
         * @param v v2
         * @return v1+v2 as new instance
         */
        public Vector3_double plus(Vector3_double v) {
            return new Vector3_double(x + v.x, y + v.y, z + v.z);
        }

        /**
         * Returns v1-v2.
         *
         * @param v v2
         * @return v1-v2 as new instance
         */
        public Vector3_double minus(Vector3_double v) {
            return new Vector3_double(x - v.x, y - v.y, z - v.z);
        }

        /**
         * Sets the coordinates to coordinates of the specified vector.
         *
         * @param vec vector
         * @return this vector
         */
        public Vector3_double set(Vector3_double vec) {
            this.x = vec.x;
            this.y = vec.y;
            this.z = vec.z;
            return this;
        }

        @Override
        public String toString() {
            return "[" + x + ", " + y + ", " + z + "]";
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final Vector3_double other = (Vector3_double) obj;
            if (this.x != other.x) {
                return false;
            }
            if (this.y != other.y) {
                return false;
            }
            if (this.z != other.z) {
                return false;
            }

            return true;
        }
    } // end Vector3
} // end SurfaceNets
