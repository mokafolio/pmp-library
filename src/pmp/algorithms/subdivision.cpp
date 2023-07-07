// Copyright 2011-2022 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/subdivision.h"
#include "pmp/algorithms/differential_geometry.h"

#include <numbers>

namespace pmp {

void catmull_clark_subdivision(SurfaceMesh& mesh,
                               BoundaryHandling boundary_handling)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto vfeature_ = mesh.get_vertex_property<bool>("v:feature");
    auto efeature_ = mesh.get_edge_property<bool>("e:feature");

    // reserve memory
    size_t nv = mesh.n_vertices();
    size_t ne = mesh.n_edges();
    size_t nf = mesh.n_faces();
    mesh.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);

    // get properties
    auto vpoint = mesh.add_vertex_property<Point>("catmull:vpoint");
    auto epoint = mesh.add_edge_property<Point>("catmull:epoint");
    auto fpoint = mesh.add_face_property<Point>("catmull:fpoint");

    // compute face vertices
    for (auto f : mesh.faces())
    {
        fpoint[f] = centroid(mesh, f);
    }

    // compute edge vertices
    for (auto e : mesh.edges())
    {
        // boundary or feature edge?
        if (mesh.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] = 0.5f * (points_[mesh.vertex(e, 0)] +
                                points_[mesh.vertex(e, 1)]);
        }

        // interior edge
        else
        {
            Point p(0, 0, 0);
            p += points_[mesh.vertex(e, 0)];
            p += points_[mesh.vertex(e, 1)];
            p += fpoint[mesh.face(e, 0)];
            p += fpoint[mesh.face(e, 1)];
            p *= 0.25f;
            epoint[e] = p;
        }
    }

    // compute new positions for old vertices
    for (auto v : mesh.vertices())
    {
        // isolated vertex?
        if (mesh.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh.is_boundary(v))
        {
            if (boundary_handling == BoundaryHandling::Preserve)
            {
                vpoint[v] = points_[v];
            }
            else
            {
                auto h1 = mesh.halfedge(v);
                auto h0 = mesh.prev_halfedge(h1);

                Point p = points_[v];
                p *= 6.0;
                p += points_[mesh.to_vertex(h1)];
                p += points_[mesh.from_vertex(h0)];
                p *= 0.125;

                vpoint[v] = p;
            }
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh.halfedges(v))
            {
                if (efeature_[mesh.edge(h)])
                {
                    p += points_[mesh.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

            const Scalar k = mesh.valence(v);
            Point p(0, 0, 0);

            for (auto vv : mesh.vertices(v))
                p += points_[vv];

            for (auto f : mesh.faces(v))
                p += fpoint[f];

            p /= (k * k);

            p += ((k - 2.0f) / k) * points_[v];

            vpoint[v] = p;
        }
    }

    // assign new positions to old vertices
    for (auto v : mesh.vertices())
    {
        points_[v] = vpoint[v];
    }

    // split edges
    for (auto e : mesh.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh.insert_vertex(e, epoint[e]);
            auto v = mesh.to_vertex(h);
            auto e0 = mesh.edge(h);
            auto e1 = mesh.edge(mesh.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    for (auto f : mesh.faces())
    {
        auto h0 = mesh.halfedge(f);
        mesh.insert_edge(h0, mesh.next_halfedge(mesh.next_halfedge(h0)));

        auto h1 = mesh.next_halfedge(h0);
        mesh.insert_vertex(mesh.edge(h1), fpoint[f]);

        auto h = mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h1)));
        while (h != h0)
        {
            mesh.insert_edge(h1, h);
            h = mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h1)));
        }
    }

    // clean-up properties
    mesh.remove_vertex_property(vpoint);
    mesh.remove_edge_property(epoint);
    mesh.remove_face_property(fpoint);
}

void loop_subdivision(SurfaceMesh& mesh, BoundaryHandling boundary_handling)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto vfeature_ = mesh.get_vertex_property<bool>("v:feature");
    auto efeature_ = mesh.get_edge_property<bool>("e:feature");

    if (!mesh.is_triangle_mesh())
    {
        auto what = std::string{__func__} + ": Not a triangle mesh.";
        throw InvalidInputException(what);
    }

    // reserve memory
    size_t nv = mesh.n_vertices();
    size_t ne = mesh.n_edges();
    size_t nf = mesh.n_faces();
    mesh.reserve(nv + ne, 2 * ne + 3 * nf, 4 * nf);

    // add properties
    auto vpoint = mesh.add_vertex_property<Point>("loop:vpoint");
    auto epoint = mesh.add_edge_property<Point>("loop:epoint");

    // compute vertex positions
    for (auto v : mesh.vertices())
    {
        // isolated vertex?
        if (mesh.is_isolated(v))
        {
            vpoint[v] = points_[v];
        }

        // boundary vertex?
        else if (mesh.is_boundary(v))
        {
            if (boundary_handling == BoundaryHandling::Preserve)
            {
                vpoint[v] = points_[v];
            }
            else
            {
                auto h1 = mesh.halfedge(v);
                auto h0 = mesh.prev_halfedge(h1);

                Point p = points_[v];
                p *= 6.0;
                p += points_[mesh.to_vertex(h1)];
                p += points_[mesh.from_vertex(h0)];
                p *= 0.125;
                vpoint[v] = p;
            }
        }

        // interior feature vertex?
        else if (vfeature_ && vfeature_[v])
        {
            Point p = points_[v];
            p *= 6.0;
            int count(0);

            for (auto h : mesh.halfedges(v))
            {
                if (efeature_[mesh.edge(h)])
                {
                    p += points_[mesh.to_vertex(h)];
                    ++count;
                }
            }

            if (count == 2) // vertex is on feature edge
            {
                p *= 0.125;
                vpoint[v] = p;
            }
            else // keep fixed
            {
                vpoint[v] = points_[v];
            }
        }

        // interior vertex
        else
        {
            Point p(0, 0, 0);
            Scalar k(0);

            for (auto vv : mesh.vertices(v))
            {
                p += points_[vv];
                ++k;
            }
            p /= k;

            Scalar beta =
                (0.625 -
                 pow(0.375 + 0.25 * std::cos(2.0 * std::numbers::pi / k), 2.0));

            vpoint[v] = points_[v] * (Scalar)(1.0 - beta) + beta * p;
        }
    }

    // compute edge positions
    for (auto e : mesh.edges())
    {
        // boundary or feature edge?
        if (mesh.is_boundary(e) || (efeature_ && efeature_[e]))
        {
            epoint[e] =
                (points_[mesh.vertex(e, 0)] + points_[mesh.vertex(e, 1)]) *
                Scalar(0.5);
        }

        // interior edge
        else
        {
            auto h0 = mesh.halfedge(e, 0);
            auto h1 = mesh.halfedge(e, 1);
            Point p = points_[mesh.to_vertex(h0)];
            p += points_[mesh.to_vertex(h1)];
            p *= 3.0;
            p += points_[mesh.to_vertex(mesh.next_halfedge(h0))];
            p += points_[mesh.to_vertex(mesh.next_halfedge(h1))];
            p *= 0.125;
            epoint[e] = p;
        }
    }

    // set new vertex positions
    for (auto v : mesh.vertices())
    {
        points_[v] = vpoint[v];
    }

    // insert new vertices on edges
    for (auto e : mesh.edges())
    {
        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto h = mesh.insert_vertex(e, epoint[e]);
            auto v = mesh.to_vertex(h);
            auto e0 = mesh.edge(h);
            auto e1 = mesh.edge(mesh.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }

        // normal edge
        else
        {
            mesh.insert_vertex(e, epoint[e]);
        }
    }

    // split faces
    Halfedge h;
    for (auto f : mesh.faces())
    {
        h = mesh.halfedge(f);
        mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
        h = mesh.next_halfedge(h);
        mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
        h = mesh.next_halfedge(h);
        mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
    }

    // clean-up properties
    mesh.remove_vertex_property(vpoint);
    mesh.remove_edge_property(epoint);
}

void quad_tri_subdivision(SurfaceMesh& mesh, BoundaryHandling boundary_handling)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto uvs_ = mesh.get_halfedge_property<TexCoord>("h:tex");

    // split each edge evenly into two parts
    for (auto e : mesh.edges())
    {
        if (uvs_)
        {
            auto ha = mesh.halfedge(e, 0);
            auto hb = mesh.halfedge(e, 1);

            auto a_from = uvs_[mesh.prev_halfedge(ha)];
            auto a_to = uvs_[ha];
            auto b_from = uvs_[mesh.prev_halfedge(hb)];
            auto b_to = uvs_[hb];
            auto nhe =
                mesh.insert_vertex(e, 0.5f * (points_[mesh.vertex(e, 0)] +
                                              points_[mesh.vertex(e, 1)]));
            auto ohe = mesh.opposite_halfedge(nhe);
            auto auv = (a_from + a_to) * 0.5f;
            auto buv = (b_from + b_to) * 0.5f;

            uvs_[nhe] = buv;
            uvs_[mesh.next_halfedge(nhe)] = b_to;
            uvs_[mesh.prev_halfedge(ohe)] = auv;
            uvs_[ohe] = a_to;
        }
        else
        {
            mesh.insert_vertex(e, 0.5f * (points_[mesh.vertex(e, 0)] +
                                          points_[mesh.vertex(e, 1)]));
        }
    }

    // subdivide faces without repositioning
    for (auto f : mesh.faces())
    {
        size_t f_val = mesh.valence(f) / 2;
        if (f_val == 3)
        {
            // face was a triangle
            auto make_next_triangle = [&](pmp::Halfedge _start_halfedge) {
                pmp::Halfedge h0 = _start_halfedge;
                pmp::Halfedge h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
                auto nh = mesh.insert_edge(h0, h1);

                if (uvs_)
                {
                    uvs_[nh] = uvs_[h1];
                    uvs_[mesh.opposite_halfedge(nh)] = uvs_[h0];
                }
            };

            auto h0 = mesh.halfedge(f);
            make_next_triangle(h0);
            h0 = mesh.next_halfedge(h0);
            make_next_triangle(h0);
            h0 = mesh.next_halfedge(h0);
            make_next_triangle(h0);
        }
        else
        {
            // quadrangulate the rest
            auto h0 = mesh.halfedge(f);
            auto h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
            //NOTE: It's important to calculate the centroid before inserting the new edge
            auto cen = pmp::centroid(mesh, f);

            // calculate the centroid UV if necessary
            pmp::TexCoord uv_cen{0.0f, 0.0f};
            if (uvs_)
            {
                pmp::Scalar n{0};
                auto hedges = mesh.halfedges(f);
                for (auto he : hedges)
                {
                    uv_cen += uvs_[he];
                    ++n;
                }
                uv_cen /= n;
            }

            //insert the new edge and centroid vertex
            auto h2 = mesh.insert_edge(h0, h1);
            auto h3 = mesh.insert_vertex(mesh.edge(h2), cen);

            //set the new UV's for the resulting new half edges
            if (uvs_)
            {
                uvs_[h2] = uv_cen;
                uvs_[mesh.opposite_halfedge(h2)] = uvs_[h0];
                uvs_[h3] = uv_cen;
                uvs_[mesh.opposite_halfedge(h3)] = uvs_[h1];
            }

            auto h =
                mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h2)));
            while (h != h0)
            {
                auto nh = mesh.insert_edge(h2, h);

                if (uvs_)
                {
                    uvs_[nh] = uvs_[h];
                    uvs_[mesh.opposite_halfedge(nh)] = uvs_[h2];
                }

                h = mesh.next_halfedge(
                    mesh.next_halfedge(mesh.next_halfedge(h2)));
            }
        }
    }

    auto new_pos =
        mesh.add_vertex_property<Point>("quad_tri:new_position", Point(0));

    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
        {
            if (boundary_handling == BoundaryHandling::Preserve)
            {
                new_pos[v] = points_[v];
            }
            else
            {
                new_pos[v] = 0.5 * points_[v];

                // add neighboring vertices on boundary
                for (auto vv : mesh.vertices(v))
                {
                    if (mesh.is_boundary(vv))
                    {
                        new_pos[v] += 0.25 * points_[vv];
                    }
                }
            }
        }
        else
        {
            // count the number of faces and quads surrounding the vertex
            int n_faces = 0;
            int n_quads = 0;
            for (auto f : mesh.faces(v))
            {
                n_faces++;
                if (mesh.valence(f) == 4)
                    n_quads++;
            }

            if (n_quads == 0)
            {
                // vertex is surrounded only by triangles
                double a =
                    2.0 *
                    pow(3.0 / 8.0 +
                            (std::cos(2.0 * std::numbers::pi / n_faces) - 1.0) /
                                4.0,
                        2.0);
                double b = (1.0 - a) / n_faces;

                new_pos[v] = a * points_[v];
                for (auto vv : mesh.vertices(v))
                {
                    new_pos[v] += b * points_[vv];
                }
            }
            else if (n_quads == n_faces)
            {
                // vertex is surrounded only by quads
                double c = (n_faces - 3.0) / n_faces;
                double d = 2.0 / pow(n_faces, 2.0);
                double e = 1.0 / pow(n_faces, 2.0);

                new_pos[v] = c * points_[v];
                for (auto h : mesh.halfedges(v))
                {
                    new_pos[v] += d * points_[mesh.to_vertex(h)];
                    new_pos[v] +=
                        e * points_[mesh.to_vertex(mesh.next_halfedge(h))];
                }
            }
            else
            {
                // vertex is surrounded by triangles and quads
                double alpha = 1.0 / (1.0 + 0.5 * n_faces + 0.25 * n_quads);
                double beta = 0.5 * alpha;
                double gamma = 0.25 * alpha;

                new_pos[v] = alpha * points_[v];
                for (auto h : mesh.halfedges(v))
                {
                    new_pos[v] += beta * points_[mesh.to_vertex(h)];
                    if (mesh.valence(mesh.face(h)) == 4)
                    {
                        new_pos[v] +=
                            gamma *
                            points_[mesh.to_vertex(mesh.next_halfedge(h))];
                    }
                }
            }
        }
    }

    // apply new positions to the mesh
    for (auto v : mesh.vertices())
    {
        points_[v] = new_pos[v];
    }

    mesh.remove_vertex_property(new_pos);
}

void linear_subdivision(SurfaceMesh& mesh)
{
    auto points_ = mesh.vertex_property<Point>("v:point");

    // linear subdivision of edges
    for (auto e : mesh.edges())
    {
        mesh.insert_vertex(e, 0.5f * (points_[mesh.vertex(e, 0)] +
                                      points_[mesh.vertex(e, 1)]));
    }

    // subdivide faces
    for (auto f : mesh.faces())
    {
        size_t f_val = mesh.valence(f) / 2;

        if (f_val == 3) // triangle
        {
            Halfedge h0 = mesh.halfedge(f);
            Halfedge h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
            mesh.insert_edge(h0, h1);

            h0 = mesh.next_halfedge(h0);
            h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
            mesh.insert_edge(h0, h1);

            h0 = mesh.next_halfedge(h0);
            h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
            mesh.insert_edge(h0, h1);
        }
        else // quadrangulate other faces
        {
            Halfedge h0 = mesh.halfedge(f);
            Halfedge h1 = mesh.next_halfedge(mesh.next_halfedge(h0));

            // NOTE: It's important to calculate the centroid before inserting the new edge
            auto cen = centroid(mesh, f);
            h1 = mesh.insert_edge(h0, h1);
            mesh.insert_vertex(mesh.edge(h1), cen);

            auto h =
                mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h1)));
            while (h != h0)
            {
                mesh.insert_edge(h1, h);
                h = mesh.next_halfedge(
                    mesh.next_halfedge(mesh.next_halfedge(h1)));
            }
        }
    }
}

} // namespace pmp
