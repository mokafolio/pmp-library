// Copyright 2011-2022 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/subdivision.h"
#include "pmp/algorithms/differential_geometry.h"

namespace pmp {

namespace {
void split_edge_uvs(SurfaceMesh& mesh, HalfedgeProperty<TexCoord>& orig_uvs,
                    HalfedgeProperty<TexCoord>& new_uvs, Edge e, Scalar factor)
{
    auto h0 = mesh.halfedge(e, 0);
    auto h1 = mesh.halfedge(e, 1);
    auto a_from = orig_uvs[mesh.prev_halfedge(h0)];
    auto a_to = orig_uvs[h0];
    auto b_from = orig_uvs[mesh.prev_halfedge(h1)];
    auto b_to = orig_uvs[h1];
    new_uvs[h0] = (a_from + a_to) * Scalar(factor);
    new_uvs[h1] = (b_from + b_to) * Scalar(factor);
}

TexCoord face_uv_centroid(SurfaceMesh& mesh, HalfedgeProperty<TexCoord>& uvs,
                          Face f)
{
    TexCoord uv_cen{0.0, 0.0};
    Scalar n{0};
    for (auto he : mesh.halfedges(f))
    {
        uv_cen += uvs[he];
        ++n;
    }
    uv_cen /= n;
    return uv_cen;
}

void triangulate_subdivided_face(SurfaceMesh& mesh,
                                 HalfedgeProperty<TexCoord>& uvs, Face f)
{
    auto make_next_triangle = [&](Halfedge he) {
        Halfedge h0 = he;
        Halfedge h1 = mesh.next_halfedge(mesh.next_halfedge(h0));
        auto nh = mesh.insert_edge(h0, h1);

        if (uvs)
        {
            uvs[nh] = uvs[h1];
            uvs[mesh.opposite_halfedge(nh)] = uvs[h0];
        }
    };

    auto h0 = mesh.halfedge(f);
    make_next_triangle(h0);
    h0 = mesh.next_halfedge(h0);
    make_next_triangle(h0);
    h0 = mesh.next_halfedge(h0);
    make_next_triangle(h0);
}

void quadrangulate_subdivided_face(SurfaceMesh& mesh,
                                   HalfedgeProperty<TexCoord>& uvs, Point ct,
                                   TexCoord uv_ct, Face f)
{
    auto h0 = mesh.halfedge(f);
    auto other_corner_he = mesh.next_halfedge(mesh.next_halfedge(h0));
    mesh.insert_edge(h0, mesh.next_halfedge(mesh.next_halfedge(h0)));

    auto h1 = mesh.next_halfedge(h0);
    auto h2 = mesh.insert_vertex(mesh.edge(h1), ct);

    //set the new UV's for the resulting new half edges
    if (uvs)
    {
        auto h3 = mesh.opposite_halfedge(h2);
        uvs[h2] = uv_ct;
        uvs[h3] = uvs[other_corner_he];
        uvs[mesh.prev_halfedge(h3)] = uv_ct;
        uvs[mesh.next_halfedge(h2)] = uvs[h0];
    }

    auto h = mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h1)));
    while (h != h0)
    {
        auto he = mesh.insert_edge(h1, h);
        if (uvs)
        {
            uvs[he] = uvs[h];
            uvs[mesh.opposite_halfedge(he)] = uvs[h1];
        }
        h = mesh.next_halfedge(mesh.next_halfedge(mesh.next_halfedge(h1)));
    }
}

Halfedge insert_vertex_and_patch_uvs(SurfaceMesh& mesh,
                                     HalfedgeProperty<TexCoord>& uvs,
                                     const HalfedgeProperty<TexCoord>& new_uvs,
                                     Edge e, const Point& v)
{
    if (uvs)
    {
        auto h0 = mesh.halfedge(e, 0);
        auto h1 = mesh.halfedge(e, 1);

        //need to cache the original he uvs before overwriting them
        auto tmpa_old = uvs[h1];
        auto tmpb_old = uvs[h0];

        auto h = mesh.insert_vertex(e, v);

        uvs[h] = new_uvs[h1];
        uvs[mesh.prev_halfedge(mesh.opposite_halfedge(h))] = new_uvs[h0];
        uvs[mesh.next_halfedge(h)] = tmpa_old;
        uvs[mesh.opposite_halfedge(h)] = tmpb_old;
        return h;
    }

    return mesh.insert_vertex(e, v);
}
} // namespace

void catmull_clark_subdivision(SurfaceMesh& mesh)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto vfeature_ = mesh.get_vertex_property<bool>("v:feature");
    auto efeature_ = mesh.get_edge_property<bool>("e:feature");
    auto uvs_ = mesh.get_halfedge_property<TexCoord>("h:tex");

    // reserve memory
    size_t nv = mesh.n_vertices();
    size_t ne = mesh.n_edges();
    size_t nf = mesh.n_faces();
    mesh.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);

    // get properties
    auto vpoint = mesh.add_vertex_property<Point>("catmull:vpoint");
    auto epoint = mesh.add_edge_property<Point>("catmull:epoint");
    auto fpoint = mesh.add_face_property<Point>("catmull:fpoint");

    HalfedgeProperty<TexCoord> new_uvs;
    FaceProperty<TexCoord> fuvs;
    if (uvs_)
    {
        new_uvs = mesh.add_halfedge_property<TexCoord>("catmull:nuvs");
        fuvs = mesh.add_face_property<TexCoord>("catmull:fuvs");
    }

    // compute face vertices
    for (auto f : mesh.faces())
    {
        fpoint[f] = centroid(mesh, f);

        // calculate the centroid UV if necessary
        if (uvs_)
        {
            fuvs[f] = face_uv_centroid(mesh, uvs_, f);
        }
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

        if (uvs_)
        {
            //TODO: Is there a better factor than 0.5 based on how loop distributes the edge points for the interior edge case?
            // Just splitting it in the middle seems to be what a lot of 3D tools such as Blender appear to be doing though.
            split_edge_uvs(mesh, uvs_, new_uvs, e, 0.5);
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
            auto h1 = mesh.halfedge(v);
            auto h0 = mesh.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh.to_vertex(h1)];
            p += points_[mesh.from_vertex(h0)];
            p *= 0.125;

            vpoint[v] = p;
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
        auto h = insert_vertex_and_patch_uvs(mesh, uvs_, new_uvs, e, epoint[e]);

        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto v = mesh.to_vertex(h);
            auto e0 = mesh.edge(h);
            auto e1 = mesh.edge(mesh.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }
    }

    // split faces
    for (auto f : mesh.faces())
    {
        quadrangulate_subdivided_face(mesh, uvs_, fpoint[f],
                                      fuvs ? fuvs[f] : TexCoord{}, f);
    }

    // clean-up properties
    mesh.remove_vertex_property(vpoint);
    mesh.remove_edge_property(epoint);
    mesh.remove_face_property(fpoint);
    if (uvs_)
    {
        mesh.remove_face_property(fuvs);
        mesh.remove_halfedge_property(new_uvs);
    }
}

void loop_subdivision(SurfaceMesh& mesh)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto vfeature_ = mesh.get_vertex_property<bool>("v:feature");
    auto efeature_ = mesh.get_edge_property<bool>("e:feature");
    auto uvs_ = mesh.get_halfedge_property<TexCoord>("h:tex");

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
            auto h1 = mesh.halfedge(v);
            auto h0 = mesh.prev_halfedge(h1);

            Point p = points_[v];
            p *= 6.0;
            p += points_[mesh.to_vertex(h1)];
            p += points_[mesh.from_vertex(h0)];
            p *= 0.125;
            vpoint[v] = p;
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
                (0.625 - pow(0.375 + 0.25 * std::cos(2.0 * M_PI / k), 2.0));

            vpoint[v] = points_[v] * (Scalar)(1.0 - beta) + beta * p;
        }
    }

    HalfedgeProperty<TexCoord> new_uvs;
    if (uvs_)
    {
        new_uvs = mesh.add_halfedge_property<TexCoord>("loop:nuvs");
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

        if (uvs_)
        {
            //TODO: Is there a better factor than 0.5 based on how loop distributes the edge points for the interior edge case?
            // Just splitting it in the middle seems to be what a lot of 3D tools such as Blender appear to be doing though.
            split_edge_uvs(mesh, uvs_, new_uvs, e, 0.5);
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
        auto h = insert_vertex_and_patch_uvs(mesh, uvs_, new_uvs, e, epoint[e]);

        // feature edge?
        if (efeature_ && efeature_[e])
        {
            auto v = mesh.to_vertex(h);
            auto e0 = mesh.edge(h);
            auto e1 = mesh.edge(mesh.next_halfedge(h));

            vfeature_[v] = true;
            efeature_[e0] = true;
            efeature_[e1] = true;
        }
    }

    // split faces
    Halfedge h;
    for (auto f : mesh.faces())
    {
        triangulate_subdivided_face(mesh, uvs_, f);
    }

    // clean-up properties
    mesh.remove_vertex_property(vpoint);
    mesh.remove_edge_property(epoint);
    if (new_uvs)
    {
        mesh.remove_halfedge_property(new_uvs);
    }
}

void quad_tri_subdivision(SurfaceMesh& mesh)
{
    auto points_ = mesh.vertex_property<Point>("v:point");
    auto uvs_ = mesh.get_halfedge_property<TexCoord>("h:tex");

    HalfedgeProperty<TexCoord> new_uvs;
    if (uvs_)
    {
        new_uvs = mesh.add_halfedge_property<TexCoord>("quad_tri:new_uvs");
    }

    // split each edge evenly into two parts
    for (auto e : mesh.edges())
    {
        if (uvs_)
        {
            split_edge_uvs(mesh, uvs_, new_uvs, e, 0.5);
        }

        insert_vertex_and_patch_uvs(
            mesh, uvs_, new_uvs, e,
            0.5f * (points_[mesh.vertex(e, 0)] + points_[mesh.vertex(e, 1)]));
    }

    // subdivide faces without repositioning
    for (auto f : mesh.faces())
    {
        size_t f_val = mesh.valence(f) / 2;
        if (f_val == 3)
        {
            // face was a triangle
            triangulate_subdivided_face(mesh, uvs_, f);
        }
        else
        {
            // quadrangulate the rest

            //face centroid
            auto cen = pmp::centroid(mesh, f);

            // calculate the centroid UV if necessary
            TexCoord uv_cen;
            if (uvs_)
            {
                uv_cen = face_uv_centroid(mesh, uvs_, f);
            }
            quadrangulate_subdivided_face(mesh, uvs_, cen, uv_cen, f);
        }
    }

    auto new_pos =
        mesh.add_vertex_property<Point>("quad_tri:new_position", Point(0));

    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
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
                    2.0 * pow(3.0 / 8.0 +
                                  (std::cos(2.0 * M_PI / n_faces) - 1.0) / 4.0,
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
    if (uvs_)
    {
        mesh.remove_halfedge_property(new_uvs);
    }
}

} // namespace pmp
