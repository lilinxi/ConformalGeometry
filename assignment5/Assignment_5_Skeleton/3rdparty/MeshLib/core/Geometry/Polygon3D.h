#ifndef _MESHLIB_POLYGON_3D_H_
#define _MESHLIB_POLYGON_3D_H_

/*!
 *      \file Polygon3D.h
 *      \brief planar polygon
 *	    \author David Gu
 *      \date 05/03/2020
 *
 */
#include <vector>

#include "Point.h"
#include "Polygon2D.h"
#include "Ray3D.h"
#include "Segment3D.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace MeshLib
{

/*!
 *	\brief spacial polygon
 *
 */
class CPolygon3D
{
  public:
    /*!
     *	CPolygon constructor
     */
    CPolygon3D(){};
    ~CPolygon3D(){};

    /*! add one segment in sequencial order */
    void add(CSegment3D& seg) { m_edges.push_back(seg); };

    /*! polygon area */
    double area();

    /*! clip by the unit disk*/
    int finite_clip(CPolygon3D& inner_poly);

    /*! the edges*/
    std::vector<CSegment3D>& edges() { return m_edges; };

    /*! the polygon is an infinite cell, clipped by the unit circle */
    int infinite_clip(CPolygon3D& inner_poly);

    /*! centeter of mass */
    CPoint mass_center();

    /*! projection to 2D*/
    CPolygon2D project();

  protected:
    std::vector<CSegment3D> m_edges;
};

/*! the total area of a planar polygon */
inline double CPolygon3D::area()
{
    double s = 0;
    for (size_t i = 1; i < m_edges.size() - 1; i++)
    {
        CPoint p = m_edges[i + 0].start() - m_edges[0].start();
        CPoint q = m_edges[i + 1].start() - m_edges[0].start();
        CPoint n = p ^ q;
        s += n.norm();
    }
    return s / 2.0;
};

/*! the finite polygon is clipped by the unit disk */
inline int CPolygon3D::finite_clip(CPolygon3D& inner_poly)
{
    double epsilon = 1e-30;

    CCircle circ(CPoint2(0, 0), 1);
    std::vector<int> on_circle_idx;

    std::vector<CPoint> inner_pts;
    std::vector<CPoint> outer_pts;

    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment3D& seg = m_edges[i];
        CPoint& pt = seg.start();

        double R = pt[0] * pt[0] + pt[1] * pt[1];

        if (R < 1 - epsilon)
        {
            inner_pts.push_back(pt);
        }
        if (R > 1 + epsilon)
        {
            outer_pts.push_back(pt);
        }
        if (R <= 1 + epsilon && R >= 1 - epsilon)
        {
            inner_pts.push_back(pt);
            outer_pts.push_back(pt);
        }
        std::vector<CPoint> intersection_points;
        bool cross = seg.intersect(circ, intersection_points);
        if (!cross)
            continue;

        for (size_t j = 0; j < intersection_points.size(); j++)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points[j]);
            outer_pts.push_back(intersection_points[j]);
        }
    }

    // completely outside
    if (inner_pts.empty())
        return -1;

    // completely inside
    if (outer_pts.empty())
        return +1;

    int n = (int) inner_pts.size();
    int on_circle_indx = 0;

    for (int i = 0; i < (int) on_circle_idx.size(); i++)
    {
        if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
        {
            on_circle_indx = on_circle_idx[i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        CSegment3D seg(inner_pts[(i + on_circle_indx) % n], inner_pts[(i + 1 + on_circle_indx) % n]);
        inner_poly.add(seg);
    }

    return 0;
}

/*! the finite polygon is clipped by the unit disk */
inline int CPolygon3D::infinite_clip(CPolygon3D& inner_poly)
{
    double epsilon = 1e-30;

    CCircle circ(CPoint2(0, 0), 1);
    std::vector<int> on_circle_idx;

    std::vector<CPoint> inner_pts;
    std::vector<CPoint> outer_pts;

    {
        CSegment3D& sg = m_edges.front();
        CPoint st = sg.start();
        CPoint dr = sg.end() - sg.start();
        double n = dr.norm();

        CRay3D ray(st, dr);
        std::vector<CPoint> intersection_points;
        bool cross = ray.intersect(circ, intersection_points);

        if (!intersection_points.empty())
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points.front());
            outer_pts.push_back(intersection_points.front());
        }
    }

    for (size_t i = 1; i < m_edges.size(); i++)
    {
        CSegment3D& seg = m_edges[i];
        CPoint& pt = seg.start();

        double R = pt[0] * pt[0] + pt[1] * pt[1];

        if (R < 1 - epsilon)
        {
            inner_pts.push_back(pt);
        }
        if (R > 1 + epsilon)
        {
            outer_pts.push_back(pt);
        }
        if (R <= 1 + epsilon && R >= 1 - epsilon)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(pt);
            outer_pts.push_back(pt);
        }

        if (i == m_edges.size() - 1)
            continue;

        std::vector<CPoint> intersection_points;
        bool cross = seg.intersect(circ, intersection_points);
        if (!cross)
            continue;

        for (size_t j = 0; j < intersection_points.size(); j++)
        {
            on_circle_idx.push_back((int) inner_pts.size());
            inner_pts.push_back(intersection_points[j]);
            outer_pts.push_back(intersection_points[j]);
        }
    }

    {
        CSegment3D& sg = m_edges.back();
        CPoint st = sg.start();
        CPoint dr = sg.end() - sg.start();
        double n = dr.norm();

        CRay3D ray(st, dr);
        std::vector<CPoint> intersection_points;
        bool cross = ray.intersect(circ, intersection_points);

        if (!intersection_points.empty())
        {
            on_circle_idx.push_back((int) inner_pts.size());
            CPoint p = intersection_points.front();
            CPoint2 q(p[0], p[1]);
            double n = q.norm();
            inner_pts.push_back(p);
            outer_pts.push_back(intersection_points.front());
        }
    }

    // completely outside
    if (inner_pts.empty())
        return -1;

    // completely inside
    if (outer_pts.empty())
        return +1;

    int n = (int) inner_pts.size();
    int on_circle_indx = 0;

    for (int i = 0; i < (int) on_circle_idx.size(); i++)
    {
        if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
        {
            on_circle_indx = on_circle_idx[i];
        }
    }

    for (int i = 0; i < n; i++)
    {
        CSegment3D seg(inner_pts[(i + on_circle_indx) % n], inner_pts[(i + 1 + on_circle_indx) % n]);
        inner_poly.add(seg);
    }

    /*
                    CSegment2D & seg = inner_poly.edges().front();
                    CPoint2 sp = seg.start();
                    CPoint2 ep = seg.end();
                    std::cout << sp.norm() << " " << ep.norm() << std::endl;
    */
    return 0;
}

inline CPoint CPolygon3D::mass_center()
{
    std::vector<CPoint> pts;
    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment3D& s = m_edges[i];
        pts.push_back(s.start());
    }

    CPoint center(0, 0, 0);
    double area = 0;

    for (int j = 1; j < pts.size() - 1; j++)
    {
        CPoint c = (pts[0] + pts[j] + pts[j + 1]) / 3;
        CPoint e = pts[j + 0] - pts[0];
        CPoint d = pts[j + 1] - pts[0];
        CPoint n = e ^ d;
        double s = n.norm() / 2;
        area += s;
        center += c * s;
    }

    center /= area;

    return center;
}

inline CPolygon2D CPolygon3D::project()
{
    CPolygon2D poly;

    for (size_t i = 0; i < m_edges.size(); i++)
    {
        CSegment3D seg = m_edges[i];
        CPoint2 st(seg.start()[0], seg.start()[1]);
        CPoint2 ed(seg.end()[0], seg.end()[1]);
        CSegment2D sg(st, ed);
        poly.add(sg);
    }

    return poly;
}

}; // namespace MeshLib
#endif
