#ifndef _MESHLIB_POLYGON_2D_H_
#define _MESHLIB_POLYGON_2D_H_

/*!
*      \file Polygon2D.h
*      \brief planar polygon
*	   \author David Gu
*      \date 05/01/2020
*
*/
#include <vector>
#include "Segment2D.h"
#include "Ray2D.h"

namespace MeshLib {


	/*!
	*	\brief planar polygon
	*
	*/
	class CPolygon2D
	{
	public:

		/*!
		*	CPolygon constructor
		*/
		CPolygon2D() {};
		~CPolygon2D() {};
		
		/*! add one segment in sequencial order */
		void add(CSegment2D& seg)
		{
			m_edges.push_back(seg);
		};

		/*! polygon area */
		double area();

		/*! clip by the unit disk*/
		int finite_clip( CPolygon2D& inner_poly );

		/*! the edges*/
		std::vector<CSegment2D>& edges() { return m_edges; };

		/*! the polygon is an infinite cell, clipped by the unit circle */
		int infinite_clip(CPolygon2D & inner_poly );

		/*! centeter of mass */
		CPoint2 mass_center();

	protected:
		std::vector<CSegment2D> m_edges;
	};

	/*! the total area of a planar polygon */

	inline double CPolygon2D::area()
	{
		double s = 0;
		for (size_t i = 0; i < m_edges.size(); i++)
		{
			CSegment2D& seg = m_edges[i];
			s += seg.start() ^ seg.end();
		}
		return s /2.0;
	};

	/*! the finite polygon is clipped by the unit disk */
	inline int CPolygon2D::finite_clip( CPolygon2D & inner_poly )
	{
		double epsilon = 1e-30;

		CCircle circ(CPoint2(0, 0), 1);
		std::vector<int> on_circle_idx;

		std::vector<CPoint2> inner_pts;
		std::vector<CPoint2> outer_pts;

		for (size_t i = 0; i < m_edges.size(); i++)
		{
			CSegment2D& seg = m_edges[i];
			CPoint2& pt = seg.start();

			double R = pt.norm2();

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
			std::vector<CPoint2> intersection_points;
			bool cross = seg.intersect(circ, intersection_points);
			if (!cross) continue;

			for (size_t j = 0; j < intersection_points.size(); j++)
			{
				on_circle_idx.push_back((int)inner_pts.size());
				inner_pts.push_back(intersection_points[j]);
				outer_pts.push_back(intersection_points[j]);
			}
		}


		//completely outside
		if ( inner_pts.empty() ) return -1;

		//completely inside
		if (outer_pts.empty() ) return +1;

		int n = (int)inner_pts.size();
		int on_circle_indx = 0;

		for (int i = 0; i < (int)on_circle_idx.size(); i++)
		{
			if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
			{
				on_circle_indx = on_circle_idx[i];
			}
		}

		for (int i = 0; i < n; i++)
		{
			CSegment2D seg(inner_pts[(i+on_circle_indx)%n], inner_pts[(i + 1+on_circle_indx) % n]);
			inner_poly.add(seg);
		}

		return 0;
	}


	/*! the finite polygon is clipped by the unit disk */
	inline int CPolygon2D::infinite_clip(CPolygon2D& inner_poly)
	{
		double epsilon = 1e-30;

		CCircle circ(CPoint2(0, 0), 1);
		std::vector<int> on_circle_idx;

		std::vector<CPoint2> inner_pts;
		std::vector<CPoint2> outer_pts;

		{
			CSegment2D& sg = m_edges.front();
			CPoint2  st = sg.start();
			CPoint2  dr = sg.end() - sg.start();
			CRay2D ray(st,dr);
			std::vector<CPoint2> intersection_points;
			bool cross = ray.intersect(circ, intersection_points);
			
			if (!intersection_points.empty())
			{
				on_circle_idx.push_back((int)inner_pts.size());
				inner_pts.push_back(intersection_points.front());
				outer_pts.push_back(intersection_points.front());
			}

		}


		for (size_t i = 1; i < m_edges.size(); i++)
		{
			CSegment2D& seg = m_edges[i];
			CPoint2& pt = seg.start();

			double R = pt.norm2();

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
				on_circle_idx.push_back((int)inner_pts.size());
				inner_pts.push_back(pt);
				outer_pts.push_back(pt);
			}

			if (i == m_edges.size() - 1) continue;

			std::vector<CPoint2> intersection_points;
			bool cross = seg.intersect(circ, intersection_points);
			if (!cross) continue;

			for (size_t j = 0; j < intersection_points.size(); j++)
			{
				on_circle_idx.push_back((int)inner_pts.size());
				inner_pts.push_back(intersection_points[j]);
				outer_pts.push_back(intersection_points[j]);
			}
		}


		{
			CSegment2D& sg = m_edges.back();
			CPoint2  st = sg.start();
			CPoint2  dr = sg.end() - sg.start();
			CRay2D ray(st, dr);
			std::vector<CPoint2> intersection_points;
			bool cross = ray.intersect(circ, intersection_points);

			if (!intersection_points.empty())
			{
				on_circle_idx.push_back((int)inner_pts.size());
				inner_pts.push_back(intersection_points.front());
				outer_pts.push_back(intersection_points.front());
			}			
		}

		//completely outside
		if (inner_pts.empty()) return -1;

		//completely inside
		if (outer_pts.empty()) return +1;

		int n = (int)inner_pts.size();
		int on_circle_indx = 0;

		for (int i = 0; i < (int)on_circle_idx.size(); i++)
		{
			if ((on_circle_idx[i] + 1) % n == on_circle_idx[(i + 1) % on_circle_idx.size()])
			{
				on_circle_indx = on_circle_idx[i];
			}
		}

		for (int i = 0; i < n; i++)
		{
			CSegment2D seg(inner_pts[(i + on_circle_indx) % n], inner_pts[(i + 1 + on_circle_indx) % n]);
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

	inline CPoint2 CPolygon2D::mass_center()
	{
		std::vector<CPoint2> pts;
		for (size_t i = 0; i < m_edges.size(); i++)
		{
			CSegment2D & s = m_edges[i];
			pts.push_back(s.start());
		}

		CPoint2 center(0,0);
		double  area = 0;

		for (int j = 1; j < pts.size() - 1; j++)
		{
			CPoint2 c = (pts[0] + pts[j] + pts[j + 1]) / 3;
			double s = (pts[j] - pts[0]) ^ (pts[j + 1] - pts[0]) / 2;
			area += s;
			center += c * s;
		}

		center /= area;

		return center;
	}
};
#endif

