#ifndef _MESHLIB_SEGMENT_3D_H_
#define _MESHLIB_SEGMENT_3D_H_

/*!
*      \file Segment 3D.h
*      \brief Planar segment
*	   \author David Gu
*      \date 05/03/2020
*
*/
#include "Circle.h"
#include "Point.h"
#include <vector>

namespace MeshLib {

	/*!
	*	\brief CSegment3D class, a finite Segment on the plane
	*
	*/
	class CSegment3D
	{
	public:

		/*!
		*	CSegment3D Constructor,
		*   \param s starting point
		*   \param e ending   point
		*/

		CSegment3D(CPoint st, CPoint ed) { m_start = st; m_end = ed; };
		CSegment3D(const CSegment3D& segment)
		{
			m_start = segment.m_start;
			m_end = segment.m_end;
		};
		CSegment3D() {};

		~CSegment3D() {};
		/*
		* intersection with the cylinder, whose section is a planar circle
		*/
		bool intersect(CCircle& circle, std::vector<CPoint>& intersection_points);

		/* start point */
		CPoint& start() { return m_start; };
		CPoint& end()   { return m_end; };

	protected:

		CPoint m_start;
		CPoint m_end;
	};

	inline bool CSegment3D::intersect(CCircle& circ, std::vector<CPoint>& intersection_point)
	{
		CPoint2 st(m_start[0], m_start[1]);
		CPoint2 ed(m_end[0], m_end[1]);
		
		double A = (st - ed).norm2();
		double B = 2 * ((st - ed) * (ed - circ.c()));
		double C = (ed - circ.c()).norm2() - circ.r() * circ.r();

		double D = B * B - 4 * A * C;
		if (D < 0) return false;

		double t = (-B - sqrt(D)) / (2 * A);

		if (t >= 0 && t <= 1)
		{
			intersection_point.push_back((m_start * t) + (m_end * (1 - t)));
		}
		t = (-B + sqrt(D)) / (2 * A);
		if (t >= 0 && t <= 1)
		{
			intersection_point.push_back((m_start * t) + (m_end * (1 - t)));
		}

		return intersection_point.size() > 0;

	};
};
#endif
