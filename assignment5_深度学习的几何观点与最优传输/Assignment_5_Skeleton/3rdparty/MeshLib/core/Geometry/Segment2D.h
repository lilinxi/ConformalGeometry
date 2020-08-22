#ifndef _MESHLIB_SEGMENT_2D_H_
#define _MESHLIB_SEGMENT_2D_H_

/*!
*      \file Segment 2D.h
*      \brief Planar segment
*	   \author David Gu
*      \date 05/01/2020
*
*/
#include "Circle.h"
#include "Point2.h"
#include <vector>

namespace MeshLib {

	/*!
	*	\brief CSegment2D class, a finite Segment on the plane
	*
	*/
	class CSegment2D
	{
	public:

		/*!
		*	CSegment2D Constructor,
		*   \param s starting point
		*   \param e ending   point
		*/

		CSegment2D(CPoint2 st, CPoint2 ed) { m_start = st; m_end = ed; };
		CSegment2D(const CSegment2D& segment)
		{
			m_start = segment.m_start;
			m_end = segment.m_end;
		};
		CSegment2D() {};

		~CSegment2D() {};

		bool intersect(CCircle& circle, std::vector<CPoint2>& intersection_points);

		/* start point */
		CPoint2 & start() { return m_start; };
		CPoint2 & end() { return m_end; };

	protected:

		CPoint2 m_start;
		CPoint2 m_end;
	};

	inline bool CSegment2D::intersect(CCircle& circ, std::vector<CPoint2>& intersection_point)
	{
		double A = (m_start - m_end).norm2();
		double B = 2 * ((m_start - m_end) * (m_end - circ.c()));
		double C = (m_end - circ.c()).norm2() - circ.r() * circ.r();

		double D = B * B - 4 * A * C;
		if (D < 0) return false;

		double t = (-B - sqrt(D)) / (2 * A);

		if (t >= 0 && t <= 1)
		{
			CPoint2 p = m_start * t;
			CPoint2 q = m_end * (1 - t);

			intersection_point.push_back(p + q);
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
