#ifndef _RAY_2D_H_
#define _RAY_2D_H_

/*!
*      \file Ray2D.h
*      \brief Planar Ray
*	   \author David Gu
*      \date 05/01/2020
*
*/
#include "Circle.h"
#include "Point2.h"
#include "Segment2D.h"
#include <vector>

namespace MeshLib {

	/*!
	*	\brief CRay2D class, Ray on the plane
	*
	*	Line on the two dimensional plane <p,n> = d
	*/
	class CRay2D : public CSegment2D
	{
	public:

		/*!
		*	CRay2D Constructor,
		*   \param n normal
		*   \param s starting point
		*/

		CRay2D(CPoint2 start, CPoint2 d) { assert(mag2(d) > 0); m_direction = d / mag(d); m_start = start; };
		CRay2D(const CRay2D& ray)
		{
			m_start = ray.m_start;
			m_direction = ray.m_direction;
		};

		~CRay2D() {};

		virtual bool intersect(CCircle& circle, std::vector<CPoint2>& intersection_points);

	protected:

		CPoint2 m_direction;
	};

	inline bool CRay2D::intersect(CCircle& circ, std::vector<CPoint2>& intersection_points)
	{
		double A = 1;
		double B = 2 * (m_direction * (m_start - circ.c()));
		double C = (m_start - circ.c()).norm2() - circ.r() * circ.r();

		double D = B * B - 4 * A * C;
		if (D < 0) return false;

		double t = (-B - sqrt(D)) / (2 * A);

		if (t >= 0)
		{
			intersection_points.push_back(m_start + m_direction * t);
		}

		t = (-B + sqrt(D)) / (2 * A);

		if (t >= 0)
		{
			intersection_points.push_back(m_start + m_direction * t);
		}

		return intersection_points.size() > 0;
	}
};
#endif