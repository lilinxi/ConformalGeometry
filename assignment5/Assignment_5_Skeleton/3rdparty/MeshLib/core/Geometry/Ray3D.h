#ifndef _RAY_3D_H_
#define _RAY_3D_H_

/*!
*      \file Ray3D.h
*      \brief Planar Ray
*	   \author David Gu
*      \date 05/03/2020
*
*/
#include "Circle.h"
#include "Point2.h"
#include "Point.h"
#include "Segment3D.h"
#include <vector>

namespace MeshLib {

	/*!
	*	\brief CRay3D class, Ray in the space
	*
	*/
	class CRay3D : public CSegment3D
	{
	public:

		/*!
		*	CRay3D Constructor,
		*   \param n normal
		*   \param s starting point
		*/

		CRay3D(CPoint start, CPoint d) { assert( d.norm() > 0); m_direction = d / d.norm(); m_start = start; };
		CRay3D(const CRay3D& ray)
		{
			m_start = ray.m_start;
			m_direction = ray.m_direction;
		};

		~CRay3D() {};

		virtual bool intersect(CCircle& circle, std::vector<CPoint>& intersection_points);

	protected:

		CPoint m_direction;
	};

	inline bool CRay3D::intersect(CCircle& circ, std::vector<CPoint>& intersection_points)
	{
		CPoint2 st(m_start[0], m_start[1]);
		CPoint2 dr(m_direction[0], m_direction[1]);
		
		//note that, dr shouldn't be normalized

		double A = dr.norm2();
		double B = 2 * (dr * (st - circ.c()));
		double C = (st - circ.c()).norm2() - circ.r() * circ.r();

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
