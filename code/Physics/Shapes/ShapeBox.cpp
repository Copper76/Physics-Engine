//
//  Shapes.cpp
//
#include "ShapeBox.h"

/*
========================================================================================================

ShapeBox

========================================================================================================
*/

/*
====================================================
ShapeBox::Build
====================================================
*/
void ShapeBox::Build( const Vec3 * pts, const int num ) {
	for (int i = 0; i < num; i++) {
		m_bounds.Expand(pts[i]);
	}

	m_points.clear();
	m_points.push_back(Vec3(m_bounds.mins.x, m_bounds.mins.y, m_bounds.mins.z));
	m_points.push_back(Vec3(m_bounds.maxs.x, m_bounds.mins.y, m_bounds.mins.z));
	m_points.push_back(Vec3(m_bounds.mins.x, m_bounds.maxs.y, m_bounds.mins.z));
	m_points.push_back(Vec3(m_bounds.mins.x, m_bounds.mins.y, m_bounds.maxs.z));

	m_points.push_back(Vec3(m_bounds.maxs.x, m_bounds.maxs.y, m_bounds.maxs.z));
	m_points.push_back(Vec3(m_bounds.mins.x, m_bounds.maxs.y, m_bounds.maxs.z));
	m_points.push_back(Vec3(m_bounds.maxs.x, m_bounds.mins.y, m_bounds.maxs.z));
	m_points.push_back(Vec3(m_bounds.maxs.x, m_bounds.maxs.y, m_bounds.mins.z));

	m_centerOfMass = (m_bounds.maxs + m_bounds.mins) * 0.5f;
}

/*
====================================================
ShapeBox::Support
====================================================
*/
Vec3 ShapeBox::Support( const Vec3 & dir, const Vec3 & pos, const Quat & orient, const float bias ) const {
	Vec3 maxPt = orient.RotatePoint(m_points[0]) + pos;
	float maxDist = dir.Dot(maxPt);
	for (int i = 1; i < m_points.size(); i++) {
		const Vec3 pt = orient.RotatePoint(m_points[i]) + pos;
		const float dist = dir.Dot(pt);
		if (dist > maxDist) {
			maxDist = dist;
			maxPt = pt;
		}
	}

	Vec3 norm = dir;
	norm.Normalize();
	norm *= bias;
	return maxPt + norm;
}

/*
====================================================
ShapeBox::InertiaTensor
====================================================
*/
Mat3 ShapeBox::InertiaTensor() const {
	//inertiaTensor around origin
	Mat3 tensor;
	const float dx = m_bounds.maxs.x - m_bounds.mins.x;
	const float dy = m_bounds.maxs.y - m_bounds.mins.y;
	const float dz = m_bounds.maxs.z - m_bounds.mins.z;

	tensor.Zero();
	tensor.rows[0][0] = (dy * dy + dz * dz) / 12.0f;
	tensor.rows[1][1] = (dx * dx + dz * dz) / 12.0f;
	tensor.rows[2][2] = (dx * dx + dy * dy) / 12.0f;

	//add modification around a different axis
	Vec3 cm = Vec3((m_bounds.maxs.x + m_bounds.mins.x) * 0.5f, (m_bounds.maxs.y + m_bounds.mins.y) * 0.5f, (m_bounds.maxs.z + m_bounds.mins.z) * 0.5f);

	const Vec3 R = Vec3(0, 0, 0) - cm;
	const float R2 = R.GetLengthSqr();

	Mat3 patTensor;
	patTensor.rows[0] = Vec3(R2 - R.x * R.x, R.x * R.y, R.x * R.z);
	patTensor.rows[1] = Vec3(R.y * R.x, R2 - R.y * R.y, R.y * R.z);
	patTensor.rows[2] = Vec3(R.z * R.x, R.z * R.y, R2 - R.z * R.z);
	tensor += patTensor;

	return tensor;
}

/*
====================================================
ShapeBox::GetBounds
====================================================
*/
Bounds ShapeBox::GetBounds( const Vec3 & pos, const Quat & orient ) const {
	Vec3 corners[8];

	corners[0] = Vec3(m_bounds.mins.x, m_bounds.mins.y, m_bounds.mins.z);
	corners[1] = Vec3(m_bounds.maxs.x, m_bounds.mins.y, m_bounds.mins.z);
	corners[2] = Vec3(m_bounds.mins.x, m_bounds.maxs.y, m_bounds.mins.z);
	corners[3] = Vec3(m_bounds.mins.x, m_bounds.mins.y, m_bounds.maxs.z);

	corners[4] = Vec3(m_bounds.maxs.x, m_bounds.maxs.y, m_bounds.maxs.z);
	corners[5] = Vec3(m_bounds.mins.x, m_bounds.maxs.y, m_bounds.maxs.z);
	corners[6] = Vec3(m_bounds.maxs.x, m_bounds.mins.y, m_bounds.maxs.z);
	corners[7] = Vec3(m_bounds.maxs.x, m_bounds.maxs.y, m_bounds.mins.z);

	Bounds bounds;
	for (int i = 0; i < 8; i++) {
		bounds.Expand(orient.RotatePoint(corners[i]) + pos);
	}

	return bounds;
}

/*
====================================================
ShapeBox::FastestLinearSpeed
====================================================
*/
float ShapeBox::FastestLinearSpeed( const Vec3 & angularVelocity, const Vec3 & dir ) const {
	float maxSpeed = 0.0f;
	for (int i = 0; i < m_points.size(); i++) {
		Vec3 r = m_points[i] - m_centerOfMass;
		Vec3 linearVelocity = angularVelocity.Cross(r);//perpendicular component of the angular speed
		float speed = dir.Dot(linearVelocity);//component of the perpendicular speed in direction of dir
		if (speed > maxSpeed) {
			maxSpeed = speed;
		}
	}

	return maxSpeed;
}