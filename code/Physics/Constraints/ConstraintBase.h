//
//	ConstraintBase.h
//
#pragma once
#include "../../Math/Vector.h"
#include "../../Math/Quat.h"
#include "../../Math/Matrix.h"
#include "../../Math/Bounds.h"
#include "../../Math/LCP.h"
#include "../Body.h"
#include <vector>

/*
====================================================
Constraint
====================================================
*/
class Constraint {
public:
	virtual void PreSolve( const float dt_sec ) {}
	virtual void Solve() {}
	virtual void PostSolve() {}

	static Mat4 Left( const Quat & q );
	static Mat4 Right( const Quat & q );

protected:
	MatMN GetInverseMassMatrix() const;
	VecN GetVelocities() const;
	void ApplyImpulses( const VecN & impulses );

public:
	Body * m_bodyA;
	Body * m_bodyB;

	Vec3 m_anchorA;		// The anchor location in bodyA's space
	Vec3 m_axisA;		// The axis direction in bodyA's space

	Vec3 m_anchorB;		// The anchor location in bodyB's space
	Vec3 m_axisB;		// The axis direction in bodyB's space
};

/*
====================================================
Constraint::GetInverseMassMatrix
====================================================
*/
inline MatMN Constraint::GetInverseMassMatrix() const {
	MatMN invMassMatrix( 12, 12 );
	invMassMatrix.Zero();

	invMassMatrix.rows[0][0] = m_bodyA->m_invMass;
	invMassMatrix.rows[1][1] = m_bodyA->m_invMass;
	invMassMatrix.rows[2][2] = m_bodyA->m_invMass;

	Mat3 invInertiaA = m_bodyA->GetInverseInertiaTensorWorldSpace();
	for (int i = 0; i < 3; i++) {
		invMassMatrix.rows[i+3][3] = invInertiaA.rows[i][0];
		invMassMatrix.rows[i+3][4] = invInertiaA.rows[i][1];
		invMassMatrix.rows[i+3][5] = invInertiaA.rows[i][2];
	}

	invMassMatrix.rows[6][6] = m_bodyB->m_invMass;
	invMassMatrix.rows[7][7] = m_bodyB->m_invMass;
	invMassMatrix.rows[8][8] = m_bodyB->m_invMass;

	Mat3 invInertiaB = m_bodyB->GetInverseInertiaTensorWorldSpace();
	for (int i = 0; i < 3; i++) {
		invMassMatrix.rows[i + 9][9] = invInertiaB.rows[i][0];
		invMassMatrix.rows[i + 9][10] = invInertiaB.rows[i][1];
		invMassMatrix.rows[i + 9][11] = invInertiaB.rows[i][2];
	}

	return invMassMatrix;
}

/*
====================================================
Constraint::GetVelocities
====================================================
*/
inline VecN Constraint::GetVelocities() const {
	VecN q_dt( 12 );

	q_dt[0] = m_bodyA->m_linearVelocity.x;
	q_dt[1] = m_bodyA->m_linearVelocity.y;
	q_dt[2] = m_bodyA->m_linearVelocity.z;

	q_dt[3] = m_bodyA->m_angularVelocity.x;
	q_dt[4] = m_bodyA->m_angularVelocity.y;
	q_dt[5] = m_bodyA->m_angularVelocity.z;

	q_dt[6] = m_bodyB->m_linearVelocity.x;
	q_dt[7] = m_bodyB->m_linearVelocity.y;
	q_dt[8] = m_bodyB->m_linearVelocity.z;

	q_dt[9] = m_bodyB->m_angularVelocity.x;
	q_dt[10] = m_bodyB->m_angularVelocity.y;
	q_dt[11] = m_bodyB->m_angularVelocity.z;

	return q_dt;
}

/*
====================================================
Constraint::ApplyImpulses
====================================================
*/
inline void Constraint::ApplyImpulses( const VecN & impulses ) {
	Vec3 tempForce(0.0f);

	//changed here to reduce amount of memory used
	tempForce[0] = impulses[0];
	tempForce[1] = impulses[1];
	tempForce[2] = impulses[2];

	m_bodyA->ApplyImpulseLinear(tempForce);

	tempForce[0] = impulses[3];
	tempForce[1] = impulses[4];
	tempForce[2] = impulses[5];

	m_bodyA->ApplyImpulseAngular(tempForce);

	tempForce[0] = impulses[6];
	tempForce[1] = impulses[7];
	tempForce[2] = impulses[8];

	m_bodyB->ApplyImpulseLinear(tempForce);

	tempForce[0] = impulses[9];
	tempForce[1] = impulses[10];
	tempForce[2] = impulses[11];

	m_bodyB->ApplyImpulseAngular(tempForce);
}

/*
====================================================
Constraint::Left
====================================================
*/
inline Mat4 Constraint::Left( const Quat & q ) {
	Mat4 L;
	L.rows[0] = Vec4(q.w, -q.x, -q.y, -q.z);
	L.rows[1] = Vec4(q.x, q.w, -q.z, q.y);
	L.rows[2] = Vec4(q.y, q.z, q.w, -q.x);
	L.rows[3] = Vec4(q.z, -q.y, q.x, q.w);

	return L.Transpose();
}

/*
====================================================
Constraint::Right
====================================================
*/
inline Mat4 Constraint::Right( const Quat & q ) {
	Mat4 R;
	R.rows[0] = Vec4(q.w, -q.x, -q.y, -q.z);
	R.rows[1] = Vec4(q.x, q.w, q.z, -q.y);
	R.rows[2] = Vec4(q.y, -q.z, q.w, q.x);
	R.rows[3] = Vec4(q.z, q.y, -q.x, q.w);

	return R.Transpose();
}