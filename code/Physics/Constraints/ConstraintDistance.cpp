//
//  ConstraintDistance.cpp
//
#include "ConstraintDistance.h"

void ConstraintDistance::PreSolve( const float dt_sec ) {
	const Vec3 worldAnchorA = m_bodyA->BodySpaceToWorldSpace(m_anchorA);
	const Vec3 worldAnchorB = m_bodyB->BodySpaceToWorldSpace(m_anchorB);

	const Vec3 r = worldAnchorB - worldAnchorA;
	const Vec3 ra = worldAnchorA - m_bodyA->GetCenterOfMassWorldSpace();
	const Vec3 rb = worldAnchorB - m_bodyB->GetCenterOfMassWorldSpace();

	m_Jacobian.Zero();

	//a = anchorA, b = anchorB
	//J = (2*(a-b),2*ra x (a-b),2*(b-a),2*rb x (b-a))
	Vec3 J1 = (worldAnchorA - worldAnchorB) * 2.0f;
	m_Jacobian.rows[0][0] = J1.x;
	m_Jacobian.rows[0][1] = J1.y;
	m_Jacobian.rows[0][2] = J1.z;

	Vec3 J2 = ra.Cross((worldAnchorA - worldAnchorB) * 2.0f);
	m_Jacobian.rows[0][3] = J2.x;
	m_Jacobian.rows[0][4] = J2.y;
	m_Jacobian.rows[0][5] = J2.z;

	Vec3 J3 = (worldAnchorB - worldAnchorA) * 2.0f;
	m_Jacobian.rows[0][6] = J3.x;
	m_Jacobian.rows[0][7] = J3.y;
	m_Jacobian.rows[0][8] = J3.z;

	Vec3 J4 = rb.Cross((worldAnchorB - worldAnchorA) * 2.0f);
	m_Jacobian.rows[0][9] = J4.x;
	m_Jacobian.rows[0][10] = J4.y;
	m_Jacobian.rows[0][11] = J4.z;

	//apply lambda from previous frame for a warm start
	ApplyImpulses(m_Jacobian.Transpose() * m_cachedLambda);

	//Calculate Baumgarte stabilisation
	float C = r.Dot(r);
	C = std::max(0.0f, C - 0.01f);//ensure non-negative, then, if the value is very close to 0, no correction needed
	const float Beta = 0.05f;//manually chosen beta factor so we don't add a lot of energy but still corrects the constraint
	m_baumgarte = Beta * C / dt_sec;
}

void ConstraintDistance::Solve() {
	const MatMN JacobianTranspose = m_Jacobian.Transpose();

	const VecN q_dt = GetVelocities();
	const MatMN invMassMatrix = GetInverseMassMatrix();
	const MatMN J_W_Jt = m_Jacobian * invMassMatrix * JacobianTranspose;//denominator
	VecN rhs = m_Jacobian * q_dt * -1.0f; //-J*v (numerator)
	rhs[0] -= m_baumgarte;

	//calculate lagrange multiplier
	const VecN lambdaN = LCP_GaussSeidel(J_W_Jt, rhs);

	//Qint = JT * lambda
	ApplyImpulses(JacobianTranspose * lambdaN);
	
	//accumulation, interesting
	m_cachedLambda += lambdaN;
}

void ConstraintDistance::PostSolve() {
	/*
	for (int i = 0; i < m_cachedLambda.N; i++) {
		if (m_cachedLambda[i] * 0.0f != m_cachedLambda[i] * 0.0f) {
			m_cachedLambda[i] = 0.0f;
		}
		const float limit = 1e2f;
		if (m_cachedLambda[i] > limit) {
			m_cachedLambda[i] = limit;
		}
		if (m_cachedLambda[i] < -limit) {
			m_cachedLambda[i] = -limit;
		}
	}
	*/
	if (m_cachedLambda[0] * 0.0f != m_cachedLambda[0] * 0.0f) {
		m_cachedLambda[0] = 0.0f;
		return;
	}
	const float limit = 1e2f;
	if (m_cachedLambda[0] > limit) {
		m_cachedLambda[0] = limit;
	}
	if (m_cachedLambda[0] < -limit) {
		m_cachedLambda[0] = -limit;
	}

}