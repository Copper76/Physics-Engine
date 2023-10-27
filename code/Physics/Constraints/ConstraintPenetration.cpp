//
//  ConstraintPenetration.cpp
//
#include "ConstraintPenetration.h"

/*
================================
ConstraintPenetration::PreSolve
================================
*/
void ConstraintPenetration::PreSolve( const float dt_sec ) {
	const Vec3 worldAnchorA = m_bodyA->BodySpaceToWorldSpace(m_anchorA);
	const Vec3 worldAnchorB = m_bodyB->BodySpaceToWorldSpace(m_anchorB);

	const Vec3 ra = worldAnchorA - m_bodyA->GetCenterOfMassWorldSpace();
	const Vec3 rb = worldAnchorB - m_bodyB->GetCenterOfMassWorldSpace();

	const float frictionA = m_bodyA->m_friction;
	const float frictionB = m_bodyB->m_friction;
	m_friction = frictionA * frictionB;

	Vec3 u;
	Vec3 v;
	m_normal.GetOrtho(u, v);//but u and v are undefined?

	Vec3 normal = m_bodyA->m_orientation.RotatePoint(m_normal);//why A?
	u = m_bodyA->m_orientation.RotatePoint(u);
	v = m_bodyA->m_orientation.RotatePoint(v);

	m_Jacobian.Zero();

	//row 1 is the primary distance constraint that holds anchor points together
	Vec3 J1 = normal * -1.0f;
	m_Jacobian.rows[0][0] = J1.x;
	m_Jacobian.rows[0][1] = J1.y;
	m_Jacobian.rows[0][2] = J1.z;

	Vec3 J2 = ra.Cross(normal * -1.0f);
	m_Jacobian.rows[0][3] = J2.x;
	m_Jacobian.rows[0][4] = J2.y;
	m_Jacobian.rows[0][5] = J2.z;

	Vec3 J3 = normal * 1.0f;
	m_Jacobian.rows[0][6] = J3.x;
	m_Jacobian.rows[0][7] = J3.y;
	m_Jacobian.rows[0][8] = J3.z;

	Vec3 J4 = rb.Cross(normal * 1.0f);
	m_Jacobian.rows[0][9] = J4.x;
	m_Jacobian.rows[0][10] = J4.y;
	m_Jacobian.rows[0][11] = J4.z;

	//Friction Jacobian
	if (m_friction > 0.0f) {
		Vec3 J1 = u * -1.0f;
		m_Jacobian.rows[1][0] = J1.x;
		m_Jacobian.rows[1][1] = J1.y;
		m_Jacobian.rows[1][2] = J1.z;

		Vec3 J2 = ra.Cross(u * -1.0f);
		m_Jacobian.rows[1][3] = J2.x;
		m_Jacobian.rows[1][4] = J2.y;
		m_Jacobian.rows[1][5] = J2.z;

		Vec3 J3 = u * 1.0f;
		m_Jacobian.rows[1][6] = J3.x;
		m_Jacobian.rows[1][7] = J3.y;
		m_Jacobian.rows[1][8] = J3.z;

		Vec3 J4 = rb.Cross(u * 1.0f);
		m_Jacobian.rows[1][9] = J4.x;
		m_Jacobian.rows[1][10] = J4.y;
		m_Jacobian.rows[1][11] = J4.z;

		J1 = v * -1.0f;
		m_Jacobian.rows[2][0] = J1.x;
		m_Jacobian.rows[2][1] = J1.y;
		m_Jacobian.rows[2][2] = J1.z;

		J2 = ra.Cross(v * -1.0f);
		m_Jacobian.rows[2][3] = J2.x;
		m_Jacobian.rows[2][4] = J2.y;
		m_Jacobian.rows[2][5] = J2.z;

		J3 = v * 1.0f;
		m_Jacobian.rows[2][6] = J3.x;
		m_Jacobian.rows[2][7] = J3.y;
		m_Jacobian.rows[2][8] = J3.z;

		J4 = rb.Cross(v * 1.0f);
		m_Jacobian.rows[2][9] = J4.x;
		m_Jacobian.rows[2][10] = J4.y;
		m_Jacobian.rows[2][11] = J4.z;
	}

	//Apply warm start (don't need post solve to resolve chain issues?)
	ApplyImpulses(m_Jacobian.Transpose() * m_cachedLambda);

	//Calculate Baumgarte stabilisation
	float C = (worldAnchorB-worldAnchorA).Dot(normal);
	C = std::min(0.0f, C + 0.02f);//ensure non-negative, then, if the value is very close to 0, no correction needed
	const float Beta = 0.25f;//manually chosen beta factor so we don't add a lot of energy but still corrects the constraint
	m_baumgarte = Beta * C / dt_sec;
}

/*
================================
ConstraintPenetration::Solve
================================
*/
void ConstraintPenetration::Solve() {
	const MatMN JacobianTranspose = m_Jacobian.Transpose();

	const VecN q_dt = GetVelocities();
	const MatMN invMassMatrix = GetInverseMassMatrix();
	const MatMN J_W_Jt = m_Jacobian * invMassMatrix * JacobianTranspose;//denominator
	VecN rhs = m_Jacobian * q_dt * -1.0f; //-J*v (numerator)
	rhs[0] -= m_baumgarte;

	//calculate lagrange multiplier
	VecN lambdaN = LCP_GaussSeidel(J_W_Jt, rhs);

	//accumulation, interesting
	VecN oldLambda = m_cachedLambda;
	m_cachedLambda += lambdaN;
	const float lambdaLimit = 0.0f;
	if (m_cachedLambda[0] < lambdaLimit) {
		m_cachedLambda[0] = lambdaLimit;//cap the lambda to be non-negative
	}
	if (m_friction > 0.0f) {
		const float umg = m_friction * 10.0f / (m_bodyA->m_invMass + m_bodyB->m_invMass);//why times by 10
		const float normalForce = fabsf(lambdaN[0] * m_friction);
		const float maxForce = (umg > normalForce) ? umg : normalForce;

		if (m_cachedLambda[1] > maxForce) {
			m_cachedLambda[1] = maxForce;
		}
		if (m_cachedLambda[1] < -maxForce) {
			m_cachedLambda[1] = -maxForce;
		}
		if (m_cachedLambda[2] > maxForce) {
			m_cachedLambda[2] = maxForce;
		}
		if (m_cachedLambda[2] < -maxForce) {
			m_cachedLambda[2] = -maxForce;
		}
	}
	lambdaN = m_cachedLambda - oldLambda;

	//Qint = JT * lambda
	ApplyImpulses(JacobianTranspose * lambdaN);
}