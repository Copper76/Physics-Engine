//
//  Body.cpp
//
#include "Body.h"

Vec3 Body::GetCenterOfMassWorldSpace() const {
	const Vec3 centerOfMass = m_shape->GetCenterOfMass();
	return m_position + m_orientation.RotatePoint(centerOfMass);
}

Vec3 Body::GetCenterOfMassModelSpace() const {
	const Vec3 centerOfMass = m_shape->GetCenterOfMass();
	return centerOfMass;
}

Vec3 Body::WorldSpaceToBodySpace(const Vec3& worldPt) const {
	Quat inverseOrient = m_orientation.Inverse();
	return inverseOrient.RotatePoint(worldPt - GetCenterOfMassWorldSpace());
}

Vec3 Body::BodySpaceToWorldSpace(const Vec3& worldPt) const {
	return GetCenterOfMassWorldSpace() + m_orientation.RotatePoint(worldPt);
}
Mat3 Body::GetInverseInertiaTensorBodySpace() const {
	Mat3 inertiaTensor = m_shape->InertiaTensor();
	return inertiaTensor.Inverse() * m_invMass;
}

//used with global vectors
Mat3 Body::GetInverseInertiaTensorWorldSpace() const {
	Mat3 inertiaTensor = m_shape->InertiaTensor();
	Mat3 invInertiaTensor = inertiaTensor.Inverse() * m_invMass;
	Mat3 orient = m_orientation.ToMat3();
	//the transposed orient converts a global vector to local space for calculation before translating it back into global space with orient
	invInertiaTensor = orient * invInertiaTensor * orient.Transpose();
	return invInertiaTensor;
}

void Body::ApplyImpulse(const Vec3& impulsePoint, const Vec3& impulse) {
	if (0.0f == m_invMass) { return; }

	ApplyImpulseLinear(impulse);

	Vec3 position = GetCenterOfMassWorldSpace();	// applying impulses must produce torques through the center of mass
	Vec3 r = impulsePoint - position;
	Vec3 dL = r.Cross(impulse);	// this is in world space
	ApplyImpulseAngular(dL);
}

void Body::ApplyImpulseLinear(const Vec3& impulse) {
	if (0.0f == m_invMass) { return; }

	//p=mv
	//dp = m*dv
	//dv = dp/m = impulse/m
	m_linearVelocity += impulse * m_invMass;
}

void Body::ApplyImpulseAngular(const Vec3& impulse) {
	if (0.0f == m_invMass) { return; }

	//change in speed = inertia*impulse
	m_angularVelocity += GetInverseInertiaTensorWorldSpace() * impulse;

	//A manual upper limit for angular speed
	const float maxAngularSpeed = 30.0f;
	if (m_angularVelocity.GetLengthSqr() > maxAngularSpeed * maxAngularSpeed) {
		m_angularVelocity.Normalize();
		m_angularVelocity *= maxAngularSpeed;
	}
}

void Body::Update(const float dt_sec) {
	m_position += m_linearVelocity * dt_sec;

	Vec3 positionCM = GetCenterOfMassWorldSpace();
	Vec3 cmToPos = m_position - positionCM;

	//Need to update the rotation torque as it needs to be always perpendicular to center of mass
	//alpha is carried over acceleration from last iteration
	Mat3 orientation = m_orientation.ToMat3();
	Mat3 inertiaTensor = orientation * m_shape->InertiaTensor() * orientation.Transpose();
	Vec3 alpha = inertiaTensor.Inverse() * (m_angularVelocity.Cross(inertiaTensor * m_angularVelocity));
	m_angularVelocity += alpha * dt_sec;

	Vec3 dAngle = m_angularVelocity * dt_sec;
	Quat dq = Quat(dAngle, dAngle.GetMagnitude());
	m_orientation = dq * m_orientation;
	m_orientation.Normalize();

	m_position = positionCM + dq.RotatePoint(cmToPos);
}