//
//  Scene.cpp
//
#include "Scene.h"
#include "Physics/Contact.h"
#include "Physics/Intersections.h"
#include "Physics/Broadphase.h"

/*
========================================================================================================

Scene

========================================================================================================
*/

/*
====================================================
Scene::~Scene
====================================================
*/
Scene::~Scene() {
	for ( int i = 0; i < m_bodies.size(); i++ ) {
		delete m_bodies[ i ].m_shape;
	}
	m_bodies.clear();
}

/*
====================================================
Scene::Reset
====================================================
*/
void Scene::Reset() {
	for ( int i = 0; i < m_bodies.size(); i++ ) {
		delete m_bodies[ i ].m_shape;
	}
	m_bodies.clear();

	for (int i = 0; i < m_constraints.size(); i++) {
		delete m_constraints[i];
	}
	m_constraints.clear();

	Initialize();
}

void AddStandardSandBox(std::vector<Body>& bodies) {
	Body body;

	body.m_position = Vec3(0, 0, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity.Zero();
	body.m_angularVelocity.Zero();
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.5f;
	body.m_shape = new ShapeBox(g_boxGround, sizeof(g_boxGround)/sizeof(Vec3));
	bodies.push_back(body);

	body.m_position = Vec3(50, 0, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity.Zero();
	body.m_angularVelocity.Zero();
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.0f;
	body.m_shape = new ShapeBox(g_boxWall0, sizeof(g_boxWall0) / sizeof(Vec3));
	bodies.push_back(body);

	body.m_position = Vec3(-50, 0, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity.Zero();
	body.m_angularVelocity.Zero();
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.0f;
	body.m_shape = new ShapeBox(g_boxWall0, sizeof(g_boxWall0) / sizeof(Vec3));
	bodies.push_back(body);

	body.m_position = Vec3(0, 25, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity.Zero();
	body.m_angularVelocity.Zero();
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.0f;
	body.m_shape = new ShapeBox(g_boxWall1, sizeof(g_boxWall1) / sizeof(Vec3));
	bodies.push_back(body);

	body.m_position = Vec3(0, -25, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity.Zero();
	body.m_angularVelocity.Zero();
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.0f;
	body.m_shape = new ShapeBox(g_boxWall1, sizeof(g_boxWall1) / sizeof(Vec3));
	bodies.push_back(body);
}

/*
====================================================
Scene::Initialize
====================================================
*/
void Scene::Initialize() {
	Body body;
	/*
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			body.m_position = Vec3(float(i-1) * 0.5f *1.5f, float(j - 1) * 0.5f * 1.5f, 20);
			body.m_orientation = Quat(0, 0, 0, 1);
			body.m_linearVelocity.Zero();
			body.m_invMass = 1.0f;
			body.m_elasticity = 0.5f;
			body.m_friction = 0.5f;
			body.m_shape = new ShapeSphere(0.5f);
			m_bodies.push_back(body);
		}
	}
	*/

	/*
	//pendulem
	const int numJoints = 3;
	body.m_position = Vec3(0, 0, 5.0f);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_invMass = 0.0f;
	body.m_shape = new ShapeBox(g_boxSmall, sizeof(g_boxSmall) / sizeof(Vec3));
	body.m_elasticity = 0.5f;
	body.m_friction = 0.5f;
	m_bodies.push_back(body);
	int anchor = m_bodies.size() - 1;
	for (int i = 0; i < numJoints; i++) {
		body.m_linearVelocity = Vec3(0, 0, 0);

		Body* bodyA = &m_bodies[anchor];

		const Vec3 jointWorldSpaceAnchor = bodyA->m_position;

		ConstraintDistance* joint = new ConstraintDistance();

		joint->m_bodyA = bodyA;
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		body.m_position = joint->m_bodyA->m_position + Vec3(i+1, 0, 0);
		body.m_orientation = Quat(0, 0, 0, 1);
		body.m_invMass = 1.0f;
		//body.m_shape = new ShapeConvex(g_diamond, sizeof(g_diamond)/sizeof(Vec3));
		//body.m_shape = new ShapeBox(g_boxSmall, sizeof(g_boxSmall) / sizeof(Vec3));
		body.m_shape = new ShapeSphere(0.5f);
		m_bodies.push_back(body);

		Body* bodyB = &m_bodies[m_bodies.size() - 1];

		joint->m_bodyB = bodyB;
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		m_constraints.push_back(joint);
	}
	*/
	/*
	//stacking box
	for (int z = 0; z < 5; z++) {
		float offset = ((z & 1) == 0) ? 0.0f : 0.15f;//offset even-numbered boxes
		float delta = 0.04f;
		float scaleHeight = 2.0f + delta;
		float deltaHeight = 1.0f + delta;
		body.m_position = Vec3(offset * scaleHeight, offset * scaleHeight, deltaHeight + (float)z * scaleHeight);
		body.m_orientation = Quat(0, 0, 0, 1);
		body.m_invMass = 1.0f;
		body.m_shape = new ShapeBox(g_boxUnit, sizeof(g_boxUnit) / sizeof(Vec3));
		body.m_elasticity = 0.5f;
		body.m_friction = 0.5f;
		m_bodies.push_back(body);
	}
	*/
	
	/*
	//
	//Ragdoll
	//

	//head
	body.m_position = Vec3(0, 0, 5.5f);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeBox(g_boxHead, sizeof(g_boxHead) / sizeof(Vec3));
	body.m_invMass = 2.0f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxHead = m_bodies.size() - 1;

	//torso
	body.m_position = Vec3(0, 0, 4);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeBox(g_boxBody, sizeof(g_boxBody) / sizeof(Vec3));
	body.m_invMass = 0.5f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxBody = m_bodies.size() - 1;

	//left arm
	body.m_position = Vec3(0, 2.0f, 4.75f);
	body.m_orientation = Quat(Vec3(0, 0, 1), -3.1415f / 2.0f);
	body.m_shape = new ShapeBox(g_boxLimb, sizeof(g_boxLimb) / sizeof(Vec3));
	body.m_invMass = 1.0f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxArmLeft = m_bodies.size() - 1;

	//right arm
	body.m_position = Vec3(0, -2.0f, 4.75f);
	body.m_orientation = Quat(Vec3(0, 0, 1), 3.1415f / 2.0f);
	body.m_shape = new ShapeBox(g_boxLimb, sizeof(g_boxLimb) / sizeof(Vec3));
	body.m_invMass = 1.0f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxArmRight = m_bodies.size() - 1;

	//left leg
	body.m_position = Vec3(0, 1.0f, 2.5f);
	body.m_orientation = Quat(Vec3(0, 1, 0), -3.1415f / 2.0f);
	body.m_shape = new ShapeBox(g_boxLimb, sizeof(g_boxLimb) / sizeof(Vec3));
	body.m_invMass = 1.0f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxLegLeft = m_bodies.size() - 1;

	//right leg
	body.m_position = Vec3(0, -1.0f, 2.5f);
	body.m_orientation = Quat(Vec3(0, 1, 0), 3.1415f / 2.0f);
	body.m_shape = new ShapeBox(g_boxLimb, sizeof(g_boxLimb) / sizeof(Vec3));
	body.m_invMass = 1.0f;
	body.m_elasticity = 1.0f;
	body.m_friction = 1.0f;
	m_bodies.push_back(body);
	const int idxLegRight = m_bodies.size() - 1;

	{
		//neck
		ConstraintHingeQuatLimited* joint = new ConstraintHingeQuatLimited();
		joint->m_bodyA = &m_bodies[idxHead];
		joint->m_bodyB = &m_bodies[idxBody];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyA->m_position + Vec3(0, 0, -0.5f);
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		//Set the initial relative orientation
		joint->m_axisA = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, 1, 0));

		m_constraints.push_back(joint);
	}

	{
		//left shoulder
		ConstraintConstantVelocityLimited* joint = new ConstraintConstantVelocityLimited();
		joint->m_bodyB = &m_bodies[idxArmLeft];
		joint->m_bodyA = &m_bodies[idxBody];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyB->m_position + Vec3(0, -1.0f, 0);
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		//Set the initial relative orientation
		joint->m_axisA = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, 1, 0));

		m_constraints.push_back(joint);
	}

	{
		//right shoulder
		ConstraintConstantVelocityLimited* joint = new ConstraintConstantVelocityLimited();
		joint->m_bodyB = &m_bodies[idxArmRight];
		joint->m_bodyA = &m_bodies[idxBody];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyB->m_position + Vec3(0, 1.0f, 0);
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		//Set the initial relative orientation
		joint->m_axisA = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, -1, 0));

		m_constraints.push_back(joint);
	}

	{
		//left hip
		ConstraintHingeQuatLimited* joint = new ConstraintHingeQuatLimited();
		joint->m_bodyB = &m_bodies[idxLegLeft];
		joint->m_bodyA = &m_bodies[idxBody];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyA->m_position + Vec3(0, 0, 0.5f);
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		//Set the initial relative orientation
		joint->m_axisA = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, 1, 0));

		m_constraints.push_back(joint);
	}

	{
		//right hip
		ConstraintHingeQuatLimited* joint = new ConstraintHingeQuatLimited();
		joint->m_bodyB = &m_bodies[idxLegRight];
		joint->m_bodyA = &m_bodies[idxBody];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyA->m_position + Vec3(0, 0, 0.5f);
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		//Set the initial relative orientation
		joint->m_axisA = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, 1, 0));

		m_constraints.push_back(joint);
	}
	*/
	/*
	//
	//Motor
	//
	body.m_position = Vec3(0, 0, 2);
	body.m_linearVelocity = Vec3(0.0f);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeBox(g_boxSmall, sizeof(g_boxSmall) / sizeof(Vec3));
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.9f;
	body.m_friction = 0.5f;
	m_bodies.push_back(body);

	body.m_position = Vec3(0, 0, 1);
	body.m_orientation = Quat(1, 0, 0, 0);
	body.m_shape = new ShapeBox(g_boxBeam, sizeof(g_boxBeam) / sizeof(Vec3));
	body.m_invMass = 0.01f;
	body.m_elasticity = 1.0f;
	body.m_friction = 0.5f;
	m_bodies.push_back(body);

	{
		ConstraintMotor* joint = new ConstraintMotor();
		joint->m_bodyA = &m_bodies[m_bodies.size() - 2];
		joint->m_bodyB = &m_bodies[m_bodies.size() - 1];

		const Vec3 jointWorldSpaceAnchor = joint->m_bodyA->m_position;
		joint->m_anchorA = joint->m_bodyA->WorldSpaceToBodySpace(jointWorldSpaceAnchor);
		joint->m_anchorB = joint->m_bodyB->WorldSpaceToBodySpace(jointWorldSpaceAnchor);

		joint->m_motorSpeed = 2.0f;
		joint->m_motorAxis = joint->m_bodyA->m_orientation.Inverse().RotatePoint(Vec3(0, 0, 1).Normalize());
		//Set the initial relative orientation
		joint->m_q0 = joint->m_bodyA->m_orientation.Inverse() * joint->m_bodyB->m_orientation;

		m_constraints.push_back(joint);
	}

	body.m_position = Vec3(2, 1, 2);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeSphere(1.0f);
	body.m_invMass = 1.0f;
	body.m_elasticity = 0.1f;
	body.m_friction = 0.9f;
	m_bodies.push_back(body);
	*/

	//
	//Mover
	//
	body.m_position = Vec3(10, 0, 5);
	body.m_linearVelocity = Vec3(0, 0, 0);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeBox(g_boxPlatform, sizeof(g_boxPlatform) / sizeof(Vec3));
	body.m_invMass = 0.0f;
	body.m_elasticity = 0.1f;
	body.m_friction = 0.9f;
	m_bodies.push_back(body);

	{
		ConstraintMoverSimple* mover = new ConstraintMoverSimple();
		mover->m_bodyA = &m_bodies[m_bodies.size() - 1];

		m_constraints.push_back(mover);
	}


	body.m_position = Vec3(10, 0, 7);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_shape = new ShapeSphere(1.0f);
	body.m_invMass = 1.0f;
	body.m_elasticity = 0.1f;
	body.m_friction = 0.9f;
	m_bodies.push_back(body);

	AddStandardSandBox(m_bodies);
}

int CompareContacts(const void* p1, const void* p2) {
	contact_t a = *(contact_t*)p1;
	contact_t b = *(contact_t*)p2;

	if (a.timeOfImpact < b.timeOfImpact) {
		return -1;
	}

	if (a.timeOfImpact == b.timeOfImpact) {
		return 0;
	}
	return 1;
}

/*
====================================================
Scene::Update
====================================================
*/
void Scene::Update( const float dt_sec ) {
	//std::vector<ConstraintPenetration> penetrationConstraints;//temporary list most constraints won't last the entire duration and there are too many
	m_manifolds.RemoveExpired();

	//gravity
	for (int i = 0; i < m_bodies.size(); i++) {
		Body* body = &m_bodies[i];
		float mass = 1.0f / body->m_invMass;
		Vec3 impulseGravity = Vec3(0, 0, -10) * mass * dt_sec;
		body->ApplyImpulseLinear(impulseGravity);
	}

	//Broadphase to retain only the possible collsions to reduce collision calculation in the narrow phase
	std::vector< collisionPair_t > collisionPairs;
	BroadPhase(m_bodies.data(), (int) m_bodies.size(), collisionPairs, dt_sec);

	//narrow phase, where general contact and collision are calculated
	int numContacts = 0;
	const int maxContacts = (int) m_bodies.size() * (int) m_bodies.size();
	contact_t* contacts = (contact_t*)alloca(sizeof(contact_t) * maxContacts);
	for (int i = 0; i < collisionPairs.size(); i++) {
		const collisionPair_t& pair = collisionPairs[i];
		Body* bodyA = &m_bodies[pair.a];
		Body* bodyB = &m_bodies[pair.b];

		//skip infinite mass pairs
		if (0.0f == bodyA->m_invMass && 0.0f == bodyB->m_invMass) {
			continue;
		}

		contact_t contact;
		if (Intersect(bodyA, bodyB, dt_sec, contact)) {
			if (contact.timeOfImpact == 0.0f) {//static, hence check penetration using constraint
				/*
				ConstraintPenetration constraint;
				constraint.m_bodyA = contact.bodyA;
				constraint.m_bodyB = contact.bodyB;

				constraint.m_anchorA = contact.ptOnA_LocalSpace;
				constraint.m_anchorB = contact.ptOnB_LocalSpace;

				constraint.m_normal = constraint.m_bodyA->m_orientation.Inverse().RotatePoint(contact.normal * -1.0f);
				constraint.m_normal.Normalize();

				penetrationConstraints.push_back(constraint);
				*/
				m_manifolds.AddContact(contact);
			}
			else {//use dynamic contact checking
				contacts[numContacts] = contact;
				numContacts++;
			}
		}
	}

	//sort the times of impact from first to last
	if (numContacts > 1) {
		qsort(contacts, numContacts, sizeof(contact_t), CompareContacts);
	}

	//
	//solve constraints for each constraint
	//
	for (int i = 0; i < m_constraints.size(); i++) {
		m_constraints[i]->PreSolve(dt_sec);
	}
	m_manifolds.PreSolve(dt_sec);

	const int numIter = 5;
	for (int j = 0; j < numIter; j++) {
		for (int i = 0; i < m_constraints.size(); i++) {
			m_constraints[i]->Solve();
		}
		m_manifolds.Solve();
	}

	for (int i = 0; i < m_constraints.size(); i++) {
		m_constraints[i]->PostSolve();
	}
	m_manifolds.PostSolve();

	float accumulatedTime = 0.0f;
	for (int i = 0; i < numContacts; i++) {
		contact_t& contact = contacts[i];
		const float dt = contact.timeOfImpact - accumulatedTime;

		//Update the position up to next collision
		for (int j = 0; j < m_bodies.size(); j++) {
			m_bodies[j].Update(dt);
		}
		ResolveContact(contact);
		accumulatedTime += dt;
	}

	//Update the position for the rest of the frame's time
	const float timeRemaining = dt_sec - accumulatedTime;
	if (timeRemaining > 0.0f) {
		for (int i = 0; i < m_bodies.size(); i++) {
			m_bodies[i].Update(timeRemaining);
		}
	}
}