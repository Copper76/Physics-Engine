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
	
	body.m_position = Vec3(0, 0, 10);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_invMass = 1.0f;
	body.m_linearVelocity.Zero();
	body.m_elasticity = 1.0f;
	body.m_friction = 0.5f;
	body.m_shape = new ShapeSphere(0.5f);
	m_bodies.push_back(body);
	
	body.m_position = Vec3(0, 0, -1000);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_invMass = 0.0f;
	body.m_linearVelocity.Zero();
	body.m_elasticity = 1.0f;
	body.m_friction = 0.5f;
	body.m_shape = new ShapeSphere(1000.0f);
	m_bodies.push_back(body);
	*/

	body.m_position = Vec3(10, 0, 3);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity = Vec3(-100, 0, 0);
	body.m_angularVelocity.Zero();
	body.m_invMass = 1.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.5f;
	body.m_shape = new ShapeSphere(0.5f);
	m_bodies.push_back(body);

	body.m_position = Vec3(-10, 0, 3);
	body.m_orientation = Quat(0, 0, 0, 1);
	body.m_linearVelocity = Vec3(100,0,0);
	body.m_angularVelocity.Zero();
	body.m_invMass = 1.0f;
	body.m_elasticity = 0.5f;
	body.m_friction = 0.5f;
	body.m_shape = new ShapeConvex(g_diamond, sizeof(g_diamond)/sizeof(Vec3));
	//body.m_shape = new ShapeBox(g_boxBody, sizeof(g_boxBody) / sizeof(Vec3));
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
			contacts[numContacts] = contact;
			numContacts++;
		}
	}

	//sort the times of impact from first to last
	if (numContacts > 1) {
		qsort(contacts, numContacts, sizeof(contact_t), CompareContacts);
	}

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