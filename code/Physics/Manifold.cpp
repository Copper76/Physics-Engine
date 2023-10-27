//
//  Manifold.cpp
//
#include "Manifold.h"


/*
================================================================================================

ManifoldCollector

================================================================================================
*/

/*
================================
ManifoldCollector::AddContact
================================
*/
void ManifoldCollector::AddContact( const contact_t & contact ) {
	int foundIdx = -1;
	for (int i = 0; i < m_manifolds.size(); i++) {
		if ((m_manifolds[i].m_bodyA == contact.bodyA && m_manifolds[i].m_bodyB == contact.bodyB) || (m_manifolds[i].m_bodyB == contact.bodyA && m_manifolds[i].m_bodyA == contact.bodyB)) {
			foundIdx = i;
			break;
		}
	}

	if (foundIdx >= 0) {
		m_manifolds[foundIdx].AddContact(contact);
	}
	else {
		Manifold manifold;
		manifold.m_bodyA = contact.bodyA;
		manifold.m_bodyB = contact.bodyB;

		manifold.AddContact(contact);
		m_manifolds.push_back(manifold);
	}
}

/*
================================
ManifoldCollector::RemoveExpired
================================
*/
void ManifoldCollector::RemoveExpired() {
	for (int i = (int)m_manifolds.size() - 1; i >= 0; i--) {
		m_manifolds[i].RemoveExpiredContacts();

		if (m_manifolds[i].m_numContacts == 0) {
			m_manifolds.erase(m_manifolds.begin() + i);
		}
	}
}

/*
================================
ManifoldCollector::PreSolve
================================
*/
void ManifoldCollector::PreSolve( const float dt_sec ) {
	for (int i = 0; i < m_manifolds.size(); i++) {
		m_manifolds[i].PreSolve(dt_sec);
	}
}

/*
================================
ManifoldCollector::Solve
================================
*/
void ManifoldCollector::Solve() {
	for (int i = 0; i < m_manifolds.size(); i++) {
		m_manifolds[i].Solve();
	}
}

/*
================================
Manifold::PostSolve
================================
*/
void ManifoldCollector::PostSolve() {
	for (int i = 0; i < m_manifolds.size(); i++) {
		m_manifolds[i].PostSolve();
	}
}

/*
================================================================================================

Manifold

================================================================================================
*/

/*
================================
Manifold::RemoveExpiredContacts
================================
*/
void Manifold::RemoveExpiredContacts() {
	for (int i = 0; i < m_numContacts; i++) {
		contact_t& contact = m_contacts[i];

		Body* bodyA = contact.bodyA;
		Body* bodyB = contact.bodyB;

		const Vec3 a = bodyA->BodySpaceToWorldSpace(contact.ptOnA_LocalSpace);
		const Vec3 b = bodyB->BodySpaceToWorldSpace(contact.ptOnB_LocalSpace);

		Vec3 normal = m_constraints[i].m_normal;//always length 1
		normal = bodyA->m_orientation.RotatePoint(normal);

		const Vec3 ab = b - a;
		float penetrationDepth = normal.Dot(ab);
		Vec3 abNormal = normal * penetrationDepth;
		Vec3 abTangent = ab - abNormal;

		const float distanceThreshold = 0.02f;//a slack on the distance. We can keep contacts that are close enough.
		if (abTangent.GetLengthSqr() < distanceThreshold * distanceThreshold && penetrationDepth <= 0) {//distance is not very far and not penetrating
			continue;
		}

		//book definition
		
		for (int j = i; j < MAX_CONTACTS - 1; j++) {
			m_constraints[j] = m_constraints[j + 1];
			m_contacts[j] = m_contacts[j + 1];
			/*
			if (j >= m_numContacts) {
				m_constraints[j].m_cachedLambda.Zero();//zero out lambdas of unused constraints? why?
			}
			*/
		}
		

		//alternative that only moves the last element
		/*
		m_constraints[i] = m_constraints[m_numContacts - 1];
		m_contacts[i] = m_contacts[m_numContacts - 1];
		*/
		m_numContacts--;
		i--;
	}
	
	for (int j = m_numContacts; j < MAX_CONTACTS - 1; j++) {
		m_constraints[j].m_cachedLambda.Zero();//zero out lambdas of unused constraints, done only once
	}
	
}

/*
================================
Manifold::AddContact
================================
*/
void Manifold::AddContact( const contact_t & contact_old ) {
	contact_t contact = contact_old;

	//fix the contact so bodies align with manifold bodies
	if (contact_old.bodyA != m_bodyA || contact_old.bodyB != m_bodyB) {//only one is necessary as long as other codes are correct
		contact.ptOnA_LocalSpace = contact_old.ptOnB_LocalSpace;
		contact.ptOnB_LocalSpace = contact_old.ptOnA_LocalSpace;
		contact.ptOnA_WorldSpace = contact_old.ptOnB_WorldSpace;
		contact.ptOnB_WorldSpace = contact_old.ptOnA_WorldSpace;

		contact.bodyA = m_bodyA;
		contact.bodyB = m_bodyB;
	}

	for (int i = 0; i < m_numContacts; i++) {
		const Body* bodyA = m_contacts[i].bodyA;
		const Body* bodyB = m_contacts[i].bodyB;

		const Vec3 oldA = bodyA->BodySpaceToWorldSpace(m_contacts[i].ptOnA_LocalSpace);
		const Vec3 oldB = bodyB->BodySpaceToWorldSpace(m_contacts[i].ptOnB_LocalSpace);

		const Vec3 newA = contact.bodyA->BodySpaceToWorldSpace(contact.ptOnA_LocalSpace);
		const Vec3 newB = contact.bodyB->BodySpaceToWorldSpace(contact.ptOnB_LocalSpace);

		const float distanceThreshold = 0.02f;
		if (((newA - oldA).GetLengthSqr() < distanceThreshold * distanceThreshold) || ((newB - oldB).GetLengthSqr() < distanceThreshold * distanceThreshold)) {
			return;
		}
	}

	//check what to discard when we are full on contacts
	int newSlot = m_numContacts;

	//book definition
	if (newSlot >= MAX_CONTACTS) {
		Vec3 avg = Vec3(0.0f);
		for (int i = 0; i < MAX_CONTACTS; i++) {
			avg += m_contacts[i].ptOnA_LocalSpace;
		}
		avg += contact.ptOnA_LocalSpace;
		avg *= 1.0f / (1+MAX_CONTACTS);

		float minDist = (avg - contact.ptOnA_LocalSpace).GetLengthSqr();
		int removeIdx = -1;

		for (int i = 0; i < MAX_CONTACTS; i++) {
			float dist2 = (avg - m_contacts[i].ptOnA_LocalSpace).GetLengthSqr();

			if (dist2 < minDist) {
				minDist = dist2;
				removeIdx = i;
			}
		}

		if (removeIdx != -1) {
			newSlot = removeIdx;
		}
		else {
			return;//don't need to add the new point
		}
	}

	/*
	if (newSlot >= MAX_CONTACTS) {//try implement a new method to find the most penetrating contact and remove the one that is furthest away
	}
	*/

	m_contacts[newSlot] = contact;

	m_constraints[newSlot].m_bodyA = contact.bodyA;
	m_constraints[newSlot].m_bodyB = contact.bodyB;
	m_constraints[newSlot].m_anchorA = contact.ptOnA_LocalSpace;
	m_constraints[newSlot].m_anchorB = contact.ptOnB_LocalSpace;

	//get the normal for constraint
	m_constraints[newSlot].m_normal = m_bodyA->m_orientation.Inverse().RotatePoint(contact.normal * -1.0f);
	m_constraints[newSlot].m_normal.Normalize();

	m_constraints[newSlot].m_cachedLambda.Zero();//replaced constraint so no cached lambda here

	if (newSlot == m_numContacts) {//it was a new contact, not a replacement
		m_numContacts++;
	}
}

/*
================================
Manifold::PreSolve
================================
*/
void Manifold::PreSolve( const float dt_sec ) {
	for (int i = 0; i < m_numContacts; i++) {
		m_constraints[i].PreSolve(dt_sec);
	}
}

/*
================================
Manifold::Solve
================================
*/
void Manifold::Solve() {
	for (int i = 0; i < m_numContacts; i++) {
		m_constraints[i].Solve();
	}
}

/*
================================
Manifold::PostSolve
================================
*/
void Manifold::PostSolve() {
	for (int i = 0; i < m_numContacts; i++) {
		m_constraints[i].PostSolve();
	}
}