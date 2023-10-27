//
//  ConstraintMover.cpp
//
#include "ConstraintMover.h"

/*
====================================================
ConstraintMoverSimple::PreSolve
====================================================
*/
void ConstraintMoverSimple::PreSolve( const float dt_sec ) {
	m_time += dt_sec;
	m_bodyA->m_linearVelocity.z = cosf(m_time * 0.25f) * 4.0f;//a platform that moves up and down in a wave fashion
}