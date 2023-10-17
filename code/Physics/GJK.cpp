//
//  GJK.cpp
//
#include "GJK.h"

struct point_t;
float EPA_Expand(const Body* bodyA, const Body* bodyB, const float bias, const point_t simplexPoints[4], Vec3& ptOnA, Vec3& ptOnB);


Vec2 SignedVolume1D(const Vec3& s1, const Vec3& s2) {
	Vec3 ab = s2 - s1;
	Vec3 ap = Vec3(0.0f) - s1;
	Vec3 p0 = s1 + ab * ab.Dot(ap)/ab.GetLengthSqr(); // project origin onto the line. Huh?

	int idx = 0;
	float mu_max = 0;
	for (int i = 0; i < 3; i++) {
		float mu = s2[i] - s1[i];
		if (mu * mu > mu_max * mu_max) {
			mu_max = mu;
			idx = i;
		}
	}

	const float a = s1[idx];
	const float b = s2[idx];
	const float p = p0[idx];

	const float C1 = p - a;
	const float C2 = b - p;

	if ((p > a && p < b) || (p<a && p>b)) {
		Vec2 lambdas;
		lambdas[0] = C2 / mu_max;
		lambdas[1] = C1 / mu_max;
		return lambdas;
	}

	if ((a <= b && p <= a) || (a>=b && p>=a)) {
		return Vec2(1.0f, 0.0f);
	}

	return Vec2(0.0f, 1.0f);
}

int CompareSigns(float a, float b) {
	if (a > 0.0f && b > 0.0f) {
		return 1;
	}

	if (a < 0.0f && b < 0.0f) {
		return 1;
	}

	return 0;
}

Vec3 SignedVolume2D(const Vec3& s1, const Vec3& s2, const Vec3& s3) {
	Vec3 normal = (s2 - s1).Cross(s3 - s1);
	Vec3 p0 = normal * s1.Dot(normal) / normal.GetLengthSqr();

	int idx = 0;
	float area_max = 0;

	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		Vec2 a = Vec2(s1[j], s1[k]);
		Vec2 b = Vec2(s2[j], s2[k]);
		Vec2 c = Vec2(s3[j], s3[k]);
		Vec2 ab = b - a;
		Vec2 ac = c - a;

		float area = ab.x * ac.y - ab.y * ac.x;
		if (area * area > area_max * area_max) {
			idx = i;
			area_max = area;
		}
	}

	int x = (idx + 1) % 3;
	int y = (idx + 2) % 3;
	Vec2 s[3];
	s[0] = Vec2(s1[x], s1[y]);
	s[1] = Vec2(s2[x], s2[y]);
	s[2] = Vec2(s3[x], s3[y]);
	Vec2 p = Vec2(p0[x], p0[y]);

	Vec3 areas;
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		Vec2 a = p;
		Vec2 b = s[j];
		Vec2 c = s[k];
		Vec2 ab = b - a;
		Vec2 ac = c - a;

		areas[i] = ab.x * ac.y - ab.y * ac.x;
	}

	if (CompareSigns(area_max, areas[0]) > 0 && CompareSigns(area_max, areas[1]) > 0 && CompareSigns(area_max, areas[2]) > 0) {
		return areas / area_max;
	}

	float dist = 1e10;//manual maximum
	Vec3 lambdas = Vec3(1, 0, 0);
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		Vec3 edgesPts[3];
		edgesPts[0] = s1;
		edgesPts[1] = s2;
		edgesPts[2] = s3;

		Vec2 lambdaEdge = SignedVolume1D(edgesPts[j], edgesPts[k]);
		Vec3 pt = edgesPts[j] * lambdaEdge[0] + edgesPts[k] * lambdaEdge[1];
		if (pt.GetLengthSqr() < dist) {
			dist = pt.GetLengthSqr();
			lambdas[i] = 0;
			lambdas[j] = lambdaEdge[0];
			lambdas[k] = lambdaEdge[1];
		}
	}

	return lambdas;
}

Vec4 SignedVolume3D(const Vec3& s1, const Vec3& s2, const Vec3& s3, const Vec3& s4) {
	Mat4 M;
	M.rows[0] = Vec4(s1.x, s2.x, s3.x, s4.x);
	M.rows[1] = Vec4(s1.y, s2.y, s3.y, s4.y);
	M.rows[2] = Vec4(s1.z, s2.z, s3.z, s4.z);
	M.rows[3] = Vec4(1.0f, 1.0f, 1.0f, 1.0f);

	Vec4 C4;
	C4[0] = M.Cofactor(3, 0);
	C4[1] = M.Cofactor(3, 1);
	C4[2] = M.Cofactor(3, 2);
	C4[3] = M.Cofactor(3, 3);

	const float detM = C4[0] + C4[1] + C4[2] + C4[3];

	if (CompareSigns(detM, C4[0]) > 0 && CompareSigns(detM, C4[1]) > 0 
		&& CompareSigns(detM, C4[2]) > 0 && CompareSigns(detM, C4[3]) > 0) {
		return C4 * (1.0f / detM);
	}

	Vec4 lambdas;
	float dist = 1e10;//manual maximum
	for (int i = 0; i < 4; i++) {
		int j = (i + 1) % 4;
		int k = (i + 2) % 4;

		Vec3 facePts[4];
		facePts[0] = s1;
		facePts[1] = s2;
		facePts[2] = s3;
		facePts[3] = s4;

		Vec3 lambdaFace = SignedVolume2D(facePts[i], facePts[j], facePts[k]);
		Vec3 pt = facePts[i] * lambdaFace[0] + facePts[j] * lambdaFace[1] + facePts[k] * lambdaFace[2];
		if (pt.GetLengthSqr() < dist) {
			dist = pt.GetLengthSqr();
			lambdas.Zero();
			lambdas[i] = lambdaFace[0];
			lambdas[j] = lambdaFace[1];
			lambdas[k] = lambdaFace[2];
		}
	}

	return lambdas;
}

void TestSignedVolumeProjection() {
	const Vec3 orgPts[4] = {
		Vec3(0,0,0),
		Vec3(1,0,0),
		Vec3(0,1,0),
		Vec3(0,0,1)
	};
	Vec3 pts[4];
	Vec4 lambdas;
	Vec3 v;

	for (int i = 0; i < 4; i++) {
		pts[i] = orgPts[i] + Vec3(1, 1, 1);
	}
	lambdas = SignedVolume3D(pts[0], pts[1], pts[2], pts[3]);
	v.Zero();
	for (int i = 0; i < 4; i++) {
		v += pts[i] * lambdas[i];
	}
	printf("lambdas: %.3f %.3f %.3f %.3f        v: %.3f %.3f %.3f\n",
		lambdas.x, lambdas.y, lambdas.z, lambdas.w,
		v.x, v.y, v.z);

	for (int i = 0; i < 4; i++) {
		pts[i] = orgPts[i] + Vec3(-1, -1, -1) * 0.25f;
	}
	lambdas = SignedVolume3D(pts[0], pts[1], pts[2], pts[3]);
	v.Zero();
	for (int i = 0; i < 4; i++) {
		v += pts[i] * lambdas[i];
	}
	printf("lambdas: %.3f %.3f %.3f %.3f        v: %.3f %.3f %.3f\n",
		lambdas.x, lambdas.y, lambdas.z, lambdas.w,
		v.x, v.y, v.z);

	for (int i = 0; i < 4; i++) {
		pts[i] = orgPts[i] + Vec3(-1, -1, -1);
	}
	lambdas = SignedVolume3D(pts[0], pts[1], pts[2], pts[3]);
	v.Zero();
	for (int i = 0; i < 4; i++) {
		v += pts[i] * lambdas[i];
	}
	printf("lambdas: %.3f %.3f %.3f %.3f        v: %.3f %.3f %.3f\n",
		lambdas.x, lambdas.y, lambdas.z, lambdas.w,
		v.x, v.y, v.z);

	for (int i = 0; i < 4; i++) {
		pts[i] = orgPts[i] + Vec3(1, 1, -0.5f);
	}
	lambdas = SignedVolume3D(pts[0], pts[1], pts[2], pts[3]);
	v.Zero();
	for (int i = 0; i < 4; i++) {
		v += pts[i] * lambdas[i];
	}
	printf("lambdas: %.3f %.3f %.3f %.3f        v: %.3f %.3f %.3f\n",
		lambdas.x, lambdas.y, lambdas.z, lambdas.w,
		v.x, v.y, v.z);

	pts[0] = Vec3(51.1996613f, 26.1989613f, 1.91339576f);
	pts[1] = Vec3(-51.0567360f, -26.0565681f, -0.436143428f);
	pts[2] = Vec3(50.8978920f, -24.1035538f, -1.04042661f);
	pts[3] = Vec3(-49.1021080f, 25.8964462f, -1.04042661f);

	lambdas = SignedVolume3D(pts[0], pts[1], pts[2], pts[3]);
	v.Zero();
	for (int i = 0; i < 4; i++) {
		v += pts[i] * lambdas[i];
	}
	printf("lambdas: %.3f %.3f %.3f %.3f        v: %.3f %.3f %.3f\n",
		lambdas.x, lambdas.y, lambdas.z, lambdas.w,
		v.x, v.y, v.z);
}

struct point_t {
	Vec3 xyz; //point on minkowski sum
	Vec3 ptA; //point on bodyA
	Vec3 ptB; //point on bodyB

	point_t() : xyz(0.0f), ptA(0.0f), ptB(0.0f) {}

	const point_t& operator= (const point_t& rhs) {
		xyz = rhs.xyz;
		ptA = rhs.ptA;
		ptB = rhs.ptB;
		return *this;
	}

	bool operator== (const point_t& rhs) const {
		return (ptA == rhs.ptA) && (ptB == rhs.ptB) && (xyz == rhs.xyz);
	}
};

point_t Support(const Body* bodyA, const Body* bodyB, Vec3 dir, const float bias) {
	dir.Normalize();//extra normalisation?

	point_t point;

	point.ptA = bodyA->m_shape->Support(dir, bodyA->m_position, bodyA->m_orientation, bias);

	dir *= -1.0f;

	point.ptB = bodyB->m_shape->Support(dir, bodyB->m_position, bodyB->m_orientation, bias);

	point.xyz = point.ptA - point.ptB;//A-B huh
	return point;
}

bool SimplexSignedVolumes(point_t* pts, const int num, Vec3& newDir, Vec4& lambdasOut) {
	const float epsilonf = 0.0001f * 0.0001f;
	lambdasOut.Zero();

	bool doesIntersect = false;
	switch (num) {
		default:
		case 2: {
			Vec2 lambdas = SignedVolume1D(pts[0].xyz, pts[1].xyz);
			Vec3 v(0.0f);
			for(int i = 0; i < 2; i++) {
				v += pts[i].xyz * lambdas[i];
			}
			newDir = v * -1.0f;
			doesIntersect = (v.GetLengthSqr() < epsilonf);
			lambdasOut[0] = lambdas[0];
			lambdasOut[1] = lambdas[1];
		} break;
		case 3: {
			Vec3 lambdas = SignedVolume2D(pts[0].xyz, pts[1].xyz, pts[2].xyz);
			Vec3 v(0.0f);
			for (int i = 0; i < 3; i++) {
				v += pts[i].xyz * lambdas[i];
			}
			newDir = v * -1.0f;
			doesIntersect = (v.GetLengthSqr() < epsilonf);
			lambdasOut[0] = lambdas[0];
			lambdasOut[1] = lambdas[1];
			lambdasOut[2] = lambdas[2];
		} break;
		case 4: {
			Vec4 lambdas = SignedVolume3D(pts[0].xyz, pts[1].xyz, pts[2].xyz, pts[3].xyz);
			Vec3 v(0.0f);
			for (int i = 0; i < 4; i++) {
				v += pts[i].xyz * lambdas[i];
			}
			newDir = v * -1.0f;
			doesIntersect = (v.GetLengthSqr() < epsilonf);
			lambdasOut[0] = lambdas[0];
			lambdasOut[1] = lambdas[1];
			lambdasOut[2] = lambdas[2];
			lambdasOut[3] = lambdas[3];
		} break;
	}

	return doesIntersect;
}

bool HasPoint(const point_t simplexPoints[4], const point_t& newPt) {
	const float precision = 1e-6f;//a good precision for approximation

	for (int i = 0; i < 4; i++) {
		Vec3 delta = simplexPoints[i].xyz - newPt.xyz;
		if (delta.GetLengthSqr() < precision * precision) {
			return true;
		}
	}
	return false;
}

void SortValids(point_t simplexPoints[4], Vec4& lambdas) {
	bool valids[4];
	//record the valid indexes
	for (int i = 0; i < 4; i++) {
		valids[i] = (lambdas[i] != 0);
	}

	Vec4 validLambdas(0.0f);
	int validCount = 0;
	point_t validPts[4];
	memset(validPts, 0, sizeof(point_t) * 4);//set the array to 0, is that the best way to do it?
	for (int i = 0; i < 4; i++) {
		if (valids[i]) {
			validPts[validCount] = simplexPoints[i];
			validLambdas[validCount] = lambdas[i];
			validCount++;
		}
	}

	for (int i = 0; i < 4; i++) {
		simplexPoints[i] = validPts[i];
		lambdas[i] = validLambdas[i];
	}
}

static int NumValids(const Vec4& lambdas) {
	int num = 0;
	for (int i = 0; i < 4; i++) {
		if (lambdas[i] != 0.0f) {
			num++;
		}
	}
	return num;
}

/*
================================
GJK_DoesIntersect
================================
*/
bool GJK_DoesIntersect(const Body* bodyA, const Body* bodyB) {
	const Vec3 origin(0.0f);

	int numPts = 1;
	point_t simplexPoints[4];
	simplexPoints[0] = Support(bodyA, bodyB, Vec3(1, 1, 1), 0.0f);

	float closestDist = 1e10f;//manual unlikely upperbound, like python's math.inf
	bool doesContainOrigin = false;
	Vec3 newDir = simplexPoints[0].xyz * -1.0f;

	//iteratively fill the tetrahedron
	do {
		point_t newPt = Support(bodyA, bodyB, newDir, 0.0f);

		//if the new point has been found before, we cannot expand further
		if (HasPoint(simplexPoints, newPt)) {
			break;
		}

		//add the new point to simplex points
		simplexPoints[numPts] = newPt;
		numPts++;

		//if the new point(furthest point) is not on the other side of the origin, then we will never include origin so break
		if (newDir.Dot(newPt.xyz - origin) < 0.0f) {
			break;
		}

		Vec4 lambdas;
		doesContainOrigin = SimplexSignedVolumes(simplexPoints, numPts, newDir, lambdas);
		if (doesContainOrigin) {
			break;
		}

		//The new projection is further from the origin point hence we cannot improve further
		float dist = newDir.GetLengthSqr();
		if (dist >= closestDist) {
			break;
		}
		closestDist = dist;

		SortValids(simplexPoints, lambdas);
		numPts = NumValids(lambdas);
		doesContainOrigin = (numPts == 4);
	} while (!doesContainOrigin);

	return doesContainOrigin;
}

/*
================================
GJK_ClosestPoints
================================
*/
void GJK_ClosestPoints(const Body* bodyA, const Body* bodyB, Vec3& ptOnA, Vec3& ptOnB) {
	const Vec3 origin(0.0f);

	float closestDist = 1e10f;
	const float bias = 0.0f;

	int numPts = 1;
	point_t simplexPoints[4];
	simplexPoints[0] = Support(bodyA, bodyB, Vec3(1, 1, 1), bias);

	Vec4 lambdas = Vec4(1, 0, 0, 0);
	Vec3 newDir = simplexPoints[0].xyz * -1.0f;

	do {
		point_t newPt = Support(bodyA, bodyB, newDir, bias);

		if (HasPoint(simplexPoints, newPt)) {
			break;
		}

		simplexPoints[numPts] = newPt;
		numPts++;

		SimplexSignedVolumes(simplexPoints, numPts, newDir, lambdas);
		SortValids(simplexPoints, lambdas);
		numPts = NumValids(lambdas);

		float dist = newDir.GetLengthSqr();
		if (dist >= closestDist) {
			break;
		}
		closestDist = dist;
	} while (numPts < 4);

	//calculate the points on A and B with projection
	ptOnA.Zero();
	ptOnB.Zero();
	for (int i = 0; i < 4; i++) {
		ptOnA += simplexPoints[i].ptA * lambdas[i];
		ptOnB += simplexPoints[i].ptB * lambdas[i];
	}
}

/*
================================
GJK_DoesIntersect
================================
*/
bool GJK_DoesIntersect(const Body* bodyA, const Body* bodyB, const float bias, Vec3& ptOnA, Vec3& ptOnB) {
	const Vec3 origin(0.0f);

	int numPts = 1;
	point_t simplexPoints[4];
	simplexPoints[0] = Support(bodyA, bodyB, Vec3(1, 1, 1), 0.0f);

	float closestDist = 1e10f;//manual unlikely upperbound, like python's math.inf
	bool doesContainOrigin = false;
	Vec3 newDir = simplexPoints[0].xyz * -1.0f;

	//iteratively fill the tetrahedron
	do {
		point_t newPt = Support(bodyA, bodyB, newDir, 0.0f);

		//if the new point has been found before, we cannot expand further
		if (HasPoint(simplexPoints, newPt)) {
			break;
		}

		//add the new point to simplex points
		simplexPoints[numPts] = newPt;
		numPts++;

		//if the new point(furthest point) is not on the other side of the origin, then we will never include origin so break
		if (newDir.Dot(newPt.xyz - origin) < 0.0f) {
			break;
		}

		Vec4 lambdas;
		doesContainOrigin = SimplexSignedVolumes(simplexPoints, numPts, newDir, lambdas);
		if (doesContainOrigin) {
			break;
		}

		//The new projection is further from the origin point hence we cannot improve further
		float dist = newDir.GetLengthSqr();
		if (dist >= closestDist) {
			break;
		}
		closestDist = dist;

		SortValids(simplexPoints, lambdas);
		numPts = NumValids(lambdas);
		doesContainOrigin = (numPts == 4);
	} while (!doesContainOrigin);

	if (!doesContainOrigin) {
		return false;
	}

	//makes sure we have a tetrahedron for EPA
	if (numPts == 1) {
		Vec3 searchDir = simplexPoints[0].xyz * -1.0f;
		simplexPoints[numPts] = Support(bodyA, bodyB, searchDir, 0.0f);
		numPts++;
	}
	if (numPts == 2) {
		Vec3 ab = simplexPoints[1].xyz - simplexPoints[0].xyz;
		Vec3 u, v;
		ab.GetOrtho(u, v);
		simplexPoints[numPts] = Support(bodyA, bodyB, u, 0.0f);
		numPts++;
	}
	if (numPts == 3) {
		Vec3 ab = simplexPoints[1].xyz - simplexPoints[0].xyz;
		Vec3 ac = simplexPoints[2].xyz - simplexPoints[0].xyz;
		Vec3 norm = ab.Cross(ac);
		simplexPoints[numPts] = Support(bodyA, bodyB, norm, 0.0f);
		numPts++;
	}

	//get the center of mass for simplex
	Vec3 avg = Vec3(0, 0, 0);
	for (int i = 0; i < 4; i++) {
		avg += simplexPoints[i].xyz;
	}
	avg *= 0.25f;

	//add bias
	for (int i = 0; i < numPts; i++) {
		point_t& pt = simplexPoints[i];

		Vec3 dir = pt.xyz - avg;
		dir.Normalize();
		pt.ptA += dir * bias;
		pt.ptB -= dir * bias;
		pt.xyz = pt.ptA - pt.ptB;
	}

	EPA_Expand(bodyA, bodyB, bias, simplexPoints, ptOnA, ptOnB);
	return true;
}

//find the center of mass for the points
Vec3 BarycentricCoordinates(Vec3 s1, Vec3 s2, Vec3 s3, const Vec3& pt) {
	s1 = s1 - pt;
	s2 = s2 - pt;
	s3 = s3 - pt;

	Vec3 normal = (s2 - s1).Cross(s3 - s1);
	Vec3 p0 = normal * s1.Dot(normal) / normal.GetLengthSqr();

	int idx = 0;
	float area_max = 0;
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		Vec2 a = Vec2(s1[j], s1[k]);
		Vec2 b = Vec2(s2[j], s2[k]);
		Vec2 c = Vec2(s3[j], s3[k]);
		Vec2 ab = b - a;
		Vec2 ac = c - a;

		float area = ab.x * ac.y - ab.y * ac.x;
		if (area * area > area_max * area_max) {
			idx = i;
			area_max = area;
		}
	}

	int x = (idx + 1) % 3;
	int y = (idx + 2) % 3;
	
	Vec2 s[3];
	s[0] = Vec2(s1[x], s1[y]);
	s[1] = Vec2(s2[x], s2[y]);
	s[2] = Vec2(s3[x], s3[y]);
	Vec2 p = Vec2(p0[x], p0[y]);

	Vec3 areas;
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;

		Vec2 a = p;
		Vec2 b = s[j];
		Vec2 c = s[k];
		Vec2 ab = b - a;
		Vec2 ac = c - a;

		areas[i] = ab.x * ac.y - ab.y * ac.x;
	}

	Vec3 lambdas = areas / area_max;
	if (!lambdas.IsValid()) {
		lambdas = Vec3(1, 0, 0);//a bit reckless here oof, but works I guess
	}
	return lambdas;
}

Vec3 NormalDirection(const tri_t& tri, const std::vector<point_t>& points) {
	const Vec3& a = points[tri.a].xyz;
	const Vec3& b = points[tri.b].xyz;
	const Vec3& c = points[tri.c].xyz;

	Vec3 ab = b - a;
	Vec3 ac = c - a;
	Vec3 normal = ab.Cross(ac);
	normal.Normalize();
	return normal;
}

float SignedDistanceToTriangle(const tri_t& tri, const Vec3& pt, const std::vector<point_t>& points) {
	const Vec3 normal = NormalDirection(tri, points);
	const Vec3& a = points[tri.a].xyz;
	const float dist = normal.Dot(pt - a);
	return dist;
}

int ClosestTriangle(const std::vector<tri_t>& triangles, const std::vector<point_t>& points) {
	float minDistSqr = 1e10;//manual unlikely upperbound, like python's math.inf

	int idx = -1;
	for (int i = 0; i < triangles.size(); i++) {
		const tri_t& tri = triangles[i];

		float dist = SignedDistanceToTriangle(tri, Vec3(0.0f), points);
		float distSqr = dist * dist;
		if (distSqr < minDistSqr) {
			idx = i;
			minDistSqr = distSqr;
		}
	}

	return idx;
}

bool HasPoint(const Vec3& w, const std::vector<tri_t> triangles, const std::vector<point_t>& points) {
	const float epsilons = 0.001f * 0.001f;
	Vec3 delta;

	for (int i = 0; i < triangles.size(); i++) {
		const tri_t& tri = triangles[i];

		delta = w - points[tri.a].xyz;
		if (delta.GetLengthSqr() < epsilons) {
			return true;
		}
		delta = w - points[tri.b].xyz;
		if (delta.GetLengthSqr() < epsilons) {
			return true;
		}
		delta = w - points[tri.c].xyz;
		if (delta.GetLengthSqr() < epsilons) {
			return true;
		}
	}
	return false;
}

int RemoveTrianglesFacingPoint(const Vec3& pt, std::vector<tri_t>& triangles, const std::vector<point_t>& points) {
	int numRemoved = 0;
	for (int i = 0; i < triangles.size(); i++) {
		const tri_t& tri = triangles[i];

		float dist = SignedDistanceToTriangle(tri, pt, points);
		if (dist > 0.0f) {//this triangle faces the point
			triangles.erase(triangles.begin() + i);
			i--;
			numRemoved++;
		}
	}
	return numRemoved;
}

void FindDanglingEdges(std::vector<edge_t>& danglingEdges, const std::vector<tri_t>& triangles) {
	danglingEdges.clear();

	for (int i = 0; i < triangles.size(); i++) {
		const tri_t tri = triangles[i];

		edge_t edges[3];
		edges[0].a = tri.a;
		edges[0].b = tri.b;

		edges[1].a = tri.b;
		edges[1].b = tri.c;

		edges[2].a = tri.c;
		edges[2].b = tri.a;

		int counts[3];
		counts[0] = 0;
		counts[1] = 0;
		counts[2] = 0;

		for (int j = 0; j < triangles.size(); j++) {
			if (j == i) {
				continue;
			}

			const tri_t tri2 = triangles[j];

			edge_t edges2[3];
			edges2[0].a = tri2.a;
			edges2[0].b = tri2.b;

			edges2[1].a = tri2.b;
			edges2[1].b = tri2.c;

			edges2[2].a = tri2.c;
			edges2[2].b = tri2.a;

			for (int k = 0; k < 3; k++) {
				if (edges[k] == edges2[0]) {
					counts[k]++;
				}
				if (edges[k] == edges2[1]) {
					counts[k]++;
				}
				if (edges[k] == edges2[2]) {
					counts[k]++;
				}
			}
		}

		for (int k = 0; k < 3; k++) {
			if (counts[k] == 0) {
				danglingEdges.push_back(edges[k]);
			}
		}
	}
}

float EPA_Expand(const Body* bodyA, const Body* bodyB, const float bias, const point_t simplexPoints[4], Vec3& ptOnA, Vec3& ptOnB) {
	std::vector<point_t> points;
	std::vector<tri_t> triangles;
	std::vector<edge_t> danglingEdges;

	Vec3 center(0.0f);
	for (int i = 0; i < 4; i++) {
		points.push_back(simplexPoints[i]);
		center += simplexPoints[i].xyz;
	}
	center *= 0.25f;

	for (int i = 0; i < 4; i++) {
		int j = (i + 1) % 4;
		int k = (i + 2) % 4;

		tri_t tri;
		tri.a = i;
		tri.b = j;
		tri.c = k;

		int unusedPt = (i + 3) % 4;
		float dist = SignedDistanceToTriangle(tri,points[unusedPt].xyz,points);

		if (dist > 0.0f) {
			std::swap(tri.a, tri.b);
		}

		triangles.push_back(tri);
	}

	while (1) {
		const int idx = ClosestTriangle(triangles, points);
		Vec3 normal = NormalDirection(triangles[idx], points);

		const point_t newPt = Support(bodyA, bodyB, normal, bias);

		//point is seen, terminate the loop
		if (HasPoint(newPt.xyz, triangles, points)) {
			break;
		}

		//check if it is on the outside of the triangle, if not terminate
		float dist = SignedDistanceToTriangle(triangles[idx], newPt.xyz, points);
		if (dist <= 0.0f) {
			break;
		}

		const int newIdx = (int)points.size();
		points.push_back(newPt);

		//remove triangles facing the new point, if there is none, we have done expanding
		if (RemoveTrianglesFacingPoint(newPt.xyz, triangles, points) == 0) {
			break;
		}

		//find all the dangling edges, if there is none, we have reached maximum expansion and terminate
		//FROM BOOK: the points should theoretically be CCW, so adding new point as a should make the triangle face away from origin
		danglingEdges.clear();
		FindDanglingEdges(danglingEdges, triangles);
		if (danglingEdges.size() == 0) {
			break;
		}

		//form new triangles using the newly added point and the dangling edges
		for (int i = 0; i < danglingEdges.size(); i++) {
			const edge_t& edge = danglingEdges[i];

			tri_t triangle;
			triangle.a = newIdx;
			triangle.b = edge.b;
			triangle.c = edge.a;

			//make sure it is facing out, ie the distance to center is negative
			float dist = SignedDistanceToTriangle(triangle, center, points);
			if (dist > 0.0f) {
				std::swap(triangle.b, triangle.c);
			}

			triangles.push_back(triangle);
		}
	}

	//get projection of the origin on the closest triangle
	const int idx = ClosestTriangle(triangles, points);
	const tri_t& tri = triangles[idx];
	Vec3 ptA_w = points[tri.a].xyz;
	Vec3 ptB_w = points[tri.b].xyz;
	Vec3 ptC_w = points[tri.c].xyz;
	//bias is necessary so point A and point B are not too close as that could cause issue calculating contact normal
	Vec3 lambdas = BarycentricCoordinates(ptA_w, ptB_w, ptC_w, Vec3(0.0f));//last term is a bias subject to change

	//get the point on body A
	Vec3 ptA_a = points[tri.a].ptA;
	Vec3 ptB_a = points[tri.b].ptA;
	Vec3 ptC_a = points[tri.c].ptA;
	ptOnA = ptA_a * lambdas[0] + ptB_a * lambdas[1] + ptC_a * lambdas[2];

	//get the point on body B
	Vec3 ptA_b = points[tri.a].ptB;
	Vec3 ptB_b = points[tri.b].ptB;
	Vec3 ptC_b = points[tri.c].ptB;
	ptOnB = ptA_b * lambdas[0] + ptB_b * lambdas[1] + ptC_b * lambdas[2];

	//return penetration distance
	Vec3 delta = ptOnB - ptOnA;
	return delta.GetMagnitude();
}