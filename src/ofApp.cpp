#include "ofApp.h"
#include <cmath>
#include <stdlib.h>
#include <cstdlib>  // for rand and RAND_MAX
#include <ctime>    // for time
ofMesh mesh;
int tickyTocky = 0;
float thetaChange = 0;
float phiChange = 0;
std::pair<float, float> getRotateCoords(float x1, float y1, float rotation) {
	//rotation formula for coords is 
// newX = oldX*cos(clockwise rotation angle) - oldY*sin(clockwise rotation angle)
// newY = oldY*cos(clockwise rotation angle) + oldX*sin(clockwise rotation angle)
	float newX = x1 * cos(rotation) - y1 * sin(rotation);
	float newY = y1 * cos(rotation) + x1 * sin(rotation);
	return std::pair<float, float>(newX, newY);
}
std::pair<float, float> getDegreesRotation(float h, float k, float x1, float y1, float rotation) {
	auto newCoords = getRotateCoords(x1 - h, y1 - k, rotation * (acos(-1) / 180));
	return std::pair<float, float>(newCoords.first + h, newCoords.second + k);
}
ofPoint sphericalToCartesian(std::vector<float> input) {
	ofPoint output;
	float rho = input[0];
	float theta = input[1];
	float phi = input[2];

	//x = rho * sin(phi) * cos(theta)
	//y = rho * sin(phi) * sin(theta)
	//z = rho * cos(phi)

	output.x = rho * sin(phi) * cos(theta);
	output.y = rho * sin(phi) * sin(theta);
	output.z = rho * cos(phi);

	return output;
}
std::vector<float> cartesianToSpherical(ofPoint input) {
	std::vector<float> output;
	float x = input.x;
	float y = input.y;
	float z = input.z;
	//rho = sqrt((x*x) + (y*y) + (z*z))
	//theta = tan-1(y/x)
	//phi = cos-1(sqrt((x*x) + (y*y) + (z*z)))

	output.push_back(sqrt((x * x) + (y * y) + (z * z)));
	output.push_back(atan(y / x));
	output.push_back(acos(z/(sqrt((x * x) + (y * y) + (z * z)))));

	if (output[1] < 0) {
		output[1] += 2 * (atan(1) * 4);
	}
	if (output[2] < 0) {
		output[2] = -output[2];
		output[1] += (atan(1) * 4);
	}
	if (output[1] > 2 * (atan(1) * 4)) {
		output[1] -= 2 * (atan(1) * 4);
	}

	return output;
}
double distance(const ofPoint p1, const ofPoint p2) {
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double dz = p2.z - p1.z;
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}
ofPoint findVector(ofPoint p1, ofPoint p2) {
	ofPoint output;
	output.x = p2.x - p1.x;
	output.y = p2.y - p1.y;
	output.z = p2.z - p1.z;
	return output;
}
ofVec3f crossProduct(ofPoint a, ofPoint b) {
	ofVec3f output(a.y * b.z - a.z * b.y,a.z * b.x - a.x * b.z,a.x * b.y - a.y * b.x);
	return output;
}
ofVec3f normalize(ofVec3f v) {
	float magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return ofVec3f(v.x / magnitude, v.y / magnitude, v.z / magnitude);
}
float angleBetween(ofVec3f v1, ofVec3f v2) {
	float dotProduct = (v1.x*v2.x)+(v1.y*v2.y)+(v1.z*v2.z);
	float magnitudeProduct = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z) * sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);
	float angle = acos(dotProduct / magnitudeProduct);
	return angle;
}
void calculateNormals() {
	mesh.clearNormals();
	//loop through every three, and add to array?
	int num = mesh.getNumVertices() - 1;
	//std::vector<int> array[1000];
	std::map<int, std::vector<int>> map;
	for (int i = 0; i < mesh.getNumIndices(); i += 3) {
		map[mesh.getIndex(i)].push_back(mesh.getIndex(i + 1));
		map[mesh.getIndex(i)].push_back(mesh.getIndex(i + 2));

		map[mesh.getIndex(i + 1)].push_back(mesh.getIndex(i));
		map[mesh.getIndex(i + 1)].push_back(mesh.getIndex(i + 2));

		map[mesh.getIndex(i + 2)].push_back(mesh.getIndex(i));
		map[mesh.getIndex(i + 2)].push_back(mesh.getIndex(i + 1));
	}
	//std::vector<ofVec3f> normalsArray[1000];
	std::map<int, std::vector<ofVec3f>> normalsMap;
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		for (int j = 0; j < map[i].size(); j+=2) {
			ofVec3f vector1 = findVector(mesh.getVertex(i), mesh.getVertex(map[i][j]));
			ofVec3f vector2 = findVector(mesh.getVertex(i), mesh.getVertex(map[i][j + 1]));
			ofVec3f normal = normalize(crossProduct(vector1, vector2));
			normalsMap[i].push_back(normal);
		}
		ofVec3f total = { 0,0,0 };
		for (int j = 0; j < normalsMap[i].size(); j++) {
			total += normalsMap[i][j];
		}
		total /= normalsMap[i].size();
		mesh.addNormal(total * 10);
	}
}
void drawTorusWithMeshSetup(float x, float y, float z, float r, float R, float sides) {
	//add vertices
	float num = sides;
	sides = 360 / sides;
	std::vector<ofPoint> plusR;
	std::vector<ofPoint> plusRr;

	//loop from 0 to num using rotation code from before to add points into plusR vector
	for (int i = 0; i < num; i++) {
		auto temp = getDegreesRotation(x, y, x, y + R, i * sides);
		//ofDrawBitmapString(i, ofPoint(temp.first, temp.second, z));
		plusR.push_back(ofPoint(temp.first, temp.second, z));
	}
	//loop from 0 to num using rotation code from before to add points into the plusR vector
	for (int i = 0; i < num; i++) {
		auto temp = getDegreesRotation(x, y, x, y + R + r, i * sides);
		//ofDrawBitmapString(i, ofPoint(temp.first, temp.second, z));
		plusRr.push_back(ofPoint(temp.first, temp.second, z));
	}
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			auto sphericalCoords = cartesianToSpherical(ofPoint(plusRr[i].x - plusR[i].x, plusRr[i].y - plusR[i].y, plusRr[i].z - plusR[i].z));
			sphericalCoords[2] += sides * (acos(-1) / 180) * j;
			auto cartesianCoords = sphericalToCartesian(sphericalCoords);
			if (plusRr[i].x < 0) {
				sphericalCoords[2] += (180) * (acos(-1) / 180);
				sphericalCoords[2] = (atan(1) * 4) - sphericalCoords[2];
				cartesianCoords = sphericalToCartesian(sphericalCoords);
			}
			mesh.addVertex(ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
			//ofDrawBitmapString(j, ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
			mesh.addTexCoord(ofVec2f(i / num, j / num));
		}
	}
	//draw vertices
	//mesh.drawVertices();
	//add indecies
	///*
	for (int j = 0; j < num-1; j++) {
		float offset = j * num;
		for (int i = 0; i < num - 1; i++) {
			mesh.addIndex(i + offset);
			mesh.addIndex(i + num + offset);
			mesh.addIndex(i + 1 + offset);
		}
		mesh.addIndex(0 + offset);
		mesh.addIndex(num - 1 + offset);
		mesh.addIndex((num - 1) + num + offset);
	}
	//*/
	//
	///*
	for (int i = 0; i < num - 1; i++) {
		mesh.addIndex(i+(num * (num - 1))+1);
		mesh.addIndex(i + 1);
		mesh.addIndex(i);
	}
	mesh.addIndex(0);
	mesh.addIndex(num-1);
	mesh.addIndex((num-1)*num);
	//*/
	//
	///*
	for (int j = 0; j < num-1; j++) {
		float offset = j * num;
		for (int i = 0; i < num - 1; i++) {
			mesh.addIndex(i + 1 + num + offset);
			mesh.addIndex(i + num + offset);
			mesh.addIndex(i + 1 + offset);
		}
		mesh.addIndex(0 + offset);
		mesh.addIndex(num + offset);
		mesh.addIndex((num - 1) + num + offset);
	}
	//*/
	///*
	for (int i = 0; i < num - 1; i++) {
		mesh.addIndex(i + 1 + (num * (num - 1)));
		mesh.addIndex(i + (num * (num - 1)));
		mesh.addIndex(i);
	}
	mesh.addIndex((num* (num - 1)) + num - 1);
	mesh.addIndex((num* (num - 1)));
	mesh.addIndex((num-1));



	//normals
	//need to add them in the same order as the vertexs
	//if i is the vertex i want to calc the normal of, use the point before and after it
	//start the loop at 1 so that it works up to num
	//then add 0 which is num - 1, 0 and 1

	//normal takes the x y z input where the values are the magnitude in each direction, and its position in the array dictates the matching vector
	//for (int i = 1; i < num-1; i++) {
		//auto vec1 = mesh.getVertex(i);
		//auto vec2 = mesh.getVertex(i + 1);
		//ofVec3f crossProduct;
		//fancier maths than anticipated
		//normalize the vector afterwards

		//need to use other circles unfortunaltly 
	//}
	for (int j = 0; j < num - 1; j++) {
		float offset = j * num;
		for (int i = 0; i < num - 1; i++) {
			/*
			auto p1 = mesh.getVertex(i + offset);
			auto p2 = mesh.getVertex(i + num + offset);
			auto p3 = mesh.getVertex(i + 1 + offset);
			auto vec1 = findVector(p1, p2);
			auto vec2 = findVector(p1, p3);
			auto crossProd = crossProduct(vec2, vec1);
			auto normalized = normalize(crossProd);
			normalized *= 10;
			mesh.addNormal(normalized);
			*/
			mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(i + offset), mesh.getVertex(i + 1 + offset)), findVector(mesh.getVertex(i + offset), mesh.getVertex(i + num + offset)))));
		}
		mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(0 + offset), mesh.getVertex((num - 1) + num + offset)),findVector(mesh.getVertex(0 + offset), mesh.getVertex(num - 1 + offset)))));
	}
	//*/
	//
	///*
	for (int i = 0; i < num - 1; i++) {
		mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(i), mesh.getVertex(i + (num * (num - 1)) + 1)), findVector(mesh.getVertex(i), mesh.getVertex(i + 1)))));
	}
	mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(0), mesh.getVertex(num - 1)), findVector(mesh.getVertex(0), mesh.getVertex((num - 1) * num)))));

	//*/
	//draw mesh
	//mesh.drawWireframe();
	//mesh.draw();

	//add texture coords
	//for (float j = 0; j < 1; j += 1 / num) {
	//	for (float i = 0; i < 1; i += 1 / num) {
	//		mesh.addTexCoord(ofVec2f(0 + i, 0 + j));
	//	}
	//}
	//im assiming this needs to happen for each vertex, is a2d coordinate of where in the texture to take image, 0-1 is the width n height of it

}
void drawSphere(float radius, float numSides, float numStacks, ofVec3f center) {
	//=Create vertices=//
	// Calculate the angle between each stack and each side of the sphere
	float phiStep = acos(-1) / numStacks;
	float thetaStep = 2 * acos(-1) / numSides;

	// Loop over each stack and side of the sphere and add the vertex
	for (int i = 0; i <= numStacks; i++) {
		for (int j = 0; j <= numSides; j++) {
			std::vector<float> vertex = { radius, j * thetaStep, i * phiStep };
			mesh.addVertex(sphericalToCartesian(vertex) + center);
		}
	}

	//=Create indices=//
	// Loop over each stack and side of the sphere (again)
	// to create the indices that define the triangles that make up the sphere's surface
	for (int i = 0; i < numStacks; i++) {
		int stackStart = i * (numSides + 1); // The starting vertex index for the current stack
		int nextStackStart = (i + 1) * (numSides + 1); // The starting vertex index for the next stack

		for (int j = 0; j < numSides; j++) {
			mesh.addIndex(stackStart + j);
			mesh.addIndex(nextStackStart + j);
			mesh.addIndex(nextStackStart + j + 1);
			//mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(stackStart + j), mesh.getVertex(nextStackStart + j)), findVector(mesh.getVertex(stackStart + j), mesh.getVertex(nextStackStart + j + 1)))));

			mesh.addIndex(stackStart + j);
			mesh.addIndex(nextStackStart + j + 1);
			mesh.addIndex(stackStart + j + 1);
			//mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(stackStart + j), mesh.getVertex(nextStackStart + j + 1)), findVector(mesh.getVertex(stackStart + j), mesh.getVertex(stackStart + j + 1)))));
		}
	}
	
	//=Create Normals=//
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		mesh.addNormal(10*normalize(findVector(center, mesh.getVertex(i))));
	}



	//=Create texture coords=//
	for (float j = 0; j < 1; j += 1 / numStacks) {
		for (float i = 0; i < 1; i += 1 / numSides) {
			mesh.addTexCoord(ofVec2f(0 + i, 0 + j));
		}
	}
}
void drawCube(ofVec3f center, float radius) {
	//add vertices
	ofVec3f bottomLeft = ofPoint(center - (radius / 2));
	mesh.addVertex(bottomLeft);
	mesh.addVertex(ofPoint(bottomLeft.x + radius, bottomLeft.y, bottomLeft.z));
	mesh.addVertex(ofPoint(bottomLeft.x + radius, bottomLeft.y + radius, bottomLeft.z));
	mesh.addVertex(ofPoint(bottomLeft.x, bottomLeft.y + radius, bottomLeft.z));
	mesh.addVertex(ofPoint(bottomLeft.x, bottomLeft.y, bottomLeft.z + radius));
	mesh.addVertex(ofPoint(bottomLeft.x + radius, bottomLeft.y, bottomLeft.z + radius));
	mesh.addVertex(ofPoint(bottomLeft.x + radius, bottomLeft.y + radius, bottomLeft.z + radius));
	mesh.addVertex(ofPoint(bottomLeft.x, bottomLeft.y + radius, bottomLeft.z + radius));
	

	//add indecies
	//0
	mesh.addIndex(0);
	mesh.addIndex(1);
	mesh.addIndex(2);
	//1
	mesh.addIndex(0);
	mesh.addIndex(2);
	mesh.addIndex(3);

	//2
	mesh.addIndex(4);
	mesh.addIndex(5);
	mesh.addIndex(6);
	//3
	mesh.addIndex(4);
	mesh.addIndex(6);
	mesh.addIndex(7);

	//4
	mesh.addIndex(6);
	mesh.addIndex(2);
	mesh.addIndex(3);
	//5
	mesh.addIndex(6);
	mesh.addIndex(7);
	mesh.addIndex(3);

	//6
	mesh.addIndex(1);
	mesh.addIndex(5);
	mesh.addIndex(2);
	//7
	mesh.addIndex(2);
	mesh.addIndex(6);
	mesh.addIndex(5);

	//8
	mesh.addIndex(0);
	mesh.addIndex(4);
	mesh.addIndex(1);
	//9
	mesh.addIndex(1);
	mesh.addIndex(5);
	mesh.addIndex(4);

	//10
	mesh.addIndex(0);
	mesh.addIndex(7);
	mesh.addIndex(3);
	//11
	mesh.addIndex(0);
	mesh.addIndex(4);
	mesh.addIndex(7);

	//add Normals
	
	//0: 1(1*3),11(11*3),8(8*3)
	//ofVec3f temp1 = normalize(crossProduct(findVector(mesh.getVertex(mesh.getIndex((1 * 3))), mesh.getVertex(mesh.getIndex((1 * 3) + 1))), findVector(mesh.getVertex(mesh.getIndex((1 * 3))), mesh.getVertex(mesh.getIndex((1 * 3) + 2)))));
	//ofVec3f temp2 = normalize(crossProduct(findVector(mesh.getVertex(mesh.getIndex((11 * 3))), mesh.getVertex(mesh.getIndex((11 * 3)+1))), findVector(mesh.getVertex(mesh.getIndex((11 * 3))), mesh.getVertex(mesh.getIndex((11 * 3)+2)))));
	//ofVec3f temp3 = normalize(crossProduct(findVector(mesh.getVertex(mesh.getIndex((8 * 3))), mesh.getVertex(mesh.getIndex((8 * 3) + 1))), findVector(mesh.getVertex(mesh.getIndex((8 * 3))), mesh.getVertex(mesh.getIndex((8 * 3) + 2)))));
	//temp1 += temp3 - temp2;
	////temp1 /= 3;
	//mesh.addNormal(-10 * temp1);

	for (int i = 0; i < mesh.getNumVertices(); i++) {
		mesh.addNormal(10 * normalize(findVector(center, mesh.getVertex(i))));
	}


}
void thetaSinShiftSphere(float frequencyMultiplier, float amplitudeMultiplier) {
	
	//for (int i = 0; i < mesh.getNumVertices(); i++) {
	//	mesh.setVertex(i, mesh.getVertex(i) + (mesh.getNormal(i)*ofRandom(10)));
	//}

	//find all the vertexs and keep in vector of original vertexs

	//if a vertex is the same as one in the old array, make it equal to the new one of the one its equal to
	/*std::vector<ofVec3f> OGVertices;
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		OGVertices.push_back(mesh.getVertex(i));
	}*/

	//for (int i = 0; i < mesh.getNumVertices(); i++) {
	//	//if mesh.getVertex(i) is in the OGVertices array excluding the value of ogVertices[i] in the array
	//	//set the value to equal mesh.getVertex(i) of the new thingy

	//	if (vectorContainsElementExcludingIndex(OGVertices, mesh.getVertex(i), i) == -1) {
	//		mesh.setVertex(i, mesh.getVertex(i) + (mesh.getNormal(i) * 10));
	//	}
	//	else {
	//		mesh.setVertex(i, mesh.getVertex(vectorContainsElementExcludingIndex(OGVertices, mesh.getVertex(i), i)));
	//	}
	//	
	//}
	
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float theta = atan2(vertex.y, vertex.x);
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier *sin(theta * frequencyMultiplier)));
	}

}
void phiSinShiftSphere(int frequencyMultiplier, int amplitudeMultiplier) {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float phi = acos(vertex.z / (sqrt((vertex.x * vertex.x) + (vertex.y * vertex.y) + (vertex.z * vertex.z)))); 
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier * sin(phi * frequencyMultiplier)));
	}
}
void thetaCosShiftSphere(int frequencyMultiplier, int amplitudeMultiplier) {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float theta = atan2(vertex.y, vertex.x);
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier*cos(theta * frequencyMultiplier)));
	}
}
void phiCosShiftSphere(int frequencyMultiplier, int amplitudeMultiplier) {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float phi = acos(vertex.z / (sqrt((vertex.x * vertex.x) + (vertex.y * vertex.y) + (vertex.z * vertex.z))));
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier *cos(phi * frequencyMultiplier)));
	}
}
void thetaTanShiftSphere(int frequencyMultiplier, int amplitudeMultiplier) {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float theta = atan2(vertex.y, vertex.x);
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier * tan(theta * frequencyMultiplier)));
	}
}
void phiTanShiftSphere(int frequencyMultiplier, int amplitudeMultiplier) {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		float phi = acos(vertex.z / (sqrt((vertex.x * vertex.x) + (vertex.y * vertex.y) + (vertex.z * vertex.z))));
		mesh.setVertex(i, vertex + (mesh.getNormal(i) * amplitudeMultiplier *tan(phi * frequencyMultiplier)));
	}
}
void createDirectionalLight(ofVec3f lightVector, ofColor color) {
	glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < mesh.getNumNormals(); i++) {
		float angle = angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4));
		if (angle < 100) {
			mesh.addColor(ofColor(color.r - (color.r * (angle / 100)), color.g - (color.g * (angle / 100)), color.b - (color.b * (angle / 100))));
		}
		else {
			mesh.addColor(ofColor(0,0,0));
		}	
	}
}
void createAmbientLight(ofColor color) {
	glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < mesh.getNumNormals(); i++) {
		mesh.addColor(color);
	}
}
void createPointLight(ofColor color, ofPoint origin) {
	glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < mesh.getNumNormals(); i++) {
		ofPoint vertex = mesh.getVertex(i);
		ofVec3f normal = mesh.getNormal(i);
		ofVec3f vectorBetween = findVector(vertex,origin);
		float angle = angleBetween(vectorBetween, normal) * (180 / (atan(1) * 4));
		if (angle < 100) {
			mesh.addColor(ofColor(color.r - (color.r * (angle / 100)), color.g - (color.g * (angle / 100)), color.b - (color.b * (angle / 100))));
		}
		else {
			mesh.addColor(ofColor(0, 0, 0));
		}
	}

	ofDrawSphere(origin, 10);

}
void drawNormals() {
	ofSetColor(ofColor(255, 0, 0));
for (int i = 0; i < mesh.getNumNormals(); i++) {
	ofDrawLine(mesh.getVertex(i), ofPoint(mesh.getVertex(i).x + mesh.getNormal(i).x, mesh.getVertex(i).y + mesh.getNormal(i).y, mesh.getVertex(i).z + mesh.getNormal(i).z));
}
ofSetColor(ofColor(255, 255, 255));
}
void drawVertexLabels() {
	for (int i = 0; i < mesh.getNumVertices(); i++) {
		ofDrawBitmapString(i, mesh.getVertex(i));
	}
}
void recalculateTorusNormals(float num) {
	mesh.clearNormals();
	for (int j = 0; j < num - 1; j++) {
		float offset = j * num;
		for (int i = 0; i < num - 1; i++) {
			mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(i + offset), mesh.getVertex(i + 1 + offset)), findVector(mesh.getVertex(i + offset), mesh.getVertex(i + num + offset)))));
		}
		mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(0 + offset), mesh.getVertex((num - 1) + num + offset)), findVector(mesh.getVertex(0 + offset), mesh.getVertex(num - 1 + offset)))));
	}
	for (int i = 0; i < num - 1; i++) {
		mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(i), mesh.getVertex(i + (num * (num - 1)) + 1)), findVector(mesh.getVertex(i), mesh.getVertex(i + 1)))));
	}
	mesh.addNormal(10 * normalize(crossProduct(findVector(mesh.getVertex(0), mesh.getVertex(num - 1)), findVector(mesh.getVertex(0), mesh.getVertex((num - 1) * num)))));
}
void recalculateSphereNormals() {
	//mesh.clearNormals();
	//for (int i = 0; i < mesh.getNumIndices(); i+=6) {
	//	//i
	//	//i+1
	//	//i+2
	//	ofVec3f temp1 = normalize(crossProduct(findVector(mesh.getVertex(mesh.getIndex(i)), mesh.getVertex(mesh.getIndex(i+1))), findVector(mesh.getVertex(mesh.getIndex(i)), mesh.getVertex(mesh.getIndex(i+2)))));


	//	//i+3
	//	//i+4
	//	//i+5
	//	ofVec3f temp2 = normalize(crossProduct(findVector(mesh.getVertex(mesh.getIndex(i+3)), mesh.getVertex(mesh.getIndex(i+4))), findVector(mesh.getVertex(mesh.getIndex(i+3)), mesh.getVertex(mesh.getIndex(i+5)))));

	//	temp1 += temp2;
	//	//temp1 /= 2;
	//	mesh.addNormal(10* temp1);
	
	//}
	calculateNormals();
}
//===============================================================================//
void ofApp::setup(){
	gui.setup();
	directionalLightGui.setup();
	ambientLightGui.setup();
	pointLightGui.setup();
	shiftGui.setup();
	sinShiftGui.setup();
	cosShiftGui.setup();
	tanShiftGui.setup();
	gui.setSize(ofGetWidth() / 2, ofGetHeight() / 2);

	
	gui.add(sphereToggle.setup("Sphere", true));
	gui.add(torusToggle.setup("Torus", false));
	gui.add(cubeToggle.setup("Cube", false));

	gui.add(radius.setup("Radius", 100, 10, 1000));
	gui.add(numSides.setup("Num sides", 50, 3, 100));
	gui.add(numStacks.setup("Num stacks", 50, 3, 100));
	gui.add(originSlider.setup("Origin", ofVec3f(0, 0, 0), ofVec3f(-1000, -1000, -1000), ofVec3f(1000, 1000, 1000)));

	gui.add(sidesAnimationToggle.setup("Number of Sides Animation", false));

	gui.add(directionalLightToggle.setup("Add Directional light", true));

	gui.add(ambientLightToggle.setup("Add Ambient light", false));

	gui.add(pointLightToggle.setup("Add Point Light", false));

	gui.add(normalsToggle.setup("Normals", false));

	gui.add(vertexLabelsToggle.setup("Vertex Labels", false));

	gui.add(textureToggle.setup("Texture", false));

	gui.add(drawWireframeToggle.setup("WireFrame", false));
	gui.add(drawToggle.setup("Draw", true));

	gui.add(shiftToggle.setup("Shifts", false));

	directionalLightGui.add(colourForDirectionalLight.setup("Colour", ofVec3f(180, 0, 100), ofVec3f(0, 0, 0), ofVec3f(255, 255, 255)));
	directionalLightGui.add(vectorForDirectionalLight.setup("Vector", ofVec3f(0, 0, 10), ofVec3f(0, 0, 0), ofVec3f(200, 200, 200)));

	ambientLightGui.add(colourForAmbientLight.setup("Colour", ofVec3f(0, 0, 0), ofVec3f(0, 0, 0), ofVec3f(255, 255, 255)));

	pointLightGui.add(colourForPointLight.setup("Colour", ofVec3f(0, 0, 0), ofVec3f(0, 0, 0), ofVec3f(255, 255, 255)));
	pointLightGui.add(originForPointLight.setup("Origin", ofVec3f(0, 0, 0), ofVec3f(-1000, -1000, -1000), ofVec3f(1000, 1000, 1000)));
	pointLightGui.add(pointLightOrbitToggle.setup("Orbit", false));

	shiftGui.add(sinToggle.setup("Sin Shifts", false));
	shiftGui.add(cosToggle.setup("Cos Shifts", false));
	shiftGui.add(tanToggle.setup("Tan Shifts", false));

	sinShiftGui.add(thetaSinShiftLabel.setup("Theta Sin Shift", ""));
	sinShiftGui.add(thetaSinShiftToggle.setup("Theta Sin Shift", false));
	sinShiftGui.add(thetaSinFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(thetaSinAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	sinShiftGui.add(thetaSinAnimation.setup("Theta Sin Shift Animation", false));

	sinShiftGui.add(phiSinShiftLabel.setup("Phi Sin Shift", ""));
	sinShiftGui.add(phiSinShiftToggle.setup("Phi Sin Shift", false));
	sinShiftGui.add(phiSinFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(phiSinAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	sinShiftGui.add(phiSinAnimation.setup("Phi Sin Shift Animation", false));

	cosShiftGui.add(thetaCosShiftLabel.setup("Theta Cos Shift", ""));
	cosShiftGui.add(thetaCosShiftToggle.setup("Theta Cos Shift", false));
	cosShiftGui.add(thetaCosFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(thetaCosAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	cosShiftGui.add(thetaCosAnimation.setup("Theta Cos Shift Animation", false));

	cosShiftGui.add(phiCosShiftLabel.setup("Phi Cos Shift", ""));
	cosShiftGui.add(phiCosShiftToggle.setup("Phi Cos Shift", false));
	cosShiftGui.add(phiCosFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(phiCosAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	cosShiftGui.add(phiCosAnimation.setup("Phi Cos Shift Animation", false));

	tanShiftGui.add(thetaTanShiftlabel.setup("Theta Tan Shift", ""));
	tanShiftGui.add(thetaTanShiftToggle.setup("Theta Tan Shift", false));
	tanShiftGui.add(thetaTanFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(thetaTanAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	tanShiftGui.add(thetaTanAnimation.setup("Theta Tan Shift Animation", false));

	tanShiftGui.add(phiTanShiftlabel.setup("Phi Tan Shift", ""));
	tanShiftGui.add(phiTanShiftToggle.setup("Phi Tan Shift", false));
	tanShiftGui.add(phiTanFrequencySlider.setup("Frequency Multiplier", 10, 0, 100));
	sinShiftGui.add(phiTanAmplitudeSlider.setup("Amplitude Multiplier", 1, 0, 100));
	tanShiftGui.add(phiTanAnimation.setup("Phi Tan Shift Animation", false));


	glPointSize(4);
	ofDisableAlphaBlending();
	ofEnableDepthTest();
	ofDisableArbTex();
	ofEnableNormalizedTexCoords();
	ofLoadImage(mTex, "sprinkles.png");
}	
//--------------------------------------------------------------
void ofApp::update(){
}

//--------------------------------------------------------------
void ofApp::draw(){
	//ofLog() << ofGetFrameNum();
	//ofLog() << ofGetFrameRate();

	ofDisableDepthTest();
	gui.draw();

	if (directionalLightToggle) {
		directionalLightGui.draw();
	}
	if (ambientLightToggle) {
		ambientLightGui.draw();
	}
	if (pointLightToggle) {
		pointLightGui.draw();
	}
	if (shiftToggle) {
		shiftGui.draw();
		if (sinToggle) {
			sinShiftGui.draw();
		}
		if (cosToggle) {
			cosShiftGui.draw();
		}
		if (tanToggle) {
			tanShiftGui.draw();
		}
	}
	ofEnableDepthTest();
	cam.begin();
	if (textureToggle) {
		mTex.bind();
	}
	ofNoFill();
	mesh.clear();
	if (sphereToggle) {
		if (sidesAnimationToggle) {
			int tens = (tickyTocky / 10) % 10; // get the tens digit
			int hundreds = (tickyTocky / 100) % 10; // get the hundreds digit
			numSides = tens + (hundreds * 10);
			numStacks = tens + (hundreds * 10);
			if (numSides < 3) {
				numSides = 3;
			}
			if (numStacks < 3) { 
				numStacks = 3;
			}
		}
		drawSphere(radius, numSides, numStacks, ofVec3f(originSlider->x, originSlider->y, originSlider->z));
	}
	if (torusToggle) {
		if (sidesAnimationToggle) {
			int tens = (tickyTocky / 10) % 10; // get the tens digit
			int hundreds = (tickyTocky / 100) % 10; // get the hundreds digit
			numSides = tens + (hundreds * 10);
			numStacks = tens + (hundreds * 10);
			if (numSides < 3) {
				numSides = 3;
			}
			if (numStacks < 3) {
				numStacks = 3;
			}
		}
		drawTorusWithMeshSetup(originSlider->x, originSlider->y, originSlider->z, radius, radius * 2, numSides);
	}
	if (cubeToggle) {
		drawCube(ofVec3f(originSlider->x, originSlider->y, originSlider->z), radius);
	}
	if (directionalLightToggle) {
		createDirectionalLight(ofVec3f(vectorForDirectionalLight->x / 10, vectorForDirectionalLight->y / 10, vectorForDirectionalLight->z / 10), ofColor(colourForDirectionalLight->x, colourForDirectionalLight->y, colourForDirectionalLight->z));
	}
	if (ambientLightToggle) {
		createAmbientLight(ofColor(colourForAmbientLight->x, colourForAmbientLight->y, colourForAmbientLight->z));
	}
	if (pointLightToggle) {
		if (pointLightOrbitToggle) {
			thetaChange += 0.01;
			phiChange += 0.05;
			auto spherical = cartesianToSpherical(ofPoint(1000, 1000, 1000));
			spherical[1] += thetaChange;
			spherical[2] += phiChange;
			auto cartesian = sphericalToCartesian(spherical);
			createPointLight(ofColor(colourForPointLight->x, colourForPointLight->y, colourForPointLight->z), cartesian);
		}
		else {
			createPointLight(ofColor(colourForPointLight->x, colourForPointLight->y, colourForPointLight->z), ofPoint(originForPointLight->x, originForPointLight->y, originForPointLight->z));
		}
	}
	if (shiftToggle) {
		if (thetaTanShiftToggle) {
			if (thetaTanAnimation) {
				thetaTanShiftSphere(thetaTanFrequencySlider + ofGetFrameNum(), thetaTanAmplitudeSlider);
			}
			else {
				thetaTanShiftSphere(thetaTanFrequencySlider, thetaTanAmplitudeSlider);
			}
		}
		if (phiTanShiftToggle) {
			if (phiTanAnimation) {
				thetaTanShiftSphere(phiTanFrequencySlider + ofGetFrameNum(), phiTanAmplitudeSlider);
			}
			else {
				thetaTanShiftSphere(phiTanFrequencySlider, phiTanAmplitudeSlider);
			}
		}
		if (thetaCosShiftToggle) {
			if (thetaCosAnimation) {
				thetaCosShiftSphere(thetaCosFrequencySlider + ofGetFrameNum(), thetaCosAmplitudeSlider);
			}
			else {
				thetaCosShiftSphere(thetaCosFrequencySlider, thetaCosAmplitudeSlider);
			}
		}
		if (phiCosShiftToggle) {
			if (phiCosAnimation) {
				phiCosShiftSphere(phiCosFrequencySlider + ofGetFrameNum(), phiCosAmplitudeSlider);
			}
			else {
				phiCosShiftSphere(phiCosFrequencySlider, phiCosAmplitudeSlider);
			}
		}
		if (thetaSinShiftToggle) {
			if (thetaSinAnimation) {
				thetaSinShiftSphere(thetaSinFrequencySlider + ofGetFrameNum(), thetaSinAmplitudeSlider);
			}
			else {
				thetaSinShiftSphere(thetaSinFrequencySlider, thetaSinAmplitudeSlider);
			}
		}
		if (phiSinShiftToggle) {
			if (phiSinAnimation) {
				phiSinShiftSphere(phiSinFrequencySlider + ofGetFrameNum(), phiSinFrequencySlider);
			}
			else {
				phiSinShiftSphere(phiSinFrequencySlider, phiSinAmplitudeSlider);
			}
		}
	}
	if (sphereToggle) {
		calculateNormals();
	}
	if (torusToggle) {
		recalculateTorusNormals(numSides);
	}
	if (drawWireframeToggle) {
		mesh.drawWireframe();
	}
	if (drawToggle) {
		mesh.draw();
	}
	if (normalsToggle) {
		drawNormals();
	}
	if (vertexLabelsToggle) {
		drawVertexLabels();
	}
	if (textureToggle) {
		mTex.unbind();
	}
	cam.end();
	tickyTocky++;
	tickyTocky++;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){	
}
//--------------------------------------------------------------
void ofApp::keyReleased(int key){
}
//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){
}
//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
}
//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
}
//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
}
//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){
}
//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){
}
//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){
}
//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
}
//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 
}