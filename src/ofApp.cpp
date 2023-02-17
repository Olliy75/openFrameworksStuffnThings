#include "ofApp.h"
#include <cmath>
#include <stdlib.h>
#include <cstdlib>  // for rand and RAND_MAX
#include <ctime>    // for time
ofMesh mesh;
void dot(float x,float y) {
	ofDrawRectangle(x, y, 1, 1);
}
void horizontalLine() {
	for (int i = 64; i < 128; i++) {
		ofDrawRectangle(i, 64, 1, 1);
	}
}
void verticalLine() {
	for (int i = 64; i < 128; i++) {
		ofDrawRectangle(64, i, 1, 1);
	}
}
void diagonalLine() {
	for (int i = 64; i < 128; i++) {
		ofDrawRectangle(i, i, 1, 1);
	}
}
void diagonalLineArbitrary(int slope) {
	int yValue = 64;
	// so say the slope is 3, then i just add another loops for the number in the slope and adds another pixel
	// this like makes a like that is horizontal, i need to mess with the y value which increases to the slope
	// need to add a dot on each x coord based on the number of slope, and need to keep increasing the y value
	for (int i = 64; i < 128; i++) {
		for (int j = 0; j <= slope; j++) {
			ofDrawRectangle(i, yValue, 1, 1);
			yValue++;
		}
		
	}
}
void straightLineBetweenPoints(int x1, int y1, int x2, int y2) {
	if (y1 == y2) {
		//means that it is horizontal
		for (int i = x1; i < x2; i++) {
			ofDrawRectangle(i, y1, 1, 1);
		}
	}
	else if (x1 == x2) {
		// means that it is vertical
		for (int i = y1; i < y2; i++) {
			ofDrawRectangle(x1, i, 1, 1);
		}
	}
}
void lineBetweenTwoPoints(int x1, int y1, int x2, int y2) {
	if ((x1 > x2)&&(y1 > y2)) {
		int tempx = x1;
		int tempy = y1;
		x1 = x2;
		y1 = y2;
		x2 = tempx;
		y2 = tempy;
	}
	int slope = ((y2 - y1) / (x2 - x1));
	if (slope == 0) {
		straightLineBetweenPoints(x1, y1, x2, y2);
	}
	if (slope > 0)
	{
		int yValue = y1;
		for (int i = x1; i < x2; i++) {
			for (int j = 0; j <= slope; j++) {
				ofDrawRectangle(i, yValue, 1, 1);
				yValue++;
			}
		}
	}
	else {
		//
		slope = slope + (slope * 2);
		int yValue = y1;
		for (int i = x1; i > x2; i--) {
			for (int j = slope; j >= 0; j--) {
				ofDrawRectangle(i, yValue, 1, 1);
				yValue--;
			}
		}
	}

}
//start again using y = mx+c, use coords and find line that goes through those two points, then plot a line using y=mx+c.
std::pair<float, float> findLineEquation(float x1, float y1, float x2, float y2) {
	float m = ((y2 - y1) / (x2 - x1));
		//y1=(m*x1)+c
		//c = y1-(m*x1)
		float c = y1 - (m * x1);
		return std::pair<float, float>(m, c);
}
void drawLineUsingFormula(float x1, float y1, float x2, float y2) {
	if (x1 == x2) {
		for (int i = y1; i < y2; i++) {
			ofDrawRectangle(x1, i, 1, 1);
		}
	}
	else {
		auto mc = findLineEquation(x1, y1, x2, y2);
		float m = mc.first;
		float c = mc.second;

		if (x1 < x2) {
			for (int i = x1; i < x2; i++) {
				ofDrawRectangle(i, ((m * i) + c), 1, 1);
			}
		}
		else {
			for (int i = x2; i < x1; i++) {
				ofDrawRectangle(i, ((m * i) + c), 1, 1);
			}
		}

		if (y1 < y2) {
			for (int i = y1; i < y2; i++) {
				ofDrawRectangle(((i - c) / m), i, 1, 1);
			}
		}
		else {
			for (int i = y2; i < y1; i++) {
				ofDrawRectangle(((i - c) / m), i, 1, 1);
			}
		}
	}
	
}
void drawSquare(float x1, float y1, float width) {
	//topLeft - topRight
	drawLineUsingFormula(x1, y1, (x1 + width), y1);
	//topLeft - bottomLeft
	drawLineUsingFormula(x1, y1, x1, (y1 + width));
	//bottomLeft - bottomRight
	drawLineUsingFormula(x1, (y1 + width), (x1 + width), (y1 + width));
	//topRight - bottomRight
	drawLineUsingFormula((x1 + width), y1, (x1 + width), (y1 + width));
}
//rotation formula for coords is 
// newX = oldX*cos(clockwise rotation angle) - oldY*sin(clockwise rotation angle)
// newY = oldY*cos(clockwise rotation angle) + oldX*sin(clockwise rotation angle)
std::pair<float, float> getRotateCoords(float x1, float y1, float rotation) {
	float newX = x1*cos(rotation) - y1*sin(rotation);
	float newY = y1*cos(rotation) + x1*sin(rotation);
	return std::pair<float, float>(newX, newY);
}
void drawSquareRotationOrigin(float x1, float y1, float width, float rotation) {
	rotation = rotation * (acos(-1) / 180);
	auto topLeft = getRotateCoords(x1, y1, rotation);
	auto topRight = getRotateCoords((x1 + width), y1, rotation);
	auto bottomLeft = getRotateCoords(x1, (y1 + width), rotation);
	auto bottomRight = getRotateCoords((x1 + width), (y1 + width), rotation);
	//topLeft - topRight
	drawLineUsingFormula(topLeft.first,topLeft.second,topRight.first,topRight.second);
	//topLeft - bottomLeft
	drawLineUsingFormula(topLeft.first, topLeft.second, bottomLeft.first, bottomLeft.second);
	//bottomLeft - bottomRight
	drawLineUsingFormula(bottomLeft.first, bottomLeft.second, bottomRight.first, bottomRight.second);
	//topRight - bottomRight
	drawLineUsingFormula(topRight.first, topRight.second, bottomRight.first, bottomRight.second);
}
void DrawSquareRotationArbitrary(float x1, float y1, float width, float rotation, float pointx, float pointy) {
	//move all points so that the arbitrary point becomes 0,0
	float translatedx = x1 - pointx;
	float translatedy = y1 - pointx;
	//throw it into the rotate about origin
	rotation = rotation * (acos(-1) / 180);
	auto topLeft = getRotateCoords(x1, y1, rotation);
	auto topRight = getRotateCoords((x1 + width), y1, rotation);
	auto bottomLeft = getRotateCoords(x1, (y1 + width), rotation);
	auto bottomRight = getRotateCoords((x1 + width), (y1 + width), rotation);
	//translate back to original state
	topLeft.first = topLeft.first + pointx;
	topLeft.second = topLeft.second + pointy;

	topRight.first = topRight.first + pointx;
	topRight.second = topRight.second + pointy;

	bottomLeft.first = bottomLeft.first + pointx;
	bottomLeft.second = bottomLeft.second + pointy;

	bottomRight.first = bottomRight.first + pointx;
	bottomRight.second = bottomRight.second + pointy;
		
	//draw
	drawLineUsingFormula(topLeft.first, topLeft.second, topRight.first, topRight.second);
	//topLeft - bottomLeft
	drawLineUsingFormula(topLeft.first, topLeft.second, bottomLeft.first, bottomLeft.second);
	//bottomLeft - bottomRight
	drawLineUsingFormula(bottomLeft.first, bottomLeft.second, bottomRight.first, bottomRight.second);
	//topRight - bottomRight
	drawLineUsingFormula(topRight.first, topRight.second, bottomRight.first, bottomRight.second);
}
void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3) {
	drawLineUsingFormula(x1, y1, x2, y2);
	drawLineUsingFormula(x1, y1, x3, y3);
	drawLineUsingFormula(x3, y3, x2, y2);
}
void drawEqualatralTriangle(float x, float y, float width) {
	//draw line from (x,y) to (x+width,y)
	//draw line from (x,y) to (x+(width/2),y-tan(cos(1/2))*(width/2))
	//draw line from (x+width,y) to (x+(width/2),y-tan(cos(1/2))*(width/2))
	//drawLineUsingFormula(x, y, x + width, y);
	//drawLineUsingFormula(x, y, x + (width / 2), y - sqrt((width * width) - ((width / 2) * (width / 2))));
	//drawLineUsingFormula(x + width, y,x + (width / 2), y - sqrt((width * width) - ((width / 2) * (width / 2))));
	float x2 = x + width;
	float y2 = y;
	float x3 = x + (width / 2);
	float y3 = y - sqrt((width * width) - ((width / 2) * (width / 2)));
	drawTriangle(x, y, x2, y2, x3, y3);
}
std::pair<float, float> midpoint(float x1,float y1,float x2,float y2) {
	return std::pair<float, float>((x1+x2)/2, (y1+y2)/2);
}
void sepinski(float x1, float y1, float x2, float y2, float x3, float y3, int n) {
	//is 3 ways to draw it seems, need to pick, maybe start by doing it with squares because i already have the code and triangles come later.

	//draw a serpinski with just 3 squares, then draw a serpinski using the last serpinski three times etc etc etc

	//first one is just three points at the bottom left of the screen
	//0,765,1,765,0.5,654
	//in reality this is (x,y)(x+width,y)(x+(width/2),y-width)
	//then set the new width to double what it was before, rice repeat?
	//float x = 0;
	//float y = 765;
	//for (float i = 1; i < 1000; i=i*2) {
	//	ofDrawRectangle(x, y, i, i);
	//	ofDrawRectangle(x+i, y, i, i);
	//	ofDrawRectangle(x+(i/2), y-i, i, i);
	//}

	//start from max size and work way back down?
	
	drawTriangle(x1, y1, x2, y2, x3, y3);

	// Seed the random number generator with the current time
	srand(time(0));

	// Generate a random integer between 1 and 10

	//int r = rand() % 10 + 1;
	//if (r <= 3) {
	//	ofSetColor(255, 0, 0);
	//}
	//else if(r <= 6){
	//	ofSetColor(0, 255, 0);
	//}
	//else if (r <= 9) {
	//	ofSetColor(0, 0, 255);
	//}
	//else {
	//	ofSetColor(0, 0, 0);
	//}

	if (n == 0) {
		return;
	}

	std::pair<float, float> m1 = midpoint(x1, y1, x2, y2);
	std::pair<float, float> m2 = midpoint(x2, y2, x3, y3);
	std::pair<float, float> m3 = midpoint(x1, y1, x3, y3);

	drawTriangle(m1.first, m1.second, m2.first, m2.second, m3.first, m3.second);
	n--;
	sepinski(m1.first, m1.second, x2, y2, m2.first, m2.second, n);
	sepinski(m3.first, m3.second, m2.first, m2.second, x3, y3, n);
	sepinski(x1, y1, m1.first, m1.second, m3.first, m3.second, n);
}
void drawCircle(float h, float k, float r) {
	//(x - h) ^ 2 + (y - k) ^ 2 = r ^ 2
	//(h,k) = center of the circle
	//r = radius of the circle

	//y = k + sqrt(r ^ 2 - (x - h) ^ 2)
	//    k - sqrt(r ^ 2 - (x - h) ^ 2)

	//x = h + sqrt(r ^ 2 - (y - k) ^ 2)
	//	  h - sqrt(r ^ 2 - (y - k) ^ 2)

	//rince and repeat plotting all whole number values of y and all of em for x?

	//min value for x and y is h - r and k - r
	//max value for x and y is h + r and k + r
	for (int x = h - r; x < h + r; x++) {
		ofDrawRectangle(x, k + sqrt((r * r) - ((x - h) * (x - h))), 1, 1);
		ofDrawRectangle(x, k - sqrt((r * r) - ((x - h) * (x - h))), 1, 1);
	}

	for (int y = k - r; y < k + r; y++) {
		ofDrawRectangle(h + sqrt((r * r) - ((y - k) * (y - k))), y, 1, 1);
		ofDrawRectangle(h - sqrt((r * r) - ((y - k) * (y - k))), y, 1, 1);
	}
}
std::pair<float, float> getDegreesRotation(float h,float k, float x1, float y1, float rotation) {
	auto newCoords = getRotateCoords(x1 - h, y1 - k, rotation * (acos(-1) / 180));
	return std::pair<float, float>(newCoords.first + h,newCoords.second + k);
}
std::pair<float, float> getDegreesRotationZ(float h,float k, float y1, float z1, float rotation) {
	rotation = rotation * (acos(-1) / 180);
	float newY = (y1-h) * cos(rotation) - (z1-k) * sin(rotation);
	float newZ = (z1-k) * cos(rotation) + (y1-h) * sin(rotation);
	return std::pair<float, float>(newY + h, newZ + k);
}
void drawCylinder(float h1, float k1, float r, float h2, float k2) {
	drawCircle(h1, k1, r);
	drawCircle(h2, k2, r);

	drawLineUsingFormula(h1, k1 + r, h2, k2 + r);
	drawLineUsingFormula(h1 + r, k1, h2 + r, k2);
	drawLineUsingFormula(h1, k1 - r, h2, k2 - r);
	drawLineUsingFormula(h1 - r, k1, h2 - r, k2);
}
void drawCone(float h, float k, float r, float x, float y) {
	drawCircle(h, k, r);

	//8 values of coordinates of the circle or something like that, need to figure out how to do that tho

	//first 4 lines
	drawLineUsingFormula(h, k + r, x, y);
	drawLineUsingFormula(h + r, k, x, y);
	drawLineUsingFormula(h, k - r, x, y);
	drawLineUsingFormula(h - r, k, x, y);

	auto points = getDegreesRotation(h, k, h, k + r, 35);
	drawLineUsingFormula(points.first, points.second, x, y);
	 points = getDegreesRotation(h, k, h + r, k,35);
	drawLineUsingFormula(points.first, points.second, x, y);
	 points = getDegreesRotation(h, k, h, k - r,35);
	drawLineUsingFormula(points.first, points.second, x, y);
	 points = getDegreesRotation(h, k, h - r, k,35);
	drawLineUsingFormula(points.first, points.second, x, y);
	
}
void drawpyramid(float x1, float y1, float width, float x2, float y2) {
	drawSquare(x1, y1, width);

	std::pair<float, float> point1(x1,y1);
	std::pair<float, float> point2(x1+width,y1);
	std::pair<float, float> point3(x1+width,y1+width);
	std::pair<float, float> point4(x1,y1+width);
	std::pair<float, float> point5(x1+(width/2),y1);
	std::pair<float, float> point6(x1+(width/2),y1+width);
	std::pair<float, float> point7(x1,y1+(width/2));
	std::pair<float, float> point8(x1+width,y1+(width/2));

	drawLineUsingFormula(point1.first, point1.second, x2, y2);
	drawLineUsingFormula(point2.first, point2.second, x2, y2);
	drawLineUsingFormula(point3.first, point3.second, x2, y2);
	drawLineUsingFormula(point4.first, point4.second, x2, y2);
	drawLineUsingFormula(point5.first, point5.second, x2, y2);
	drawLineUsingFormula(point6.first, point6.second, x2, y2);
	drawLineUsingFormula(point7.first, point7.second, x2, y2);
	drawLineUsingFormula(point8.first, point8.second, x2, y2);

}
void drawCylinderWithStackedCircles(float h, float k, float r, float height, string direction) {
	if (direction == "horizontal") {
		//draw circles with incrementing y
		for (int y = k; y < k+height; y=y+10) {
			drawCircle(h, y, r);
		}
	}
	else {
		//draw cirlces with incrementing x
		for (int x = h; x < h + height; x=x+10) {
			drawCircle(x, k, r);
		}
	}
}
void drawMeshCube(float x, float y, float z, float size) {
	//so first you need to create all of the vertexs in such a way where the order is known
	//add them as vertexs in an order where three at a time create a triangle
	ofMesh mesh;
	//a: first side bottom left
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y + size, z));
	mesh.addVertex(ofPoint(x + size, y + size, z));
	//a: first side top right
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x + size, y, z));
	mesh.addVertex(ofPoint(x + size, y + size, z));

	//b:second side bottom left
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x + size, y, z));
	mesh.addVertex(ofPoint(x + size, y, z + size));
	//b:second side top right
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y, z + size));
	mesh.addVertex(ofPoint(x + size, y, z + size));

	//c:third side bottom left
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y + size, z));
	mesh.addVertex(ofPoint(x, y + size, z + size));
	//c:third side top right
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y, z + size));
	mesh.addVertex(ofPoint(x, y + size, z + size));

	//switching to other side of the cube
	x += size;
	y += size;
	z += size;

	//a: first side bottom left OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y - size, z));
	mesh.addVertex(ofPoint(x - size, y - size, z));
	//a: first side top right OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x - size, y, z));
	mesh.addVertex(ofPoint(x - size, y - size, z));

	//b:second side bottom left OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x - size, y, z));
	mesh.addVertex(ofPoint(x - size, y, z - size));
	//b:second side top right OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y, z - size));
	mesh.addVertex(ofPoint(x - size, y, z - size));

	//c:third side bottom left OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y - size, z));
	mesh.addVertex(ofPoint(x, y - size, z - size));
	//c:third side top right OPPOSITE
	mesh.addVertex(ofPoint(x, y, z));
	mesh.addVertex(ofPoint(x, y, z - size));
	mesh.addVertex(ofPoint(x, y - size, z - size));

	for (int i = 0; i < 36; i++) {
		mesh.addIndex(i);
	}
	mesh.drawWireframe();
}
void cleanDrawMeshCube(float x, float y, float z, float size) {

}
void drawLine3Dimensionally(ofPoint a, ofPoint b) {

}
std::pair<std::vector<ofPoint>, std::vector<float>> drawCircleWithMesh(float x, float y, float z, float r, float sides) {
	int counter = 0;
	std::vector<ofPoint> vertex;
	std::vector<float> index;
	sides = 360 / sides;
	//add first point
	ofMesh mesh;
	mesh.addVertex(ofPoint(x,y,z));
	vertex.push_back(ofPoint(x,y,z));
	counter++;
	//centre is now set to vertex 0
	for (int i = 0; i < sides; i++) {
		auto temp = getDegreesRotation(x, y, x, y + r, i * sides);
		mesh.addVertex(ofPoint(temp.first,temp.second,z));
		vertex.push_back(ofPoint(temp.first, temp.second, z));
		counter++;
	}

	for (int i = 0; i < counter -1; i++) {
		mesh.addIndex(0);
		index.push_back(0);
		mesh.addIndex(i);
		index.push_back(i);
		mesh.addIndex(i + 1);
		index.push_back(i + 1);
	}
	mesh.drawWireframe();
	return std::pair<std::vector<ofPoint>, std::vector<float>>(vertex,index);
}
void drawCylinderWithMesh(float x, float y, float z,float r,float sides,float width) {
	std::vector<ofPoint> vertex1;
	std::vector<float> index1;
	std::vector<ofPoint> vertex2;
	std::vector<float> index2;
	std::pair<std::vector<ofPoint>, std::vector<float>> pair;

	pair = drawCircleWithMesh(x, y, z + width / 2, r, sides);
	vertex1 = pair.first;
	index1 = pair.second;

	pair = drawCircleWithMesh(x, y, z - width / 2, r, sides);
	vertex2 = pair.first;
	index2 = pair.second;
	ofMesh mesh;
	for (int i = 1; i < vertex1.size(); i++) {
		//need a draw line function that works 3 dimensionally
		//function(vertex1(i),vertex2(i));
		mesh.addVertex(vertex1[i]);
		mesh.addVertex(vertex2[i]);
	}
	for (int i = 0; i < vertex1.size()-1; i++) {
		mesh.addVertex(vertex1[i]);
		mesh.addVertex(vertex1[i + 1]);
		mesh.addVertex(vertex2[i]);
		mesh.addVertex(vertex2[i + 1]);

		mesh.addIndex(i);
		mesh.addIndex(i + 1);
		mesh.addIndex(i + 2);

		mesh.addIndex(i + 1);
		mesh.addIndex(i + 2);
		mesh.addIndex(i + 3);
	}
	mesh.drawWireframe();
}
void drawSphereWithMesh(float x, float y, float z, float r, float sides) {
	float zPrime = z;
	float num = sides;
	sides = 360 / sides;
	float h = x;
	float k = y;
	int counter = 0;

	ofMesh mesh;

	//for (int i = 0; i < sides; i++) {
		//auto temp = getDegreesRotationZ(h, zPrime, y + r, z, i * sides);
		//mesh.addVertex(ofPoint(z,temp.first, temp.second));
	//}

	for (int i = 0; i < num; i++) {
		auto temp = getDegreesRotation(h, k, x, y+r, i * sides);
		mesh.addVertex(ofPoint(temp.first, temp.second, z));
		ofDrawBitmapString((i * sides), mesh.getVertex(mesh.getNumVertices() - 1));
		//std::cout << temp.first;
		//std::cout << ",";
		//std::cout << temp.second;
		//std::cout << ",";
		//std::cout << z;
		//std::cout << "\n";
		counter++;
	}
	for (int i = 1; i < ceil(num/4)+1; i++) {
		auto newLayer = getDegreesRotationZ(h, zPrime, y + r, z, i*sides);

		for (int j = 0; j < num; j++) {
			auto temp = getDegreesRotation(h, k, x, newLayer.first, j * sides);
			mesh.addVertex(ofPoint(temp.first, temp.second, newLayer.second));
			ofDrawBitmapString((j * sides), mesh.getVertex(mesh.getNumVertices() - 1));
			//std::cout << temp.first;
			//std::cout << ",";
			//std::cout << temp.second;
			//std::cout << ",";
			//std::cout << newLayer.second;
			//std::cout << "\n";
			counter++;
		}
	}
	for (int i = 1; i < ceil(num/4)+1; i++) {
		auto newLayer = getDegreesRotationZ(h, zPrime, y + r, z, (i * sides)*-1);

		for (int j = 0; j < num; j++) {
			auto temp = getDegreesRotation(h, k, x, newLayer.first, j * sides);
			mesh.addVertex(ofPoint(temp.first, temp.second, newLayer.second));
			ofDrawBitmapString((j * sides) * -1, mesh.getVertex(mesh.getNumVertices()-1));
			//std::cout << temp.first;
			//std::cout << ",";
			//std::cout << temp.second;
			//std::cout << ",";
			//std::cout << newLayer.second;
			//std::cout << "\n";
			counter++;
		}
	}
	mesh.drawVertices();
	//adding index is driving me crazy

	//first triangle
	
	for (int j = 0; j < ceil(num / 4); j++) {
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
	
	
	for (int i = 0; i < num - 1; i++) {
		mesh.addIndex(i);
		mesh.addIndex(i + (ceil(num / 4) + 1)*num);
		mesh.addIndex(i + 1);
	}
	mesh.addIndex(0);
	mesh.addIndex(num-1 + ((ceil(num / 4) + 1) * num));
	mesh.addIndex(num-1);
	
	
	for (int j = ceil(num / 4) + 1; j < ceil(num / 2); j++) {
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
	
	//second triangle
	
	for (int j = 0; j < ceil(num / 4); j++) {
		float offset = j * num;
		for (int i = 0; i < num - 1; i++) {
			mesh.addIndex(i + 1+ num +offset);
			mesh.addIndex(i + num + offset);
			mesh.addIndex(i + 1 + offset);
		}
		mesh.addIndex(0 + offset);
		mesh.addIndex(num + offset);
		mesh.addIndex((num - 1) + num + offset);
	}
	

	for (int i = 0; i < num - 1; i++) {
		mesh.addIndex(i + 1 + ((ceil(num / 4) + 1) * num));
		mesh.addIndex(i + ((ceil(num / 4) + 1) * num));
		mesh.addIndex(i + 1);
	}
	mesh.addIndex(0);
	mesh.addIndex(num - 1 + ((ceil(num / 4) + 1) * num));
	mesh.addIndex(((ceil(num / 4) + 1)* num));

	for (int j = ceil(num / 4) + 1; j < ceil(num / 2); j++) {
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

	//mesh.addIndex(0);
	//mesh.addIndex(1);
	//mesh.addIndex(30);
	//((num/2) -1)/2
	//mesh.addIndex(0);
	//mesh.addIndex(1);
	//mesh.addIndex(num);


		//for (int i = num; i < (num*2) - 1; i++) {
		//	mesh.addIndex(i + offset);
		//	mesh.addIndex(i + num + offset);
		//	mesh.addIndex(i + 1 + offset);
		//}
		//mesh.addIndex(0 + offset+ num);
		//mesh.addIndex(num - 1 + offset + num);
		//mesh.addIndex((num - 1) + num + offset + num);
	//mesh.addIndex(num);
	//mesh.addIndex(num + 1);
	//mesh.addIndex(num + num);
	for (int i = 0; i < mesh.getNumIndices(); i+=3)
	{
		glEnable(GL_DEPTH_TEST);
		//mesh.addColor(ofColor(ofRandom(255), ofRandom(255), ofRandom(255)));
		//mesh.addColor(ofColor(ofRandom(255), ofRandom(255), ofRandom(255)));
		//mesh.addColor(ofColor(0, 0, 255));
	}
	
	mesh.drawWireframe();
	
	//ofDrawBitmapString("test word",20,20,20);
	
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
			//auto sphericalCoords = cartesianToSpherical(ofPoint(plusRr[i].x - plusR[i].x, plusRr[i].y - plusR[i].y, plusRr[i].z - plusR[i].z));

			//something here needs to change so that the rotation is always in the same direction and starts at the same place
			//has something to do with calcing the difference and ending up with the wrong point but the point is still technically in the circle
			//might have something to do with the spherical coords have the same value twice or something along those lines
			// dfference is biggest - smallest

			//check if next point is in the array, if not, rotate by 180 and add it, then change the roation to be in the opposite direction
			//so if j = 0, check if it contains the point, and blah blah blah go again
			//-2-2
			//		-1,-1
			// 
			// if x y z are less than 0 then *-1, if xyz are less than plusRr[i], do the flip stuff, else do normal stuff
			//ofPoint tempPoint = sphericalToCartesian(sphericalCoords);
			//tempPoint.x + plusR[i].x;
			//tempPoint.y + plusR[i].y;
			//tempPoint.z + plusR[i].z;
				//if (distance(tempPoint,ofPoint(0,0,0))<distance(plusRr[i],ofPoint(0,0,0))){
					//sphericalCoords[2] += (((sides * j) + 180) * (acos(-1) / 180)*-1);
					//auto cartesianCoords = sphericalToCartesian(sphericalCoords);
					//mesh.addVertex(ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
				//}
				//else {
					//ofDrawBitmapString(sphericalCoords, 0, 0, 0);
					//sphericalCoords[2] += sides * (acos(-1) / 180) * j;
					//ofDrawBitmapString(sphericalCoords[2], 0, 0, 0);
					//auto cartesianCoords = sphericalToCartesian(sphericalCoords);
					//ofDrawBitmapString(cartesianCoords,10,10,10);
					//mesh.addVertex(ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
					//ofDrawBitmapString(sphericalCoords[2], ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
				//}
/*
	for (int i = 0; i < num; i++) {
		auto sphericalCoords = cartesianToSpherical(ofPoint(plusRr[i].x - plusR[i].x, plusRr[i].y - plusR[i].y, plusRr[i].z - plusR[i].z));
		auto cartesianCoords = sphericalToCartesian(sphericalCoords);
		if (plusRr[i].x < 0) {
			sphericalCoords[2] += (((sides * 0) + 180) * (acos(-1) / 180) * -1);
			cartesianCoords = sphericalToCartesian(sphericalCoords);
		}
		mesh.addVertex(ofPoint(cartesianCoords.x + plusR[i].x, cartesianCoords.y + plusR[i].y, cartesianCoords.z + plusR[i].z));
	}
	*/
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

}
void drawSphereWithMeshSetup(float x, float y, float z, float r, float sides) {
	////add vertices
	ofPoint center;
	center.x = x;
	center.y = y;
	center.z = z;
	float num = sides;
	sides = 360 / sides;
	z = z - r;
	x = x + r;
	int counter = 0;
	//180
	//360
	for (int phi = 0; phi < 180; phi += sides) {
		for (int theta = 0; theta < 360; theta += sides){
			auto sphericalCoords = cartesianToSpherical(ofPoint(x, y, z));
			sphericalCoords[1] += (phi) * (acos(-1) / 180);
			sphericalCoords[2] += (theta) * (acos(-1) / 180);
			auto cartesianCoords = sphericalToCartesian(sphericalCoords);
			/*if (cartesianCoords.x < 0) {
				sphericalCoords[1] += (180) * (acos(-1) / 180);
				sphericalCoords[1] = (atan(1) * 4) - sphericalCoords[1];
				cartesianCoords = sphericalToCartesian(sphericalCoords);
			}*/
			mesh.addVertex(ofPoint(cartesianCoords.x,cartesianCoords.y,cartesianCoords.z));
			//ofDrawBitmapString(counter, cartesianCoords);
			counter++;
		}
	}
	//(180 / sides)-1
	//
	for (int j = 0; j < (180 / sides) - 1; j++) {
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
	//
	//
	/*for (int i = 0; i < num - 1; i++) {
		mesh.addIndex(((180 / sides)*num-1)+i);
		mesh.addIndex(i + 1);
		mesh.addIndex(i);
	}
	mesh.addIndex(0);
	mesh.addIndex(num - 1);
	mesh.addIndex((num - 1) * num);*/
	//
	/*
	mesh.addIndex(60);
	mesh.addIndex(71);
	mesh.addIndex(4);
	//a-11,a,b+4 ((num/2)-1)-1

	mesh.addIndex(71);
	mesh.addIndex(70);
	mesh.addIndex(5);
	//a,a-1,b+5 ((num/2)-1)
	mesh.addIndex(70);
	mesh.addIndex(69);
	mesh.addIndex(6);
	//a-1,a-2,b+6 ((num/2)-1)+1
	mesh.addIndex(69);
	mesh.addIndex(68);
	mesh.addIndex(7);
	//a-2,a-3,b+7 ((num/2)-1)+2
	mesh.addIndex(68);
	mesh.addIndex(67);
	mesh.addIndex(8);
	//a-3,a-4,b+8 ((num/2)-1)+3
	mesh.addIndex(67);
	mesh.addIndex(66);
	mesh.addIndex(9);
	//a-4,a-5,b+9 ((num/2)-1)+4
	mesh.addIndex(66);
	mesh.addIndex(65);
	mesh.addIndex(10);
	//a-5,a-6,b+10 ((num/2)-1)+5
	mesh.addIndex(65);
	mesh.addIndex(64);
	mesh.addIndex(11);
	//a-6,a-7,b+11 ((num/2)-1)+6
	mesh.addIndex(64);
	mesh.addIndex(63);
	mesh.addIndex(0);
	//a-7,a-8,b 0
	mesh.addIndex(63);
	mesh.addIndex(62);
	mesh.addIndex(1);
	//a-8,a-9,b+1 1
	mesh.addIndex(62);
	mesh.addIndex(61);
	mesh.addIndex(2);
	//a-10,a-11,b+2 2
	mesh.addIndex(61);
	mesh.addIndex(60);
	mesh.addIndex(3);
	//a-11,a-12,b+3 3
	*/

	//=====================PROBLEM BELOW - BAD CODE - DOESNT WORK=========================//

	//((180/sides)*num)-1
	float max = ((180 / sides) * num) - 1;
	float smallCounter = 0;
	for (int i = 0; i < num -1; i++) {
		float firstRing = (num / 2) - 1 + i;
		if (firstRing > num-1) {
			firstRing = smallCounter;
			smallCounter++;
		}
		mesh.addIndex(max - i);
		mesh.addIndex(max - i - 1);
		mesh.addIndex(firstRing);
	}
	mesh.addIndex(max - num + 1);
	mesh.addIndex(max);
	mesh.addIndex(smallCounter);
	//a-11,a,b+3
	
	//========================PROBLEM ABOVE - BAD CODE - DOESNT WORK==============================//

	 
	for (int j = 0; j < (180 / sides) - 1; j++) {
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
	
/*
	//4,5,71 b+4,b+5,a
	mesh.addIndex(4);
	mesh.addIndex(5);
	mesh.addIndex(71);
	//5,6,70 b+5,b+6,a-1
	mesh.addIndex(5);
	mesh.addIndex(6);
	mesh.addIndex(70);
	//6,7,69 b+6,b+7,a-2
	mesh.addIndex(6);
	mesh.addIndex(7);
	mesh.addIndex(69);
	//7,8,68 b+7,b+8,a-3
	mesh.addIndex(7);
	mesh.addIndex(8);
	mesh.addIndex(68);
	//8,9,67 b+8,b+9,a-4
	mesh.addIndex(8);
	mesh.addIndex(9);
	mesh.addIndex(67);
	//9,10,66 b+9,b+10,a-5
	mesh.addIndex(9);
	mesh.addIndex(10);
	mesh.addIndex(66);
	//10,11,65 b+9,b+10,a-6
	mesh.addIndex(10);
	mesh.addIndex(11);
	mesh.addIndex(65);


	//11,0,64 b+10,b,a-7
	mesh.addIndex(11);
	mesh.addIndex(0);
	mesh.addIndex(64);


	//0,1,63 b,b+1,a-8
	mesh.addIndex(0);
	mesh.addIndex(1);
	mesh.addIndex(63);
	//1,2,62 b+1,b+2,a-9
	mesh.addIndex(1);
	mesh.addIndex(2);
	mesh.addIndex(62);
	//2,3,61 b+2,b+3,a-10
	mesh.addIndex(2);
	mesh.addIndex(3);
	mesh.addIndex(61);
	//3,4,60 b+3,b+4,a-11
	mesh.addIndex(3);
	mesh.addIndex(4);
	mesh.addIndex(60);
	*/

	//======================PROBLEM BELOW - BAD CODE - DOESNT WORK=========================//

	max = ((180 / sides) * num) - 1;
	smallCounter = 0;
	for (int i = 0; i <= (num/2); i++) {
		float firstRing = (num / 2) - 1 + i;
			mesh.addIndex(max - i);
			mesh.addIndex(firstRing-1);
			mesh.addIndex(firstRing);
		}

	mesh.addIndex(0);
	mesh.addIndex(max - (num/2) - 1);
	mesh.addIndex(num - 1);

	float reset = 1;
	for (int i = (num / 2)+2; i < num; i++) {
		mesh.addIndex(max - i);
		mesh.addIndex(reset - 1);
		mesh.addIndex(reset);
		reset++;
	}

	//========================PROBLEM ABOVE - BAD CODE - DOESNT WORK==============================//


	for (int i = 0; i < mesh.getNumVertices(); i++) {
		mesh.addNormal(normalize(findVector(center, mesh.getVertex(i))));
	}



	// 
	//mesh.drawVertices();
	//mesh.draw();
}
void createDirectionalLight(ofVec3f lightVector, ofColor color) {
	glEnable(GL_DEPTH_TEST);
	//random colour test
	/*for (int i = 0; i < mesh.getNumVertices() / 3; i++) {
		mesh.addColor(ofColor(ofRandom(255), ofRandom(255), ofRandom(255)));
		mesh.addColor(ofColor(ofRandom(255), ofRandom(255), ofRandom(255)));
		mesh.addColor(ofColor(ofRandom(255), ofRandom(255), ofRandom(255)));

	}*/

	//possibly change it so that as long as its less than 90 then you reduce the colour values by that percentage
	//for (int i = 0; i < mesh.getNumNormals(); i++) {
	//	if (angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4)) < 10) {
	//		mesh.addColor(ofColor(color.r, color.g, color.b));
	//	}
	//	else if (angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4)) < 30) {
	//		mesh.addColor(ofColor((color.r/5)*4, (color.g/5)*4, (color.b/5)*4));
	//	}
	//	else if (angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4)) < 50) {
	//		mesh.addColor(ofColor((color.r / 5) * 3, (color.g / 5) * 3, (color.b / 5) * 3));
	//	}
	//	else if (angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4)) < 70) {
	//		mesh.addColor(ofColor((color.r / 5) * 2, (color.g / 5) * 2, (color.b / 5) * 2));
	//	}
	//	else if (angleBetween(lightVector, mesh.getNormal(i)) * (180 / (atan(1) * 4)) < 90) {
	//		mesh.addColor(ofColor(color.r / 5, color.g / 5, color.b / 5));
	//	}
	//	else {
	//		mesh.addColor(ofColor(0, 0, 0));
	//	}
	//}
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
void diamondStep(float BL, float TR, float BR, float TL, float resolution, float maxChange) {
	float num = floor((TR - TL) / 2);
	float middle = (num * resolution) + BL + num;
	ofPoint temp = mesh.getVertex(middle);
	temp.z = ((mesh.getVertex(BL).z + mesh.getVertex(TR).z + mesh.getVertex(BR).z + mesh.getVertex(TL).z) / 4) + ofRandom(maxChange);
	mesh.setVertex(middle, temp);
}
void squareStep(float BL, float TR, float BR, float TL, float resolution, float maxChange) {
	float num = floor((TR - TL) / 2);
	float LM = BL + (num * resolution);
	float RM = BR + (num * resolution);
	float TM = TL + num;
	float BM = BL + num;
	//middle left
	ofPoint temp = mesh.getVertex(LM);
	temp.z = ((mesh.getVertex(BL).z + mesh.getVertex(TR).z + mesh.getVertex(BR).z + mesh.getVertex(TL).z) / 4) + ofRandom(maxChange);
	mesh.setVertex(LM, temp);
	//middle right
	temp = mesh.getVertex(RM);
	temp.z = ((mesh.getVertex(BL).z + mesh.getVertex(TR).z + mesh.getVertex(BR).z + mesh.getVertex(TL).z) / 4) + ofRandom(maxChange);
	mesh.setVertex(RM, temp);
	//top middle
	temp = mesh.getVertex(TM);
	temp.z = ((mesh.getVertex(BL).z + mesh.getVertex(TR).z + mesh.getVertex(BR).z + mesh.getVertex(TL).z) / 4) + ofRandom(maxChange);
	mesh.setVertex(TM, temp);
	//bottom middle
	temp = mesh.getVertex(BM);
	temp.z = ((mesh.getVertex(BL).z + mesh.getVertex(TR).z + mesh.getVertex(BR).z + mesh.getVertex(TL).z) / 4) + ofRandom(maxChange);
	mesh.setVertex(BM, temp);
}
void diamondSquare(float BL, float TR, float BR, float TL, float resolution, float maxChange) {
	if (TR - TL < 2) {
		return;
	}
	diamondStep(BL, TR, BR, TL, resolution, maxChange);
	squareStep(BL, TR, BR, TL, resolution, maxChange);
	//recursive call function for each small square
	//BL,BM,LM,middle
	//TL,TM,LM,middle
	//TM,TR,middle,RM
	//middle,MR,BM,BR
	float num = floor((TR - TL) / 2);
	float middle = (num * resolution) + BL + num;
	float LM = BL + (num * resolution);
	float RM = BR + (num * resolution);
	float TM = TL + num;
	float BM = BL + num;
	diamondSquare(BL, middle, BM, LM, resolution, maxChange * 0.8);
	diamondSquare(middle, TR, RM, TM, resolution, maxChange * 0.8);
	diamondSquare(BM, RM, BR, middle, resolution, maxChange * 0.8);
	diamondSquare(LM, TM, middle, TL, resolution, maxChange * 0.8);
}
void meshTerrain(float x, float y, float z, float width, float length, float maxChange, float resolution) {
	//need to add a grid of points
	float counter = 0;
	for (int j = 0; j < length; j += length / resolution) {
		for (int i = 0; i < width; i += width / resolution) {
			mesh.addVertex(ofPoint(x + i, y + j, z));
			//ofDrawBitmapString(counter, ofPoint(x + i, y + j, z + ofRandom(maxChange)));
			counter++;
		}
	}
	for (int j = 0; j < resolution - 1; j++) {
		float offset = j * resolution;
		for (int i = 0; i < resolution - 1; i++) {
			mesh.addIndex(i + offset);
			mesh.addIndex(i + offset + 1);
			mesh.addIndex(i + offset + resolution);
		}
	}

	for (int j = 0; j < resolution - 1; j++) {
		float offset = j * resolution;
		for (int i = 0; i < resolution - 1; i++) {
			mesh.addIndex(i + offset + 1);
			mesh.addIndex(i + offset + resolution + 1);
			mesh.addIndex(i + offset + resolution);
		}
	}

	//diamondSquare shenanigans
	//maxChange = 100;
	//setting first point to random height

	float BL = 0;
	float TR = mesh.getNumVertices() - 1;
	float BR = resolution - 1;
	float TL = mesh.getNumVertices() - resolution;
	mesh.setVertex(BL, ofPoint(mesh.getVertex(BL).x,mesh.getVertex(BL).y,mesh.getVertex(BL).z+ofRandom(maxChange)));
	//setting the last point to a random height
	mesh.setVertex(TR, ofPoint(mesh.getVertex(TR).x, mesh.getVertex(TR).y, mesh.getVertex(TR).z + ofRandom(maxChange)));
	//setting the bottom right point to a random height
	mesh.setVertex(BR, ofPoint(mesh.getVertex(BR).x, mesh.getVertex(BR).y, mesh.getVertex(BR).z + ofRandom(maxChange)));
	//setting the top left to a random height
	mesh.setVertex(TL, ofPoint(mesh.getVertex(TL).x, mesh.getVertex(TL).y, mesh.getVertex(TL).z + ofRandom(maxChange)));

	diamondSquare(BL,TR,BR,TL,resolution,maxChange);

	for (int j = 0; j < resolution -1; j++) {
		float offset = j * resolution;
		for (int i = 0; i < resolution -1; i++) {
			mesh.addNormal(20 * normalize(crossProduct(findVector(mesh.getVertex(i + offset), mesh.getVertex(i + offset + 1)), findVector(mesh.getVertex(i + offset), mesh.getVertex(i + offset + resolution)))));
		}
		mesh.addNormal(20 * normalize(crossProduct(findVector(mesh.getVertex(resolution + offset -1), mesh.getVertex(resolution + resolution + offset - 2)), findVector(mesh.getVertex(resolution + offset - 1), mesh.getVertex(resolution + offset - 2)))));
	}

	for (int i = (resolution -1)*resolution; i < (resolution*resolution)-1; i++) {
		mesh.addNormal(20 * normalize(crossProduct(findVector(mesh.getVertex(i), mesh.getVertex(i-resolution)), findVector(mesh.getVertex(i), mesh.getVertex(i-resolution+1)))));
	}
	mesh.addNormal(20 * normalize(crossProduct(findVector(mesh.getVertex((resolution * resolution) - 1), mesh.getVertex((resolution * resolution) - 1 -1)), findVector(mesh.getVertex((resolution * resolution) - 1), mesh.getVertex((resolution * resolution) - 1)))));

}
ofPoint midpoint3d(float x1,float y1,float z1,float x2,float y2,float z2) {
	ofPoint output;

	//(x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2
	output.x = (x1 + x2) / 2;
	output.y = (y1 + y2) / 2;
	output.z = (z1 + z2) / 2;

	return output;
}
void fractalMountain(float x1,float y1,float z1, float x2, float y2, float z2, float x3, float y3, float z3, float resolution, float maxHeight) {
	if (resolution == 0) {
	/*	for (int i = 0; i < mesh.getNumVertices(); i++) {
			float random1 = ofRandom(maxHeight);
			if (random1 > (maxHeight / 2)) {
				random1 -= maxHeight;
			}
			mesh.setVertex(i, ofVec3f(mesh.getVertex(i).x, mesh.getVertex(i).y, mesh.getVertex(i).z + ofRandom(random1)));
		}*/
		return;
	}
	//draws the triangle//
	ofColor randCol;
	randCol.r = ofRandom(255);
	randCol.g = ofRandom(255);
	randCol.b = ofRandom(255);
	mesh.addVertex(ofPoint(x1, y1, z1));
	mesh.addIndex(mesh.getNumVertices() - 1);
	mesh.addColor(randCol);
	mesh.addVertex(ofPoint(x2, y2, z2));
	mesh.addIndex(mesh.getNumVertices() - 1);
	mesh.addColor(randCol);
	mesh.addVertex(ofPoint(x3, y3, z3));
	mesh.addIndex(mesh.getNumVertices() - 1);
	mesh.addColor(randCol);
	/////////////////////////
	//if (mesh.getNumVertices() > 3) {
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addColor(randCol);
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addColor(randCol);
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addIndex(mesh.getNumVertices());
	//	mesh.addColor(randCol);
	//}
	resolution--;
	ofPoint mid1 = midpoint3d(x1, y1, z1, x2, y2, z2);
	//mesh.addVertex(mid1);
	//mesh.addIndex(mesh.getNumVertices() - 6);
	//mesh.addIndex(mesh.getNumVertices() - 5);
	//mesh.addIndex(mesh.getNumVertices() - 1);
	//mesh.addColor(randCol);
	ofPoint mid2 = midpoint3d(x2, y2, z2, x3, y3, z3);
	//mesh.addVertex(mid2);
	//mesh.addIndex(mesh.getNumVertices() - 5);
	//mesh.addIndex(mesh.getNumVertices() - 6);
	//mesh.addIndex(mesh.getNumVertices() - 1);
	//mesh.addColor(randCol);
	ofPoint mid3 = midpoint3d(x1, y1, z1, x3, y3, z3);
	//mesh.addVertex(mid3);
	//mesh.addIndex(mesh.getNumVertices() - 6);
	//mesh.addIndex(mesh.getNumVertices() - 4);
	//mesh.addIndex(mesh.getNumVertices() - 1);
	//mesh.addColor(randCol);
	////////////=============////////////////==============//////////
	//add midpoint to each point around it as a index
	////////////===========//////////////////==============//////////

	//float random1 = ofRandom(maxHeight);
	//if (random1 > (maxHeight/2)) {
	//	random1 -= maxHeight;
	//}
	//mid1.z += random1;

	//float random2 = ofRandom(maxHeight);
	//if (random2 > (maxHeight / 2)) {
	//	random2 -= maxHeight;
	//}
	//mid2.z += random2;

	//float random3 = ofRandom(maxHeight);
	//if (random3 > (maxHeight / 2)) {
	//	random3 -= maxHeight;
	//}
	//mid3.z += random3;

	maxHeight = maxHeight * 0.5;
	fractalMountain(mid1.x, mid1.y, mid1.z, mid2.x, mid2.y, mid2.z, mid3.x, mid3.y, mid3.z, resolution, maxHeight);
	fractalMountain(mid1.x, mid1.y, mid1.z, x2,y2,z2, mid2.x, mid2.y, mid2.z, resolution, maxHeight);
	fractalMountain(x1, y1, z1, mid1.x, mid1.y, mid1.z, mid3.x, mid3.y, mid3.z, resolution, maxHeight);
	fractalMountain(mid3.x, mid3.y, mid3.z, mid2.x, mid2.y, mid2.z, x3, y3, z3, resolution, maxHeight);

	//plan to do all of the randomness after?
	//also need to change code to do one whole thing at once then recurse i think maybe? would be very hard 
}
void fractalMountain2ElectricBoogaloo(float x1, float y1, float z1,float width, float length, float resolution, float maxHeight) {
	glEnable(GL_DEPTH_TEST);
	for (int j = 0; j < length; j += length / resolution) {
		x1 += (length / resolution);
		for (int i = 0; i < width; i += width / resolution) {
			mesh.addVertex(ofPoint(x1 + i, y1 + j, z1));
		}
	}

	for (int j = 0; j < resolution - 1; j++) {
		float offset = j * resolution;
		for (int i = 0; i < resolution - 1; i++) {
			mesh.addIndex(i + offset);
			mesh.addIndex(i + offset + 1);
			mesh.addIndex(i + offset + resolution);
		}
	}

	for (int j = 0; j < resolution - 1; j++) {
		float offset = j * resolution;
		for (int i = 0; i < resolution - 1; i++) {
			mesh.addIndex(i + offset + 1);
			mesh.addIndex(i + offset + resolution + 1);
			mesh.addIndex(i + offset + resolution);
		}
	}
	float temp = ofRandom(mesh.getNumVertices());
	mesh.setVertex(((resolution/4)*(resolution/4)), ofVec3f(mesh.getVertex(((resolution / 4) * (resolution / 4))).x, mesh.getVertex(((resolution / 4) * (resolution / 4))).y, mesh.getVertex(((resolution / 4) * (resolution / 4))).z + maxHeight+3));

	temp = ofRandom(mesh.getNumVertices());
	mesh.setVertex(((resolution / 2) * (resolution / 2)), ofVec3f(mesh.getVertex(((resolution / 2) * (resolution / 2))).x, mesh.getVertex(((resolution / 2) * (resolution / 2))).y, mesh.getVertex(((resolution / 2) * (resolution / 2))).z - maxHeight*1.5));

		for (int i = 1; i < mesh.getNumVertices()-resolution; i++) {
		float random1 = ofRandom(maxHeight);
		if (random1 > (maxHeight / 2)) {
			random1 -= maxHeight/2;
		}
		//maxHeight = maxHeight * 0.9;
		mesh.setVertex(i, ofVec3f(mesh.getVertex(i).x, mesh.getVertex(i).y, (mesh.getVertex(i-1).z+mesh.getVertex(i+resolution).z+mesh.getVertex(i).z)/3 + ofRandom(random1)));
	}

		for (int i = mesh.getNumVertices() - resolution; i < mesh.getNumVertices(); i++) {
			float random1 = ofRandom(maxHeight);
			if (random1 > (maxHeight / 2)) {
				random1 -= maxHeight / 2;
			}
			mesh.setVertex(i, ofVec3f(mesh.getVertex(i).x, mesh.getVertex(i).y, (mesh.getVertex(i - 1).z + mesh.getVertex(i - resolution).z + mesh.getVertex(i).z) / 3 + ofRandom(random1)));
		}

		for (int i = 0; i < mesh.getNumVertices()-2; i++) {
			ofColor randCol;
			randCol.r = ofRandom(255);
			randCol.g = ofRandom(255);
			randCol.b = ofRandom(255);
			mesh.addColor(randCol);
			mesh.addColor(randCol);
			mesh.addColor(randCol);
		}
		

}
//--------------------------------------------------------------
void ofApp::setup(){
	//drawSphereWithMesh(0, 0, 0, 100, 18);
	glPointSize(4);
	//drawSphereWithMeshSetup(0, 0, 0, 100, 12);
	//drawTorusWithMeshSetup(0, 0, 0, 100, 200, 240);
	//createDirectionalLight(ofVec3f(0,1,1),ofColor(255, 180, 128));
	/*
	ofLog() << "vertices:" << mesh.getNumVertices(); 
	for (int i = 0; i < mesh.getNumIndices(); i++) {
		ofLog() << mesh.getIndex(i);
	}
	*/
	//meshTerrain(-50, -50, 0, 100, 100, 5, 20);

	//-1000,-500,0
	//0,500,0
	//1000,-500,0
	//fractalMountain(-1000, -500, 0, 0, 500, 0, 1000, -500, 0, 4, 1000);
	//fractalMountain2ElectricBoogaloo(-1000, -500, 0, 2000, 1000, 40, 200);
	//meshTerrain(-2450, -2450, 0, 4900, 4900, 500, 49);
	meshTerrain(-450, -450, 0, 900, 900, 100, 9);
	createDirectionalLight(ofVec3f(0,1,1),ofColor(178, 190, 181));
}	

//--------------------------------------------------------------
void ofApp::update(){
}

//--------------------------------------------------------------
void ofApp::draw(){
	cam.begin();
	ofNoFill();	
	//ofLight light;
	//light.setPosition(0, 0, 0);
	//light.enable();
	//drawLineUsingFormula(64, 64, 128, 256);
	//drawSquareRotationOrigin(300, 300, 64, 0);
	//drawSquareRotationOrigin(600, 300, 64, 0);
	//drawSquareRotationOrigin(600, 300, 64, 45);
	//300,300 364,300
	//DrawSquareRotationArbitrary(600, 300, 64, 45,64,64);
	//765
	//dot(0, 765);
	//dot(1, 765);
	//dot(0.5, 654);

	//drawEqualatralTriangle(128, 128, 60);
	//drawLineUsingFormula(64, 64, 74, 64);
	//drawLineUsingFormula(64, 64, 69, 64 - 8.66);
	//drawLineUsingFormula(74, 64, 69, 64 - 8.66);

	//drawLineUsingFormula(128, 128, 178, 128);
	//drawLineUsingFormula(128, 128, 153, 128 - 43.3);
	//drawLineUsingFormula(178, 128, 153, 128 - 43.3);

	//sepinski(0,765,765,765, 382.5, 765 - sqrt((765 * 765) - ((765 / 2) * (765 / 2))),3);
	//sepinski(10,750,700,750,300,500,10);
	//float value = 2220;
	//float offset = 500;
	//sepinski(offset, value - (offset/2), value + offset, value - (offset / 2), value / 2 + offset, value - sqrt((value * value) - ((value / 2) * (value / 2))) - (offset / 2), 8);

	//drawCircle(300,300,50);
	
	//drawCylinder(300, 300, 50, 300,500);

	//drawCone(300, 300, 50, 500, 400);
	
	//drawpyramid(300, 300, 50, 500, 400);
	//testing3D();
	//ofNoFill();
	//ofDrawSphere(64);
	//ofDrawCircle(0,0,72);
	//ofRotateXDeg(rotate);

	//drawCylinderWithStackedCircles(64, 64, 100, 400, "horizontal");

	
	//drawMeshCube(0, 0, 0, 100);

	//drawCircleWithMesh(0,0,0,100,20);

	//drawCylinderWithMesh(0, 0, 0, 100, 15, 300);
	
	//drawSphereWithMesh(0, 0, 0, 100,24);
	//drawTorusWithMeshSetup(0, 0, 0, 100, 200, 24);
		
	//mesh.drawWireframe();

	/*mesh.drawVertices();*/

	mesh.draw();
	//ofSetColor(ofColor(255, 0, 0));
	//for (int i = 0; i < mesh.getNumVertices(); i++) {
	//	ofDrawLine(mesh.getVertex(i), ofPoint(mesh.getVertex(i).x + mesh.getNormal(i).x, mesh.getVertex(i).y + mesh.getNormal(i).y, mesh.getVertex(i).z + mesh.getNormal(i).z));
	//}
	//ofSetColor(ofColor(255, 255, 255));
	//for (int i = 0; i < mesh.getNumVertices(); i++) {
	//	ofDrawBitmapString(i, mesh.getVertex(i));
	//}
	//drawSphereWithMeshSetup(0, 0, 0, 100, 12);
	//meshTerrain(0, 0, 0, 100, 100, 20, 10);
	cam.end();
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
