#include"rendering_API.h"
#include"maths.h"
#include<vector>
#include<cstring>
#include<iostream>
#include<fstream>
#define IMGUI_DEFINE_MATH_OPERATORS
std::vector<int> ReadStartFromTxt(const std::string& filePath) {
	std::vector<int> numbers(0);
	std::ifstream file(filePath);

	if (!file.is_open()) {
		std::cerr << "Error: Unable to open file " << filePath << std::endl;
		return numbers;
	}

	int number;
	while (file >> number) {
		numbers.push_back(number * 4);
	}

	file.close();
	return numbers;
}
std::vector<ks::vec3> ReadVertFromTxt(const std::string& filePath) {
	std::vector<ks::vec3> verts(0);
	std::ifstream file(filePath);

	if (!file.is_open()) {
		std::cerr << "Error: Unable to open file " << filePath << std::endl;
		return verts;
	}

	float x, y, z;
	while (file >> x >> y >> z) {
		verts.push_back(ks::vec3(x, y, z) * 30);
	}

	file.close();
	return verts;
}

struct yarnFile
{

	struct BCCheader
	{
		char sign[3];
		unsigned char byteCount;
		char curveType[2];
		char dimensions;
		char upDimension;
		uint64_t curveCount;
		uint64_t totalControlPointCount;
		char fileInfo[40];
	};

	int LoadFromFile(const char* filename)
	{
		struct BCCheader header_;
		FILE* pFile = fopen(filename, "rb");
		fread(&header_, sizeof(BCCheader), 1, pFile);

		if (header_.sign[0] != 'B')
			return -1; // Invalid file signature
		if (header_.sign[1] != 'C')
			return -1; // Invalid file signature
		if (header_.sign[2] != 'C')
			return -1; // Invalid file signature
		if (header_.byteCount != 0x44)
			return -1; // Only supporting 4-byte integers and floats

		if (header_.curveType[0] != 'C')
			return -1; // Not a Catmull-Rom curve
		if (header_.curveType[1] != '0')
			return -1; // Not uniform parameterization
		if (header_.dimensions != 3)
			return -1; // Only curves in

		std::vector<ks::vec3> controlPoints;
		controlPoints.clear();
		std::vector<int> firstControlPoint(header_.curveCount + 1);
		std::vector<char> isCurveLoop(header_.curveCount);
		//float *cp = (float *)controlPoints.data();
		int prevCP = 0;
		int sum_p = 0;
		for (uint64_t i = 0; i < header_.curveCount; i++) {
			int curveControlPointCount;
			fread(&curveControlPointCount, sizeof(int), 1, pFile);
			isCurveLoop[i] = curveControlPointCount < 0;
			if (curveControlPointCount < 0)
				curveControlPointCount = -curveControlPointCount;
			//std::cout << curveControlPointCount << std::endl;  
			sum_p += curveControlPointCount;
			for (int j = 0; j < curveControlPointCount; j++) {
				float x, y, z;
				fread(&x, sizeof(float), 1, pFile);
				fread(&y, sizeof(float), 1, pFile);
				fread(&z, sizeof(float), 1, pFile);
				//std::cout << x << " " << y << " " <<z <<std::endl;
				controlPoints.push_back(ks::vec3(x / 3, y / 3, z / 3));
			}
			//fread(cp, sizeof(float), curveControlPointCount*3, pFile);
			//cp += curveControlPointCount;
			firstControlPoint[i] = prevCP;
			prevCP += curveControlPointCount;
		}
		std::cout << sum_p << " is supposed to be " << header_.totalControlPointCount << std::endl;
		std::cout << " curvecount is " << header_.curveCount << std::endl;
		firstControlPoint[header_.curveCount] = prevCP;
		controlPoints_.assign(controlPoints.begin(), controlPoints.end());
		firstControlPoint_.assign(firstControlPoint.begin(), firstControlPoint.end());
		isCurveLoop_.assign(isCurveLoop.begin(), isCurveLoop.end());
		header = header_;
		fclose(pFile);
		return 0;
	}
	std::vector<ks::vec3> controlPoints_;
	std::vector<int> firstControlPoint_;
	std::vector<char> isCurveLoop_;
	struct BCCheader header;
};

std::string WORK_PATH;
int main() {
	std::vector<int> yarn_start = ReadStartFromTxt("./rib/yarn_start.txt");
	std::vector<ks::vec3> yarn_ctrP = ReadVertFromTxt("./rib/spline_points1.txt");

	//yarnFile Yarn_file;
	//Yarn_file.LoadFromFile("D:/StitchMesh4.0/glove.bcc");

	//std::string filename = "D:/StitchMesh4.0/teapot_tiny.poly";
	//std::ifstream file(filename);
	//std::vector<ks::vec3> points;
	//if (!file.is_open()) {
	//	std::cerr << "Failed to open file: " << filename << std::endl;
	//}

	//std::string line;
	//// Skip the header "POINTS"
	//std::getline(file, line);

	//std::vector<ks::vec3> verts(0);
	//std::vector<bool>  vert_edge(0);
	//while (std::getline(file, line)) {
	//	std::istringstream iss(line);
	//	ks::vec3 point;
	//	char colon; // To capture the colon after ID
	//	int id;
	//	if (iss >> id >> colon >> point.x() >> point.z() >> point.y()) {
	//		verts.push_back(point / 2);
	//		vert_edge.push_back(true);
	//	}
	//	else break;
	//}

	//int new_one = true;
	//std::vector<int> yarn_start(0);
	//std::vector<ks::vec3> yarn_ctrP(0);
	//while (std::getline(file, line)) {
	//	std::istringstream iss(line);
	//	char colon;
	//	int id;
	//	int vertex1, vertex2;
	//	if (iss >> id >> colon >> vertex1 >> vertex2) {
	//		if ((vert_edge[vertex1 - 1]) || (vert_edge[vertex2 - 1])) {
	//			if (vert_edge[vertex1 - 1] == false) {
	//				vert_edge[vertex2 - 1] == false;
	//				yarn_ctrP.push_back(verts[vertex2 - 1]);
	//			}
	//			else if (vert_edge[vertex2 - 1] == false) {
	//				vert_edge[vertex1 - 1] == false;
	//				yarn_ctrP.push_back(verts[vertex1 - 1]);
	//			}
	//			else {
	//				yarn_start.push_back(yarn_ctrP.size());
	//				yarn_ctrP.push_back(verts[vertex1 - 1]);
	//				yarn_ctrP.push_back(verts[vertex2 - 1]);
	//			}
	//			vert_edge[vertex1 - 1] = false;
	//			vert_edge[vertex2 - 1] = false;
	//		}
	//		else {
	//			yarn_ctrP.push_back(verts[vertex2 - 1]);
	//		}
	//	}
	//}
	//yarn_start.push_back(yarn_ctrP.size());
	//file.close();

	GLWidget glwidget;
	glwidget.initializeGL();
	glwidget.render_task((float*)&yarn_ctrP[0], (int*)&yarn_start[0], yarn_start.size());
	//glwidget.render_task((float*)&Yarn_file.controlPoints_[0], (int*)&Yarn_file.firstControlPoint_[0], Yarn_file.firstControlPoint_.size());

	return 0;
}