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
        numbers.push_back(number);
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

std::string WORK_PATH;
int main() {
    std::vector<int> yarn_start = ReadStartFromTxt("D:/Rendering/data_yarn/yarn_start.txt");
    std::vector<ks::vec3> yarn_ctrP = ReadVertFromTxt("D:/Rendering/data_yarn/spline_points.txt");
    GLWidget glwidget;
    glwidget.initializeGL();
    glwidget.render_task((float*)&yarn_ctrP[0], (int*)&yarn_start[0], yarn_start.size());

    return 0;
}