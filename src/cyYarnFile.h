#pragma once
#include "maths.h"
using namespace ks;
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

    int LoadFromFile(const char *filename)
    {
        struct BCCheader header_;
        FILE *pFile = fopen(filename, "rb");
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
                controlPoints.push_back(ks::vec3(x/3 , y/3, z/3));
            }
            //fread(cp, sizeof(float), curveControlPointCount*3, pFile);
            //cp += curveControlPointCount;
            firstControlPoint[i] = prevCP;
            prevCP += curveControlPointCount;
        }
        std::cout << sum_p << " is supposed to be " << header_.totalControlPointCount<< std::endl;   
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
