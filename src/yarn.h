//
//  yarn.hpp
//  KnittingCapturer
//
//  Created by kui wu on 11/21/15.
//  Copyright © 2015 kui wu. All rights reserved.
//

#ifndef yarn_hpp
#define yarn_hpp

#include <stdio.h>
#include <iostream>
#include "cyPoint.h"

class yarn {
public:
    cyPoint3f** _ctrlPts;   // always 4
    cyPoint3f* _bezierPts;
    cyPoint4f* _renderPts;
    
    int _ctrlPtsNum;
    int _segNum;
    int _subSegNum;
    int _bezierPtsNum;
    int _renderPtsNum;
    
    yarn():_ctrlPtsNum(0),_segNum(0),_subSegNum(11) {};
    ~yarn() {};
    
    void InitBezierCurve(float *ctrlpts, int ctrlPtsNum, int segNum, float offset);

    void GenerateBezierPoints();
};

#endif /* yarn_hpp */
