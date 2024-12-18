//
//  yarn.cpp
//  KnittingCapturer
//
//  Created by kui wu on 11/21/15.
//  Copyright © 2015 kui wu. All rights reserved.
//

#include "yarn.h"

cyPoint3f v3max(cyPoint3f a, cyPoint3f b) {
	cyPoint3f result;
	result.x = fmax(a.x, b.x);
	result.y = fmax(a.y, b.y);
	result.z = fmax(a.z, b.z);
	return result;
}

cyPoint3f v3min(cyPoint3f a, cyPoint3f b) {
	cyPoint3f result;
	result.x = fmin(a.x, b.x);
	result.y = fmin(a.y, b.y);
	result.z = fmin(a.z, b.z);
	return result;
}

void yarn::InitBezierCurve(float *ctrlpts, int ctrlPtsNum, int segNum, float offset) {
    _ctrlPtsNum = ctrlPtsNum;
    _segNum = segNum;
    
    if (_segNum < 1) {
        std::cerr << "segment number cannot be less than 1\n";
    }
    
    float initialCtrlPoints[4][3];
    for (int i = 0; i < 12; i++) {
        int j = i/3;
        int k = i%3;
        initialCtrlPoints[j][k] = ctrlpts[i];
    }
    
    int stitchNumPerRow = segNum;
    _ctrlPts = new cyPoint3f*[stitchNumPerRow * 2];
    
    // first two segment
    cyPoint3f* ctrlpoint0 = new cyPoint3f[4];
    cyPoint3f* ctrlpoint1 = new cyPoint3f[4];
    
    _ctrlPts[0] = ctrlpoint0;
    _ctrlPts[1] = ctrlpoint1;
    
    for (int i = 0; i < 4; i++) {
        _ctrlPts[0][i].Set(initialCtrlPoints[i][0], initialCtrlPoints[i][1], initialCtrlPoints[i][2]);
        _ctrlPts[1][i] = ctrlpoint0[i];
    }
    
    // reflection
    for (int i = 0; i < 3; i++)
        _ctrlPts[1][i][0] = 2*_ctrlPts[1][3][0] - _ctrlPts[1][i][0];
    
    // reorder
    cyPoint3f tmp = _ctrlPts[1][0];
    _ctrlPts[1][0] = _ctrlPts[1][3];
    _ctrlPts[1][3] = tmp;
    tmp = _ctrlPts[1][1];
    _ctrlPts[1][1] = _ctrlPts[1][2];
    _ctrlPts[1][2] = tmp;
    
    for (int sn = 1; sn < segNum; sn++) {
        ctrlpoint0 = new cyPoint3f[4];
        ctrlpoint1 = new cyPoint3f[4];
        _ctrlPts[sn*2] = ctrlpoint0;
        _ctrlPts[sn*2+1] = ctrlpoint1;
        cyPoint3f rotateOrgin = _ctrlPts[sn*2-1][3];
        
        for (int i = 0; i < 4; i++) {
            _ctrlPts[sn*2+1][i] = _ctrlPts[sn*2-2][3-i] - rotateOrgin;
            _ctrlPts[sn*2][i] = _ctrlPts[sn*2-1][3-i] - rotateOrgin;
        }
        
        for (int i = 0; i < 4; i++) {
            _ctrlPts[sn*2+1][i][0] = -_ctrlPts[sn*2+1][i][0];
            _ctrlPts[sn*2+1][i][1] = -_ctrlPts[sn*2+1][i][1];
            _ctrlPts[sn*2][i][0] = -_ctrlPts[sn*2][i][0];
            _ctrlPts[sn*2][i][1] = -_ctrlPts[sn*2][i][1];
        }
        
        for (int i = 0; i < 4; i++) {
            _ctrlPts[sn*2+1][i] = _ctrlPts[sn*2+1][i] + rotateOrgin;
            _ctrlPts[sn*2][i] = _ctrlPts[sn*2][i] + rotateOrgin;
        }
    }
    
    // move by offset
    for (int i = 0; i < 2*segNum; i++) {
        for (int j = 0; j < 4; j++) {
            _ctrlPts[i][j][1] = _ctrlPts[i][j][1] + offset;
        }
    }
    
    GenerateBezierPoints();
    
    // generate render pts
    _renderPtsNum = 0;
    for (int i = 0; i < 2*_segNum; i++) {
        _renderPtsNum += 4;
    }
    
    _renderPts = new cyPoint4f[_renderPtsNum];
    int pn = 0;
	float sumLength = 0;
	int subSegNum = 11;
    for (int i = 0; i < 2*_segNum; i++) {

		float localLength = 0.0;
		for (int j = 0; j < subSegNum - 1; j++) {
			float t0 = j / float(subSegNum - 1);
			float t1 = (j + 1) / float(subSegNum - 1);
			cyPoint3f p0 = pow(1 - t0, 3) * _ctrlPts[i][0] +
				3 * pow(1 - t0, 2) * t0 *_ctrlPts[i][1] +
				3 * pow(t0, 2) * (1 - t0) *_ctrlPts[i][2] +
				pow(t0, 3) * _ctrlPts[i][3];

			cyPoint3f p1 = pow(1 - t1, 3) * _ctrlPts[i][0] +
				3 * pow(1 - t1, 2) * t1 * _ctrlPts[i][1] +
				3 * pow(t1, 2) * (1 - t1) * _ctrlPts[i][2] +
				pow(t1, 3) * _ctrlPts[i][3];

			cyPoint3f diff = p0 - p1;

			localLength += diff.Length();
		}


		_renderPts[pn++].Set(_ctrlPts[i][0][0], _ctrlPts[i][0][1], _ctrlPts[i][0][2], sumLength);
		_renderPts[pn++].Set(_ctrlPts[i][1][0], _ctrlPts[i][1][1], _ctrlPts[i][1][2], sumLength);
		_renderPts[pn++].Set(_ctrlPts[i][2][0], _ctrlPts[i][2][1], _ctrlPts[i][2][2], sumLength);
		_renderPts[pn++].Set(_ctrlPts[i][3][0], _ctrlPts[i][3][1], _ctrlPts[i][3][2], sumLength + localLength);

		sumLength += localLength;
    }

}

void yarn::GenerateBezierPoints() {
    _bezierPtsNum = (_subSegNum-1)*2*_segNum+1;
    _bezierPts = new cyPoint3f[_bezierPtsNum];
    
    int i, j, ptsCount = 0;
    float t;
    for (i = 0; i < 2*_segNum; i++) {
        for (j = 0; j < _subSegNum-1; j++) {
            t = j/float(_subSegNum-1);
            _bezierPts[ptsCount++] = pow(1-t, 3) * _ctrlPts[i][0] +
                                        3 * pow(1-t, 2) * t * _ctrlPts[i][1] +
                                        3 * pow(t, 2) * (1-t) * _ctrlPts[i][2] +
                                        pow(t, 3) * _ctrlPts[i][3];
        }
        
    }
    t = 1.0;
    i = 2*_segNum - 1;
    _bezierPts[ptsCount] = pow(1-t, 3) * _ctrlPts[i][0] +
                            3 * pow(1-t, 2) * t * _ctrlPts[i][1] +
                            3 * pow(t, 2) * (1-t) * _ctrlPts[i][2] +
                            pow(t, 3) * _ctrlPts[i][3];
}