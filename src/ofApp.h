#pragma once

#include "ofMain.h"
#include "ofxGui.h"

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

		ofEasyCam cam;
		ofSpherePrimitive sphere;
		ofTexture mTex;
		ofLight light;

		ofxPanel gui;
		ofxPanel torusGui;
		ofxPanel sphereGui;
		ofxPanel directionalLightGui;
		ofxPanel ambientLightGui;
		ofxPanel shiftGui;
		ofxPanel sinShiftGui;
		ofxPanel cosShiftGui;
		ofxPanel tanShiftGui;

		ofxIntSlider thetaSinAmplitudeSlider;
		ofxIntSlider phiSinAmplitudeSlider;
		ofxIntSlider thetaCosAmplitudeSlider;
		ofxIntSlider phiCosAmplitudeSlider;
		ofxIntSlider thetaTanAmplitudeSlider;
		ofxIntSlider phiTanAmplitudeSlider;
		ofxIntSlider thetaSinFrequencySlider;
		ofxIntSlider phiSinFrequencySlider;
		ofxIntSlider thetaCosFrequencySlider;
		ofxIntSlider phiCosFrequencySlider;
		ofxIntSlider thetaTanFrequencySlider;
		ofxIntSlider phiTanFrequencySlider;

		ofxToggle sphereToggle;
		ofxToggle torusToggle;

		ofxToggle sidesAnimationToggle;

		ofxToggle drawWireframeToggle;
		ofxToggle drawToggle;
		ofxToggle normalsToggle;

		ofxToggle shiftToggle;

		ofxToggle thetaSinShiftToggle;
		ofxToggle phiSinShiftToggle;
		ofxToggle thetaSinAnimation;
		ofxToggle phiSinAnimation;

		ofxToggle thetaCosShiftToggle;
		ofxToggle phiCosShiftToggle;
		ofxToggle thetaCosAnimation;
		ofxToggle phiCosAnimation;

		ofxToggle thetaTanShiftToggle;
		ofxToggle phiTanShiftToggle;
		ofxToggle thetaTanAnimation;
		ofxToggle phiTanAnimation;

		ofxToggle directionalLightToggle;
		ofxToggle ambientLightToggle;

		ofxToggle sinToggle;
		ofxToggle cosToggle;
		ofxToggle tanToggle;


		ofxLabel normalsLabel;
		ofxLabel shiftLabel;
		ofxLabel thetaSinShiftLabel;
		ofxLabel phiSinShiftLabel;
		ofxLabel thetaCosShiftLabel;
		ofxLabel phiCosShiftLabel;
		ofxLabel thetaTanShiftlabel;
		ofxLabel phiTanShiftlabel;
		ofxLabel directionalLightLabel;
		ofxLabel ambientLightLabel;

		ofxIntField radius;
		ofxIntField numSides;
		ofxIntField numStacks;

		ofxVec3Slider originSlider;
		ofxVec3Slider colourForDirectionalLight;
		ofxVec3Slider vectorForDirectionalLight;
		ofxVec3Slider colourForAmbientLight;
};
