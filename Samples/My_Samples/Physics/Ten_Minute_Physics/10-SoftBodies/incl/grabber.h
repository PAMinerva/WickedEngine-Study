#pragma once

#include "softBody.h"
#include "wiScene.h"
#include <DirectXMath.h>

class Grabber
{
public:
    // void start(wi::primitive::Ray ray, wi::scene::PickResult pick, SoftBodyInstance* sbInstance);
    void start(wi::scene::PickResult pick, SoftBody* softBody);
    void move(wi::primitive::Ray ray, wi::scene::PickResult pick);
	// void move(wi::scene::PickResult pick);
    void end();

	void increaseTime(float dt);

    wi::scene::PickResult pickInfo {};
    // SoftBodyInstance* sbInstance = nullptr;
	SoftBody* sbInstance = nullptr;
    std::vector<float> prevPos = {0, 0, 0};
    std::vector<float> vel = {0, 0, 0};
    float time = 0.0f;
    bool isGrabbing = false;
	float startDistance = 0.0f;
};
