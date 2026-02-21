#include <vector>
#include <wiPrimitive.h>
#include "grabber.h"
#include "cloth.h"
#include "simulation_utils.h"

void Grabber::start(wi::scene::PickResult pick, ClothMesh* clothMesh)
{
    if (isGrabbing || clothMesh == nullptr)
        return;

    pickInfo = pick;
    this->clothMesh = clothMesh;

	this->clothMesh->StartGrab(pick);

	std::vector<float> pos {pick.position.x, pick.position.y, pick.position.z};
    prevPos = pos;
    vel = {0, 0, 0};
    time = 0.0f;
	startDistance = pick.distance;

    isGrabbing = true;

	if (gPhysicsScene.paused)
		gPhysicsScene.Run();
}

void Grabber::move(wi::primitive::Ray ray, wi::scene::PickResult pick)
{
    if (!isGrabbing || clothMesh == nullptr)
        return;

	XMVECTOR origin = XMLoadFloat3(&ray.origin);
	XMVECTOR direction = XMLoadFloat3(&ray.direction);
	XMVECTOR posVec = XMVectorAdd(origin, XMVectorScale(direction, startDistance));
	XMFLOAT3 posF;
	XMStoreFloat3(&posF, posVec);
	std::vector<float> pos = {posF.x, posF.y, posF.z};

	vel = {pos};
	vel[0] -= prevPos[0];
	vel[1] -= prevPos[1];
	vel[2] -= prevPos[2];

	if (time > 0.0f)
	{
		vel[0] /= time;
		vel[1] /= time;
		vel[2] /= time;
	}
	else
		vel = {0, 0, 0};

    prevPos = {pos};
    time = 0.0f;

    // Move the grabbed vertex
    this->clothMesh->MoveGrabbed(pos, vel);
}

void Grabber::end()
{
    if (isGrabbing && clothMesh != nullptr)
    {
        std::vector<float> posVec = {prevPos};
        std::vector<float> velVec = {vel};
        clothMesh->EndGrab(posVec, velVec);
        clothMesh = nullptr;
        isGrabbing = false;
    }
}

void Grabber::increaseTime(float dt) { time += dt; }
