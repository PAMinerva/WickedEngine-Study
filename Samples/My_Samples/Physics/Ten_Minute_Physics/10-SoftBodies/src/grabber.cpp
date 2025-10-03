#include <vector>
#include <wiPrimitive.h>
#include "grabber.h"
#include "softBody.h"
#include "simulation_utils.h"

// void Grabber::start(wi::primitive::Ray ray, wi::scene::PickResult pick, SoftBodyInstance* sbInstance)
void Grabber::start(wi::scene::PickResult pick, SoftBody* sbInstance)
{
    if (isGrabbing || sbInstance == nullptr)
        return;

    pickInfo = pick;
    this->sbInstance = sbInstance;
	// this->softBody = softBody;

	// XMVECTOR origin = XMLoadFloat3(&ray.origin);
	// XMVECTOR direction = XMLoadFloat3(&ray.direction);
	// XMVECTOR posVec = XMVectorAdd(origin, XMVectorScale(direction, pick.distance));
	// XMFLOAT3 posF;
	// XMStoreFloat3(&posF, posVec);
	// std::vector<float> pos = {posF.x, posF.y, posF.z};
	//
	//    physicsObject->StartGrab(pos);
	

	// auto objtrans = wi::scene::GetScene().transforms.GetComponent(sbInstance->objectEntity);
	// XMMATRIX objMatrix = objtrans->GetWorldMatrix();
	// XMMATRIX invObjMatrix = XMMatrixInverse(nullptr, objMatrix);
	//
	// XMVECTOR pickPosWorld = XMLoadFloat3(&pick.position);
	// XMVECTOR pickPosObject = XMVector3Transform(pickPosWorld, invObjMatrix);
	//
	// XMFLOAT3 pickPosObj;
	// XMStoreFloat3(&pickPosObj, pickPosObject);
	//
	// std::vector<float> pickPosObjVec = { pickPosObj.x, pickPosObj.y, pickPosObj.z };
	// sbInstance->softBody->StartGrab(pickPosObjVec);
	
	sbInstance->StartGrab(pick);

	std::vector<float> pos {pick.position.x, pick.position.y, pick.position.z};
    prevPos = pos;
	// prevPos = {pickPosObjVec};
    vel = {0, 0, 0};
    time = 0.0f;
	startDistance = pick.distance;
	// XMVECTOR originVec = XMLoadFloat3(&ray.origin);
	// XMVECTOR pickPosVec = XMVectorSet(
	// 	pickPosObjVec[0],
	// 	pickPosObjVec[1],
	// 	pickPosObjVec[2],
	// 	0.0f
	// );
	// XMVECTOR diff = XMVectorSubtract(originVec, pickPosVec);
	// XMVECTOR distVec = XMVector3Length(diff);
	// startDistance = XMVectorGetX(distVec);
    isGrabbing = true;

	if (gPhysicsScene.paused)
		gPhysicsScene.Run();
}

void Grabber::move(wi::primitive::Ray ray, wi::scene::PickResult pick)
// void Grabber::move(wi::scene::PickResult pick)
{
    if (!isGrabbing || sbInstance == nullptr)
        return;

	XMVECTOR origin = XMLoadFloat3(&ray.origin);
	XMVECTOR direction = XMLoadFloat3(&ray.direction);
	XMVECTOR posVec = XMVectorAdd(origin, XMVectorScale(direction, startDistance));
	XMFLOAT3 posF;
	XMStoreFloat3(&posF, posVec);
	std::vector<float> pos = {posF.x, posF.y, posF.z};
	//
	// vel = {pos};
	//
	// XMFLOAT3 velF = XMFLOAT3(vel[0], vel[1], vel[2]);
	// XMFLOAT3 prevPosF = XMFLOAT3(prevPos[0], prevPos[1], prevPos[2]);
	// XMVECTOR velVec = XMLoadFloat3(&velF);
	// XMVECTOR prevPosVec = XMLoadFloat3(&prevPosF);
	// velVec = XMVectorSubtract(velVec, prevPosVec);
	//
	// if (time > 0)
	//    {
	// 	velVec = XMVectorScale(velVec, 1.0f / time);
	//    }
	// else
	//        vel = {0, 0, 0};

	// std::vector<float> pos {pick.position.x, pick.position.y, pick.position.z};
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
	// sbInstance->softBody->MoveGrabbed(pos, vel);
    sbInstance->MoveGrabbed(pos, vel);
}

void Grabber::end()
{
    if (isGrabbing && sbInstance != nullptr)
    {
        std::vector<float> posVec = {prevPos};
        std::vector<float> velVec = {vel};
        sbInstance->EndGrab(posVec, velVec);
        sbInstance = nullptr;
        isGrabbing = false;
    }
}

void Grabber::increaseTime(float dt) { time += dt; }
