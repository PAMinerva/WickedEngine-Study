#pragma once

#include "wiRenderPath3D.h"

class MyRender3D : public wi::RenderPath3D
{
    //wi::ecs::Entity triangleMeshEntity;

public:
    MyRender3D()
    {
        wi::scene::Scene& scene = wi::scene::GetScene();

        auto materialEntity = scene.Entity_CreateMaterial("default");
        auto& material = *scene.materials.GetComponent(materialEntity);
        material.shaderType = wi::scene::MaterialComponent::SHADERTYPE_UNLIT;
        material.SetUseVertexColors(true);
        material.SetDoubleSided(true);

        auto triangleMeshEntity = scene.Entity_CreateMesh("triangle");

        auto& objectTriangle = scene.objects.Create(triangleMeshEntity);
        objectTriangle.meshID = triangleMeshEntity;

        auto& meshTriangle = *scene.meshes.GetComponent(triangleMeshEntity);
        meshTriangle.subsets.push_back(wi::scene::MeshComponent::MeshSubset());
        meshTriangle.subsets.back().materialID = materialEntity;
        meshTriangle.indices.resize(3);
        meshTriangle.subsets.back().indexOffset = (uint32_t)0;
        meshTriangle.subsets.back().indexCount = (uint32_t)3;
        meshTriangle.indices[0] = 0;
        meshTriangle.indices[1] = 2;
        meshTriangle.indices[2] = 1;
        meshTriangle.vertex_positions.resize(3);
        meshTriangle.vertex_positions[0] = XMFLOAT3(-1.0f, -0.5f, 2.0f);
        meshTriangle.vertex_positions[1] = XMFLOAT3(1.0f, -0.5f, 2.0f);
        meshTriangle.vertex_positions[2] = XMFLOAT3(0.0f, 1.0f, 2.0f);
        meshTriangle.vertex_colors.resize(3);

        const XMFLOAT4 red = XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f);
        const XMFLOAT4 green = XMFLOAT4(0.0f, 1.0f, 0.0f, 1.0f);
        const XMFLOAT4 blue = XMFLOAT4(0.0f, 0.0f, 1.0f, 1.0f);
        uint32_t rgbaRed = wi::math::CompressColor(red);
        uint32_t rgbaGreen = wi::math::CompressColor(green);
        uint32_t rgbaBlue = wi::math::CompressColor(blue);
        meshTriangle.vertex_colors[0] = rgbaRed;
        meshTriangle.vertex_colors[1] = rgbaGreen;
        meshTriangle.vertex_colors[2] = rgbaBlue;

        meshTriangle.CreateRenderData();

        scene.transforms.Create(triangleMeshEntity);
        auto meshtrans = scene.transforms.GetComponent(triangleMeshEntity);
        meshtrans->Translate(XMFLOAT3(0, -0.2f, 0));
        meshtrans->Scale(XMFLOAT3(1, 1, 1));
    }

    /*void Update(float dt) override
    {
        using namespace wi::scene;
        static float speed = 1.5f;
        static float direction = 1.0f;
        static float delta_x = 0.0f;

        TransformComponent* transform = GetScene().transforms.GetComponent(triangleMeshEntity);

        if (transform != nullptr)
        {
            float pos_x = transform->GetPosition().x;
            delta_x = speed * direction * dt;

            if (pos_x > 1.25f)
            {
                direction = -1.0f;
            }
            else if (pos_x < -1.25f)
            {
                direction = 1.0f;
            }

            transform->Translate(XMFLOAT3(delta_x, 0, 0));
        }

        RenderPath3D::Update(dt);
    }*/
};

 class Sample
{
public:
    MyRender3D render3D;
    const char *title {"01 - Hello Triangle"};

    wi::RenderPath* GetRenderPath3D() {return &render3D;};
    const char* GetTitle() {return title;};
};