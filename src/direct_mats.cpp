#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectMaterialSampling : public Integrator
{
public:
    DirectMaterialSampling(const PropertyList &props){ }

    ///Preprocess() method can be implemented, useful for importance sampling lights??

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        Color3f Lo(0.);

        Intersection its;
        if(!scene->rayIntersect(ray, its))
        {
            return scene->getBackground(ray);                   //TODO keep this? i guess yes
        }

        if(its.mesh->isEmitter())
        {
            //We only build an EmitterQueryRecord to get the normals right when evaluating
            EmitterQueryRecord emitterRecord(ray.o);
            its.mesh->getEmitter()->sample(emitterRecord, sampler->next2D(), 0.);   
            Lo += its.mesh->getEmitter()->eval(emitterRecord);      //return instead of Lo+= ???
        }

        //if intersects and is not emitter
        // we sample the material
        BSDFQueryRecord materialRecord(its.toLocal(-ray.d), its.uv);    //only wi and uv at intersection point are known
        Color3f fr_sample_weighted = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

        //now we cast another ray to try to reach a light source
        //3 cases are possible now:
        //      1) Don't intersect: Le = background illumination (i.e envmap)
        //      2) Intersect non-emitter: Le = 0, product with fr_sample_weighted will be ~0
        //      3) Intersect an emitter: Le = eval(emitter), multiply with fr_sample_weighted

        Color3f Le(Epsilon);
        Ray3f sray(its.p, its.toWorld(materialRecord.wo));
        Intersection it_sray;
        if(scene->rayIntersect(sray, it_sray))
        {
            if(it_sray.mesh->isEmitter())
            {
                //We only build an EmitterQueryRecord to get the normals right when evaluating
                EmitterQueryRecord emitterRecord2(it_sray.mesh->getEmitter(), its.p, it_sray.p, its.shFrame.n, it_sray.uv);       //ref, p, n, uv
                it_sray.mesh->getEmitter()->sample(emitterRecord2, sampler->next2D(), 0.);
                Le = it_sray.mesh->getEmitter()->eval(emitterRecord2);
            }
            //else{ Le ~= 0,0,0 }
        }
        else
        {
            Le = scene->getBackground(sray);
        }

        Lo += Le * fr_sample_weighted;

        return Lo;
    }

    std::string toString() const
    {
        return "Direct Material Sampling Integrator[]";
    }

private:
    Color3f getCameraDirectIllumination(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Color3f(0., 0., 0.);
    }
};

NORI_REGISTER_CLASS(DirectMaterialSampling, "direct_mats");
NORI_NAMESPACE_END