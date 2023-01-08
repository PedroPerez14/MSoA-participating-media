#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathTracing : public Integrator
{
public:
    PathTracing(const PropertyList &props){ }

    ///Preprocess() method can be implemented, useful for importance sampling lights??

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {   
        Ray3f path_ray = ray;
        Color3f throughput(1.0f, 1.0f, 1.0f);
        Intersection its;
        bool intersected = scene->rayIntersect(path_ray, its);
        int n_bounces = 0;
        EmitterQueryRecord emitterRecord(its.p);
        BSDFQueryRecord materialRecord(its.toLocal(-path_ray.d), its.uv);
        float prob_continuepath = 0.95f;
        
        while(intersected && throughput.getLuminance() > 0.f)
        {
            emitterRecord = EmitterQueryRecord(its.p);
            materialRecord = BSDFQueryRecord(its.toLocal(-path_ray.d), its.uv);
            
            if(its.mesh->isEmitter())
            {
                emitterRecord = EmitterQueryRecord(its.mesh->getEmitter(), path_ray.o, its.p, its.shFrame.n, its.uv);  //TODO CONSTRUIR BIEN JODER COJONES MIERDA HOSTIA
                return throughput * its.mesh->getEmitter()->eval(emitterRecord);
            }
            else
            {   
                Color3f iter_weight = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());
                prob_continuepath = std::min(throughput.maxCoeff() * pow(materialRecord.eta, 2.f), 0.99f);
                if(n_bounces > 3 && sampler->next1D() > prob_continuepath/*russianRoulette(sampler->next1D(), iter_weight, prob_continuepath)*/)
                {
                    return Color3f(0.f);
                }
                //else
                if(n_bounces > 3)
                {
                    throughput *= iter_weight / prob_continuepath;
                }
                else
                {
                    throughput *= iter_weight;
                }
                n_bounces++; 
            }

            path_ray.o = its.p;
            path_ray.d = its.toWorld(materialRecord.wo);
            intersected = scene->rayIntersect(Ray3f(its.p, its.toWorld(materialRecord.wo)), its);
        }
        return scene->getBackground(path_ray) * throughput;
    }

    std::string toString() const
    {
        return "Naive Path Tracing Sampling Integrator[]";
    }

private:
    Color3f getCameraDirectIllumination(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Color3f(0., 0., 0.);
    }

    bool russianRoulette(const float& sample, const Color3f& iter_weight, float &prob_alive) const
    {
        //if this condition is fulfilled, kill the path.
        //This happens when it's carrying too little energy (random number > energy carried)
        float prob_kill = std::max(std::max(iter_weight.x(), iter_weight.y()), iter_weight.z());
        prob_alive = 1.f - prob_kill;
        return (sample < prob_kill);
    }
};

NORI_REGISTER_CLASS(PathTracing, "path");
NORI_NAMESPACE_END