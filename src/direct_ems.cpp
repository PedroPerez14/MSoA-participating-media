#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectEmitterSampling : public Integrator
{
public:
    DirectEmitterSampling(const PropertyList &props){ }

    ///Preprocess() method can be implemented, useful for importance sampling lights??

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        Color3f Lo(0.);

        Intersection its;
        if(!scene->rayIntersect(ray, its))
        {
            return scene->getBackground(ray);
        }

        float pdflight;
        if(its.mesh->isEmitter())
        {
            EmitterQueryRecord emitterRecord(ray.o);
            its.mesh->getEmitter()->sample(emitterRecord, sampler->next2D(), 0.);   //We only call this to get the normals right
            Lo += its.mesh->getEmitter()->eval(emitterRecord);    //It's the same as the line above, but the other option looks better
        }

        EmitterQueryRecord emitterRecord(its.p);

        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.);

        float V(1.0f);   //Visibility term
        Ray3f sray(its.p, emitterRecord.wi);
        Intersection it_shadow;
        if(scene->rayIntersect(sray, it_shadow))
        {
            if(it_shadow.t > (emitterRecord.dist - Epsilon))
            {
                V = 1.0f;
            }
            else
            {
                V = 0.0f;
            }
        }
        else
        {
            V = 1.0f;
        }
        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d),
                                its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
        Lo += Le * V * its.shFrame.n.dot(emitterRecord.wi) * 
                                its.mesh->getBSDF()->eval(bsdfRecord) / pdflight;
        return Lo;
    }

    std::string toString() const
    {
        return "Direct Emitter Sampling Integrator[]";
    }

private:
    Color3f getCameraDirectIllumination(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Color3f(0., 0., 0.);
    }
};

NORI_REGISTER_CLASS(DirectEmitterSampling, "direct_ems");
NORI_NAMESPACE_END