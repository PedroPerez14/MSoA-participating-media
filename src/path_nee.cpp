#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathTracingNEE : public Integrator
{
public:
    PathTracingNEE(const PropertyList &props){ }
    
    //Failed attempt of making it iterative, now working
    Color3f Li_iter(const Scene* scene, Sampler* sampler, Ray3f ray) const
    {
        Intersection its;
        Color3f indirect(1.f);
        Color3f direct(0.f);
        Color3f contrib_indirect(1.f);
        int n_bounces = 0;
        float prob_alive = 0.95f;
        BSDFQueryRecord materialRecord(Vector3f(0.f), Vector2f(0.f));

        while(scene->rayIntersect(ray, its) && indirect.getLuminance() > 0.f)
        {
            if(its.mesh->isEmitter())
            {
                if(n_bounces == 0 || materialRecord.measure == EDiscrete)
                {
                    EmitterQueryRecord emitterRecord(its.mesh->getEmitter(), 
                                        ray.o, its.p, its.shFrame.n, its.uv);
                    return direct + indirect * its.mesh->getEmitter()->eval(emitterRecord);
                }//else
                return direct;
            }
            else
            {
                float pdflight(1.f);
                EmitterQueryRecord emitterRecord(its.p);
                const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);                    //Uniform light sampling
                //const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight, emitterRecord);      //Importance light sampling
                Intersection shray_its = Intersection();

                materialRecord = BSDFQueryRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                Color3f light_spectrum = em->sample(emitterRecord, sampler->next2D(), 0.f);
                contrib_indirect = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());
                Ray3f shadowray = Ray3f(its.p, emitterRecord.wi);

                if(emitterShRayIntersectFree(scene, shadowray, emitterRecord, shray_its) && materialRecord.measure != EDiscrete)
                {
                    BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                    Color3f contrib_direct = indirect * light_spectrum * its.mesh->getBSDF()->eval(bsdfRecord) * abs(its.shFrame.n.dot(emitterRecord.wi)) / pdflight;
                    direct += contrib_direct;
                }

                prob_alive = std::min(indirect.maxCoeff() * pow(materialRecord.eta, 2.f), 0.99f);
                if(n_bounces > 3 && sampler->next1D() > prob_alive)
                {
                    return direct;
                }

                if(n_bounces > 3)
                    indirect *= (contrib_indirect / prob_alive);
                else
                    indirect *= (contrib_indirect);
                
                n_bounces++;
                ray = Ray3f(its.p, its.toWorld(materialRecord.wo));
                its = Intersection();
            }
        }
        return direct + indirect * scene->getBackground(ray);
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Li_iter(scene, sampler, ray);
        //BSDFQueryRecord bsdf_rec(Vector3f(0.f), Vector2f(0.f));
        //Color3f indirect(1.f);
        //return LiRec(scene, sampler, ray, indirect, 0, bsdf_rec);
    }
    
    //Recursive version of our NEE path tracer. Works perfectly, as well as the iterative version
    Color3f LiRec(const Scene* scene, Sampler* sampler, const Ray3f& ray, Color3f indirect, int n_bounces, BSDFQueryRecord& materialRecord) const
    {
        Intersection its, its_sray;
        float pdflight = 1.f;
        float prob_alive_rr = 0.95f;            // Recalculated at each bounce according to https://wjakob.github.io/nori/ (Assignment 5 section)

        if(scene->rayIntersect(ray, its) && indirect.getLuminance() > 0.f)
        {
            if(its.mesh->isEmitter())
            {
                if(n_bounces == 0 || materialRecord.measure == EDiscrete)
                {
                    EmitterQueryRecord emitterRecord(its.mesh->getEmitter(), ray.o, its.p, its.shFrame.n, its.uv);  
                    return its.mesh->getEmitter()->eval(emitterRecord);
                }
                else
                {
                    return Color3f(0.f);
                }
            }
            //else

            Color3f direct(0.f);
            EmitterQueryRecord emitterRecord(its.p);
            const Emitter* em = scene->sampleEmitter(sampler, pdflight, emitterRecord);             //Importance emitter sampling
            Color3f light_spectrum = em->sample(emitterRecord, sampler->next2D(), 0.);
            Ray3f sray(its.p, emitterRecord.wi);

            //Indirect light
            materialRecord = BSDFQueryRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
            Color3f bsdf_spectrum = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

            prob_alive_rr = std::min(indirect.maxCoeff() * pow(materialRecord.eta, 2.f), 0.99f);
            //rr
            if(n_bounces > 3 && sampler->next1D() > prob_alive_rr)
            {
                return Color3f(0.f);
            }//else, rr kills the path. Check number of bounces to apply prob. correctly


            indirect *= LiRec(scene, sampler, Ray3f(its.p, its.toWorld(materialRecord.wo)), indirect, ++n_bounces, materialRecord) * bsdf_spectrum;

            //Direct light
            if(emitterShRayIntersectFree(scene, sray, emitterRecord, its_sray) && materialRecord.measure != EDiscrete)
            {   
                BSDFQueryRecord bsdfRecord = BSDFQueryRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                direct = (light_spectrum * abs(its.shFrame.n.dot(emitterRecord.wi)) * its.mesh->getBSDF()->eval(bsdfRecord) / pdflight);
            }

            if(n_bounces > 3)
            {
                return (direct + indirect) / prob_alive_rr;
            }
            //else
            return (direct + indirect);
        }
        else
        {
            return scene->getBackground(ray);
        }
    }
    
    std::string toString() const
    {
        return "Next Event Estimation Path Tracing Sampling Integrator[]";
    }

private:
    Color3f getCameraDirectIllumination(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Color3f(0., 0., 0.);
    }

    bool emitterShRayIntersectFree(const Scene *scene, Ray3f sray, const EmitterQueryRecord& emitterRecord, Intersection it_shadow) const
    {
        if(!scene->rayIntersect(sray, it_shadow))
        {
            return true;
        }
        else
        {
            return (it_shadow.t > (emitterRecord.dist - Epsilon));
        }
    }
};

NORI_REGISTER_CLASS(PathTracingNEE, "path_nee");
NORI_NAMESPACE_END