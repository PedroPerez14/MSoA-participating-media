#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathTracingMIS : public Integrator
{
public:
    PathTracingMIS(const PropertyList &props){ }
    
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Li_iter(scene, sampler, ray);
        //BSDFQueryRecord bsdf_rec(Vector3f(0.f), Vector2f(0.f));
        //return LiRec(scene, sampler, ray, Color3f(1.f), 0, bsdf_rec);
    }
    
    //Iterative version of the path tracer
    Color3f Li_iter(const Scene* scene, Sampler* sampler, Ray3f ray) const
    {
        Intersection its;
        Color3f indirect(1.f);
        Color3f direct(0.f);
        Color3f contrib_indirect(1.f);
        int n_bounces = 0;
        float prob_alive = 0.95f;
        BSDFQueryRecord materialRecord(Vector3f(0.f), Vector2f(0.f));
        float pdir_wdir, pdir_wbsdf, pbsdf_wbsdf, pbsdf_wdir;
        float w_mis_indir, w_mis_dir;

        while(scene->rayIntersect(ray, its) && indirect.getLuminance() > 0.f)
        {
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord emitterRecord(its.mesh->getEmitter(), 
                                ray.o, its.p, its.shFrame.n, its.uv);
                /*cout << "------------------------------------------------" << endl;
                cout << "N_BOUNCES: " << n_bounces << endl;
                cout << "DIRECT: " << direct.toString() << endl;
                cout << "INDIRECT: " << indirect.toString() << endl;
                cout << "INDIRECT * EMITEREVAL(): " << Color3f(indirect * its.mesh->getEmitter()->eval(emitterRecord)).toString() << endl;
                cout << "------------------------------------------------" << endl;*/
                
                return direct + indirect * its.mesh->getEmitter()->eval(emitterRecord);
            }
            else
            {
                float pdflight(0.f);
                const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
                EmitterQueryRecord emitterRecord(its.p);
                Intersection shray_its = Intersection();

                materialRecord = BSDFQueryRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                Color3f light_spectrum = em->sample(emitterRecord, sampler->next2D(), 0.f);
                contrib_indirect = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

                Ray3f shadowray = Ray3f(its.p, emitterRecord.wi);

                if(emitterShRayIntersectFree(scene, shadowray, emitterRecord, shray_its) && materialRecord.measure != EDiscrete)
                {
                    
                    BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                    //its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());
                    Color3f contrib_direct = indirect * light_spectrum * its.mesh->getBSDF()->eval(bsdfRecord) * abs(its.shFrame.n.dot(emitterRecord.wi)) / pdflight;

                    //MIS for direct light
                    pdir_wdir = em->pdf(emitterRecord);
                    pbsdf_wdir = its.mesh->getBSDF()->pdf(bsdfRecord);
                    
                    
                    if(!its.mesh->getBSDF()->isDiffuse())
                        pbsdf_wdir = pbsdf_wdir = emitterRecord.pdf; 
                    w_mis_dir = balanceHeuristic(pdir_wdir, pbsdf_wdir);
                    

                    //cout << "W_MIS_DIR " << w_mis_dir << " " << n_bounces << endl;
                    
                    if(!isnan(w_mis_dir))
                        contrib_direct *= w_mis_dir;
                    direct += contrib_direct;
                }

                prob_alive = std::min(indirect.maxCoeff() * pow(materialRecord.eta, 2.f), 0.99f);
                if(n_bounces > 3 && sampler->next1D() > prob_alive)
                {
                    return direct;
                }                

                Intersection itsaux;
                bool intersects_aux = scene->rayIntersect(Ray3f(its.p, its.toWorld(materialRecord.wo)), itsaux);
                if((intersects_aux && itsaux.mesh->isEmitter()))
                {
                    //MIS for indirect light
                    pbsdf_wbsdf = its.mesh->getBSDF()->pdf(materialRecord);
                    EmitterQueryRecord emitterRecordAux(itsaux.mesh->getEmitter(), 
                                    its.p, itsaux.p, its.shFrame.n, itsaux.uv);

                    //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                    pdir_wbsdf = (pow(emitterRecordAux.dist, 2.0f) / abs(emitterRecordAux.n.dot(emitterRecordAux.wi)));//emitterRecordAux.emitter->pdf(emitterRecordAux) * scene->pdfEmitter(emitterRecordAux.emitter);
                    
                    if(materialRecord.measure == EDiscrete)
                    {
                        pbsdf_wbsdf = 1.f;
                    }
                    w_mis_indir = balanceHeuristic(pbsdf_wbsdf, pdir_wbsdf);
                    

                    if(!isnan(w_mis_indir))
                        contrib_indirect *= w_mis_indir;
                    //cout << "W_MIS_INDIR " << w_mis_indir << " " << n_bounces <<endl;
                }
                if(!intersects_aux && scene->getEnvironmentalEmitter())
                {
                    //MIS for indirect light over the env emitter
                    pbsdf_wbsdf = its.mesh->getBSDF()->pdf(materialRecord);
                    
                    EmitterQueryRecord emitterRecordAux(its.p);
                    
                    emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                    scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                    pdir_wbsdf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                    
                    if(materialRecord.measure == EDiscrete)
                    {
                        pbsdf_wbsdf = 1.f;
                    }
                    w_mis_indir = balanceHeuristic(pbsdf_wbsdf, pdir_wbsdf);
                    

                    //cout << "W_MIS_INDIR " << w_mis_indir << " " << n_bounces << endl;
                    contrib_indirect *= w_mis_indir;
                }

                

                if(n_bounces > 3)
                    indirect *= (contrib_indirect / prob_alive);
                else
                    indirect *= (contrib_indirect);

                ++n_bounces;
                ray = Ray3f(its.p, its.toWorld(materialRecord.wo));
                its = Intersection();
            }
        }
        return direct + indirect * scene->getBackground(ray);
    }
    
    std::string toString() const
    {
        return "Path Tracing Sampling Integrator with Multiple Importance Sampling[]";
    }

private:
    Color3f getCameraDirectIllumination(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return Color3f(0., 0., 0.);
    }

    float balanceHeuristic(float px_wx, float py_wx) const
    {
        return (px_wx / (px_wx + py_wx));
    }

    float powerHeuristic(float px_wx, float py_wx, float beta) const
    {
        return (pow(px_wx, beta) / pow(px_wx + py_wx, beta));
    }

    bool emitterShRayIntersectFree(const Scene *scene, Ray3f &sray, const EmitterQueryRecord& emitterRecord, Intersection &it_shadow) const
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

NORI_REGISTER_CLASS(PathTracingMIS, "path_mis");
NORI_NAMESPACE_END