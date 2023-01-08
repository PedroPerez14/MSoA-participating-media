#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectMIS : public Integrator
{
public:
    DirectMIS(const PropertyList &props){ }

    ///Preprocess() method can be implemented, useful for importance sampling lights??

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        Color3f Lo(Epsilon);

        Intersection its;
        if(!scene->rayIntersect(ray, its))
        {
            return scene->getBackground(ray);
        }

        if(its.mesh->isEmitter())
        {
            //We only build an EmitterQueryRecord to get the normals right when evaluating
            EmitterQueryRecord emitterRecord(ray.o);
            its.mesh->getEmitter()->sample(emitterRecord, sampler->next2D(), 0.);   
            Lo += its.mesh->getEmitter()->eval(emitterRecord);          //return instead of Lo+= ???
        }

        float pem_wem, pmat_wmat, pem_wmat, pmat_wem;
        EmitterQueryRecord emitterRecord(its.p);                        //sample = emitterRecord.wi;
        BSDFQueryRecord materialRecord(its.toLocal(-ray.d), its.uv);    //only wi and uv at intersection point are known              

        Color3f mis_emitterSampling = EmitterSampling(scene, sampler, ray, its, pem_wem, emitterRecord);
        float lightpdf(0.f);
        Color3f mis_bsdfSampling = BSDFSampling(scene, sampler, ray, its, pmat_wmat, materialRecord, lightpdf);

        Vector3f wh = (materialRecord.wi + materialRecord.wo).normalized();
        EmitterQueryRecord eqr_aux(emitterRecord.emitter, its.p, its.p + its.toWorld(materialRecord.wo) , its.toWorld(wh), materialRecord.uv);  //huele a mierda
        BSDFQueryRecord mqr_aux(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);


        pmat_wem = its.mesh->getBSDF()->pdf(mqr_aux);   //puta mierda casi seguro
        pem_wmat = scene->pdfEmitter(eqr_aux.emitter) * lightpdf;                           //LE FALTA(ba) PESO A MATS       
        //pem_wmat = eqr_aux.emitter->pdf(eqr_aux) / scene->pdfEmitter(eqr_aux.emitter); //puta mierda seguro 100% //SI ECUENTRO LUZ, 1/N * PROB DE SAMPLEAR EL PUNTO QUE HE PILLADO DE LA LUZ
        //pem_wmat = (emitterRecord.pdf * (pow((eqr_aux.p - eqr_aux.ref).norm(), 2.0f) / eqr_aux.n.dot(eqr_aux.wi))) / scene->pdfEmitter(eqr_aux.emitter); //puta mierda seguro 100% //SI ECUENTRO LUZ, 1/N * PROB DE SAMPLEAR EL PUNTO QUE HE PILLADO DE LA LUZ
        
        
        float w_emitterSampling(0.);
        if(pem_wem != 0.)
        {
           w_emitterSampling = pem_wem / (pem_wem + pmat_wem);
        } 
        
        float w_bsdfSampling(0.);
        if(pmat_wmat != 0.)
        {
            w_bsdfSampling = pmat_wmat / (pem_wmat + pmat_wmat);
        }        
        return Lo + (w_emitterSampling * mis_emitterSampling) + (w_bsdfSampling * mis_bsdfSampling);
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

    Color3f EmitterSampling(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& its, float& _pdflight, EmitterQueryRecord& emitterRecord) const
    {
        float pdflight(0.);
        //const Emitter* em = scene->sampleEmitter(sampler, pdflight, EmitterQueryRecord(its.p));   //Importance emitter sampling   //TODO can use it if i change the 1/n for the correct term
        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);                      //Uniform light sampling (use here for now)
        emitterRecord.emitter = em;
        Color3f Li = em->sample(emitterRecord, sampler->next2D(), 0.);

        _pdflight = emitterRecord.pdf;
        
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
        //else V = (1,1,1)

        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d),
                                its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
        return Li * V * its.shFrame.n.dot(emitterRecord.wi) * its.mesh->getBSDF()->eval(bsdfRecord) / pdflight;
    }

    Color3f BSDFSampling(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& its, float& pdflight, BSDFQueryRecord& materialRecord, float& lightpdf) const
    {
        Color3f fr_sample_weighted = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());
        pdflight = its.mesh->getBSDF()->pdf(materialRecord);

        //now we cast another ray to try to reach a light source
        //3 cases are possible now:
        //      1) Don't intersect: Le = background illumination (i.e envmap)
        //      2) Intersect non-emitter: Le = 0, product with fr_sample_weighted will be ~0
        //      3) Intersect an emitter: Le = eval(emitter), multiply with fr_sample_weighted

        Color3f Le(0.f);
        Ray3f sray(its.p, its.toWorld(materialRecord.wo));
        Intersection it_sray;
        lightpdf = 0.f;
        if(scene->rayIntersect(sray, it_sray))
        {
            if(it_sray.mesh->isEmitter())
            {
                //We only build an EmitterQueryRecord to get the normals right when evaluating
                EmitterQueryRecord emitterRecord2(it_sray.mesh->getEmitter(), its.p, it_sray.p, its.shFrame.n, it_sray.uv);       //ref, p, n, uv
                it_sray.mesh->getEmitter()->sample(emitterRecord2, sampler->next2D(), 0.);
                lightpdf = emitterRecord2.emitter->pdf(emitterRecord2);
                Le = it_sray.mesh->getEmitter()->eval(emitterRecord2);
            }
            //else{ Le ~= 0,0,0 }
        }
        else
        {
            Le = scene->getBackground(sray);
        }
        return Le * fr_sample_weighted;
    }
};

NORI_REGISTER_CLASS(DirectMIS, "direct_mis");
NORI_NAMESPACE_END