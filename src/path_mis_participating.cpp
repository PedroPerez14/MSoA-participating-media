#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <exception>

NORI_NAMESPACE_BEGIN

class PathTracingMISParticipatingMedia : public Integrator
{
public:
    PathTracingMISParticipatingMedia(const PropertyList &props){ }
    
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        /// TODO: We assume that the medium the camera is inside is the global one
        std::shared_ptr<Volume> currentVolumeMedium = scene->getEnviromentalVolumeMedium();
        return LiRec(scene, sampler, ray, currentVolumeMedium, 0);
    }

    Color3f LiRec(const Scene* scene, Sampler* sampler, const Ray3f& ray, std::shared_ptr<Volume> currentVolumeMedium, int depth) const
    {
        std::cout << "---------------------------------------------------" << std::endl;
        // if(depth > 0)
        //     return Color3f(1.f);
        /// TODO: ver qué valor darle a max_depth y no ponerlo feo con un const int
        //const int max_depth = 128;
        Intersection its;
        //Color3f indirect(1.f);
        //Color3f direct(0.f);
        //Color3f contrib_indirect(1.f);
        //int n_bounces = 0;
        //float prob_alive = 0.95f;
        //BSDFQueryRecord materialRecord(Vector3f(0.f), Vector2f(0.f));
        //float pdir_wdir, pdir_wbsdf, pbsdf_wbsdf, pbsdf_wdir;
        //float w_mis_indir, w_mis_dir;

        
        //std::vector<unsigned int> volumes_through = std::vector<unsigned int>();
        //volumes_through.reserve(scene->getVolumes().size());    //We shouldn't be able to traverse a volume twice at the same time

        //Check if the camera is inside any volume at start of execution
        //isCameraInsideAnyVolume(scene, volumes_through);
        //getCurrentVolume(volumes_through, scene->getVolumes(), scene);


        //Color3f Lems(0.f), Lmat(0.f);
        Color3f L(0.f);
        float pdf_pathstep, pdf_succ, pdf_fail;
        
        //Boundary of the medium

        Point3f xt(0.f);
        xt = currentVolumeMedium->samplePathStep(ray.o, Point3f(10000.f), ray.d, sampler->next2D(), pdf_succ, pdf_fail);
        std::vector<VolumetricSegmentRecord> volumeSegments;
        //std::cout << "PRE INTERSECT" << std::endl;
        //std::cout << "SAMPLED XT: " << xt.x() << " " << xt.y() << " " << xt.z() << std::endl;
        bool intersected = scene->rayIntersectThroughVolumes(sampler, ray, xt, its, currentVolumeMedium, volumeSegments);
        //std::cout << "POST INTERSECT 1" << std::endl;
        //xz = its.p
        
        //Inscattering
        //Point3f xt(0.f, 0.f, 0.f);
        
        
        float z = abs(Vector3f(its.p - ray.o).norm());
        float t = abs(Vector3f(xt - ray.o).norm());
        
        if(intersected && t >= z)
        {
            //We hit the surface!
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord er(its.mesh->getEmitter(), ray.o, its.p, its.shFrame.n, its.uv);
                L = volumetric_transmittance(volumeSegments) * its.mesh->getEmitter()->eval(er);
                pdf_pathstep = currentVolumeMedium->pdfFail(its.p, z, sampler->next2D());//;pdf_fail;               // fail: 1 - cdf()      ///////////// xt??
            }
            else
            {
                std::shared_ptr<Volume> xz_volume = volumeSegments[volumeSegments.size() - 1].segment_vol;
                L = volumetric_transmittance(volumeSegments) * directLight(scene, sampler, currentVolumeMedium, xz_volume, its, ray.d, depth);    //TODO
                pdf_pathstep = currentVolumeMedium->pdfFail(its.p, z, sampler->next2D());//;pdf_fail;               // fail: 1 - cdf()       ///////////// xt???
            }
        }
        else
        {
            std::shared_ptr<Volume> xt_volume = volumeSegments[volumeSegments.size() - 1].segment_vol;
            if(!intersected)
            {
                //  Yeet into oblivion
                L = volumetric_transmittance(volumeSegments) * scene->getBackground(ray);
                pdf_pathstep = 1.f;
            }
            else   // t<z, 
            {
                //  Medium interaction
                L = volumetric_transmittance(volumeSegments) * inscattering(scene, sampler, currentVolumeMedium, xt_volume, xt, ray.d, depth);
                pdf_pathstep = (volumeSegments[0].segment_vol->transmittance(volumeSegments[0].x0, volumeSegments[0].xs)).sum() / 3.f;//pdf_succ;
            }
            //pdf_pathstep = pdf_succ;
        }
        
        return (L / pdf_pathstep);
    }
    
    std::string toString() const
    {
        return "Path Tracing Sampling Integrator for Single Scattering with Multiple Importance Sampling[]";
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

    bool emitterShRayIntersectFree(const Scene *scene, Sampler* sampler, const Ray3f &sray, const EmitterQueryRecord& emitterRecord, std::shared_ptr<Volume> currVolMedium, std::vector<VolumetricSegmentRecord>& segs_shadow) const
    {
        float t;
        if(!scene->shadowRayThroughVolumes(sampler, sray, currVolMedium, segs_shadow, t))
        {
            //std::cout << "emitterShRayIntersectFree NO INTERSECTA" << std::endl;
            return true;
        }
        else
        {
           //std::cout << "emitterShRayIntersectFree SI INTERSECTA" << std::endl;
            return (t > (emitterRecord.dist - Epsilon));            
        }
    }

    Color3f emitterSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w, float& w_mis_dir) const
    {
        Color3f Lems(Epsilon);
        float pdflight(1.f);
        float pdir_wdir(1.f), pbsdf_wdir(0.f);

        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(its.p);
        BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
        its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

        Intersection shray_its = Intersection();
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray(its.p, emitterRecord.wi);
        //std::cout << "EmitterSampling: Pre-ShadowRay" << std::endl;
        std::vector<VolumetricSegmentRecord> shadow_vsr;
        if(emitterShRayIntersectFree(scene, sampler, shadowray, emitterRecord, currentVolumeMedium, shadow_vsr) && materialRecord.measure != EDiscrete)
        {
            //std::cout << "EmitterSampling: IN-ShadowRay" << std::endl;
            BSDFQueryRecord bsdfRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * volumetric_transmittance(shadow_vsr) * its.mesh->getBSDF()->eval(bsdfRecord) * abs(its.shFrame.n.dot(emitterRecord.wi)) / pdflight;

            //MIS for emitter sampling
            pdir_wdir = em->pdf(emitterRecord);
            pbsdf_wdir = its.mesh->getBSDF()->pdf(bsdfRecord);
            
            if(!its.mesh->getBSDF()->isDiffuse())
                pbsdf_wdir = emitterRecord.pdf; 
            w_mis_dir = balanceHeuristic(pdir_wdir, pbsdf_wdir);
        }
        //std::cout << "EmitterSampling: Post-ShadowRay" << std::endl;
        return Lems;
        
    }

    Color3f volumetric_transmittance(const std::vector<VolumetricSegmentRecord>& vsr) const
    {
        std::cout << "VOL_TRANS: [";
        Color3f throughput(1.f);

        for(size_t i = 0; i < vsr.size(); i++)
        {
            std::cout << "[" << vsr[i].segment_vol->isHeterogeneous() << " " << vsr[i].x0.toString() << " " << vsr[i].xs.toString() << "]";
        }
        std::cout << "]" << std::endl;
        
        for(size_t i = 0; i < vsr.size(); i++)
        {
            Color3f a = vsr[i].segment_vol->transmittance(vsr[i].x0, vsr[i].xs);
            float pdf = vsr[i].vol_pdf;                                               /////// QUITARRR
            throughput *= (a / pdf);
        }
        
        return throughput;
    }


    /// Emitter Sampling, but the weights for MIS are calculated using phase function instead of BSDF
    Color3f emitterSamplingPF(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Point3f& xt, const Vector3f& w, float& w_mis_dir) const
    {
        Color3f Lems(Epsilon);
        float pdflight(1.f);
        float pdir_wdir(1.f), ppf_wdir(0.f);
        
        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(xt);
        Intersection shray_its;
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray = Ray3f(xt, emitterRecord.wi);

        std::vector<VolumetricSegmentRecord> shadow_vsr;
        if(emitterShRayIntersectFree(scene, sampler, shadowray, emitterRecord, currentVolumeMedium, shadow_vsr))
        {
            PFQueryRecord pfRecord(-w, emitterRecord.wi);
            Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * volumetric_transmittance(shadow_vsr) * currentVolumeMedium->getPhaseFunction()->eval(pfRecord) * mu_s / pdflight;

            //MIS for emitter sampling
            pdir_wdir = em->pdf(emitterRecord);
            ppf_wdir = currentVolumeMedium->getPhaseFunction()->pdf(pfRecord);
            
            w_mis_dir = balanceHeuristic(pdir_wdir, ppf_wdir);
                        
        }
        return Lems;
    }


    Color3f brdfSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w, float& w_mis_indir, Vector3f& wo_brdf, Color3f& fs) const
    {
        Color3f Lmat(Epsilon);
        //Second parameter, its wo, is random, since it will get properly calculated with sample()
        BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-w), its.toLocal(w), its.uv, ESolidAngle);
        fs = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());
        Intersection itsaux;
        Ray3f ray(its.p, its.toWorld(materialRecord.wo));

        std::vector<VolumetricSegmentRecord> vsr_shadow;
        float inf = std::numeric_limits<float>::infinity();
        bool intersects_aux = scene->rayIntersectThroughVolumes(sampler, ray, itsaux, currentVolumeMedium, vsr_shadow);
        // std::cout << "POST INTERSECT 2" << std::endl;
        if((intersects_aux && itsaux.mesh->isEmitter()) || (!intersects_aux && scene->getEnvironmentalEmitter()))
        {
            float pbsdf_wbsdf = its.mesh->getBSDF()->pdf(materialRecord);
            float pdir_wbsdf(0.f);
            Color3f emit, transmittance;
            if(intersects_aux && itsaux.mesh->isEmitter())
            {
                // Normal light emitters case
                EmitterQueryRecord emitterRecordAux(itsaux.mesh->getEmitter(), its.p, itsaux.p, itsaux.shFrame.n, itsaux.uv);
                pdir_wbsdf = (pow(emitterRecordAux.dist, 2.0f) / abs(emitterRecordAux.n.dot(emitterRecordAux.wi)));//emitterRecordAux.emitter->pdf(emitterRecordAux) * scene->pdfEmitter(emitterRecordAux.emitter);
                emit = emitterRecordAux.emitter->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(vsr_shadow);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(its.p);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wbsdf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(vsr_shadow);
            }

            if(materialRecord.measure == EDiscrete)
            {
                pbsdf_wbsdf = 1.f;
            }
            w_mis_indir = balanceHeuristic(pbsdf_wbsdf, pdir_wbsdf);
            
            Lmat = fs * emit * transmittance;
        }
        wo_brdf = its.toWorld(materialRecord.wo);
        //std::cout << "WO_BRDF: " << wo_brdf.x() << " " << wo_brdf.y() << " " << wo_brdf.z() << std::endl;
        return Lmat;
    }

    Color3f phaseFunctionSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Point3f& xt, const Vector3f& w, float& w_mis_indir, Vector3f& wo_pf, float& fs) const
    {
        Color3f Lpf(Epsilon);

        PFQueryRecord pfRecord(-w);
        fs = currentVolumeMedium->getPhaseFunction()->sample(pfRecord, sampler->next2D());
        wo_pf = pfRecord.wo;
        Intersection itsaux;

        Ray3f ray(xt, pfRecord.wo);
        std::vector<VolumetricSegmentRecord> vsr_shadow;
        float inf = std::numeric_limits<float>::infinity();
        bool intersects_aux = scene->rayIntersectThroughVolumes(sampler, ray, itsaux, currentVolumeMedium, vsr_shadow);
        // std::cout << "POST INTERSECT 3" << std::endl;
        if((intersects_aux && itsaux.mesh->isEmitter()) || (!intersects_aux && scene->getEnvironmentalEmitter()))
        {
            float ppf_wpf = pfRecord.m_pdf;
            float pdir_wpf(0.f);
            Color3f emit, transmittance;
            if(intersects_aux && itsaux.mesh->isEmitter())
            {
                // Normal light emitters case
                EmitterQueryRecord emitterRecordAux(itsaux.mesh->getEmitter(), xt, itsaux.p, itsaux.geoFrame.n, itsaux.uv);
                emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                emitterRecordAux = EmitterQueryRecord(itsaux.mesh->getEmitter(), xt, itsaux.p, itsaux.geoFrame.n, itsaux.uv);
                pdir_wpf = (pow(emitterRecordAux.dist, 2.0f) / abs(emitterRecordAux.n.dot(emitterRecordAux.wi)));//emitterRecordAux.emitter->pdf(emitterRecordAux) * scene->pdfEmitter(emitterRecordAux.emitter);
                emit = emitterRecordAux.emitter->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(vsr_shadow);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(xt);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wpf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(vsr_shadow);
            }

            w_mis_indir = balanceHeuristic(ppf_wpf, pdir_wpf);

            Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lpf = fs * emit * transmittance * mu_s;
        }
        return Lpf;
    }


    Color3f directLight(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, std::shared_ptr<Volume> xzVolumeMedium, const Intersection& its, const Vector3f& w, int n_bounces) const
    {
        Color3f Lems(0.f), Lmat(0.f), fs_bsdf(0.f);
        float w_mis_dir(0.f), w_mis_indir(0.f);
        Vector3f wo_brdf;

        Lems = emitterSampling(scene, sampler, xzVolumeMedium, its, w, w_mis_dir);
        if(!isnan(w_mis_dir))
                Lems *= w_mis_dir;
        //std::cout << "directLight postEmitterSampling" << std::endl;
        Lmat = brdfSampling(scene, sampler, xzVolumeMedium, its, w, w_mis_indir, wo_brdf, fs_bsdf);
        //std::cout << "directLight postBRDFSampling" << std::endl;
        if(!isnan(w_mis_indir))
                Lmat *= w_mis_indir;

        Color3f Lnee = Lems + Lmat;

        float prob_alive = 0.95f;
        if(n_bounces > 2 && sampler->next1D() > prob_alive)
        {
            return Lnee;
        }    

        Ray3f newRay(its.p, wo_brdf);
        if(n_bounces > 2)
            return Lnee * (fs_bsdf * LiRec(scene, sampler, Ray3f(its.p, wo_brdf), xzVolumeMedium, n_bounces + 1) / (prob_alive));
        else
            return Lnee * (fs_bsdf * LiRec(scene, sampler, Ray3f(its.p, wo_brdf), xzVolumeMedium, n_bounces + 1));
    }

    Color3f inscattering(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, std::shared_ptr<Volume> xtVolumeMedium, const Point3f& xt, const Vector3f& w, int n_bounces) const
    {
        
        Color3f Lems(0.f), Lpf(0.f);
        float w_mis_dir(0.f), w_mis_pf(0.f), fs_pf(Epsilon);
        Vector3f wo_pf;

        Lems = emitterSamplingPF(scene, sampler, xtVolumeMedium, xt, w, w_mis_dir);
        if(!isnan(w_mis_dir))
                Lems *= w_mis_dir;

        //std::cout << "inscattering postEmitterSamplingPF" << std::endl;
        Lpf = phaseFunctionSampling(scene, sampler, xtVolumeMedium, xt, w, w_mis_pf, wo_pf, fs_pf);
        if(!isnan(w_mis_pf))
                Lpf *= w_mis_pf;
        //std::cout << "inscattering postPhaseFunctionSampling" << std::endl;
        Color3f Lnee = Lems + Lpf;

        float prob_alive = 0.95f;
        if(n_bounces > 2 && sampler->next1D() > prob_alive)
        {
            return Lnee;
        }    

        /// HERE
        /// @brief 
        /// @param scene 
        /// @param sampler 
        /// @param currentVolumeMedium 
        /// @param xtVolumeMedium 
        /// @param xt 
        /// @param w 
        /// @param n_bounces 
        /// @return 
        fs_pf = 1.0f;
        if(n_bounces > 2)
            return Lnee * (fs_pf * LiRec(scene, sampler, Ray3f(xt, wo_pf), xtVolumeMedium, n_bounces + 1) / prob_alive);
        else
            return Lnee * (fs_pf * LiRec(scene, sampler, Ray3f(xt, wo_pf), xtVolumeMedium, n_bounces + 1));
    }
};

NORI_REGISTER_CLASS(PathTracingMISParticipatingMedia, "path_mis_participating_media");
NORI_NAMESPACE_END