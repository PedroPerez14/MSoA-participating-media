#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/volume.h>
#include <nori/phasefunction.h>
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

    Color3f LiRec(const Scene* scene, Sampler* sampler, Ray3f ray, std::shared_ptr<Volume>& currentVolumeMedium, int depth) const
    {
        /// TODO: Change this for testing and faster rendering, I guess
        const int maxDepth = 9999999999;

        bool specularBounce = false;
        Color3f L(0.f);
        Color3f beta(1.f, 1.f, 1.f);     //Same notation as PBR Book, this would be Tr() / p(t) or Tr() / 1-cdf() depending on the interaction

        for(int bounces = 0; ; ++bounces)
        {
            //std::cout << "---------------------------------------------------" << std::endl;
            Intersection its;
            Point3f xt;
            bool foundIntersection = scene->rayIntersect(ray, its); 
            std::shared_ptr<Volume> nextVolumeMedium;
            bool sampledMedium = false;

            /// If we intersect anything and our current ray comes from any medium
            /// Sample the participating medium, if present
            if(foundIntersection && currentVolumeMedium)
            {
                // The intersection distance will be stored in its.t
                // So we sample an interaction
                // And later we will check if it's < its.t or >= its.t to get a medium interaction or geometry intersection
                Color3f _beta(1.f);
                xt = currentVolumeMedium->samplePathStep(ray, its, sampler, _beta, nextVolumeMedium, currentVolumeMedium, sampledMedium);
                beta *= _beta;
            }

            if(sqrt(beta.abs2().sum()) < Epsilon) // isBlack check
            {
                break;
            }

            /// If it's a medium interaction
            if(sampledMedium)
            {
                if(bounces >= maxDepth) break;              //check this
                L += beta * inscattering(scene, sampler, currentVolumeMedium, xt, ray.d);
                PFQueryRecord pfqr(-ray.d);
                currentVolumeMedium->getPhaseFunction()->sample(pfqr, sampler->next2D());
                ray = Ray3f(xt, pfqr.wo);
                specularBounce = false;
            }
            else
            {

                if(bounces == 0 || specularBounce)
                {
                    if(foundIntersection)
                    {
                        if(its.mesh->isEmitter())
                        {
                            EmitterQueryRecord er(its.mesh->getEmitter(), ray.o, its.p, its.shFrame.n, its.uv);
                            L += beta * its.mesh->getEmitter()->eval(er);
                            return L;           /// TODO: ojo! a lo mejor este return sobra!
                        }
                    }
                    else
                    {
                        L += beta * scene->getBackground(ray);
                        return L;               /// TODO: ojo! a lo mejor este return sobra!
                    }
                }

                /// Terminate path if no intersection was found or maxDepth was reached
                if(!foundIntersection || bounces >= maxDepth)
                    break;  /// TODO: Revisar, si se va a la nada no tengo que ponerle la luz de fondo?
                
                /// Skip over medium boundaries
                if(its.mesh->isVolume())
                {
                    ray = Ray3f(its.p, ray.d);  /// TODO: Sumo epsilon o no?
                    bounces--;
                    
                    if(!sampledMedium && currentVolumeMedium == nextVolumeMedium)
                    {
                        currentVolumeMedium = scene->getEnviromentalVolumeMedium();
                    }
                    else
                    {
                        currentVolumeMedium = nextVolumeMedium;
                    }
                    continue;
                }

                /// Sample illumination from lights to find attenuated path contribution
                L += beta * directLight(scene, sampler, currentVolumeMedium, its, ray.d);

                /// TODO: Creo que esto lo puedo hacer directamente en el UniformSampleOneLight (que mejor le dejo el nombre original o ké)
                /// Sample BSDF to get new path direction
                BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-ray.d), its.uv); //Esto me lo están sacando de los sample_direct, sample_material, sample_pf aquí para poder reusar el mismo método en todas partes
                Color3f fs = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

                /// TODO: revisar
                if(sqrt(fs.abs2().sum()) < Epsilon || isnan(fs.x()) || isnan(fs.y()) || isnan(fs.z()))     // If isBlack() or pdf == 0 (which produces NaNs in fs)
                {
                    break;
                }
                beta *= fs;     //TODO: acomodar luego a lo que me vaya pidiendo el codiguín
                specularBounce = !its.mesh->getBSDF()->isDiffuse();
                ray = Ray3f(its.p + ray.d * Epsilon, its.toWorld(materialRecord.wo));
            }

            /// Possibly terminate the path with Russian Roulette
            if(bounces > 3)
            {
                float rr = std::max(0.01f, 1 - beta.y());
                if(sampler->next1D() < rr)
                {
                    break;
                }
                beta /= (1.f - rr);
            }
        }
        return L;
            


        ///
        /// OLD VERSION, RECURSIVE, DOES NOT WORK PROPERLY
        ///


        /*
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
        */
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
        Color3f Lems(0.f);
        float pdflight(1.f);
        float pdir_wdir(1.f), pbsdf_wdir(0.f);

        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(its.p);
        BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
        its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());

        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray(its.p, emitterRecord.wi);
        //std::cout << "EmitterSampling: Pre-ShadowRay" << std::endl;
        std::vector<VolumetricSegmentRecord> shadow_vsr;
        if(emitterShRayIntersectFree(scene, sampler, shadowray, emitterRecord, currentVolumeMedium, shadow_vsr) && materialRecord.measure != EDiscrete)
        {
            //std::cout << "EmitterSampling: IN-ShadowRay" << std::endl;
            BSDFQueryRecord bsdfRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * volumetric_transmittance(sampler, shadow_vsr) * its.mesh->getBSDF()->eval(bsdfRecord) * abs(its.shFrame.n.dot(emitterRecord.wi)) / pdflight;

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

    Color3f volumetric_transmittance(Sampler*& sampler, const std::vector<VolumetricSegmentRecord>& vsr) const
    {
        //std::cout << "VOL_TRANS: [";
        Color3f throughput(1.f);

        // for(size_t i = 0; i < vsr.size(); i++)
        // {
        //     std::cout << "[" << vsr[i].segment_vol->isHeterogeneous() << " " << vsr[i].x0.toString() << " " << vsr[i].xs.toString() << "]";
        // }
        // std::cout << "]" << std::endl;
        
        for(size_t i = 0; i < vsr.size(); i++)
        {
            Color3f a = vsr[i].segment_vol->transmittance(sampler, vsr[i].x0, vsr[i].xs);
            float pdf = 1.0f; //vsr[i].vol_pdf;                                               /////// TODO: QUITAR
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
        
        /// Sample light source with Multiple Importance Sampling
        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(xt);
        Intersection shray_its;
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray = Ray3f(xt, emitterRecord.wi);

        std::vector<VolumetricSegmentRecord> shadow_vsr;
        if(emitterShRayIntersectFree(scene, sampler, shadowray, emitterRecord, currentVolumeMedium, shadow_vsr))
        {
            /// Compute Phase Function value using Emitter Sampling sampled direction
            PFQueryRecord pfRecord(-w, emitterRecord.wi);
            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * volumetric_transmittance(sampler, shadow_vsr) * currentVolumeMedium->getPhaseFunction()->eval(pfRecord) / pdflight;      /// TODO: He quitado el * mu_s en las reformas a iterativo según PBRBook

            //MIS weight for emitter sampling
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
                transmittance = volumetric_transmittance(sampler, vsr_shadow);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(its.p);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wbsdf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(sampler, vsr_shadow);
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
        Color3f Lpf(0.f);

        PFQueryRecord pfRecord(-w);
        fs = currentVolumeMedium->getPhaseFunction()->sample(pfRecord, sampler->next2D());
        wo_pf = pfRecord.wo;
        Intersection itsaux;

        Ray3f ray(xt, pfRecord.wo);
        std::vector<VolumetricSegmentRecord> vsr_shadow;
        bool intersects_aux = scene->rayIntersectThroughVolumes(sampler, ray, itsaux, currentVolumeMedium, vsr_shadow);
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
                transmittance = volumetric_transmittance(sampler, vsr_shadow);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(xt);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wpf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = volumetric_transmittance(sampler, vsr_shadow);
            }

            w_mis_indir = balanceHeuristic(ppf_wpf, pdir_wpf);

            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lpf = fs * emit * transmittance;   /// TODO: He quitado el * mu_s en las reformas a iterativo según PBRBook
        }
        return Lpf;
    }


    Color3f directLight(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w) const
    {
        Color3f Lems(0.f), Lmat(0.f), fs_bsdf(0.f);
        float w_mis_dir(0.f), w_mis_indir(0.f);
        Vector3f wo_brdf;

        Lems = emitterSampling(scene, sampler, currentVolumeMedium, its, w, w_mis_dir);
        if(!isnan(w_mis_dir))
                Lems *= w_mis_dir;
        //std::cout << "directLight postEmitterSampling" << std::endl;
        Lmat = brdfSampling(scene, sampler, currentVolumeMedium, its, w, w_mis_indir, wo_brdf, fs_bsdf);
        //std::cout << "directLight postBRDFSampling" << std::endl;
        if(!isnan(w_mis_indir))
                Lmat *= w_mis_indir;

        Color3f Lnee = Lems + Lmat;
        return Lnee;
        
        /*
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
        */
    }

    Color3f inscattering(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Point3f& xt, const Vector3f& w) const
    {
        
        Color3f Lems(0.f), Lpf(0.f);
        float w_mis_dir(0.f), w_mis_pf(0.f), fs_pf(Epsilon);
        Vector3f wo_pf;

        Lems = emitterSamplingPF(scene, sampler, currentVolumeMedium, xt, w, w_mis_dir);
        if(!isnan(w_mis_dir))
                Lems *= w_mis_dir;

        Lpf = phaseFunctionSampling(scene, sampler, currentVolumeMedium, xt, w, w_mis_pf, wo_pf, fs_pf);
        if(!isnan(w_mis_pf))
                Lpf *= w_mis_pf;

        Color3f Lnee = Lems + Lpf;
        
        return Lnee * fs_pf;    /// fs_pf???
        /*
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
        */
    }
};

NORI_REGISTER_CLASS(PathTracingMISParticipatingMedia, "path_mis_participating_media");
NORI_NAMESPACE_END