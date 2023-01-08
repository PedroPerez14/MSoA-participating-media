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
        std::shared_ptr<Volume> currentVolumeMedium = scene->getEnviromentalVolumeMedium();//getCurrentVolume(volumes_through, scene->getVolumes(), scene);


        /// TODO: Voy a asumir por el bien de mi salud mental, que DE MOMENTO solo hay un volumen ambiental infinito
        //Color3f Lems(0.f), Lmat(0.f);
        Color3f L(0.f);
        float pdf_pathstep, pdf_succ, pdf_fail;
        
        //Boundary of the medium
        bool intersected = scene->rayIntersect(ray, its);      //xz = its.p
        
        //Inscattering
        //Point3f xt(0.f, 0.f, 0.f);
        Point3f xt = currentVolumeMedium->samplePathStep(ray.o, its.p, ray.d, sampler->next2D(), pdf_succ, pdf_fail);
        float z = Vector3f(its.p - ray.o).norm();
        float t = Vector3f(xt - ray.o).norm();
        
        if(intersected && t >= z)
        {
            //We hit the surface!
            if(its.mesh->isEmitter())
            {
                EmitterQueryRecord er(its.mesh->getEmitter(), ray.o, its.p, its.shFrame.n, its.uv);
                L = currentVolumeMedium->transmittance(ray.o, its.p) * its.mesh->getEmitter()->eval(er);
            }
            else
            {
                L = currentVolumeMedium->transmittance(ray.o, its.p) * directLight(scene, sampler, currentVolumeMedium, its, ray.d);
            }
            pdf_pathstep = pdf_fail;
        }
        else
        {
            //Medium interaction
            L = currentVolumeMedium->transmittance(ray.o, xt) * inscattering(scene, sampler, currentVolumeMedium, its, xt, ray.d);
            pdf_pathstep = pdf_succ;
        }
        return L / pdf_pathstep;
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

    /// TODO: WARNING: PROBLEM WITH POINT LIGHTS AND TRANSMITTANCE (NO INTERSECTION, SO ITS_SHRAY WILL NOT BE INITIALIZED, CAREFUL!!!!!!!)
    bool emitterShRayIntersectFree(const Scene *scene, Ray3f &sray, const EmitterQueryRecord& emitterRecord, Intersection &it_shadow) const
    {
        if(!scene->rayIntersect(sray, it_shadow))
        {
            return true;
        }
        else
        {
            /// TOOD: De momento lo modifico para que salga del propio volumen donde está
            /// TODO: Modificar pata volúmenes anidados
            //if(it_shadow.mesh->isVolume())
            //{
                /// TODO:  Voy a tener que implementar una cola o algo así
                // Aunque para SS puedo suponer solo el scatterer global
                // Y aparcadísimo
            //}

            return (it_shadow.t > (emitterRecord.dist - Epsilon));
        }
    }

    Color3f emitterSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w, float& w_mis_dir) const
    {
        Color3f Lems(Epsilon);
        float pdflight(1.f);
        float pdir_wdir(0.f), pbsdf_wdir(0.f);

        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(its.p);
        BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
        Intersection shray_its;
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray(its.p, emitterRecord.wi);
        

        ///TODO: We are ignoring intersections with other volumes, be careful!
        if(emitterShRayIntersectFree(scene, shadowray, emitterRecord, shray_its) && materialRecord.measure != EDiscrete)
        {
            BSDFQueryRecord bsdfRecord(its.toLocal(-w), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * currentVolumeMedium->transmittance(its.p, emitterRecord.p) * its.mesh->getBSDF()->eval(bsdfRecord) * abs(its.shFrame.n.dot(emitterRecord.wi)) / pdflight;

            //MIS for emitter sampling
            pdir_wdir = em->pdf(emitterRecord);
            pbsdf_wdir = its.mesh->getBSDF()->pdf(bsdfRecord);
            
            if(!its.mesh->getBSDF()->isDiffuse())
                pbsdf_wdir = emitterRecord.pdf; 
            w_mis_dir = balanceHeuristic(pdir_wdir, pbsdf_wdir);
        }
        return Lems;
        /*
        Ray3f _ray(ray);
        return Li_iter(scene, sampler, _ray);
        //BSDFQueryRecord bsdf_rec(Vector3f(0.f), Vector2f(0.f));
        //return LiRec(scene, sampler, ray, Color3f(1.f), 0, bsdf_rec);
        */
        
    }

    /// Emitter Sampling, but the weights for MIS are calculated using phase function instead of BSDF
    Color3f emitterSamplingPF(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Point3f& xt, const Vector3f& w, float& w_mis_dir) const
    {
        Color3f Lems(Epsilon);
        float pdflight(1.f);
        float pdir_wdir(0.f), ppf_wdir(0.f);
        
        const Emitter* em = scene->sampleEmitter(sampler->next1D(), pdflight);
        EmitterQueryRecord emitterRecord(xt);
        Intersection shray_its;
        Color3f Le = em->sample(emitterRecord, sampler->next2D(), 0.f);
        
        Ray3f shadowray = Ray3f(xt, emitterRecord.wi);

        ///TODO: We are ignoring intersections with other volumes, be careful!
        if(emitterShRayIntersectFree(scene, shadowray, emitterRecord, shray_its))
        {
            PFQueryRecord pfRecord(-w, emitterRecord.wi);
            Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lems = Le * currentVolumeMedium->transmittance(xt, emitterRecord.p) * currentVolumeMedium->getPhaseFunction()->eval(pfRecord) * mu_s / pdflight;

            //MIS for emitter sampling
            pdir_wdir = em->pdf(emitterRecord);
            ppf_wdir = currentVolumeMedium->getPhaseFunction()->pdf(pfRecord);
            
            w_mis_dir = balanceHeuristic(pdir_wdir, ppf_wdir);
                        
        }
        return Lems;
    }


    Color3f brdfSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w, float& w_mis_indir) const
    {
        Color3f Lmat(Epsilon);
        //Second parameter, its wo, is random, since it will get properly calculated with sample()
        BSDFQueryRecord materialRecord = BSDFQueryRecord(its.toLocal(-w), its.toLocal(w), its.uv, ESolidAngle);
        Color3f fs = its.mesh->getBSDF()->sample(materialRecord, sampler->next2D());
        Intersection itsaux;
        Ray3f ray(its.p, its.toWorld(materialRecord.wo));
        bool intersects_aux = scene->rayIntersect(ray, itsaux);
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
                transmittance = currentVolumeMedium->transmittance(its.p, itsaux.p);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(its.p);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wbsdf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = currentVolumeMedium->transmittance(its.p, emitterRecordAux.p);
            }

            if(materialRecord.measure == EDiscrete)
            {
                pbsdf_wbsdf = 1.f;
            }
            w_mis_indir = balanceHeuristic(pbsdf_wbsdf, pdir_wbsdf);
            
            //Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lmat = fs * emit * transmittance;
        }
        return Lmat;
    }

    Color3f phaseFunctionSampling(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Point3f& xt, const Vector3f& w, float& w_mis_indir) const
    {
        Color3f Lpf(0.f);

        PFQueryRecord pfRecord(-w);
        Color3f fs = currentVolumeMedium->getPhaseFunction()->sample(pfRecord, sampler->next2D());
        Intersection itsaux;

        bool intersects_aux = scene->rayIntersect(Ray3f(xt, pfRecord.wo), itsaux);
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
                transmittance = currentVolumeMedium->transmittance(xt, itsaux.p);
            }
            else
            {
                //Enviromental emitters case
                EmitterQueryRecord emitterRecordAux(xt);
                //emitterRecordAux.emitter->sample(emitterRecordAux, sampler->next2D(), 0.);
                scene->getEnvironmentalEmitter()->sample(emitterRecordAux, sampler->next2D(), 0.f);
                pdir_wpf = scene->getEnvironmentalEmitter()->pdf(emitterRecordAux);
                emit = scene->getEnvironmentalEmitter()->eval(emitterRecordAux);
                transmittance = currentVolumeMedium->transmittance(xt, emitterRecordAux.p);
            }

            w_mis_indir = balanceHeuristic(ppf_wpf, pdir_wpf);

            Color3f mu_s = currentVolumeMedium->getPhaseFunction()->get_mu_s();
            Lpf = fs * emit * transmittance * mu_s;
        }
        return Lpf;
    }


    Color3f directLight(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Vector3f& w) const
    {
        Color3f Lems(0.f), Lmat(0.f);
        float w_mis_dir(0.f), w_mis_indir(0.f);

        Lems = emitterSampling(scene, sampler, currentVolumeMedium, its, w, w_mis_dir);
        if(!isnan(w_mis_dir) && w_mis_dir != 0.f)
                Lems *= w_mis_dir;

        Lmat = brdfSampling(scene, sampler, currentVolumeMedium, its, w, w_mis_indir);
        if(!isnan(w_mis_indir) && w_mis_indir != 0.f)
                Lmat *= w_mis_indir;

        return Lems + Lmat;
    }

    Color3f inscattering(const Scene* scene, Sampler* sampler, std::shared_ptr<Volume> currentVolumeMedium, const Intersection& its, const Point3f& xt, const Vector3f& w) const
    {
        Color3f Lems(0.f), Lpf(0.f);
        float w_mis_dir(0.f), w_mis_pf(0.f);

        Lems = emitterSamplingPF(scene, sampler, currentVolumeMedium, xt, w, w_mis_dir);
        if(!isnan(w_mis_dir) && w_mis_dir != 0.f)
                Lems *= w_mis_dir;

        Lpf = phaseFunctionSampling(scene, sampler, currentVolumeMedium, xt, w, w_mis_pf);
        if(!isnan(w_mis_pf) && w_mis_pf != 0.f)
                Lpf *= w_mis_pf;


        return Lems + Lpf;
    }

    std::shared_ptr<Volume> getCurrentVolume(const std::vector<unsigned int>& vols_through
            , const std::vector<std::shared_ptr<Volume>> &vols, const Scene* scene) const
    {
        if(vols_through.size() == 0)
        {
            return scene->getEnviromentalVolumeMedium();
        }
        // else
        //always return "deepest" volume if nested
        return vols[vols_through[vols_through.size() - 1]];
    }

    /// TODO: Check, I don't trust this one
    bool isVolumeBeingTraversed(const Volume * vol, const std::vector<unsigned int>& vols_through,
                                const std::vector<Volume *>& vols, int& _index) const
    {
        for(size_t i = 0; i < vols.size(); i++)
        {
            if(vol == vols[vols_through[i]])
            {
                _index = vols_through[i];
                return true;
            }
        }
        _index = -1;
        return false;
    }

    void markVolumeAsTraversed(const Volume * vol, std::vector<unsigned int>& vols_through,
                                    const std::vector<Volume *>& vols, int& _index) const
    {
        for(int i = 0; i < static_cast<int>(vols.size()); i++)
        {
            if(vols[i] == vol)
            {
                vols_through.push_back(static_cast<unsigned int>(i));
                _index = i;
            }
        }
    }

    //Also reassings currentVolumeMedium pointer in case that was the medium we just got out of
    void removeVolumeFromTraversed(std::vector<unsigned int>& vols_through, int& _index, std::shared_ptr<Volume>& currentVol,
                                             const std::vector<std::shared_ptr<Volume>> &scene_volumes, const Scene* scene) const
    {
        for(size_t i = 0; i < vols_through.size(); i++)
        {
            int vols_through_i = vols_through[i];
            if(vols_through_i == _index)
            {
                vols_through.erase(vols_through.begin() + i);
                if(currentVol == scene_volumes[vols_through_i])
                {
                    if(vols_through.size() == 0)
                    {
                        currentVol = scene->getEnviromentalVolumeMedium();
                    }
                    else
                    {
                        currentVol = scene_volumes[vols_through[(i - 1) % vols_through.size()]];
                    }
                }
            }
        }
    }

    bool isCameraInsideVolume(const Camera* cam, const std::shared_ptr<Volume> vol) const
    {
        /// TODO: Hacer, tal vez lanzando un rayo hacia arriba y guardando en orden inverso a lo que encuentre 
        ///        + checkear prod. vect para saber si entro o salgo
        /*if(true)
            return true;
        return false;*/
        return false;
    }

    void isCameraInsideAnyVolume(const Scene *scene, std::vector<unsigned int> &vols_through) const
    {
        std::vector<std::shared_ptr<Volume>> volumes = scene->getVolumes();
        const Camera * cam = scene->getCamera();
        for(unsigned int i = 0; i < volumes.size(); i++)
        {
            if(isCameraInsideVolume(cam, volumes[i]))
                vols_through.push_back(i);
        }
    }
};

NORI_REGISTER_CLASS(PathTracingMISParticipatingMedia, "path_mis_participating_media");
NORI_NAMESPACE_END