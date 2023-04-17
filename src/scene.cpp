/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    v1 - Dec 01 2020
    Copyright (c) 2020 by Adrian Jarabo

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_accel = new Accel();
    m_enviromentalEmitter = nullptr;
    m_enviromentalVolumeMedium = nullptr;
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
}

void Scene::activate() {

    // Check if there's emitters attached to meshes, and
    // add them to the scene. 
    int n_emitters(0);
    for(unsigned int i=0; i<m_meshes.size(); ++i )
        if (m_meshes[i]->isEmitter())
            m_emitters.push_back(m_meshes[i]->getEmitter());

    for(unsigned int i=0; i<m_meshes.size(); ++i )
        if (m_meshes[i]->isVolume())
            m_volumes.push_back(m_meshes[i]->getVolume());
            


    m_accel->build();
    m_pdf.clear();
    m_pdf.reserve(n_emitters);

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }
    if(!m_enviromentalVolumeMedium)
    {
        m_enviromentalVolumeMedium = std::shared_ptr<Volume>(static_cast<Volume*>(NoriObjectFactory::createInstance("volumevdb", PropertyList())));
        m_enviromentalVolumeMedium->addChild(NoriObjectFactory::createInstance("henyey-greenstein", PropertyList()));
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

/// Sample emitter
const Emitter * Scene::sampleEmitter(float rnd, float &pdf) const {
	auto const & n = m_emitters.size();
	size_t index = std::min(static_cast<size_t>(std::floor(n*rnd)), n - 1);
	pdf = 1. / float(n);
	return m_emitters[index];
}

const std::shared_ptr<Volume> Scene::sampleVolume(Sampler* sampler, std::vector<std::shared_ptr<Volume>> vols, float& pdf) const
{
    if(vols.size() == 0)
    {
        pdf = 1.f;
        return m_enviromentalVolumeMedium;
    }
    //else
    auto const & n = vols.size();
	size_t index = std::min(static_cast<size_t>(std::floor(n*sampler->next1D())), n - 1);
	pdf = 1. / float(n);
	return vols[index];
}

bool Scene::rayIntersectThroughVolumes(Sampler* sampler, const Ray3f &ray, Intersection& its_out, std::shared_ptr<Volume> startVol, std::vector<VolumetricSegmentRecord>& segs) const
{
    bool done = false;
    std::vector<std::shared_ptr<Volume>> volumes_traversed;
    volumes_traversed.reserve(m_volumes.size());
    volumes_traversed.push_back(startVol);
    Ray3f _ray(ray);
    int n_bounces = 0;
    while(true)
    {
        if(n_bounces > 100)
        {
            std::cout << "rayIntersectThroughVolumes(no xt): Too many bounces. Is everything alright? Aborting..." << std::endl;
            return false;
        }
        n_bounces++;


        Intersection its;
        bool intersects = rayIntersect(_ray, its);
        // If we haven´t found the final point of our ray (either xt or its.p)
        if(intersects)
        {
            // float t = (xt - _ray.o).norm();
            // float z = (its.p - _ray.o).norm();
            if(its.mesh->isVolume())
            {
                /// Now we choose one volume from the ones we are currently traversing
                ///     This is done uniformly (for now), if none traversed, return envVolume
                float volume_pdf;
                Vector3f dir = _ray.d;
                bool delete_vol = false;
                _ray = Ray3f(its.p + (Epsilon * dir), dir);
                std::shared_ptr<Volume> sampled_vol;
                
                // Add the volume to our "currently traversed" list if we aren't already traversing it
                // If we are, remove it from the list, because we have just got out of it
                auto iter = std::find(volumes_traversed.begin(), volumes_traversed.end(), its.mesh->getVolume());
                if(iter != volumes_traversed.end())
                {
                    //If on the list, delete it because we are getting outside of it
                    sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);
                    volumes_traversed.erase(iter);
                }
                else
                {
                    //Otherwise, it's a new volume (we are entering it) and we add it to the list
                    sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);
                    volumes_traversed.push_back(its.mesh->getVolume());
                }

                
                VolumetricSegmentRecord vsr(_ray.o, its.p, sampled_vol, volume_pdf);
                //ISNAN ?
                //std::cout << "_NO_XT_INTER_VOL: " << vsr.xs << std::endl;
                segs.push_back(vsr);

                continue;
            }
            else
            {
                float volume_pdf;
                std::shared_ptr<Volume> sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);
                VolumetricSegmentRecord vsr(_ray.o, its.p, sampled_vol, volume_pdf);
                //ISNAN
                //std::cout << "_NO_XT_INTER_NOVOL: " << vsr.xs << std::endl;
                segs.push_back(vsr);
                its_out = its;
                return true;
            }
        }
        else
        {
            float volume_pdf;
            std::shared_ptr<Volume> sampled_vol = m_enviromentalVolumeMedium;//sampleVolume(sampler, volumes_traversed, volume_pdf);

            VolumetricSegmentRecord vsr(_ray.o, _ray.o + ray.d * 1.f, sampled_vol, volume_pdf);
            //ISNAN
            //std::cout << "KAMEHAMEHA: " << vsr.xs << std::endl;
            //std::cout << "_NO_XT_NOINTER: " << vsr.xs << std::endl;
            segs.push_back(vsr);
            return false;
        }
    }
}


bool Scene::shadowRayThroughVolumes(Sampler* sampler, const Ray3f& sray, std::shared_ptr<Volume> startVol, std::vector<VolumetricSegmentRecord>& segs_shadow, float& t) const
{
    t = 0.0f;
    std::vector<std::shared_ptr<Volume>> volumes_traversed;
    volumes_traversed.reserve(m_volumes.size());
    volumes_traversed.push_back(startVol);
    Ray3f ray(sray);
    int n_bounces = 0;
    while(true)
    {
        if(n_bounces > 100)
        {
            //std::cout << "rayIntersectThroughVolumes(no xt): Too many bounces. Is everything alright? Aborting..." << std::endl;
            return false;
        }
        n_bounces++;        //I could put this here or before every continue, this one's easier though

        Intersection its;
        bool intersects = rayIntersect(ray, its);
        if(!intersects)
        {
            float volume_pdf;
            std::shared_ptr<Volume> sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);

            VolumetricSegmentRecord vsr(ray.o, ray.o + ray.d * 1.f, sampled_vol, volume_pdf);
            //std::cout << "_SHRAY_NO_INTER: " << vsr.xs << std::endl;
            segs_shadow.push_back(vsr);
            t += 1000.f;            //TODO distancia si se me va al infinito un rayo
            return false;
        }
        else
        {
            if(its.mesh->isVolume())
            {
                float volume_pdf;
                std::shared_ptr<Volume> sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);
                VolumetricSegmentRecord vsr(ray.o, its.p, sampled_vol, volume_pdf);
                t += its.t;
                ray = Ray3f(its.p + (ray.d * Epsilon), ray.d);
                //continue;
            }
            else
            {
                float volume_pdf;
                std::shared_ptr<Volume> sampled_vol = sampleVolume(sampler, volumes_traversed, volume_pdf);
                VolumetricSegmentRecord vsr(ray.o, its.p, sampled_vol, volume_pdf);
                t += its.t;
                return true;
            }
        }
    }
}




/// Sample emitter with importance sampling
const Emitter * Scene::sampleEmitter(Sampler *sampler, float &pdf, EmitterQueryRecord lRec) const {
	
    DiscretePDF m_pdf = this->m_pdf;
    auto const & n = m_emitters.size();
    if(n == 1)  //If only one emitter, no need to do calculations and samples
    {
        pdf = 1.0f;
        return m_emitters[0];
    }
    //else
    for(size_t i = 0; i < n; i++)
    {
        Color3f rad = m_emitters[i]->sample(lRec, sampler->next2D(), 0.f);
        float maxRadianceCoeff = std::max(std::max(rad.x(), rad.y()), rad.z());
        m_pdf.append(maxRadianceCoeff);
    }
    m_pdf.normalize();
    size_t idx = m_pdf.sample(sampler->next1D(), pdf);
    m_pdf.clear();
    pdf = 1. / std::max(Epsilon,pdf);
    return m_emitters[idx];
}

float Scene::pdfEmitter(const Emitter *em) const {
    return 1. / float(m_emitters.size());
}


void Scene::addChild(NoriObject *obj, const std::string& name) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
            }
            break;
        
        case EEmitter: {
				Emitter *emitter = static_cast<Emitter *>(obj);
				if (emitter->getEmitterType() == EmitterType::EMITTER_ENVIRONMENT)
				{
					if (m_enviromentalEmitter)
						throw NoriException("There can only be one enviromental emitter per scene!");
					m_enviromentalEmitter = emitter;
				}
				
                m_emitters.push_back(emitter);
			}
            break;

        case ESampler: {
                if (m_sampler)
                    throw NoriException("There can only be one sampler per scene!");
                m_sampler = static_cast<Sampler *>(obj);
            }
            break;

        case ECamera:
            {
                if (m_camera)
                throw NoriException("There can only be one camera per scene!");
                m_camera = static_cast<Camera *>(obj);
            }
            break;
        
        case EIntegrator: {
                if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
                m_integrator = static_cast<Integrator *>(obj);
            }
            break;

        case EVolume: {
                if (m_enviromentalVolumeMedium)
                {
                    throw NoriException("There can only be one enviromental volume medium per scene!");
                }
                Volume * vol_ptr = static_cast<Volume *>(obj);
                m_enviromentalVolumeMedium = std::shared_ptr<Volume>(vol_ptr);
            }
            break;

        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

Color3f Scene::getBackground(const Ray3f& ray) const
{
    if (!m_enviromentalEmitter)
        return Color3f(0);

    EmitterQueryRecord lRec(m_enviromentalEmitter, ray.o, ray.o + ray.d, Normal3f(0, 0, 1), Vector2f());
	return m_enviromentalEmitter->eval(lRec);
}


std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    std::string volumes;
    for (size_t i = 0; i < m_volumes.size(); ++i) {
		volumes += std::string("  ") + indent(m_volumes[i]->toString(), 2);
		if (i + 1 < m_volumes.size())
			volumes += ",";
		volumes += "\n";
	}

	std::string lights;
	for (size_t i = 0; i < m_emitters.size(); ++i) {
		lights += std::string("  ") + indent(m_emitters[i]->toString(), 2);
		if (i + 1 < m_emitters.size())
			lights += ",";
		lights += "\n";
	}


    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "  volumes = {\n"
        "  %s  }\n"
		"  emitters = {\n"
		"  %s  }\n"
        "  envVolume = %s\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2),
        indent(volumes, 2),
		indent(lights, 2),
        indent(m_enviromentalVolumeMedium->toString())
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
