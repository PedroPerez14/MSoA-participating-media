/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#include <nori/volume.h>
#include <nori/mesh.h>
#include <nori/volumedatabase.h>

NORI_NAMESPACE_BEGIN

class VolumeVDB : public Volume
{
public:

    VolumeVDB()
    {
        mu_t = 0.05f;
        mu_a = 0.f;
        /*m_vdb_file_name = "null";
        m_mu_t_vdb_name = "density";
        m_mu_a_vdb_name = "temperature";*/
    }

    VolumeVDB(const PropertyList& props)
    {
        mu_t = props.getColor("mu_t", Color3f(0.f));
        mu_a = props.getColor("mu_a", Color3f(0.f));        // TODO: not supported (yet?), but we can load it even though it does nothing
        m_vdb_file_name = props.getString("vdb_filename", "null");
        m_mu_t_vdb_name = props.getString("vdb_mu_t_name", "density");
        m_heterogeneous = false;
        //m_mu_a_vdb_name = props.getString("vdb_mu_a_name", "temperature");*/
        /// TODO: Call the constructor for mu_{t,a}_structure_to_be_determined_dont_use_this when I implement it
        if(m_vdb_file_name != "null")
        {
            Transform trafo = props.getTransform("toWorld", Transform());
            std::cout << "TRANSFORM: " << trafo.toString() << std::endl;
            m_volumegrid_mu_t = std::make_shared<Volumedatabase>(m_vdb_file_name, trafo);
            m_heterogeneous = true;
        }
    }

    /// TODO: Free memory here if needed
    ~VolumeVDB()
    {
        std::cout << "DESTRUCTOR VDB" << std::endl;
        /*if(m_phase_function)
            delete m_phase_function;

        /// TODO: Change this in a future version when i finally load vdb files and know how to store them
        if(mu_t_structure_to_be_determined_dont_use_this)
            delete mu_t_structure_to_be_determined_dont_use_this;
        if(mu_a_structure_to_be_determined_dont_use_this)
            delete mu_a_structure_to_be_determined_dont_use_this;
            */
    }

    float pdfFail(const Point3f& xz, const float& z, const Vector2f& sample) const
    {
        float pdf_failure = 1.f;
        float spectrum_mu_t = sample_mu_t(xz).x();
        float tmp = exp(-spectrum_mu_t * z);
        pdf_failure += tmp;

        spectrum_mu_t = sample_mu_t(xz).y();
        tmp = exp(-spectrum_mu_t * z);
        pdf_failure += tmp;

        spectrum_mu_t = sample_mu_t(xz).z();
        tmp = exp(-spectrum_mu_t * z);
        pdf_failure += tmp;

        pdf_failure /= 3.f;

        return pdf_failure;
    }


    Point3f samplePathStep(const Ray3f& ray, const Intersection& its, Sampler*& sampler, Color3f& _beta, std::shared_ptr<Volume>& nextVolume, std::shared_ptr<Volume>& currentVolume, bool& sampledMedium) const
    {
        //If our medium is homogeneous, it's easier
        if(!m_heterogeneous)
        {
            Point3f sampled_point;
            int channel = std::min(static_cast<int>(std::floor(3 * sampler->next1D())), 2);

            /// We compute a t and check if it is medium interaction or surface interaction
            //Since this is a homogeneous medium the point we pass here is irrelevant
            float _mu_t = sample_mu_t(Point3f(0.f))[channel];
            float dist = -std::log(1.f - sampler->next1D()) / _mu_t;
            float t = std::min(dist * ray.d.norm(), its.t);
            sampledMedium = t < its.t;
            nextVolume = currentVolume;

            if(!sampledMedium && its.mesh->isVolume())
            {
                nextVolume = its.mesh->getVolume();
            }

            sampled_point = ray(t);

            /// Compute transmittance and sample density
            Color3f Tr = std::exp(-_mu_t * t * ray.d.norm());   // Hace falta poner ray.d.norm() si es 1 siempre? (creo)????

            Color3f density = sampledMedium ? (_mu_t * Tr) : Tr;
            float pdf = 0.f;
            for(int i = 0; i < 3; i++)
            {
                pdf += density[i];
            }
            pdf *= 1.f / 3.f;

            /// If we sampled within the medium, use the right pdf for the weighting, otherwise, 1-cdf
            if(sampledMedium)
            {
                _beta = Tr * currentVolume->getPhaseFunction()->get_mu_s() / pdf;
            }
            else
            {
                _beta = Tr / pdf;
            }
            return sampled_point;
        }
        //else, heterogeneous volume

        Color3f mu_s = m_phase_function->get_mu_s();
        Point3f sampled_point = m_volumegrid_mu_t->samplePathStep(ray, its, sampler, mu_s, mu_t, _beta, sampledMedium);

        if(sampledMedium)
        {
            nextVolume = currentVolume;
        }
        else if(its.mesh->isVolume())
        {
            nextVolume = its.mesh->getVolume();
        }

        return sampled_point;
    }

    Color3f transmittance(Sampler*& sampler, const Point3f& x0, const Point3f& xz) const
    {
        //if this is a homogeneous volume
        if(!m_heterogeneous)
        {
            float norm = std::min((xz-x0).norm(), std::numeric_limits<float>::infinity());
            Color3f expterm = sample_mu_t(xz) * -1.f;
            return exp(expterm * norm);
        }
        //else, we have a heterogeneous one
        /// Perform ratio tracking

        /// TODO: monocanal + revisar para mi queridísimo renderizador 2.0
        return m_volumegrid_mu_t->ratioTracking(sampler, x0, xz, (mu_t.sum() / 3.f));      
    }

    Color3f sample_mu_t(const Point3f& p_world) const
    {
        //If we have a heterogeneous volume
        //Look up for the values
        if(!m_heterogeneous)
        {
            return mu_t;
        }
        //else, we have an heterogeneous volume, return mu_t
        /// TODO: We have to check the grid
        return m_volumegrid_mu_t->sample_density(p_world);

        //float rand = (static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX)) * 2.f;
        //return Color3f(rand * mu_t);
        
        /*if(mu_t_structure_to_be_determined_dont_use_this)
        {
            return mu_t * 1.f; // change the 1.f for mu_t_str.sample() or eval();   /// TODO: Future version
        }*/
        
    }

    Color3f sample_mu_a(const Point3f& p_world) const
    {
        return mu_a * 1.f;
        /*if(mu_a_structure_to_be_determined_dont_use_this)
        {
            return mu_a * 1.f; // change the 1.f for mu_a_str.sample() or eval();   /// TODO: Future version
        }
        */
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case EPhaseFunction:
            {
                if(m_phase_function)
                    throw NoriException("Volume: Tried to register multiple PF instances!");
                PhaseFunction* _pf = static_cast<PhaseFunction*>(obj);
                m_phase_function = std::shared_ptr<PhaseFunction>(_pf);
            }
            break;
        default:
            throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
            break;
        }
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "VolumeVDB[\n"
            "  mu_t = %s\n"
            "  mu_a = %s\n"
            "  mu_t_channel = %s\n"
            "  mu_a_channel = %s\n"
            "  vdb_file = %s\n"
            "  phase_function = {\n"
            "  %s  }\n"
            "]", mu_t, mu_a, m_mu_t_vdb_name, mu_a, m_vdb_file_name, m_phase_function->toString());
    }

private:
    std::string     m_vdb_file_name;
    std::string     m_mu_t_vdb_name;
    std::shared_ptr<Volumedatabase> m_volumegrid_mu_t;

    /// TODO: Tengo que decidir si lo controlo todo en función de albedo y mu_t o de mu_t y mu_a
    /// TODO: Mira el chat de Discord con Nestor para esto!!
    //Volumedatabase  *m_volumegrid_albedo;
};


NORI_REGISTER_CLASS(VolumeVDB, "volumevdb")
NORI_NAMESPACE_END