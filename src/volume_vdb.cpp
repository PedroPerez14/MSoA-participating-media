/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#include <nori/volume.h>
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


    Point3f samplePathStep(const Point3f& x, const Point3f& xz, const Vector3f& dir, 
                const Vector2f& sample, float& pdf_success, float& pdf_failure) const
    {
        if(!m_heterogeneous)
        {
            float _mu_t = sample_mu_t(x).sum() / 3.f;
            float t = -log(1.f - sample.y()) / _mu_t;
            float spectrum_mu_t = sample_mu_t(x).x();
            float tmp = exp(-spectrum_mu_t * t);
            pdf_failure += tmp;
            pdf_success += spectrum_mu_t * tmp;

            spectrum_mu_t = sample_mu_t(x).y();
            tmp = exp(-spectrum_mu_t * t);
            pdf_failure += tmp;
            pdf_success += spectrum_mu_t * tmp;

            spectrum_mu_t = sample_mu_t(x).z();
            tmp = exp(-spectrum_mu_t * t);
            pdf_failure += tmp;
            pdf_success += spectrum_mu_t * tmp;

            pdf_failure /= 3.f;
            pdf_success /= 3.f;

            if(isnan(x.z()))
                std::cout << "X.Z NAN" << std::endl;
            if(isnan(t))
                std::cout << "T NAN" << std::endl;
            if(isnan(dir.z()))
                std::cout << "DIR.Z NAN" << std::endl;

            return x + t * dir;
        }
        //else, heterogeneous volume
        // srand (static_cast <unsigned> (time(0)));

        // float tMax = Vector3f(xz - x0).norm();
        // Vector3f dir = (xz - x0).normalized();
        // float tr = 1.f;
        // float t = 0.f;
        // int maxSteps = 256;
        // int nSteps = 0;
        // while (true) {
        //     float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        //     t -= std::log(1.0f - r) / m_max / mu_t;
        //     //std::cout << "t_ratio: " << t << " " << tMax << std::endl;
        //     if(t >= tMax || nSteps >= maxSteps)
        //         return Color3f(tr, tr, tr);
        //     std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        //     std::cout << "DIR: " << dir.toString() << std::endl;
        //     std::cout << "T: " << t << std::endl;
        //     std::cout << "TMAX: " << tMax << std::endl;
        //     std::cout << "XZ: " << xz.toString() << std::endl;
        //     std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        //     float density = sample_density(x0 + t*dir).x();         /// TODO: We assume only 1 color for now, change it later(?)
        //     tr *= 1 - std::max((float)0, density / m_max);
        //     nSteps++;
        // }
        // return Color3f(tr, tr, tr);
        
    }

    Color3f transmittance(const Point3f& x0, const Point3f& xz) const
    {
        //if this is a homogeneous volume
        if(!m_heterogeneous)
        {
            float norm = (x0-xz).norm();
            Color3f expterm = sample_mu_t(xz) * -1.f;
            return exp(expterm * norm);
        }
        //else, we have a heterogeneous one

        /// TODO: Implementar uno de estos 3?:
        //deltaTracking()
        //ratioTracking()
        //monteCarloRayMarching()

        /// TODO: monocanal
        return m_volumegrid_mu_t->ratioTracking(x0, xz, mu_t.x());      
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