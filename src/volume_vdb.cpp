/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García and Adolfo Ber San Agustín
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#include <nori/volume.h>

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
        mu_t = props.getColor("mu_t", Color3f(.15f));
        mu_a = props.getColor("mu_a", Color3f(0.f));        // TODO: not supported (yet?), but we can load it even though it does nothing
        /*m_vdb_file_name = props.getString("vdb_filename", "null");
        m_mu_a_vdb_name = props.getString("vdb_mu_t_name", "density");
        m_mu_a_vdb_name = props.getString("vdb_mu_a_name", "temperature");*/
        /// TODO: Call the constructor for mu_{t,a}_structure_to_be_determined_dont_use_this when I implement it
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

    Point3f samplePathStep(const Point3f& x, const Point3f& xz, const 
        Vector3f& dir, const Vector2f& sample, float& pdf_success, float& pdf_failure) const
    {
        float _mu_t;
        if(sample.x() <= (1.f/3.f))
            _mu_t = sample_mu_t(x).x();
        else if(sample.x() > (1.f/3.f) && sample.x() <= (2.f/3.f))
            _mu_t = sample_mu_t(x).y();
        else
            _mu_t = sample_mu_t(x).z();    //else

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
        return Point3f(x + t * dir);
    }

    Color3f transmittance(const Point3f x0, const Point3f xz) const
    {
        return exp(-sample_mu_t(x0) * Vector3f(x0-xz).norm());
    }

    Color3f sample_mu_t(const Point3f& p_world) const
    {
        return mu_t * 1.f;
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
            "]", mu_t, mu_a, mu_a, mu_a, mu_a, m_phase_function->toString());
    }

private:
    /*std::string     m_vdb_file_name;
    std::string     m_mu_t_vdb_name;
    std::string     m_mu_a_vdb_name;*/
    ///TODO: Don't use the 2 variables below, at least for now until i support (or quit attempting to) .vdb files
    //void            *mu_t_structure_to_be_determined_dont_use_this;
    //void            *mu_a_structure_to_be_determined_dont_use_this;
};


NORI_REGISTER_CLASS(VolumeVDB, "volumevdb")
NORI_NAMESPACE_END