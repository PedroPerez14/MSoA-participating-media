/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García and Adolfo Ber San Agustín
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#include <nori/phasefunction.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 h* \brief Rayleign Phase Function model
 */

class Rayleigh : public PhaseFunction {
public:
    Rayleigh()
    {
        wavelengths = Color3f(
            0.000000450f,
            0.000000532f,
            0.000000620f
        );
        m_diam = 0.000000001f;
        m_ior = 1.00029f;
        m_density = 1.0f;           ///TODO: Creo que esto varía exponencialmente con la altura
    }

    Rayleigh(const PropertyList &propList) {
        //m_albedo = new ConstantSpectrumTexture(propList.getColor("albedo", Color3f(0.5f)));
        wavelengths = propList.getColor("lambda",Color3f(0.000000450f,0.000000532f,0.000000620f));
        m_diam = float(propList.getFloat("diam", 0.000000001f));
        m_ior = float(propList.getFloat("ior", 1.00029f));
        m_density = float(propList.getFloat("ro", 1.0f));
    }

    ~Rayleigh()
    {
        std::cout << "DESTRUCTOR RS" << std::endl;
    }

    float eval(const PFQueryRecord &pRec) const
    {
        //WARNING: This assumes that you have already sampled a direction !!!!
        //          Call sample() before calling this !!!!1!
        float cos_th = pRec.wi.dot(pRec.wo); 
        return (3.f / 16.f) * INV_PI * (1.f + (cos_th * cos_th));
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(PFQueryRecord &pRec) const {
        /// TODO: hacer esto
        pRec.m_pdf = Warp::squareToRayleighPdf(pRec.wo);
        return pRec.m_pdf;
    }

    /// Draw a a sample from the BRDF model
    float sample(PFQueryRecord &pRec, const Point2f &sample) const
    {
        //TODO ME HE QUEDADO AQUÍ
        pRec.wo = Warp::squareToRayleigh(sample);
        float retVal = eval(pRec);
        pRec.m_pdf = pdf(pRec);
        //std::cout << "WO: " << pRec.wo.toString() << std::endl;
        return retVal / pRec.m_pdf;
    }

    const Color3f get_mu_s() const
    {
        float ior_2 = m_ior * m_ior;
        float ior_term = (ior_2 - 1.f) / (ior_2 + 2.f);
        Color3f mu_s = m_density * ((2.f * pow(M_PI, 5.f) * pow(m_diam, 6.f) / (3.f * pow(wavelengths, 4.f)) ) * pow(ior_term, 2.f));
        //std::cout << "mu_s: " << mu_s << std::endl;
        return mu_s;
    }


    bool isDiffuse() const {
        return false;       //TODO creo que puedo borrar esta función por completo
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Rayleigh\n"
            "  m_diam = %s\n"
            "  m_ior = %s\n"
            "  m_density = %s\n"
            "  wavelengths used = %s\n"
            "]", m_diam, m_ior, m_density, wavelengths);
    }

    EClassType getClassType() const { return EPhaseFunction; }
private:
    Color3f wavelengths = Color3f();         //Wavelengths for each color channel (m)
    float m_diam; //= 0.000000001f;         //diameter of our particles (m)
    float m_ior;                            //Index of refraction (try 1 or 1.00029 iirc)
    float m_density;                        //Density of the medium (m^⁻3)
};

NORI_REGISTER_CLASS(Rayleigh, "rayleigh");
NORI_NAMESPACE_END