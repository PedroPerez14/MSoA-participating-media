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
 * \brief Henyey-Greenstein Phase Function model
 */

class HenyeyGreenstein : public PhaseFunction {
public:
    HenyeyGreenstein()
    {
        g = 0.f;
        mu_s = Color3f(0.05f);
    }

    HenyeyGreenstein(const PropertyList &propList) {
        //m_albedo = new ConstantSpectrumTexture(propList.getColor("albedo", Color3f(0.5f)));
        g = float(propList.getFloat("g", 0.f));
        mu_s = Color3f(propList.getColor("mu_s", Color3f(0.05f)));
    }

    ~HenyeyGreenstein()
    {
        std::cout << "DESTRUCTOR HG" << std::endl;
    }

    /// Evaluate the BRDF model
    float eval(const PFQueryRecord &pRec) const
    {
        //WARNING: This assumes that you have already sampled a direction !!!!
        //          Call sample() before calling this !!!!1!
        float cos_th = pRec.wi.dot(pRec.wo); 
        float denom_term = sqrt(pow((1.f + pow(g, 2.f) - (2.f * g) * cos_th), 3.f));
        return (0.25f * INV_PI) * ((1.f - pow(g, 2.f)) / (denom_term));   // (1/4pi) * (... / ...)
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(PFQueryRecord &pRec) const {
        pRec.m_pdf = Warp::squareToHenyeyGreensteinPdf(pRec.wi.dot(pRec.wo), get_g());
        //std::cout << "PDF: " << pRec.m_pdf << std::endl;
        return pRec.m_pdf;
    }

    /// Draw a a sample from the BRDF model
    float sample(PFQueryRecord &pRec, const Point2f &sample) const
    {
        //TODO ME HE QUEDADO AQUÍ
        pRec.wo = Warp::squareToHenyeyGreenstein(sample, get_g());
        float retVal = eval(pRec);
        pRec.m_pdf = pdf(pRec);
        //std::cout << "WO: " << pRec.wo.toString() << std::endl;
        return retVal / pRec.m_pdf;
    }

    float get_g() const
    {
        return g;
    }

    /// TODO: Is this always constant?????
    const Color3f get_mu_s() const
    {
        return mu_s;
    }


    bool isDiffuse() const {
        return false;       //TODO creo que puedo borrar esta función por completo
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Henyey-Greenstein[\n"
            "  g = %s\n"
            "  mu_s = %s\n"
            "]", g, mu_s);
    }

    /*
    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "albedo")
            {
                delete m_albedo;
                m_albedo = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("Diffuse::addChild(<%s>,%s) is not supported!",
                classTypeName(obj->getClassType()), name);
            break;

        default:
            throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }
    */

    EClassType getClassType() const { return EPhaseFunction; }
private:
    float g;                        //Should this be calculated using wi and wo? Who knows
    Color3f mu_s;                   //We might want to scatter colors differently, i.e Rayleigh Scattering
};

NORI_REGISTER_CLASS(HenyeyGreenstein, "henyey-greenstein");
NORI_NAMESPACE_END