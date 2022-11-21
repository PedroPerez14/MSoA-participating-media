/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
	
	v1 - Dec 01 2020
    v2 - Oct 30 2021
	Copyright (c) 2021 by Adrian Jarabo

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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/reflectance.h>
#include <nori/texture.h>
#include <nori/dpdf.h>

NORI_NAMESPACE_BEGIN

#define KS_THRES 0.

class RoughConductor : public BSDF {
public:
    RoughConductor(const PropertyList& propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Reflectance at direction of normal incidence.
           To be used when defining the Fresnel term using the Schlick's approximation*/
        m_R0 = new ConstantSpectrumTexture(propList.getColor("R0", Color3f(0.5f)));
    }


    float D(const Vector3f& wh, const Vector2f& tex_uv) const   //Beckmann-Spizzichino normal distribution
    {
        return Reflectance::BeckmannNDF(wh, m_alpha->eval(tex_uv).getLuminance());
    }

    Color3f F(const float& cos_th_i, const Vector2f& tex_uv) const        //Schlick approximation (1993)
    {
        return Reflectance::fresnel(cos_th_i, m_R0->eval(tex_uv));
    }

    Color3f G(const Vector3f& wi, const Vector3f& wo, const Vector3f& wh, const Vector2f& uv_tex) const
    {
        float alpha = m_alpha->eval(uv_tex).getLuminance();
        return Reflectance::G1(wi, wh, alpha) * Reflectance::G1(wo, wh, alpha);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const
    {

        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Color3f f_r(0,0,0);
        Vector3f wi = bRec.wi, wo = bRec.wo;   // Si explota cambiar por      Vector3f wi = bRec.wo, wo = -bRec.wi       tal y como estaba
        Vector3f wh = (wi + wo).normalized();
        
        float cos_th_i = Frame::cosTheta(wi), cos_th_o = Frame::cosTheta(wo);

        f_r = (D(wh, bRec.uv) * F(cos_th_i, bRec.uv) * G(wi, wo, wh, bRec.uv)) / (4.0f * cos_th_i * cos_th_o);
        return f_r;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        Vector3f wh = bRec.wi + bRec.wo;
        wh.normalize();
        return Warp::squareToBeckmannPdf(wh, m_alpha->eval(bRec.uv).getLuminance());
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;


        Vector3f wh = Warp::squareToBeckmann(_sample, m_alpha->eval(bRec.uv).getLuminance());
        bRec.wo = (-bRec.wi) - 2.0f * ((-bRec.wi).dot(wh)) * wh;                                //Reflect wi with respect to the sampled wh
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "R0")
            {
                delete m_R0;
                m_R0 = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("RoughConductor::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughConductor::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughConductor[\n"
            "  alpha = %f,\n"
            "  R0 = %s,\n"
            "]",
            m_alpha->toString(),
            m_R0->toString()
        );
    }
private:
    Texture* m_alpha;
    Texture* m_R0;
};


class RoughDielectric : public BSDF {
public:
    RoughDielectric(const PropertyList& propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Tint of the glass, modeling its color */
        m_ka = new ConstantSpectrumTexture(propList.getColor("ka", Color3f(1.f)));
    }


    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BSDF -- return zero if the measure is wrong */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);


        throw NoriException("RoughDielectric::eval() is not yet implemented!");
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BSDF -- return zero if the measure is wrong */
        if (bRec.measure != ESolidAngle)
            return 0.0f;

        throw NoriException("RoughDielectric::pdf() is not yet implemented!");
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        bRec.measure = ESolidAngle;

        throw NoriException("RoughDielectric::sample() is not yet implemented!");
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "m_ka")
            {
                delete m_ka;
                m_ka = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("RoughDielectric::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughDielectric::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  ka = %s,\n"
            "]",
            m_alpha->toString(),
            m_intIOR,
            m_extIOR,
            m_ka->toString()
        );
    }
private:
    float m_intIOR, m_extIOR;
    Texture* m_alpha;
    Texture* m_ka;
};



class RoughSubstrate : public BSDF {
public:
    RoughSubstrate(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = new ConstantSpectrumTexture(propList.getColor("kd", Color3f(0.5f)));
    }

    float D(const Vector3f& wh, const Vector2f& tex_uv) const   //Beckmann-Spizzichino normal distribution
    {
        return Reflectance::BeckmannNDF(wh, m_alpha->eval(tex_uv).getLuminance());
    }

    float F(const float& cos_th_i, const float& extIOR, const float& intIOR) const
    {
        return Reflectance::fresnel(cos_th_i, extIOR, intIOR);
    }

    Color3f G(const Vector3f& wi, const Vector3f& wo, const Vector3f& wh, const Vector2f& uv_tex) const
    {
        float alpha = m_alpha->eval(uv_tex).getLuminance();
        return Reflectance::G1(wi, wh, alpha) * Reflectance::G1(wo, wh, alpha);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Color3f f_r(0.);
        Color3f f_diff(0.), f_mf(0.);
        Vector3f wi = bRec.wi, wo = bRec.wo;
        Vector3f wh = (wi + wo).normalized();
        Vector2f uv_tex = bRec.uv;


        float cos_th_i = Frame::cosTheta(wi), cos_th_o = Frame::cosTheta(wo);
        f_diff = (28.0f / 23.0f) * (m_kd->eval(uv_tex) * INV_PI) 
                * (1. - pow(((m_extIOR - m_intIOR) / (m_extIOR + m_intIOR)), 2.0f)) 
                * (1. - pow(1. - (0.5f * cos_th_i), 5.)) 
                * (1. - pow(1. - (0.5f * cos_th_o), 5.));

        f_mf = (D(wh, bRec.uv) * F(cos_th_i, m_extIOR, m_intIOR) * G(wi, wo, wh, bRec.uv)) / (4.0f * cos_th_i * cos_th_o);
        
        f_r = f_diff + f_mf;

        return f_r;
	}

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
       is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // Versión ponderada de los pesos
        Vector3f wh = bRec.wi + bRec.wo;
        wh.normalize();
        float prob_mf = F((Frame::cosTheta(bRec.wi)), m_extIOR, m_intIOR);
        float prob_diff = 1. - prob_mf;

        float weighted_beckmann_pdf = Warp::squareToBeckmannPdf(wh, m_alpha->eval(bRec.uv).getLuminance()) * prob_mf;
        float weighted_cosine_pdf = Warp::squareToCosineHemispherePdf(bRec.wo) * prob_diff;

        return weighted_beckmann_pdf + weighted_cosine_pdf;
    }

    // A second version for the pdf method, this time return only one pdf, no weighting
    //TODO probably will need to do 1. / this
    float pdf(const BSDFQueryRecord& bRec, const bool& sampleCosine) const
    {
        if(sampleCosine)
        {
            return Warp::squareToCosineHemispherePdf(bRec.wo);
        }
        //else, sample was generated using Beckmann, get its pdf
        Vector3f wh = bRec.wi + bRec.wo;
        return Warp::squareToBeckmannPdf(wh, m_alpha->eval(bRec.uv).getLuminance());
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);
        bRec.measure = ESolidAngle;

        Color3f col_sample(Epsilon);
        //Russian roulette , using p(f_ml) = F((n . wi), Next, Nint) and p(f_diff) = 1 - p(f_ml)
        Vector3f wi = bRec.wi;
        float prob_mf = F((Frame::cosTheta(wi)), m_extIOR, m_intIOR);
        float prob_diff = 1. - prob_mf;

        DiscretePDF rr_pdf(2);
        float sampleReused = _sample.x();
        float pdfval;
        rr_pdf.append(prob_mf);     // First, mf
        rr_pdf.append(prob_diff);   // Second, diff
        //rr_pdf.normalize();
        size_t index = rr_pdf.sampleReuse(sampleReused, pdfval);
        if(index == 0)              //Dielectric, beckmann sampling and pdf
        {
            //cout << "dielectric" << endl;
            Vector3f wh = Warp::squareToBeckmann(Point2f(sampleReused, _sample.y()), m_alpha->eval(bRec.uv).getLuminance()).normalized();
            bRec.wo = (-bRec.wi) - 2.0f * ((-bRec.wi).dot(wh)) * wh;                                //Reflect wi with respect to the sampled wh
            col_sample = eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        }
        else                        //Diffuse, sample direction by cosine
        {
            //cout << "diffuse" << endl;
            bRec.wo = Warp::squareToCosineHemisphere(Point2f(sampleReused, _sample.y()));
            col_sample = eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        }

        return col_sample;
	}

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "kd")
            {
                delete m_kd;
                m_kd = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else 
                throw NoriException("RoughSubstrate::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughSubstrate::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughSubstrate[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "]",
            m_alpha->toString(),
            m_intIOR,
            m_extIOR,
            m_kd->toString()
        );
    }
private:
    float m_intIOR, m_extIOR;
    Texture* m_alpha;
    Texture* m_kd;
};

NORI_REGISTER_CLASS(RoughConductor, "roughconductor");
NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_REGISTER_CLASS(RoughSubstrate, "roughsubstrate");

NORI_NAMESPACE_END
