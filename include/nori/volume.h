/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García and Adolfo Ber San Agustín
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#pragma once

#include <nori/object.h>
#include <nori/phasefunction.h>
#include <nori/frame.h>             //TODO: Quitar? no sé si me hará falta para cuando tenga los .vdb
#include <nori/warp.h>              //TODO: Quitar? ...
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass for volumes (TODO: openVDB / 3d perlin noise textures support)
 */
class Volume : public NoriObject
{
public:

    virtual Point3f samplePathStep(const Point3f& x, const Point3f& xz, const 
        Vector3f& dir, const Vector2f& sample, float& pdf_failure, float& pdf_success) const = 0;

    virtual Color3f transmittance(const Point3f x0, const Point3f xz) const = 0;

    virtual Color3f sample_mu_t(const Point3f& p_world) const = 0;

    virtual Color3f sample_mu_a(const Point3f& p_world) const = 0;

    const std::shared_ptr<PhaseFunction> getPhaseFunction() const
    {
        /// WARNING: Might be nullptr, programmer has to check it
        /// TODO: He comprobado y aunque pueda ser nullptr, yo no lo leo mal nunca y siempre está inicializado
        return m_phase_function;
    }

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
    * provided by this instance
    * */
    EClassType getClassType() const { return EVolume; }

protected:
    Color3f                         mu_t;
    Color3f                         mu_a;
    std::shared_ptr<PhaseFunction>  m_phase_function;
};

NORI_NAMESPACE_END