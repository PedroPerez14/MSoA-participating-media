/* 
    Date: 2-1-2023
    Authors: Pedro José Pérez García 
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#pragma once

#include <memory>
#include <nori/phasefunction.h>
#include <nori/intersection.h>
#include <nori/object.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass for volumes (TODO: openVDB / 3d perlin noise textures support)
 */
class Volume : public NoriObject
{
public:

    virtual Point3f samplePathStep(const Ray3f& ray, const Intersection& its, Sampler*& sampler, Color3f& _beta, std::shared_ptr<Volume>& nextVolume, std::shared_ptr<Volume>& currentVolume, bool& sampledMedium) const = 0;

    virtual float pdfFail(const Point3f& xz, const float& z, const Vector2f& sample) const = 0;

    virtual Color3f transmittance(Sampler*& sampler, const Point3f& x0, const Point3f& xz) const = 0;

    virtual Color3f sample_mu_t(const Point3f& p_world) const = 0;

    virtual Color3f sample_mu_a(const Point3f& p_world) const = 0;

    const std::shared_ptr<PhaseFunction> getPhaseFunction() const
    {
        /// WARNING: Might be nullptr, programmer has to check it
        /// TODO: He comprobado y aunque pueda ser nullptr, yo no lo leo mal nunca y siempre está inicializado
        return m_phase_function;
    }

    const bool isHeterogeneous() const
    {
        return m_heterogeneous;
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
    bool                            m_heterogeneous = false;
};



/// @brief We might have several volumes through a scene, maybe even overlapping among themselves
///             This record will split the "path" a ray makes in different sections.
///             For example, we have {ray.origin - vol1 - vol1&vol2 - vol1,2,3 - vol1,2 - vol2 - xt} 
///             We will create a vector storing 6 VolumetricSegmentRecords to correctly calcualte and 
///             accumulate transmittance. 
///             Idea from: https://graphics.cg.uni-saarland.de/courses/ris-2019/slides/RIS14_VolumeRendering.pdf
struct VolumetricSegmentRecord {
    public:
        /// Point where the volumetric segment starts
        Point3f x0;

        /// Point where the volumetric segment ends
        Point3f xs;

        /// We choose only one volume for every segment
        ///     Using some importance sampling technique
        std::shared_ptr<Volume> segment_vol;

        /// This makes it have an associated pdf
        float vol_pdf = 1.f;

        VolumetricSegmentRecord() {};

        /// Create a new record
        VolumetricSegmentRecord(const Point3f& x0, const Point3f& xs,
            const std::shared_ptr<Volume>& segment_vol, const float& vol_pdf)
            : x0(x0), xs(xs), segment_vol(segment_vol), vol_pdf(vol_pdf) 
            { }
};

NORI_NAMESPACE_END