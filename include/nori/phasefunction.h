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

NORI_NAMESPACE_BEGIN

/**
 * 
 *  //TODO: Me estoy dejando algo???
 *  //TODO: Always assuming ESolidAngle, might change it later?
 * 
 */
struct PFQueryRecord {
    public:
        /// Incident direction (in the global frame)
        Vector3f wi;

        /// Outgoing direction (in the global frame)
        Vector3f wo;

        /// Probability associated to that sample
        float m_pdf;

        /// Create a new record for sampling the PF
        PFQueryRecord(const Vector3f &wi)
            : wi(wi) { }

        /// Create a new record for querying the PF
        PFQueryRecord(const Vector3f &wi,
                const Vector3f &wo)
            : wi(wi), wo(wo), m_pdf(0.f) { }
};


/**
 * \brief Superclass of all Phase Functions
 */
class PhaseFunction : public NoriObject {
public:
    /**
     * \brief Sample the PF and return the importance weight (i.e. the
     * value of the PF * cos(theta_o) divided by the probability density
     * of the sample with respect to solid angles).
     *
     * \param pRec    A PF query record
     * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
     *
     * \return The PF value divided by the probability density of the sample
     *         sample. The returned value also includes the cosine
     *         foreshortening factor associated with the outgoing direction,
     *         when this is appropriate. A zero value means that sampling
     *         failed.
     */
    virtual float sample(PFQueryRecord &pRec, const Point2f &sample) const = 0;

    /**
     * \brief Evaluate the PF for a pair of point and a direction
     * specified in \code pRec
     *
     * \param pRec
     *     A record with detailed information on the PF query
     * \return
     *     The PF value, evaluated for each color channel
     */
    virtual float eval(const PFQueryRecord &pRec) const = 0;

    /**
     * \brief Compute the probability of sampling \c pRec.wo
     * (conditioned on \c pRec.wi).
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param pRec
     *     A record with detailed information on the PF query
     *
     * \return
     *     A probability/density value expressed with respect
     *     to the specified measure
     */

    virtual float pdf(PFQueryRecord &pRec) const = 0;

    virtual const Color3f get_mu_s() const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EPhaseFunction; }

    /**
     * \brief Return whether or not this BRDF is diffuse. This
     * is primarily used by photon mapping to decide whether
     * or not to store photons on a surface
     */
    //No sé si puedo quitar esto o no, la verdad, pero lo dejo por si acaso
    virtual bool isDiffuse() const { return false; }
};

NORI_NAMESPACE_END
