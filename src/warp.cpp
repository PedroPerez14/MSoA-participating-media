/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <cmath>
#include <algorithm>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) {
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) 
{
    float r1 = sample.x() * 2.0f - 1.0f;    //from 0..1 to -1..1
    float r2 = sample.y() * 2.0f - 1.0f;
    float r;
    float phi;
    if(r1 > -r2)
    {
        if(r1 > r2)
        {
            r = r1;
            phi = (M_PI / 4.0f) * (r2 / r1);
        }
        else
        {
            r = r2;
            phi = (M_PI / 4.0f) * (2.0f - (r1 / r2));
        }      
    }
    else //if(r1 <= -r2)
    {
        if(r1 < r2)
        {
            r = -r1;
            phi = (M_PI / 4.0f) * (4.0f + (r2 / r1));
        }
        else
        {
            r = -r2;
            phi = (M_PI / 4.0f) * (6.0f - (r1 / r2));
        }
    }
    return Point2f(r * cos(phi), r * sin(phi));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) 
{
    if(p.norm() < (1.0f + Epsilon))
    {
        return 1.0f / M_PI;
    }
    return 0.0f;
}

Point2f Warp::squareToUniformTriangle(const Point2f& sample) 
{
    float r1 = sample.x();
    float r2 = sample.y();

    float alpha = 1 - sqrt(r1);
    float beta = (1 - r2) * sqrt(r1);
    //float gamma = r2 * sqrt(r1);
    //gamma can be recovered with 1 - alpha - beta
    return Point2f(alpha, beta);
}

float Warp::squareToUniformTrianglePdf(const Point2f& p)
{
    float alpha = p.x();
    float beta = p.y();
    float gamma = 1 - alpha - beta;
    //point p is not i triangle if alfa, beta or gamma are not in range (0,1)
    if(!((alpha > 0.0f && alpha < 1.0f) && (beta > 0.0f && beta < 1.0f) && (gamma > 0.0f && gamma < 1.0f)))
    {
        return 0.0f;
    }
    //else{}

    // A = (base * height) / 2; base and height ==1 , so A=1/2
    // Our PDF is 1/A, so 1/1/2, 2/1 == 2.
    return 2.0f;
}


Vector3f Warp::squareToUniformSphere(const Point2f &sample)
{
    float th = 2 * M_PI * sample.x();
    float phi = acos(2 * sample.y() - 1.0f);

    return Point3f(cos(th) * sin(phi), sin(th) * sin(phi), cos(phi));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v)
{
    if(v.norm() > (1.0f - Epsilon) && v.norm() < (1.0f + Epsilon))
    {
        return 1.0 / (4.0f * M_PI);
    }
    //else
    return 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample)
{
    float r1 = sample.x();
    float r2 = sample.y();
    float phi = 2.0f * M_PI * r1;
    //float th = acos(r2);          //Unused

    float x = cos(phi) * sqrt(1.0f - (r2 * r2));
    float y = sin(phi) * sqrt(1.0f - (r2 * r2));
    float z = r2;

    return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v)
{
    if(v.z() > 0.0f && v.norm() > (1.0f - Epsilon) && v.norm() < (1.0f + Epsilon))
    {
        return 1.0f / (2.0f * M_PI);
    }
    //else
    return 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample)
{
    float r1 = sample.x();
    float r2 = sample.y();
    float phi = 2.0f * M_PI * r1;
    //float th = acos(r2);          //Unused

    float x = cos(phi) * sqrt(1.0f - r2);
    float y = sin(phi) * sqrt(1.0f - r2);
    float z = sqrt(r2);

    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v)
{
    if(v.z() > 0.0f && v.norm() > (1.0f - Epsilon) && v.norm() < (1.0f + Epsilon))
    {
        float cos_angle = ((v.dot(Vector3f(0.0f,0.0f,1.0f))) / (v.norm()));
        return cos_angle / M_PI;
    }
    //else
    return 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha)
{
    float phi = 2.0f * M_PI * sample.x();
    float logSample = log(1-sample.y());
    if(isinf(logSample)){ logSample = 0.0f;}
    float tan2th = -alpha * alpha * logSample;
    float cosTheta = 1.0f / sqrt(1.0f + tan2th);
    float sinTheta = sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    Vector3f ret = Vector3f(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    if(ret.z() < 0)
    {
        ret = -ret;
    }
    return ret;
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) 
{
    if(m.z() > 0.0f && m.norm() > (1.0f - Epsilon) && m.norm() < (1.0f + Epsilon))
    {
        float cos_th = m.dot(Vector3f(0,0,1)) / m.norm();
        float sin_th = sqrt(1.0f - (cos_th * cos_th));
        float tan2 = pow(sin_th / cos_th, 2.0f);
        float alpha2 = alpha * alpha;
        return ((exp(-tan2/alpha2)) / (M_PI * alpha2 * pow(cos_th, 3.0f)));
    }
    return 0.0f;
    
}


/// Añadido para el trabajo final
Vector3f Warp::squareToHenyeyGreenstein(const Point2f& sample, float g)
{
    float cos_th = 1.f - 2.f * sample.x();
    if(g != 0.f)
    {
        cos_th = (1.f / (2.f * g)) * (1.f + pow(g, 2.f) - pow(((1 - pow(g, 2.f)) / (1.f - g + 2.f *  g * sample.x())), 2.f));
    }
    float sin_th = std::max(0.f, sqrt(1.f - pow(cos_th, 2.f)));
    float phi = 2.f * M_PI * sample.y();

    return Vector3f(sin_th * cos(phi), sin_th * sin(phi), cos_th);
}

/// Añadido para el trabajo final
float Warp::squareToHenyeyGreensteinPdf(const float &cos_th, float g)
{
    float denom_term = sqrt(pow((1.f + pow(g, 2.f) - (2.f * g) * cos_th), 3.f));
    return (0.25f * INV_PI) * ((1.f - pow(g, 2.f)) / (denom_term));   // (1/4pi) * (... / ...)
}

NORI_NAMESPACE_END
