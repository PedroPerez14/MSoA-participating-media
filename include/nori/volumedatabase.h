/* 
    Date: 19-1-2023
    Authors: Pedro José Pérez García 
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

#pragma once

#include <fstream>
#include <nori/bbox.h>
#include <nori/common.h>
#include <string>

NORI_NAMESPACE_BEGIN

class Volumedatabase {
    public:

        Volumedatabase(const std::string& filename, const Transform& trafo);

        Color3f sample_density(const Point3f& pos_world);

        Color3f ratioTracking(const Point3f& x0, const Point3f& xz, const float& mu_t);

    private:

        //loadVOLdata loads the file into memory
        //readVOLdata stores the volume info in an appropiate way to use it later
        void loadVOLfile(const std::string& filename);

        template <typename T> T read(std::ifstream &f) {
            T v;
            f.read(reinterpret_cast<char *>(&v), sizeof(v));
            return v;
        }

        /// TODO: Unused
        /// Adapted from https://github.com/mitsuba-renderer/mitsuba2-vdb-converter
        /// This work is protected under the BSD 3-Clause License
        /// Copyright (c) 2022, Delio Vicini, All rights reserved.
        void convertVDBtoVOL(const std::string& filename, const std::string& gridname);

        template<typename T>
        void writeGeneric(std::ofstream &f, T data)
        {
            f.write(reinterpret_cast<const char *>(&data), sizeof(data));
        }

        std::string m_volfilename;
        std::string m_volgridname;

        /// Relevant info about the .VOL
        int32_t m_cellsX;
        int32_t m_cellsY;
        int32_t m_cellsZ;
        int32_t m_numChannels;

        Transform m_transform;

        //In local coordinates (IMPORTANT!)
        float m_bb_xmin;
        float m_bb_xmax;
        float m_bb_ymin;
        float m_bb_ymax;
        float m_bb_zmin;
        float m_bb_zmax;
        BoundingBox3f m_bbox;                /// Bounding box of the volume

        double m_mean;
        float m_max;

        std::ifstream VOL_stream;
        int m_dataCount;
        float* VOL_data;
};

NORI_NAMESPACE_END
