/* 
    Date: 19-1-2023
    Authors: Pedro José Pérez García 
    Comms: Part of our nori extensions for the final assignment
            of the course "Modelling and Simullation of Appearance",
            Master in Robotics, Graphics and Computer Vision (MRGCV),
            Universidad de Zaragoza, course 2022-2023
*/

/// TODO: remove openvdb imports, no va ni a tiros
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>

#include <nori/color.h>
#include <nori/mesh.h>
#include <nori/vector.h>
#include <nori/transform.h>
#include <nori/volumedatabase.h>

#include <fstream>
#include <iostream>
#include <exception>
#include <algorithm>

NORI_NAMESPACE_BEGIN

Volumedatabase::Volumedatabase(const std::string& filename, const Transform& trafo)
{
    m_volfilename = filename;
    // m_volgridname = gridname;
    m_transform = trafo;
    std::string extension = filename.substr(filename.find_last_of(".") + 1);
    if(extension == "vdb" || extension == "VDB")
    {
        std::cout << "\nOpenVDB format is not supported yet" << std::endl;
        std::cout << "Please convert your .vdb files to .vol using the following Mitsuba 2 converter" << std::endl;
        std::cout << "https://github.com/mitsuba-renderer/mitsuba2-vdb-converter \n" << std::endl;
        exit(1);
    }
    else if(extension == "vol" || extension == "VOL")
    {
        loadVOLfile(m_volfilename);          //Initializes data according to the file's content
    }
    else    //Unknown file type, throw a NoriException
    {
        std::cout << "\nUnrecognized volume database file format! Aborting..." << std::endl;
        exit(1);
    }
}

void Volumedatabase::loadVOLfile(const std::string& filename)
{
    VOL_stream.open(filename, std::ifstream::binary);
    if(!VOL_stream.is_open())
    {
        std::cout << "Error trying to read .VOL file! Aborting... " << filename << std::endl;
        exit(1);
    } //else
    
    int offset = 0;

    VOL_stream.seekg(0, VOL_stream.end);
    size_t length = VOL_stream.tellg();
    
    VOL_stream.seekg(0, VOL_stream.beg);

    char header[3];
    uint8_t version, encoding;

    VOL_stream.read(header, sizeof(char) * 3);
    if(header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
    {
        std::cout << "Wrong file format! Aborting..." << std::endl;
        exit(1);
    }
    version = read<uint8_t>(VOL_stream);
    if((int)version != 3)
    {
        std::cout << "Wrong .VOL version! Aborting..." << std::endl;
        exit(1);
    }
    encoding = read<int32_t>(VOL_stream);
    if((int)encoding != 1)
    {
        std::cout << "Unknown encoding! Only '1' (Float32) is supported currently. Aborting..." << std::endl;
        exit(1);
    }
    m_cellsX = read<int32_t>(VOL_stream);
    m_cellsY = read<int32_t>(VOL_stream);
    m_cellsZ = read<int32_t>(VOL_stream);
    m_numChannels = read<int32_t>(VOL_stream);

    m_dataCount = m_cellsX * m_cellsY * m_cellsZ * m_numChannels;

    // This worked properly last time I checked 
    // (it would be so fun to read this in a week after completely breaking the code)
    //
    // std::cout << "VOL: " << header << std::endl;
    // std::cout << "format: " << version << std::endl;
    // std::cout << "Encoding: " << encoding << std::endl;
    // std::cout << "_Cells X: " << m_cellsX << std::endl;
    // std::cout << "_Cells Y: " << m_cellsY << std::endl;
    // std::cout << "_Cells Z: " << m_cellsZ << std::endl;
    // std::cout << "_Number of Channels: " << m_numChannels << std::endl;

    m_bb_xmin = read<float>(VOL_stream);        
    m_bb_ymin = read<float>(VOL_stream);
    m_bb_zmin = read<float>(VOL_stream);
    m_bb_xmax = read<float>(VOL_stream);
    m_bb_ymax = read<float>(VOL_stream);
    m_bb_zmax = read<float>(VOL_stream);
    std:: cout << "BOUNDING BOX: " << m_bb_xmin << " " << m_bb_ymin << " " << m_bb_zmin << " " << m_bb_xmax << " " << m_bb_ymax<< " " << m_bb_zmax << std::endl;

    m_bbox = BoundingBox3f(
            m_transform * Point3f(m_bb_xmin, m_bb_ymin, m_bb_zmin)
        ,   m_transform * Point3f(m_bb_xmax, m_bb_ymax, m_bb_zmax
        ));

    m_bb_xmin = m_bbox.min.x();
    m_bb_ymin = m_bbox.min.y();
    m_bb_zmin = m_bbox.min.z();
    m_bb_xmax = m_bbox.max.x();
    m_bb_ymax = m_bbox.max.y();
    m_bb_zmax = m_bbox.max.z();


    std:: cout << "BOUNDING BOX TRAFO: " << m_bbox.min.x() << " " << m_bbox.min.y() << " " << m_bbox.min.z() << " " << m_bbox.max.x() << " " << m_bbox.max.y() << " " << m_bbox.max.z() << std::endl;

    m_mean = 0.f;
    m_max = -std::numeric_limits<float>::infinity();

    VOL_data = new float[m_dataCount];    //Assign size for the data buffer 

    // Read the float values one by one, for every grid in the volume
    for(int i = 0; i < m_dataCount; i++)
    {
        float val = read<float>(VOL_stream);
        m_mean += (double)val;
        m_max = std::max(m_max, val);
        VOL_data[i] = val;
    }

    m_mean /= (double)m_dataCount;

    std::cout << "Loaded grid volume data with dimensions: (" << m_cellsX << ", " 
    << m_cellsY << ", " << m_cellsZ << ") , mean " << m_mean << " and max " << m_max << std::endl;
    VOL_stream.close();
}

Color3f Volumedatabase::sample_density(const Point3f& pos_world)
{
    Point3f pos_local = Point3f((float)clamp(pos_world.x(), std::min(m_bb_xmin, m_bb_xmax), std::max(m_bb_xmin, m_bb_xmax))
                            , (float)clamp(pos_world.y(), std::min(m_bb_ymin, m_bb_ymax), std::max(m_bb_ymin, m_bb_ymax))
                            , (float)clamp(pos_world.z(), std::min(m_bb_zmin, m_bb_zmax), std::max(m_bb_zmin, m_bb_zmax)));
    // std::cout << "................................................" << std::endl;
    // std::cout << "POINT: " << pos_local.toString() << "\nMIN: " << m_bbox.min.toString() << "\nMAX: " << m_bbox.max.toString() << std::endl;
    // We use (for now) the convention that min = 0,0,0 and max = m_cells{X,Y,Z}
    float t_x = abs(pos_local.x() - m_bb_xmin) / abs(m_bb_xmax - m_bb_xmin);
    float t_y = abs(pos_local.y() - m_bb_ymin) / abs(m_bb_ymax - m_bb_ymin);
    float t_z = abs(pos_local.z() - m_bb_zmin) / abs(m_bb_zmax - m_bb_zmin);
    // std::cout << "T_X: " << t_x << " T_Y: " << t_y << " T_Z: " << t_z << std::endl;
    int idx_x = (int)floor(lerp(t_x, 0, m_cellsX));
    int idx_y = (int)floor(lerp(t_y, 0, m_cellsY));
    int idx_z = (int)floor(lerp(t_z, 0, m_cellsZ));
    // std::cout << "X: " << idx_x << " Y: " << idx_y << " Z: " << idx_z << std::endl;
    // std::cout << "CELLS_X: " << m_cellsX << " CELLS_Y: " << m_cellsY << " CELLS_Z: " << m_cellsZ << std::endl;
    // std::cout << "................................................" << std::endl;
    // Finally, access the data structure
    if(m_numChannels == 3)
    {
        /// TODO: Trilinear interpolation!!!!

        float mu_t_r = VOL_data[((idx_z*m_cellsY + idx_y)*m_cellsX + idx_x)*m_numChannels + 0];
        float mu_t_g = VOL_data[((idx_z*m_cellsY + idx_y)*m_cellsX + idx_x)*m_numChannels + 1];
        float mu_t_b = VOL_data[((idx_z*m_cellsY + idx_y)*m_cellsX + idx_x)*m_numChannels + 2];
        return Color3f(mu_t_r, mu_t_g, mu_t_b);
    }
    // Otherwise we will assume only 1 channel
    if(idx_z > m_cellsZ || idx_x > m_cellsX || idx_y > m_cellsY)
    {
        return Color3f(m_mean);
    }

    float mu_t = VOL_data[((idx_z*m_cellsY + idx_y)*m_cellsX + idx_x)*m_numChannels + 0];
    //std::cout << "mu_t: " << mu_t << std::endl;
    return Color3f(mu_t);
}


// Unlike PBBook, I think I don't need to calculate tMin and tMax
// and my intersections will be from 0 to its.p? (i think)
Point3f Volumedatabase::samplePathStep(const Ray3f& ray, const Intersection& its, Sampler* sampler, const Color3f& mu_s, const Color3f& mu_t, Color3f& _beta, bool& sampledMedium)
{
    float tMax = its.t;
    float t = 0.f;
    while (true) {
        t -= std::log(1.0f - sampler->next1D()) / m_max / (mu_t.sum() / 3.f);
        //std::cout << "t_ratio: " << t << " " << tMax << std::endl;
        if(t >= tMax)
        {
            sampledMedium = false;
            _beta = Color3f(1.f);
            return ray(tMax);
        }

        float density = sample_density(ray(t)).sum() / 3.f;         /// TODO: We assume only 1 channel for now, change it later(?)
        
        // Check if we sample an interaction with the medium
        if(density / m_max > sampler->next1D())
        {
            sampledMedium = true;
            _beta = Color3f(mu_s / mu_t);
            return ray(t);
        }
    }
    return ray(0.f);
}

Color3f Volumedatabase::ratioTracking(Sampler* sampler, const Point3f& x0, const Point3f& xz, const float& mu_t)
{
    float tMax = Vector3f(xz - x0).norm();
    Vector3f dir = (xz - x0).normalized();

    float tr = 1.f;
    float t = 0.f;

    while (true) {
        t -= std::log(1.0f - sampler->next1D()) / m_max / mu_t;
        //std::cout << "t_ratio: " << t << " " << tMax << std::endl;
        if(t >= tMax)
            return Color3f(tr, tr, tr);
        // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        // std::cout << "DIR: " << dir.toString() << std::endl;
        // std::cout << "T: " << t << std::endl;
        // std::cout << "TMAX: " << tMax << std::endl;
        // std::cout << "XZ: " << xz.toString() << std::endl;
        // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        float density = sample_density(x0 + t*dir).sum() / 3.f;         /// TODO: We assume only 1 color for now, change it later(?)
        tr *= (1 - std::max((float)0, density / m_max));
    }
    return Color3f(tr, tr, tr);
}



/// TODO: Below code is unused, for now we convert from vdb to vol with external tools and load the vol directly
// Adapted from https://github.com/mitsuba-renderer/mitsuba2-vdb-converter
/// This work is protected under the BSD 3-Clause License
/// Copyright (c) 2022, Delio Vicini, All rights reserved.
void Volumedatabase::convertVDBtoVOL(const std::string& filename, const std::string& gridname)
{
    std::string outputformat = "binary";
    ///////////////////////////////
    ///                       ||||
    /// TODO: Core dump here! vvvv
    ///////////////////////////////
    openvdb::initialize();
    
    openvdb::io::File file(filename);
    file.open();
    // Loop over all grids in the file and retrieve a shared pointer
    // to the one named "LevelSetSphere".  (This can also be done
    // more simply by calling file.readGrid("LevelSetSphere").)
    openvdb::GridBase::Ptr base_grid;
    for (openvdb::io::File::NameIterator name_iter = file.beginName();
        name_iter != file.endName(); ++name_iter)
    {
        // Read in only the grid we are interested in.
        if (gridname == "" || name_iter.gridName() == gridname) {
            base_grid = file.readGrid(name_iter.gridName());
            if (gridname == "")
                break;
        } else {
            std::cout << "skipping grid " << name_iter.gridName() << std::endl;
        }
    }
    openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);

    auto bbox_dim = grid->evalActiveVoxelDim();
    auto bbox = grid->evalActiveVoxelBoundingBox();

    auto ws_min = grid->indexToWorld(bbox.min());
    auto ws_max = grid->indexToWorld(bbox.max() - openvdb::Vec3R(1, 1, 1));
    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*grid);
    // Compute the value of the grid at fractional coordinates in index space.

    std::vector<float> values;
    for (int k = bbox.min().z(); k < bbox.max().z(); ++k) {
        for (int j = bbox.min().y(); j < bbox.max().y(); ++j) {
            for (int i = bbox.min().x(); i < bbox.max().x(); ++i) {                
                float value = sampler.isSample(openvdb::Vec3R(i, j, k));
                values.push_back(value);
            }
        }
    }
    std::string output_filename = filename.substr(0, filename.find_last_of('.')) + ".vol";
    if (outputformat == "ascii") {
        std::ofstream output_file(output_filename);
        if (!output_file.is_open()) {
            std::cout << "Could not open output file!\n";
            throw NoriException("Could not open .vol output file!");
        }
        for (size_t i = 0; i < values.size(); ++i) {
            if (i > 0) {
                output_file << ", ";
            }
            output_file << std::to_string(values[i]);
        }
        output_file.close();
    } 
    else //if (outputformat == "binary")
    {   
        std::ofstream output_file(output_filename);
        if (!output_file.is_open()) {
            std::cout << "Could not open output file!\n";
            throw NoriException("Could not open .vol output file!");
        }
        output_file << 'V';
        output_file << 'O';
        output_file << 'L';
        uint8_t version = 3;
        writeGeneric(output_file, version);

        writeGeneric(output_file, (int32_t) 1); // type
        writeGeneric(output_file, (int32_t)bbox_dim.x() - 1);
        writeGeneric(output_file, (int32_t)bbox_dim.y() - 1);
        writeGeneric(output_file, (int32_t)bbox_dim.z() - 1);
        writeGeneric(output_file, (int32_t) 1); // #n channels

        float xmin = ws_min.x();
        float ymin = ws_min.y();
        float zmin = ws_min.z();
        float xmax = ws_max.x();
        float ymax = ws_max.y();
        float zmax = ws_max.z();
        writeGeneric(output_file, xmin);
        writeGeneric(output_file, ymin);
        writeGeneric(output_file, zmin);
        writeGeneric(output_file, xmax);
        writeGeneric(output_file, ymax);
        writeGeneric(output_file, zmax);

        for (size_t i = 0; i < values.size(); ++i) {
            writeGeneric(output_file, values[i]);
        }
        output_file.close();
    } 
    file.close();
}

NORI_NAMESPACE_END