/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    v1 - Dec 2020
    Copyright (c) 2020 by Adrian Jarabo

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

#include <nori/bitmap.h>
#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfVersion.h>
#include <ImfIO.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>


NORI_NAMESPACE_BEGIN

Bitmap::Bitmap(const std::string &filename) {
    Imf::InputFile file(filename.c_str());
    const Imf::Header &header = file.header();
    const Imf::ChannelList &channels = header.channels();

    Imath::Box2i dw = file.header().dataWindow();
    resize(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);

    cout << "Reading a " << cols() << "x" << rows() << " OpenEXR file from \""
         << filename << "\"" << endl;

    const char *ch_r = nullptr, *ch_g = nullptr, *ch_b = nullptr;
    for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it) {
        std::string name = toLower(it.name());

        if (it.channel().xSampling != 1 || it.channel().ySampling != 1) {
            /* Sub-sampled layers are not supported */
            continue;
        }

        if (!ch_r && (name == "r" || name == "red" ||
                endsWith(name, ".r") || endsWith(name, ".red"))) {
            ch_r = it.name();
        } else if (!ch_g && (name == "g" || name == "green" ||
                endsWith(name, ".g") || endsWith(name, ".green"))) {
            ch_g = it.name();
        } else if (!ch_b && (name == "b" || name == "blue" ||
                endsWith(name, ".b") || endsWith(name, ".blue"))) {
            ch_b = it.name();
        }
    }

    if (!ch_r || !ch_g || !ch_b)
        throw NoriException("This is not a standard RGB OpenEXR file!");

    size_t compStride = sizeof(float),
           pixelStride = 3 * compStride,
           rowStride = pixelStride * cols();

    char *ptr = reinterpret_cast<char *>(data());

    Imf::FrameBuffer frameBuffer;
    frameBuffer.insert(ch_r, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert(ch_g, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert(ch_b, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride));
    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);
}

void Bitmap::saveEXR(const std::string &filename) {
    cout << "Writing a " << cols() << "x" << rows()
         << " OpenEXR file to \"" << filename << "\"" << endl;

    std::string path = filename + ".exr";

    Imf::Header header((int) cols(), (int) rows());
    header.insert("comments", Imf::StringAttribute("Generated by Nori"));

    Imf::ChannelList &channels = header.channels();
    channels.insert("R", Imf::Channel(Imf::FLOAT));
    channels.insert("G", Imf::Channel(Imf::FLOAT));
    channels.insert("B", Imf::Channel(Imf::FLOAT));

    Imf::FrameBuffer frameBuffer;
    size_t compStride = sizeof(float),
           pixelStride = 3 * compStride,
           rowStride = pixelStride * cols();

    char *ptr = reinterpret_cast<char *>(data());
    frameBuffer.insert("R", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert("G", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
    frameBuffer.insert("B", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride));

    Imf::OutputFile file(path.c_str(), header);
    file.setFrameBuffer(frameBuffer);
    file.writePixels((int) rows());
}

void Bitmap::savePNG(const std::string &filename) {
    cout << "Writing a " << cols() << "x" << rows()
         << " PNG file to \"" << filename << "\"" << endl;

    std::string path = filename + ".png";

    uint8_t *rgb8 = new uint8_t[3 * cols() * rows()];
    uint8_t *dst = rgb8;
    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            Color3f tonemapped = coeffRef(i, j).toSRGB();
            dst[0] = (uint8_t) clamp(255.f * tonemapped[0], 0.f, 255.f);
            dst[1] = (uint8_t) clamp(255.f * tonemapped[1], 0.f, 255.f);
            dst[2] = (uint8_t) clamp(255.f * tonemapped[2], 0.f, 255.f);
            dst += 3;
        }
    }
    int ret = stbi_write_png(path.c_str(), (int) cols(), (int) rows(), 3, rgb8, 3 * (int) cols());
    if (ret == 0) {
        cout << "Bitmap::savePNG(): Could not save PNG file \"" << path << "%s\"" << endl;
    }

    delete[] rgb8;
}

Color3f Bitmap::eval(const Point2f& uv) const
{
    float x = (1.f - uv[0]) * cols();;
    float y = (1.f - uv[1]) * rows();

    int ix = x, iy = y;
    float wx = x - ix, wy = y - iy;

    // Only warp suported is repeat
    if (ix >= cols() || ix < 0) ix = ix % cols();
    if (iy >= rows() || iy < 0) iy = iy % rows();

    int ix1 = ix + 1, iy1 = iy + 1;
    // Only warp suported is repeat
    if (ix1 >= cols() || ix1 < 0) ix1 = ix1 % cols();
    if (iy1 >= rows() || iy1 < 0) iy1 = iy1 % rows();

    
    Color3f color = ((1.f - wx) * (1.f - wy)) * coeff(iy, ix) + (wx * (1.f - wy)) * coeff(iy, ix1)
        +((1.f - wx) * wy) * coeff(iy1, ix) + (wx * wy) * coeff(iy1, ix1);

    return color / 255.f;
}


LDRBitmap::LDRBitmap(const std::string& filename)
{
    m_use_linear_interpolation = true;
    int x,y,n;
    unsigned char *im = stbi_load(filename.c_str(), &x, &y, &n, 3);
    
    if(!im)
        throw NoriException("Failed opening Bitmap");

    resize(y, x);

    cout << "Reading a " << cols() << "x" << rows() << " LDR file from \""
        << filename << "\"" << endl;

    memcpy((void*)this->data(), (void*)im, y * x * 3 * sizeof(uint8_t));
    
    stbi_image_free(im);
}

LDRBitmap::LDRBitmap(const std::string& filename, const bool& use_linear_interpolation)
{
    m_use_linear_interpolation = use_linear_interpolation;
    int x,y,n;
    unsigned char *im = stbi_load(filename.c_str(), &x, &y, &n, 3);
    
    if(!im)
        throw NoriException("Failed opening Bitmap");

    resize(y, x);

    cout << "Reading a " << cols() << "x" << rows() << " LDR file from \""
        << filename << "\"" << endl;

    memcpy((void*)this->data(), (void*)im, y * x * 3 * sizeof(uint8_t));
    
    stbi_image_free(im);
}


Color3f LDRBitmap::eval_bilinear_interp(const Point2f& uv) const
{
    float x = (1.f - uv[0]) * cols();
    float y = (1.f - uv[1]) * rows();

    int ix = x, iy = y;
    float wx = x - ix, wy = y - iy;

    // Only warp suported is repeat
    if (ix >= cols() || ix < 0) ix = ix % cols();
    if (iy >= rows() || iy < 0) iy = iy % rows();

    int ix1 = ix + 1, iy1 = iy + 1;
    // Only warp suported is repeat
    if (ix1 >= cols() || ix1 < 0) ix1 = ix1 % cols();
    if (iy1 >= rows() || iy1 < 0) iy1 = iy1 % rows();


    Color3f color = ((1.f - wx) * (1.f - wy)) * coeff(iy, ix).cast<float>() + (wx * (1.f - wy)) * coeff(iy, ix1).cast<float>()
                  + ((1.f - wx) * wy) * coeff(iy1, ix).cast<float>() + (wx * wy) * coeff(iy1, ix1).cast<float>();
}

Color3f LDRBitmap::eval_closest_interp(const Point2f& uv) const
{
    float x = (1.f - uv[0]) * cols();
    float y = (1.f - uv[1]) * rows();

    int ix = x, iy = y;

    // Only warp suported is repeat
    if (ix >= cols() || ix < 0) ix = ix % cols();
    if (iy >= rows() || iy < 0) iy = iy % rows();

    Color3f color = coeff(iy, ix).cast<float>();
    
    return color / 255.f;
}

Color3f LDRBitmap::eval(const Point2f& uv) const
{
    Color3f color(0.f);
    if(m_use_linear_interpolation)
    {
        return eval_bilinear_interp(uv);
    }
    else
    {
        return eval_closest_interp(uv);
    }
    //Maybe an enum can help here instead of a bool, but anyway
}



NORI_NAMESPACE_END
