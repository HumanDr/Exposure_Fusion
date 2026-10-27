#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "exposure_fusion.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

Image loadImage(const std::string &path) {
  int w = 0;
  int h = 0;
  int channels = 0;
  unsigned char *pixels = stbi_load(path.c_str(), &w, &h, &channels, 3);
  if (pixels == nullptr) {
    throw std::runtime_error("failed to load image: " + path + " (" + stbi_failure_reason() + ")");
  }
  Image image;
  image.w = w;
  image.h = h;
  image.data.resize(static_cast<size_t>(w) * h * 3);
  for (size_t i = 0; i < image.data.size(); ++i) {
    image.data[i] = pixels[i] / 255.0f;
  }
  stbi_image_free(pixels);
  return image;
}

void saveImage(const std::string &path, const Image &image) {
  std::vector<unsigned char> bytes(image.data.size());
  for (size_t i = 0; i < image.data.size(); ++i) {
    bytes[i] = static_cast<unsigned char>(std::min(1.0f, std::max(0.0f, image.data[i])) * 255.0f + 0.5f);
  }
  std::filesystem::path out(path);
  if (out.has_parent_path()) {
    std::filesystem::create_directories(out.parent_path());
  }
  const int ok = stbi_write_png(path.c_str(), image.w, image.h, 3, bytes.data(), image.w * 3);
  if (ok == 0) {
    throw std::runtime_error("failed to write image: " + path);
  }
}

void savePlane(const std::string &path, const Plane &plane, float scale) {
  std::vector<unsigned char> bytes(static_cast<size_t>(plane.w) * plane.h);
  for (size_t i = 0; i < bytes.size(); ++i) {
    bytes[i] = static_cast<unsigned char>(std::min(1.0f, std::max(0.0f, plane.data[i] * scale)) * 255.0f + 0.5f);
  }
  std::filesystem::path out(path);
  if (out.has_parent_path()) {
    std::filesystem::create_directories(out.parent_path());
  }
  const int ok = stbi_write_png(path.c_str(), plane.w, plane.h, 1, bytes.data(), plane.w);
  if (ok == 0) {
    throw std::runtime_error("failed to write plane: " + path);
  }
}

Plane makePlane(int w, int h, float value) {
  Plane plane;
  plane.w = w;
  plane.h = h;
  plane.data.assign(static_cast<size_t>(w) * h, value);
  return plane;
}

Plane toGray(const Image &image) {
  Plane gray = makePlane(image.w, image.h);
  for (int y = 0; y < image.h; ++y) {
    for (int x = 0; x < image.w; ++x) {
      const float *p = image.px(x, y);
      gray.at(x, y) = 0.2126f * p[0] + 0.7152f * p[1] + 0.0722f * p[2];
    }
  }
  return gray;
}

Plane laplacianMagnitude(const Plane &gray) {
  static const float kernel[9] = {0, 1, 0, 1, -4, 1, 0, 1, 0};
  Plane out = makePlane(gray.w, gray.h);
  for (int y = 1; y < gray.h - 1; ++y) {
    for (int x = 1; x < gray.w - 1; ++x) {
      float sum = 0.0f;
      for (int ky = -1; ky <= 1; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          sum += gray.at(x + kx, y + ky) * kernel[(ky + 1) * 3 + (kx + 1)];
        }
      }
      out.at(x, y) = std::fabs(sum);
    }
  }
  return out;
}
