#pragma once

#include <cstdint>
#include <string>
#include <vector>

/** 8-bit-ish RGB image stored as float in [0,1] (linear-ish sRGB values). */
struct Image {
  int w = 0;
  int h = 0;
  std::vector<float> data; // w*h*3

  float *px(int x, int y) { return &data[(static_cast<size_t>(y) * w + x) * 3]; }
  const float *px(int x, int y) const { return &data[(static_cast<size_t>(y) * w + x) * 3]; }
};

/** Single-channel float image used for weights and pyramid levels. */
struct Plane {
  int w = 0;
  int h = 0;
  std::vector<float> data;

  float &at(int x, int y) { return data[static_cast<size_t>(y) * w + x]; }
  float at(int x, int y) const { return data[static_cast<size_t>(y) * w + x]; }
};

struct FusionOptions {
  int levels = 8;
  float contrastPower = 1.0f;
  float saturationPower = 1.0f;
  float exposurePower = 1.0f;
  float wellExposedSigma = 0.2f;
  bool savePyramids = true;
  int weightBlur = 6;
  float contrastFloor = 0.02f;
  float saturationFloor = 0.02f;
};

struct FusionResult {
  Image fused;
  Image naive;
  std::vector<Plane> weights;
  std::vector<Plane> contrast;
  std::vector<Plane> saturation;
  std::vector<Plane> wellExposed;
  int levels = 0;
  double milliseconds = 0.0;
};

Image loadImage(const std::string &path);
void saveImage(const std::string &path, const Image &image);
void savePlane(const std::string &path, const Plane &plane, float scale = 1.0f);

Plane makePlane(int w, int h, float value = 0.0f);
Plane toGray(const Image &image);
Plane laplacianMagnitude(const Plane &gray);

std::vector<Plane> gaussianPyramid(const Plane &src, int levels);
std::vector<Plane> laplacianPyramid(const Plane &src, int levels);
Plane reconstructFromLaplacian(const std::vector<Plane> &pyramid);

FusionResult fuseExposures(const std::vector<Image> &images, const FusionOptions &options);
