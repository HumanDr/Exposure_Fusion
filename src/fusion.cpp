#include "exposure_fusion.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace {

inline int clampInt(int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); }

/** Separable 5-tap binomial blur, clamped borders. */
Plane blur5(const Plane &src) {
  static const float k[5] = {1.0f / 16, 4.0f / 16, 6.0f / 16, 4.0f / 16, 1.0f / 16};
  Plane tmp = makePlane(src.w, src.h);
  Plane out = makePlane(src.w, src.h);
  for (int y = 0; y < src.h; ++y) {
    for (int x = 0; x < src.w; ++x) {
      float sum = 0.0f;
      for (int i = -2; i <= 2; ++i) {
        sum += src.at(clampInt(x + i, 0, src.w - 1), y) * k[i + 2];
      }
      tmp.at(x, y) = sum;
    }
  }
  for (int y = 0; y < src.h; ++y) {
    for (int x = 0; x < src.w; ++x) {
      float sum = 0.0f;
      for (int i = -2; i <= 2; ++i) {
        sum += tmp.at(x, clampInt(y + i, 0, src.h - 1)) * k[i + 2];
      }
      out.at(x, y) = sum;
    }
  }
  return out;
}

Plane downsample2(const Plane &src) {
  const int w = (src.w + 1) / 2;
  const int h = (src.h + 1) / 2;
  Plane out = makePlane(w, h);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      out.at(x, y) = src.at(std::min(src.w - 1, x * 2), std::min(src.h - 1, y * 2));
    }
  }
  return out;
}

Plane upsampleTo(const Plane &src, int w, int h) {
  Plane out = makePlane(w, h);
  const float sx = static_cast<float>(src.w) / w;
  const float sy = static_cast<float>(src.h) / h;
  for (int y = 0; y < h; ++y) {
    const float fy = (y + 0.5f) * sy - 0.5f;
    const int y0 = clampInt(static_cast<int>(std::floor(fy)), 0, src.h - 1);
    const int y1 = clampInt(y0 + 1, 0, src.h - 1);
    const float ty = std::max(0.0f, std::min(1.0f, fy - y0));
    for (int x = 0; x < w; ++x) {
      const float fx = (x + 0.5f) * sx - 0.5f;
      const int x0 = clampInt(static_cast<int>(std::floor(fx)), 0, src.w - 1);
      const int x1 = clampInt(x0 + 1, 0, src.w - 1);
      const float tx = std::max(0.0f, std::min(1.0f, fx - x0));
      const float a = src.at(x0, y0) * (1 - tx) + src.at(x1, y0) * tx;
      const float b = src.at(x0, y1) * (1 - tx) + src.at(x1, y1) * tx;
      out.at(x, y) = a * (1 - ty) + b * ty;
    }
  }
  return out;
}

Plane channelToPlane(const Image &image, int channel) {
  Plane plane = makePlane(image.w, image.h);
  for (int y = 0; y < image.h; ++y) {
    for (int x = 0; x < image.w; ++x) {
      plane.at(x, y) = image.px(x, y)[channel];
    }
  }
  return plane;
}

Image mergePlanes(const std::vector<Plane> &planes) {
  Image image;
  image.w = planes[0].w;
  image.h = planes[0].h;
  image.data.assign(static_cast<size_t>(image.w) * image.h * 3, 0.0f);
  for (int c = 0; c < 3; ++c) {
    for (int y = 0; y < image.h; ++y) {
      for (int x = 0; x < image.w; ++x) {
        image.px(x, y)[c] = std::max(0.0f, std::min(1.0f, planes[c].at(x, y)));
      }
    }
  }
  return image;
}

} // namespace

std::vector<Plane> gaussianPyramid(const Plane &src, int levels) {
  std::vector<Plane> pyramid;
  pyramid.reserve(levels);
  pyramid.push_back(src);
  for (int i = 1; i < levels; ++i) {
    pyramid.push_back(downsample2(blur5(pyramid.back())));
  }
  return pyramid;
}

std::vector<Plane> laplacianPyramid(const Plane &src, int levels) {
  const std::vector<Plane> gaussian = gaussianPyramid(src, levels);
  std::vector<Plane> laplacian(levels);
  for (int i = 0; i < levels; ++i) {
    if (i == levels - 1) {
      laplacian[i] = gaussian[i];
      continue;
    }
    const Plane up = upsampleTo(gaussian[i + 1], gaussian[i].w, gaussian[i].h);
    Plane diff = makePlane(gaussian[i].w, gaussian[i].h);
    for (size_t p = 0; p < diff.data.size(); ++p) {
      diff.data[p] = gaussian[i].data[p] - up.data[p];
    }
    laplacian[i] = std::move(diff);
  }
  return laplacian;
}

Plane reconstructFromLaplacian(const std::vector<Plane> &pyramid) {
  Plane current = pyramid.back();
  for (int i = static_cast<int>(pyramid.size()) - 2; i >= 0; --i) {
    const Plane up = upsampleTo(current, pyramid[i].w, pyramid[i].h);
    Plane sum = makePlane(pyramid[i].w, pyramid[i].h);
    for (size_t p = 0; p < sum.data.size(); ++p) {
      sum.data[p] = pyramid[i].data[p] + up.data[p];
    }
    current = std::move(sum);
  }
  return current;
}

FusionResult fuseExposures(const std::vector<Image> &images, const FusionOptions &options) {
  if (images.empty()) {
    throw std::runtime_error("no input images");
  }
  const auto start = std::chrono::steady_clock::now();
  const int n = static_cast<int>(images.size());
  const int w = images[0].w;
  const int h = images[0].h;
  for (const Image &image : images) {
    if (image.w != w || image.h != h) {
      throw std::runtime_error("all inputs must share the same resolution");
    }
  }

  // ---- quality measures -------------------------------------------------
  std::vector<Plane> contrast(n);
  std::vector<Plane> saturation(n);
  std::vector<Plane> wellExposed(n);
  std::vector<Plane> weights(n);
  for (int i = 0; i < n; ++i) {
    contrast[i] = laplacianMagnitude(toGray(images[i]));
    saturation[i] = makePlane(w, h);
    wellExposed[i] = makePlane(w, h, 1.0f);
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        const float *p = images[i].px(x, y);
        const float mean = (p[0] + p[1] + p[2]) / 3.0f;
        const float var =
            ((p[0] - mean) * (p[0] - mean) + (p[1] - mean) * (p[1] - mean) +
             (p[2] - mean) * (p[2] - mean)) /
            3.0f;
        saturation[i].at(x, y) = std::sqrt(std::max(0.0f, var));
        float well = 1.0f;
        for (int c = 0; c < 3; ++c) {
          const float d = p[c] - 0.5f;
          well *= std::exp(-(d * d) / (2.0f * options.wellExposedSigma * options.wellExposedSigma));
        }
        wellExposed[i].at(x, y) = well;
      }
    }
  }

  // ---- weight maps ------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    weights[i] = makePlane(w, h);
    for (size_t p = 0; p < weights[i].data.size(); ++p) {
      // In flat regions C and S are ~0, so the ratio between exposures would be
      // pure noise. Small floors let well-exposedness decide there instead.
      const float c = std::pow(std::max(options.contrastFloor, contrast[i].data[p]), options.contrastPower);
      const float s = std::pow(std::max(options.saturationFloor, saturation[i].data[p]), options.saturationPower);
      const float e = std::pow(wellExposed[i].data[p], options.exposurePower);
      weights[i].data[p] = c * s * e + 1e-6f;
    }
  }
  // A light blur removes pixel-level flip-flopping between exposures (8-bit
  // noise in flat areas otherwise turns the naive blend into salt-and-pepper).
  for (int i = 0; i < n; ++i) {
    for (int b = 0; b < options.weightBlur; ++b) {
      weights[i] = blur5(weights[i]);
    }
  }

  Plane weightSum = makePlane(w, h, 0.0f);
  for (int i = 0; i < n; ++i) {
    for (size_t p = 0; p < weightSum.data.size(); ++p) {
      weightSum.data[p] += weights[i].data[p];
    }
  }
  for (int i = 0; i < n; ++i) {
    for (size_t p = 0; p < weightSum.data.size(); ++p) {
      weights[i].data[p] /= weightSum.data[p];
    }
  }

  // ---- naive weighted average (seam-prone baseline) ---------------------
  Image naive;
  naive.w = w;
  naive.h = h;
  naive.data.assign(static_cast<size_t>(w) * h * 3, 0.0f);
  for (int i = 0; i < n; ++i) {
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        const float weight = weights[i].at(x, y);
        const float *src = images[i].px(x, y);
        float *dst = naive.px(x, y);
        for (int c = 0; c < 3; ++c) {
          dst[c] += src[c] * weight;
        }
      }
    }
  }

  // ---- multiresolution fusion ------------------------------------------
  int levels = std::max(1, options.levels);
  const int minDim = std::min(w, h);
  while (levels > 1 && (1 << (levels - 1)) >= minDim) {
    --levels;
  }

  std::vector<std::vector<Plane>> weightPyramids(n);
  for (int i = 0; i < n; ++i) {
    weightPyramids[i] = gaussianPyramid(weights[i], levels);
  }

  std::vector<Plane> fusedChannels(3);
  for (int c = 0; c < 3; ++c) {
    std::vector<Plane> fusedPyramid(levels);
    for (int l = 0; l < levels; ++l) {
      const int lw = weightPyramids[0][l].w;
      const int lh = weightPyramids[0][l].h;
      fusedPyramid[l] = makePlane(lw, lh, 0.0f);
    }
    for (int i = 0; i < n; ++i) {
      const std::vector<Plane> lap = laplacianPyramid(channelToPlane(images[i], c), levels);
      for (int l = 0; l < levels; ++l) {
        for (size_t p = 0; p < fusedPyramid[l].data.size(); ++p) {
          fusedPyramid[l].data[p] += weightPyramids[i][l].data[p] * lap[l].data[p];
        }
      }
    }
    fusedChannels[c] = reconstructFromLaplacian(fusedPyramid);
  }

  FusionResult result;
  result.fused = mergePlanes(fusedChannels);
  result.naive = naive;
  result.weights = std::move(weights);
  result.contrast = std::move(contrast);
  result.saturation = std::move(saturation);
  result.wellExposed = std::move(wellExposed);
  result.levels = levels;
  const auto end = std::chrono::steady_clock::now();
  result.milliseconds = std::chrono::duration<double, std::milli>(end - start).count();
  return result;
}
