#include "exposure_fusion.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split(const std::string &text, char separator) {
  std::vector<std::string> parts;
  size_t start = 0;
  while (start <= text.size()) {
    const size_t end = text.find(separator, start);
    const std::string part = text.substr(start, end == std::string::npos ? std::string::npos : end - start);
    if (!part.empty()) parts.push_back(part);
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return parts;
}

std::string baseName(const std::string &path) {
  return std::filesystem::path(path).stem().string();
}

/** Laplacian coefficients are signed; visualise them around 0.5 grey. */
void saveLaplacianImage(const std::string &path, const Plane &plane, float gain = 2.5f) {
  Image image;
  image.w = plane.w;
  image.h = plane.h;
  image.data.resize(static_cast<size_t>(plane.w) * plane.h * 3);
  for (size_t i = 0; i < plane.data.size(); ++i) {
    const float v = 0.5f + plane.data[i] * gain;
    for (int c = 0; c < 3; ++c) image.data[i * 3 + c] = std::max(0.0f, std::min(1.0f, v));
  }
  saveImage(path, image);
}

void usage() {
  std::cout
      << "Usage:\n"
      << "  exposure_fusion --inputs dark.png,bright.png[,mid.png] [--outdir DIR]\n"
      << "                  [--levels 8] [--sigma 0.2] [--no-pyramids]\n\n"
      << "Outputs: fused.png, naive.png, weights/, measures/, pyramid/, stats.json\n";
}

} // namespace

int main(int argc, char **argv) {
  std::vector<std::string> inputs;
  std::string outdir = "data/output";
  FusionOptions options;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&](const char *name) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << name << "\n";
        std::exit(1);
      }
      return argv[++i];
    };
    if (arg == "--inputs" || arg == "-i") {
      for (const std::string &part : split(next("--inputs"), ',')) inputs.push_back(part);
    } else if (arg == "--outdir" || arg == "-o") {
      outdir = next("--outdir");
    } else if (arg == "--levels") {
      options.levels = std::stoi(next("--levels"));
    } else if (arg == "--sigma") {
      options.wellExposedSigma = std::stof(next("--sigma"));
    } else if (arg == "--weight-blur") {
      options.weightBlur = std::stoi(next("--weight-blur"));
    } else if (arg == "--no-pyramids") {
      options.savePyramids = false;
    } else if (arg == "--help" || arg == "-h") {
      usage();
      return 0;
    } else if (!arg.empty() && arg[0] != '-') {
      inputs.push_back(arg);
    } else {
      std::cerr << "unknown argument: " << arg << "\n";
      return 1;
    }
  }

  if (inputs.size() < 2) {
    usage();
    return 1;
  }

  std::vector<Image> images;
  for (const std::string &path : inputs) {
    std::cout << "loading " << path << "\n";
    images.push_back(loadImage(path));
  }

  const FusionResult result = fuseExposures(images, options);

  std::filesystem::create_directories(outdir);
  saveImage(outdir + "/fused.png", result.fused);
  saveImage(outdir + "/naive.png", result.naive);

  std::filesystem::create_directories(outdir + "/weights");
  std::filesystem::create_directories(outdir + "/measures");
  for (int i = 0; i < static_cast<int>(images.size()); ++i) {
    savePlane(outdir + "/weights/weight_" + std::to_string(i + 1) + ".png", result.weights[i]);
    // measures are scaled for display only
    float cmax = 1e-6f, smax = 1e-6f, emax = 1e-6f;
    for (size_t p = 0; p < result.contrast[i].data.size(); ++p) {
      cmax = std::max(cmax, result.contrast[i].data[p]);
      smax = std::max(smax, result.saturation[i].data[p]);
      emax = std::max(emax, result.wellExposed[i].data[p]);
    }
    savePlane(outdir + "/measures/contrast_" + std::to_string(i + 1) + ".png", result.contrast[i], 1.0f / cmax);
    savePlane(outdir + "/measures/saturation_" + std::to_string(i + 1) + ".png", result.saturation[i], 1.0f / smax);
    savePlane(outdir + "/measures/well_exposed_" + std::to_string(i + 1) + ".png", result.wellExposed[i], 1.0f / emax);
  }

  if (options.savePyramids) {
    std::filesystem::create_directories(outdir + "/pyramid");
    // Gaussian pyramids of the weight maps
    for (int i = 0; i < static_cast<int>(images.size()); ++i) {
      const std::vector<Plane> gp = gaussianPyramid(result.weights[i], result.levels);
      for (int l = 0; l < result.levels; ++l) {
        savePlane(outdir + "/pyramid/weight_" + std::to_string(i + 1) + "_g" + std::to_string(l) + ".png", gp[l]);
      }
    }
    // Laplacian pyramid of the fused image (channel 0 for display)
    Image fusedGrayImage = result.fused;
    const Plane gray = toGray(fusedGrayImage);
    const std::vector<Plane> lp = laplacianPyramid(gray, result.levels);
    for (int l = 0; l < result.levels; ++l) {
      saveLaplacianImage(outdir + "/pyramid/fused_l" + std::to_string(l) + ".png", lp[l]);
    }
  }

  std::ofstream stats(outdir + "/stats.json");
  stats << "{\n";
  stats << "  \"inputs\": [";
  for (size_t i = 0; i < inputs.size(); ++i) {
    stats << "\"" << inputs[i] << "\"" << (i + 1 < inputs.size() ? ", " : "");
  }
  stats << "],\n";
  stats << "  \"width\": " << images[0].w << ",\n";
  stats << "  \"height\": " << images[0].h << ",\n";
  stats << "  \"levels\": " << result.levels << ",\n";
  stats << "  \"milliseconds\": " << result.milliseconds << "\n";
  stats << "}\n";
  stats.close();

  std::cout << "\ninputs : " << inputs.size() << "\n";
  std::cout << "size   : " << images[0].w << "x" << images[0].h << "\n";
  std::cout << "levels : " << result.levels << "\n";
  std::cout << "time   : " << result.milliseconds << " ms\n";
  std::cout << "output : " << outdir << "\n";
  return 0;
}
