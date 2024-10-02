#include "ParticleData.H"

namespace spades::particles {

ParticleContainerInfo::ParticleContainerInfo(std::string basename)
    : m_basename(std::move(basename))
{}

ParticleContainerInfo::~ParticleContainerInfo() = default;

} // namespace spades::particles
