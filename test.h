#pragma once

#include "CubicalComplexRefinement.h"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

void printGridInfo(
    std::pair<const nsCCR::Coordinate, nsCCR::CubicalComplexRefinement::CellData> simplex);
void printNodeRef(nsCCR::Coordinate coord);

// Simple in-memory OctreeSF used to test the leaf-iterator input interface.
class TestOctreeSF final : public nsSF::OctreeSF<float>
{
public:
    using NodeInfo = nsSF::NodeInfo<float>;

    class TestLeafIterator final : public LeafIterator
    {
    public:
        explicit TestLeafIterator(const std::vector<NodeInfo>& leaves);
        bool next(NodeInfo& nodeInfo) override;

    private:
        const std::vector<NodeInfo>& m_leaves;
        size_t m_nextLeaf = 0;
    };

    explicit TestOctreeSF(const std::array<int, 3>& rootSize);
    void addLeaf(const std::array<uint16_t, 3>& location, uint16_t gridSize);
    std::unique_ptr<LeafIterator> createLeafIterator() const override;
    void rootNodeSize(int* size) const override;
    float getValue(const uint16_t* location) const override;

private:
    std::array<int, 3> m_rootSize;
    std::vector<NodeInfo> m_leaves;
};

// Description: 
// Generate cubical complex from 3 dimensional digital images. The image is assumed to be 
// gray style and presented in an array with pixel(x,y,z)=image(x*100+y*10+z).
// The corresponding octree have pixels as vertices of a cube, and all nodes are leaf nodes.
nsCCR::CubicalComplexRefinement generateFrom3dImg(
    const std::string& filename,
    const std::array<size_t, 3>& dimensions = {10, 10, 10},
    uint16_t prime = 2);

nsCCR::CubicalComplexRefinement ccrExample(int n);

void testBuildVFilt();
int isApparentPair(nsCCR::CubicalComplexRefinement* ccr,
                   nsCCR::Coordinate coord1,
                   nsCCR::Coordinate coord2);
int testHaveApparentPair(nsCCR::CubicalComplexRefinement* ccr, 
    nsCCR::CubicalComplexRefinement::CellOrder splx);
void testcomputeApntPair();
void testcomputePersPair();
