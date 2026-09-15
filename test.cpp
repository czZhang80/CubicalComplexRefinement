#include "test.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

void printGridInfo(
    std::pair<const nsCCR::Coordinate, nsCCR::CubicalComplexRefinement::CellData> simplex)
{
    nsCCR::Coordinate coord = simplex.first;

    std::cout << "coordinate: " << std::bitset<64>(coord.data()) << std::endl;
    std::cout << "dimension: " << coord.dim() << std::endl;

    int extendDir[3];
    coord.getExtendDir(extendDir);
    std::cout << "extend direction: ("
              << extendDir[0] << ", " << extendDir[1] << ", " << extendDir[2] << ")"
              << std::endl;

    std::cout << "grid size: " << coord.gridSize() << std::endl;

    uint16_t location[3];
    coord.location(location);
    std::cout << "location: ("
              << location[0] << ", " << location[1] << ", " << location[2] << ")"
              << std::endl;
    //std::cout << "Filtration value " << std::to_string(simplex.second.) << std::endl;
}

void printNodeRef(nsCCR::Coordinate coord)
{
    uint16_t loc[3];
    int extDir[3];
    coord.location(loc);
    coord.getExtendDir(extDir);

    std::cout << "node: "
              << "dimension: " << std::to_string(coord.dim()) + "\n"
              << std::to_string(coord.gridSize()) << std::endl
              << "location: (" + std::to_string(loc[0]) + ", " + std::to_string(loc[1]) + ", "
                     + std::to_string(loc[2]) + ")"
              << std::endl
              << "extend direction: (" + std::to_string(extDir[0]) + ", "
                     + std::to_string(extDir[1]) + ", " + std::to_string(extDir[2]) + ")\n";
}

TestOctreeSF::TestLeafIterator::TestLeafIterator(const std::vector<NodeInfo>& leaves)
    : m_leaves(leaves)
{
}

bool TestOctreeSF::TestLeafIterator::next(NodeInfo& nodeInfo)
{
    if(m_nextLeaf >= m_leaves.size())
    {
        return false;
    }

    nodeInfo = m_leaves[m_nextLeaf++];
    return true;
}

TestOctreeSF::TestOctreeSF(const std::array<int, 3>& rootSize)
    : m_rootSize(rootSize)
{
}

void TestOctreeSF::addLeaf(const std::array<uint16_t, 3>& location, uint16_t gridSize)
{
    NodeInfo leaf{};
    int numSubDiv = 0;
    for(uint16_t currentGridSize = static_cast<uint16_t>(1u << nsCCR::g_secondSubdivNum);
        currentGridSize > gridSize;
        currentGridSize >>= 1)
    {
        ++numSubDiv;
    }
    leaf.m_numSubDiv = numSubDiv;
    leaf.m_location[0] = location[0];
    leaf.m_location[1] = location[1];
    leaf.m_location[2] = location[2];

    m_leaves.push_back(leaf);
}

std::unique_ptr<TestOctreeSF::LeafIterator> TestOctreeSF::createLeafIterator() const
{
    return std::make_unique<TestLeafIterator>(m_leaves);
}

void TestOctreeSF::rootNodeSize(int* size) const
{
    size[0] = m_rootSize[0];
    size[1] = m_rootSize[1];
    size[2] = m_rootSize[2];
}

float TestOctreeSF::getValue(const uint16_t* location) const
{
    return static_cast<float>(location[0] + location[1] + location[2]);
}

class ImageOctreeSF final : public nsSF::OctreeSF<float>
{
public:
    using NodeInfo = nsSF::NodeInfo<float>;

    class ImageLeafIterator final : public LeafIterator
    {
    public:
        explicit ImageLeafIterator(const std::array<size_t, 3>& dimensions)
            : m_dimensions(dimensions)
        {
        }

        bool next(NodeInfo& nodeInfo) override
        {
            if(m_x + 1 >= m_dimensions[0])
            {
                return false;
            }

            constexpr uint16_t gridSize =
                static_cast<uint16_t>(1u << nsCCR::g_secondSubdivNum);
            nodeInfo.m_numSubDiv = 0;
            nodeInfo.m_location[0] = static_cast<uint16_t>(m_x * gridSize);
            nodeInfo.m_location[1] = static_cast<uint16_t>(m_y * gridSize);
            nodeInfo.m_location[2] = static_cast<uint16_t>(m_z * gridSize);

            if(++m_z + 1 >= m_dimensions[2])
            {
                m_z = 0;
                if(++m_y + 1 >= m_dimensions[1])
                {
                    m_y = 0;
                    ++m_x;
                }
            }
            return true;
        }

    private:
        std::array<size_t, 3> m_dimensions;
        size_t m_x = 0;
        size_t m_y = 0;
        size_t m_z = 0;
    };

    ImageOctreeSF(std::array<size_t, 3> dimensions, std::vector<float> image)
        : m_dimensions(dimensions), m_image(std::move(image))
    {
    }

    std::unique_ptr<LeafIterator> createLeafIterator() const override
    {
        return std::make_unique<ImageLeafIterator>(m_dimensions);
    }

    void rootNodeSize(int* size) const override
    {
        for(size_t axis = 0; axis < 3; ++axis)
        {
            size[axis] = static_cast<int>(m_dimensions[axis] - 1);
        }
    }

    float getValue(const uint16_t* location) const override
    {
        constexpr size_t gridSize = size_t{1} << nsCCR::g_secondSubdivNum;
        std::array<size_t, 3> pixel{};
        for(size_t axis = 0; axis < 3; ++axis)
        {
            if(location[axis] % gridSize != 0)
            {
                throw std::out_of_range("Image value requested away from a pixel vertex.");
            }
            pixel[axis] = location[axis] / gridSize;
            if(pixel[axis] >= m_dimensions[axis])
            {
                throw std::out_of_range("Image pixel coordinate is outside the volume.");
            }
        }

        const size_t index =
            pixel[0] * m_dimensions[1] * m_dimensions[2]
            + pixel[1] * m_dimensions[2]
            + pixel[2];
        return m_image[index];
    }

private:
    std::array<size_t, 3> m_dimensions;
    std::vector<float> m_image;
};

std::vector<float> readGrayscaleVolume(const std::string& filename, size_t expectedValues)
{
    std::ifstream input(filename);
    if(!input)
    {
        throw std::runtime_error("Cannot open 3D image file: " + filename);
    }

    std::vector<float> image;
    image.reserve(expectedValues);
    std::string line;
    while(std::getline(input, line))
    {
        const size_t comment = line.find('#');
        if(comment != std::string::npos)
        {
            line.erase(comment);
        }

        std::istringstream values(line);
        int value;
        while(values >> value)
        {
            if(value < 0 || value > 255)
            {
                throw std::runtime_error("Grayscale values must be in the range [0, 255].");
            }
            image.push_back(static_cast<float>(value));
        }
        if(!values.eof())
        {
            throw std::runtime_error("Invalid token in 3D image file: " + filename);
        }
    }

    if(image.size() != expectedValues)
    {
        throw std::runtime_error(
            "Expected " + std::to_string(expectedValues) + " grayscale values, but read "
            + std::to_string(image.size()) + ".");
    }
    return image;
}


nsCCR::CubicalComplexRefinement generateFrom3dImg(
    const std::string& filename,
    const std::array<size_t, 3>& dimensions,
    uint16_t prime)
{
    size_t valueCount = 1;
    constexpr size_t coordinateScale = size_t{1} << nsCCR::g_secondSubdivNum;
    for(size_t dimension : dimensions)
    {
        if(dimension < 2)
        {
            throw std::invalid_argument("Every image dimension must contain at least two pixels.");
        }
        if(dimension - 1 > std::numeric_limits<uint16_t>::max() / coordinateScale)
        {
            throw std::invalid_argument("Image dimensions exceed the coordinate encoding capacity.");
        }
        if(valueCount > std::numeric_limits<size_t>::max() / dimension)
        {
            // protect the reserving vector process
            throw std::invalid_argument("Image dimensions are too large.");
        }
        valueCount *= dimension;
    }

    ImageOctreeSF scalarField(dimensions, readGrayscaleVolume(filename, valueCount));
    return nsCCR::CubicalComplexRefinement(&scalarField, prime);
}

void testBuildVFilt()
{
    constexpr uint16_t rootGridSize = static_cast<uint16_t>(1u << nsCCR::g_secondSubdivNum);
    constexpr uint16_t childGridSize = static_cast<uint16_t>(rootGridSize / 2);

    // Build a 2 x 1 x 1 root grid. The first root cube is subdivided once;
    // the second root cube remains a single leaf.
    TestOctreeSF scalarField({2, 1, 1});
    for(int x = 0; x < 2; ++x)
    {
        for(int y = 0; y < 2; ++y)
        {
            for(int z = 0; z < 2; ++z)
            {
                scalarField.addLeaf({
                    static_cast<uint16_t>(x * childGridSize),
                    static_cast<uint16_t>(y * childGridSize),
                    static_cast<uint16_t>(z * childGridSize)}, childGridSize);
            }
        }
    }
    scalarField.addLeaf({rootGridSize, 0, 0}, rootGridSize);

    nsCCR::CubicalComplexRefinement ccr;
    ccr.buildVFilt(&scalarField);
    std::cout << ccr.m_splxFilt.size();
    for(auto splx : ccr.m_splxFilt)
    {
        printGridInfo(splx);
    }
}

// Return:
// -1: invalid input
//  0: not an apparent pair
//  1/2: zero-persistence apparent pair
//  3/4: nonzero-persistence apparent pair
int isApparentPair(nsCCR::CubicalComplexRefinement* ccr,
                   nsCCR::Coordinate coord1,
                   nsCCR::Coordinate coord2)
{
    if(ccr == nullptr)
    {
        return -1;
    }
    if(ccr->m_splxFilt.find(coord1) == ccr->m_splxFilt.end()
       || ccr->m_splxFilt.find(coord2) == ccr->m_splxFilt.end())
    {
        return -1;
    }

    using CellOrder = nsCCR::CubicalComplexRefinement::CellOrder;
    CellOrder cell1, cell2;
    int returnVal;
    if(coord1.dim() - coord2.dim() == 1)
    {
        cell1 = CellOrder(ccr->m_splxFilt[coord2].filtration, coord2);
        cell2 = CellOrder(ccr->m_splxFilt[coord1].filtration, coord1);
        returnVal = 2;
    }
    else if(coord1.dim() - coord2.dim() == -1)
    {
        cell1 = CellOrder(ccr->m_splxFilt[coord1].filtration, coord1);
        cell2 = CellOrder(ccr->m_splxFilt[coord2].filtration, coord2);
        returnVal = 1;
    }
    else
    {
        return 0;
    }

    bool incidence = false;
    for(int i = 0; i < (3 - cell1.coord.dim()) * 2; ++i)
    {
        nsCCR::Coordinate cofacetI;
        if(ccr->getCofacetI(i, cell1.coord, cofacetI) == 0)
        {
            continue;
        }
        if(cell2 > CellOrder(ccr->m_splxFilt[cofacetI].filtration, cofacetI))
        {
            return 0;
        }
        if(cell2 == CellOrder(ccr->m_splxFilt[cofacetI].filtration, cofacetI))
        {
            incidence = true;
        }
    }
    if(!incidence)
    {
        return 0;
    }
    for(int i = 0; i < cell2.coord.dim() * 2; ++i)
    {
        std::vector<nsCCR::Coordinate> facetI;
        if(ccr->getFacetI(i, cell2.coord, facetI) <= 0)
        {
            throw std::logic_error("Failed to enumerate a facet of a stored cell.");
        }
        for(nsCCR::Coordinate f : facetI)
        {
            if(cell1 < CellOrder(ccr->m_splxFilt[f].filtration, f))
            {
                return 0;
            }
        }
    }
    if(std::abs(cell1.filtration - cell2.filtration) > nsCCR::g_zero)
    {
        returnVal += 2;
    }
    return returnVal;
}

void testProcessLowerStar()
{
    constexpr uint16_t rootGridSize = static_cast<uint16_t>(1u << nsCCR::g_secondSubdivNum);
    constexpr uint16_t childGridSize = static_cast<uint16_t>(rootGridSize / 2);

    TestOctreeSF scalarField({2, 1, 1});
    for(int x = 0; x < 2; ++x)
    {
        for(int y = 0; y < 2; ++y)
        {
            for(int z = 0; z < 2; ++z)
            {
                scalarField.addLeaf({
                    static_cast<uint16_t>(x * childGridSize),
                    static_cast<uint16_t>(y * childGridSize),
                    static_cast<uint16_t>(z * childGridSize)}, childGridSize);
            }
        }
    }
    scalarField.addLeaf({rootGridSize, 0, 0}, rootGridSize);

    nsCCR::CubicalComplexRefinement ccr;
    ccr.buildVFilt(&scalarField);
    ccr.computeGradPair();
    for(auto node : ccr.m_splxFilt)
    {
        if(node.second.gradPair == nullptr)
        {
            printNodeRef(node.first);
            const auto& column = ccr.m_columnToReduce[node.first.dim()];
            const bool contains = std::find(
                column.begin(), column.end(),
                nsCCR::CubicalComplexRefinement::CellOrder(node.second.filtration, node.first))
                != column.end();
            std::cout << (contains ? "true" : "false") << std::endl;
        }
        else if(isApparentPair(&ccr, node.first, *node.second.gradPair) <= 0)
        {
            throw std::range_error("wrong");
        }
    }
}

void testcomputePersPair()
{
    // constexpr uint16_t rootGridSize = static_cast<uint16_t>(1u << nsCCR::g_secondSubdivNum);
    // constexpr uint16_t childGridSize = static_cast<uint16_t>(rootGridSize / 2);

    // TestOctreeSF scalarField({2, 1, 1});
    // for(int x = 0; x < 2; ++x)
    // {
    //     for(int y = 0; y < 2; ++y)
    //     {
    //         for(int z = 0; z < 2; ++z)
    //         {
    //             scalarField.addLeaf({
    //                 static_cast<uint16_t>(x * childGridSize),
    //                 static_cast<uint16_t>(y * childGridSize),
    //                 static_cast<uint16_t>(z * childGridSize)}, childGridSize);
    //         }
    //     }
    // }
    // scalarField.addLeaf({rootGridSize, 0, 0}, rootGridSize);

    

    nsCCR::CubicalComplexRefinement ccr=
    generateFrom3dImg("D:\\projects\\Cubical Complex Refinement\\CubicalComplexRefinement\\example\\image_10x10x10_spherical_shell.txt",
    {10,10,10},5);
    std::vector<std::pair<nsCCR::CubicalComplexRefinement::CellOrder,
                          nsCCR::CubicalComplexRefinement::CellOrder>> persPair;
    ccr.computeGradPair();
    for(int i = 0; i < 3; ++i)
    {
        persPair.clear();
        std::cout << "Persistent pair in dim " + std::to_string(i) << std::endl;
        ccr.computePersPair(i, persPair);
        for(auto pp : persPair)
        {
            std::string death=(pp.second.filtration-nsCCR::g_infty)>-nsCCR::g_zero?"inf":std::to_string(pp.second.filtration);

            std::cout << "[" + std::to_string(pp.first.filtration) + ","
                             + death + ")"
                      << std::endl;
        }
    }
}
