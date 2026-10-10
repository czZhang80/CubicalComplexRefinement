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
    std::pair<const nsCCR::Coordinate, nsCCR::CubicalCplxSubd::CellData> simplex)
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


nsCCR::CubicalCplxSubd generateFrom3dImg(
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
    return nsCCR::CubicalCplxSubd(&scalarField, prime);
}

nsCCR::CubicalCplxSubd ccrExample(int n)
{
    switch(n)
    {
        case 0:
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
            return nsCCR::CubicalCplxSubd(&scalarField, 5);
        }
        case 1:
            return generateFrom3dImg("D:\\projects\\CubicalComplexRefinement\\example\\image_10x10x10_spherical_shell.txt",
            {10,10,10},5);
        case 2:
            return generateFrom3dImg("D:\\projects\\CubicalComplexRefinement\\example\\image_10x10x10.txt",
            {10,10,10},5);
        case 3:
        {
            // 3 x 2 x 2 vertices: two cubes along the x direction, over F_5.
            // Layout: image[x * 4 + y * 2 + z], with z varying fastest.
            // After computePH(), moveLeft(1, 4, 0) takes the else branch
            // at k = 2, then the swap branch at k = 0 (k = 3, 1 are skipped).
            ImageOctreeSF scalarField({3, 2, 2}, {
                 8,  4,  2,  0, // x = 0
                10,  3,  1, 11, // x = 1
                 5,  7,  6,  9  // x = 2
            });
            return nsCCR::CubicalCplxSubd(&scalarField, 5);
        }

        default:
            throw std::invalid_argument("Invalid example number.");
    }
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

    nsCCR::CubicalCplxSubd ccr;
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
int isApparentPair(nsCCR::CubicalCplxSubd* ccr,
                   nsCCR::Coordinate coord1,
                   nsCCR::Coordinate coord2)
{
    if(ccr == nullptr)
    {
        return -1;
    }
    if(!ccr->contain(coord1) || !ccr->contain(coord2))
    {
        return -1;
    }

    using CellOrder = nsCCR::CubicalCplxSubd::CellOrder;
    CellOrder cell1, cell2;
    int returnVal;
    if(coord1.dim() - coord2.dim() == 1)
    {
        cell1 = ccr->getCellOrder(coord2);
        cell2 = ccr->getCellOrder(coord1);
        returnVal = 2;
    }
    else if(coord1.dim() - coord2.dim() == -1)
    {
        cell1 = ccr->getCellOrder(coord1);
        cell2 = ccr->getCellOrder(coord2);
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
        const CellOrder cofacetOrder = ccr->getCellOrder(cofacetI);
        if(cell2 > cofacetOrder)
        {
            return 0;
        }
        if(cell2 == cofacetOrder)
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
            if(cell1 < ccr->getCellOrder(f))
            {
                return 0;
            }
        }
    }
    if(cell1.filtration < cell2.filtration)
    {
        returnVal += 2;
    }
    return returnVal;
}

int testHaveApparentPair(nsCCR::CubicalCplxSubd *ccr, nsCCR::CubicalCplxSubd::CellOrder splx)
{
    if(ccr == nullptr || !ccr->contain(splx.coord))
    {
        return -1;
    }

    const nsCCR::Coordinate& coord = splx.coord;

    // Check whether this cell is the lower-dimensional member of an
    // apparent pair.
    if(coord.dim() < 3)
    {
        for(int i = 0; i < (3 - coord.dim()) * 2; ++i)
        {
            nsCCR::Coordinate cofacet;
            if(ccr->getCofacetI(i, coord, cofacet) != 0)
            {
                const int pairType = isApparentPair(ccr, coord, cofacet);
                if(pairType > 0)
                {
                    return pairType;
                }
            }
        }
    }

    // Check whether this cell is the higher-dimensional member of an
    // apparent pair. A geometric facet may consist of several stored cells.
    for(int i = 0; i < coord.dim() * 2; ++i)
    {
        std::vector<nsCCR::Coordinate> facets;
        if(ccr->getFacetI(i, coord, facets) < 0)
        {
            return -1;
        }

        for(const nsCCR::Coordinate& facet : facets)
        {
            const int pairType = isApparentPair(ccr, coord, facet);
            if(pairType > 0)
            {
                return pairType;
            }
        }
    }

    return 0;
}

void testcomputeApntPair()
{
    nsCCR::CubicalCplxSubd ccr = ccrExample(1);
    ccr.computeApntPair();

    int apparentPairCount = 0, criticalSplxCount = 0;
    // step1: diirect check: every paired cell are apparent cells, and unpaired cells 
    // added into column to reduce
    for(auto vert:ccr.m_vertFilt)
    {
        if(!vert.second.gradPair.isNull())
        {
            const auto partnerIt = ccr.m_splxFilt.find(vert.second.gradPair);
            if(partnerIt == ccr.m_splxFilt.end()||
                partnerIt->second.gradPair.data()!=vert.first.data())
            {
                throw std::range_error("unclosed apparent pair");
            }
            if(isApparentPair(&ccr, vert.first, vert.second.gradPair) <= 0)
            {
                throw std::range_error("Invalid apparent pair");
            }
            ++apparentPairCount;
        }
        else
        {            
            const auto& column = ccr.m_columnToReduce[0];
            if(std::find(column.begin(), column.end(), ccr.getCellOrder(vert.first))
               == column.end())
            {
                throw std::range_error("critical pair not in column to reduce");
            }
            criticalSplxCount++;
        }
    }
    for(auto splx:ccr.m_splxFilt)
    {
        if(!splx.second.gradPair.isNull())
        {
            nsCCR::Coordinate reversePair;
            if(splx.second.gradPair.dim()==0)
            {
                const auto partnerIt = ccr.m_vertFilt.find(splx.second.gradPair);
                if(partnerIt==ccr.m_vertFilt.end())
                {
                    throw std::range_error("invalid coordinate in gradient pair");
                }
                reversePair=partnerIt->second.gradPair;
            }
            else
            {
                const auto partnerIt = ccr.m_splxFilt.find(splx.second.gradPair);
                                if(partnerIt==ccr.m_splxFilt.end())
                {
                    throw std::range_error("invalid coordinate in gradient pair");
                }

                reversePair=partnerIt->second.gradPair;
            }
            if(reversePair.data()!=splx.first.data())
            {
                throw std::range_error("unclosed apparent pair");
            }
            if(isApparentPair(&ccr, splx.first, splx.second.gradPair) <= 0)
            {
                throw std::range_error("Invalid apparent pair");
            }
            ++apparentPairCount;
        }
        else
        {
            const auto& column = ccr.m_columnToReduce[splx.first.dim()];
            if(std::find(column.begin(), column.end(), ccr.getCellOrder(splx.first))
               == column.end())
            {
                throw std::range_error("critical pair not in column to reduce");
            }
            criticalSplxCount++;
        }
    }
    std::cout << "Paired cell count: " << apparentPairCount << std::endl;
    std::cout << "Critical simplex count: " << criticalSplxCount << std::endl;
    // step 2: reverse check: every simplex in column to reduce sre not apparent
    for(int i=0;i<ccr.m_columnToReduce.size();i++)
    {
        auto& currectCol=ccr.m_columnToReduce[i];
        for(auto splx:currectCol)
        {
            const int pairType = testHaveApparentPair(&ccr, splx);
            if(pairType < 0)
            {
                throw std::range_error("Invalid cell in column to reduce");
            }
            if(pairType > 0)
            {
                throw std::range_error("Apparent pair was missed");
            }
        }
    }
}

void testcomputePersPair()
{
    nsCCR::CubicalCplxSubd ccr=ccrExample(2);
    std::vector<std::pair<nsCCR::CubicalCplxSubd::real,
                          nsCCR::CubicalCplxSubd::real>> persPair;
    ccr.computeApntPair();
    for(int i=0;i<ccr.m_Vk.size();i++)
    {
        ccr.m_Vk[i].resize(ccr.m_columnToReduce[i].size());
        ccr.m_pivotCol[i].resize(ccr.m_columnToReduce[i].size());
    }

    for(int i=0;i<ccr.m_columnToReduce.size();i++)
    {
        persPair.clear();
        std::cout << "Persistent pair in dim" << i << std::endl;
        ccr.computePersPair(i,persPair);
    
        for(auto pp : persPair)
        {
            if(pp.first==pp.second)
            {
                continue;
            }
            std::string death=(pp.second-nsCCR::g_infty)>-nsCCR::g_zero?"inf":std::to_string(pp.second);
            std::cout << "[" + std::to_string(pp.first) + "," + death + ")" << std::endl;
        }
    }
}

void testMoveLeft()
{
    using CellOrder = nsCCR::CubicalCplxSubd::CellOrder;
    using nsCCR::Coefficient;
    using nsCCR::Coordinate;

    // Synthetic matrix example, not a geometric cubical complex.
    // The coordinates below are distinct row labels only.
    nsCCR::CubicalCplxSubd ccr;
    ccr.m_modulus = 5;
    nsCCR::multInvrVect(ccr.m_modulus, ccr.m_multInvr);
    const auto makeCell = [](int dim, uint16_t x, float filtration,
                             Coefficient coefficient = 0)
    {
        int extendDir[3] = {1, dim == 2 ? 1 : 0, 0};
        uint16_t location[3] = {x, 0, 0};
        return CellOrder(filtration, Coordinate(0, 0, nullptr, location),
                         Coordinate(dim, 1, extendDir, location), coefficient);
    };

    const CellOrder a = makeCell(1, 0, 4);
    const CellOrder b = makeCell(1, 2, 3);
    const CellOrder c = makeCell(1, 4, 2);
    const CellOrder d = makeCell(1, 6, 1);
    const CellOrder pa = makeCell(2, 0, 10, 4);
    const CellOrder pb = makeCell(2, 2, 40, 1);
    const CellOrder pc = makeCell(2, 4, 30, 2);
    const CellOrder pd = makeCell(2, 6, 20, 3);

    // V = [a, b, c, 2a + 3c + d]; diagonal entries are implicit.
    // A compatible D has columns [4pa, pb, 2pc, 2pa + 4pc + 3pd],
    // giving R = DV = [4pa, pb, 2pc, 3pd] over F_5.
    // Pivot order is pa < pd < pc < pb.
    ccr.m_columnToReduce[1] = {a, b, c, d};
    ccr.m_pivotCol[1] = {pa, pb, pc, pd};
    auto& V = ccr.m_Vk[1];
    V.resize(4);
    CellOrder entry = a;
    entry.coef = 2;
    V.addEntry(3, entry, ccr.m_modulus);
    entry = c;
    entry.coef = 3;
    V.addEntry(3, entry, ccr.m_modulus);

    // k = 2: eliminate; k = 1: skip; k = 0: swap and rescale.
    ccr.moveLeft(1, 3, 0);

    const std::vector<CellOrder> expectedCells = {d, a, b, c};
    const std::array<CellOrder, 4> expectedPivots = {pa, pd, pb, pc};
    const std::array<Coefficient, 4> expectedCoefficients = {2, 4, 1, 2};
    if(ccr.m_columnToReduce[1] != expectedCells)
    {
        throw std::runtime_error("moveLeft: incorrect cell order.");
    }
    // V = [d, a + 3d, b, c]; only the entry (d, column 1) is stored.
    for(size_t column = 0; column < expectedCells.size(); ++column)
    {
        const auto& stored = V.column(column);
        if((column == 1 && (stored.size() != 1 || !(stored[0] == d) ||
                            stored[0].coef != 3)) ||
           (column != 1 && !stored.empty()))
        {
            throw std::runtime_error("moveLeft: incorrect stored V column.");
        }
        const auto& pivot = ccr.m_pivotCol[1][column];
        if(!(pivot == expectedPivots[column]) ||
           pivot.coef != expectedCoefficients[column])
        {
            throw std::runtime_error("moveLeft: incorrect pivot or coefficient.");
        }
    }
    std::cout << "moveLeft four-column example passed." << std::endl;
}

void testMoveRight()
{
    using CellOrder = nsCCR::CubicalCplxSubd::CellOrder;
    using nsCCR::Coefficient;
    using nsCCR::Coordinate;

    // Synthetic matrix example over F_5; coordinates are distinct labels.
    nsCCR::CubicalCplxSubd ccr;
    ccr.m_modulus = 5;
    nsCCR::multInvrVect(ccr.m_modulus, ccr.m_multInvr);
    const auto makeCell = [](int dim, uint16_t x, float filtration,
                             Coefficient coefficient = 0)
    {
        int extendDir[3] = {1, dim == 2 ? 1 : 0, 0};
        uint16_t location[3] = {x, 0, 0};
        return CellOrder(filtration, Coordinate(0, 0, nullptr, location),
                         Coordinate(dim, 1, extendDir, location), coefficient);
    };

    const CellOrder a = makeCell(1, 0, 5);
    const CellOrder b = makeCell(1, 2, 4);
    const CellOrder c = makeCell(1, 4, 3);
    const CellOrder d = makeCell(1, 6, 2);
    const CellOrder e = makeCell(1, 8, 1);
    const CellOrder pa = makeCell(2, 0, 20, 4);
    const CellOrder pb = makeCell(2, 2, 10, 1);
    const CellOrder pc = makeCell(2, 4, 50, 2);
    const CellOrder pd = makeCell(2, 6, 40, 3);
    const CellOrder pe = makeCell(2, 8, 30, 4);

    // V = [a, 2a + b, c, 3a + d, 4a + e].
    // D = [4pa, 2pa + pb, 2pc, 3pa + 3pd, 4pa + 4pe], so
    // R = DV = [4pa, pb, 2pc, 3pd, 4pe], with pb < pa < pe < pd < pc.
    ccr.m_columnToReduce[1] = {a, b, c, d, e};
    ccr.m_pivotCol[1] = {pa, pb, pc, pd, pe};
    auto& V = ccr.m_Vk[1];
    V.resize(5);
    CellOrder entry = a;
    entry.coef = 2;
    V.addEntry(1, entry, ccr.m_modulus);
    entry.coef = 3;
    V.addEntry(3, entry, ccr.m_modulus);
    entry.coef = 4;
    V.addEntry(4, entry, ccr.m_modulus);

    // K = {1, 3, 4}: column 2 is skipped.
    // k = 1: eliminate with donor a, retaining pivot pb.
    // k = 3: eliminate and replace donor with the old column 3a + d.
    // k = 4: factor = -4/3 = 2; target becomes 2d + e.
    // Normalize the final donor by 1/3 = 2 to obtain a + 2d.
    ccr.moveRight(1, 0, 4);

    const std::vector<CellOrder> expectedCells = {b, c, d, e, a};
    const std::array<CellOrder, 5> expectedPivots = {pb, pc, pa, pe, pd};
    const std::array<Coefficient, 5> expectedCoefficients = {1, 2, 3, 4, 1};
    if(ccr.m_columnToReduce[1] != expectedCells)
    {
        throw std::runtime_error("moveRight: incorrect cell order.");
    }
    // V = [b, c, d, 2d + e, a + 2d], with implicit diagonal entries.
    for(size_t column = 0; column < expectedCells.size(); ++column)
    {
        const auto& stored = V.column(column);
        const Coefficient expectedEntry = 2;
        if((column < 3 && !stored.empty()) ||
           (column >= 3 && (stored.size() != 1 || !(stored[0] == d) ||
                             stored[0].coef != expectedEntry)))
        {
            throw std::runtime_error("moveRight: incorrect stored V column.");
        }
        const auto& pivot = ccr.m_pivotCol[1][column];
        if(!(pivot == expectedPivots[column]) ||
           pivot.coef != expectedCoefficients[column])
        {
            throw std::runtime_error("moveRight: incorrect pivot or coefficient.");
        }
    }
    std::cout << "moveRight five-column example passed." << std::endl;
}

namespace
{
using RowMoveMatrix = std::array<std::array<nsCCR::Coefficient, 4>, 4>;

struct RowMoveExample
{
    nsCCR::CubicalCplxSubd ccr;
    std::array<nsCCR::Coordinate, 3> cubes; // A, B, C along the x-axis
};

RowMoveExample makeRowMoveExample()
{
    constexpr uint16_t size = 1u << nsCCR::g_secondSubdivNum;
    TestOctreeSF scalarField({3, 1, 1});
    for(uint16_t x = 0; x < 3; ++x)
    {
        scalarField.addLeaf({static_cast<uint16_t>(x * size), 0, 0}, size);
    }
    RowMoveExample example;
    auto& ccr = example.ccr;
    ccr = nsCCR::CubicalCplxSubd(&scalarField, 5);
    int cubeAxes[3] = {1, 1, 1};
    int faceAxes[3] = {0, 1, 1};
    for(uint16_t x = 0; x < 3; ++x)
    {
        uint16_t location[3] = {static_cast<uint16_t>(x * size), 0, 0};
        example.cubes[x] = nsCCR::Coordinate(3, size, cubeAxes, location);
    }
    // Select four yz faces in descending x order: f3, f2, f1, f0.
    // Their actual coboundaries are C, B-C, A-B, -A.
    for(int x = 3; x >= 0; --x)
    {
        uint16_t location[3] = {static_cast<uint16_t>(x * size), 0, 0};
        ccr.m_columnToReduce[2].push_back(
            ccr.getCellOrder(nsCCR::Coordinate(2, size, faceAxes, location)));
    }
    for(int x = 2; x >= 0; --x)
    {
        ccr.m_columnToReduce[3].push_back(ccr.getCellOrder(example.cubes[x]));
        ccr.m_pivotCol[2].push_back(ccr.getCellOrder(example.cubes[x], 1));
    }
    ccr.m_pivotCol[2].push_back(nsCCR::CubicalCplxSubd::CellOrder{});
    ccr.m_Vk[2].resize(4);
    // V = [f3, f2, f1, f0+f1+f2+f3], so R = [C, B-C, A-B, 0].
    for(size_t column = 0; column < 3; ++column)
    {
        auto entry = ccr.m_columnToReduce[2][column];
        entry.coef = 1;
        ccr.m_Vk[2].addEntry(3, entry, ccr.m_modulus);
    }
    return example;
}

void checkRowMove(RowMoveExample& example, const RowMoveMatrix& expectedV,
                  const std::string& name)
{
    auto& ccr = example.ccr;
    const auto& cells = ccr.m_columnToReduce[2];
    // Independent, explicit D: rows A,B,C; columns f3,f2,f1,f0.
    constexpr int D[3][4] = {{0, 0, 1, 4}, {0, 1, 4, 0}, {1, 4, 0, 0}};
    for(size_t column = 0; column < 4; ++column)
    {
        std::array<nsCCR::Coefficient, 4> fullV{};
        fullV[column] = 1;
        for(const auto& entry : ccr.m_Vk[2].column(column))
        {
            const auto row = std::find(cells.begin(), cells.end(), entry);
            if(row == cells.end() || static_cast<size_t>(row - cells.begin()) >= column)
            {
                throw std::runtime_error(name + ": V is not strictly upper triangular in storage.");
            }
            auto& value = fullV[row - cells.begin()];
            value = static_cast<nsCCR::Coefficient>((value + entry.coef) % 5);
        }
        if(fullV != expectedV[column])
        {
            throw std::runtime_error(name + ": unexpected V column " + std::to_string(column));
        }
        nsCCR::CubicalCplxSubd::CellOrder expectedPivot{};
        for(size_t row = 0; row < 3; ++row)
        {
            int coefficient = 0;
            for(size_t k = 0; k < 4; ++k) coefficient += D[row][k] * fullV[k];
            coefficient %= 5;
            if(coefficient == 0) continue;
            const auto candidate = ccr.getCellOrder(example.cubes[row],
                static_cast<nsCCR::Coefficient>(coefficient));
            if(expectedPivot.coord.isNull() || candidate < expectedPivot) expectedPivot = candidate;
        }
        const auto& actual = ccr.m_pivotCol[2][column];
        if(!(actual == expectedPivot) || actual.coef != expectedPivot.coef)
        {
            throw std::runtime_error(name + ": pivot mismatch in column " +
                std::to_string(column) + " (expected coefficient " +
                std::to_string(expectedPivot.coef) + ", got " + std::to_string(actual.coef) + ").");
        }
        if(!expectedPivot.coord.isNull() && actual.filtration != expectedPivot.filtration)
        {
            throw std::runtime_error(name + ": stale pivot filtration in column " + std::to_string(column));
        }
        for(size_t earlier = 0; earlier < column; ++earlier)
        {
            if(!actual.coord.isNull() && actual == ccr.m_pivotCol[2][earlier])
            {
                throw std::runtime_error(name + ": duplicate nonzero pivot.");
            }
        }
    }
}

constexpr RowMoveMatrix initialRowMoveV = {{{1, 0, 0, 0}, {0, 1, 0, 0},
                                           {0, 0, 1, 0}, {1, 1, 1, 1}}};
}

void testMoveUp()
{
    auto example = makeRowMoveExample();
    checkRowMove(example, initialRowMoveV, "MoveUp initial state");
    auto& ccr = example.ccr;
    // Isolate the row update: change A's filtration, retaining the old row
    // vector [C,B,A], as required by MoveUp's calling convention.
    ccr.m_vertFilt.at(ccr.m_splxFilt.at(example.cubes[0]).domiVert).filtration = 1536;
    // A moves from row 2 to row 0. Column 2 collides first with B, then C:
    // (A-B) + (B-C) + C = A. The zero column must remain zero.
    ccr.MoveUp(3, 2, 0);
    const RowMoveMatrix expectedV = {{{1, 0, 0, 0}, {0, 1, 0, 0},
                                     {1, 1, 1, 0}, {1, 1, 1, 1}}};
    checkRowMove(example, expectedV, "MoveUp");
    std::cout << "MoveUp three-cube example passed." << std::endl;
}

void testMoveDown()
{
    auto example = makeRowMoveExample();
    checkRowMove(example, initialRowMoveV, "moveDown initial state");
    auto& ccr = example.ccr;
    // C moves from row 0 to row 2; retain old row/pivot records during repair.
    ccr.m_vertFilt.at(ccr.m_splxFilt.at(example.cubes[2]).domiVert).filtration = 640;
    // K = {0,1}: eliminate C from column 1, obtaining (B-C) + C = B.
    // Expected R = [C, B, A-B, 0], with pivot coefficients [1,1,1,0].
    ccr.moveDown(3, 0, 2);
    const RowMoveMatrix expectedV = {{{1, 0, 0, 0}, {1, 1, 0, 0},
                                     {0, 0, 1, 0}, {1, 1, 1, 1}}};
    checkRowMove(example, expectedV, "moveDown");
    std::cout << "moveDown three-cube example passed." << std::endl;
}
