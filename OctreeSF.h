// Octree scalar field class, used as an input for CCR class
#pragma once

#include <array>
#include <cstdint>
#include <memory>

namespace nsSF
{
template<typename real>
struct NodeInfo
{
    // // conbinatorical location of the lower-left corner
    // std::array<uint16_t, 3> location;
    // uint16_t gridSize;
    // // In x-y-z order, i.e. vertex (x,y,z) have valuevertVal(x*4+y*2+z)
    // std::array<real, 8> vertVal;
    
    //number of subdivision from root node to this node
    int m_numSubDiv;

    // conbinatorical location: got from recursive subdivsion
    // first 8 nodes encedes the combinatorial location of the root node
    // later 8 nodes can be obtained from subdividion information
    // e.g. first node of subdiviison corr. to coordinate (0,0,0) in its depth
    uint16_t m_location[3];
};

template<typename real>
class OctreeSF
{
public:
    class LeafIterator
    {
    public:
        virtual ~LeafIterator() = default;

        // Write the next leaf description to nodeInfo and return true.
        // Return false once every leaf has been visited.
        virtual bool next(NodeInfo<real>& nodeInfo) = 0;
    };

    virtual ~OctreeSF() = default;

    virtual std::unique_ptr<LeafIterator> createLeafIterator() const = 0;

    // Write the three root-grid dimensions to size[0], size[1], and size[2].
    virtual void rootNodeSize(int* size) const = 0;

    // Query the scalar-field value at a point in combinatorial coordinates.
    virtual real getValue(const uint16_t* location) const = 0;
};

}
