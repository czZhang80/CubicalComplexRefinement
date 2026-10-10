// CubicalComplexRefinement.h: 
// class implementing cubical complex refinement data structure

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <queue>
#include <set>
#include <stack>
#include <map>
#include "OctreeSF.h"

namespace nsCCR
{

const int g_maxDim=5;
const int g_maxDepth=16;
const int g_secondSubdivNum=8;
const int g_firstSubdivNum=8;//resolution of initial grid is 2^8
inline constexpr float g_zero = std::numeric_limits<float>::epsilon();
inline constexpr float g_infty = std::numeric_limits<float>::max();

typedef uint16_t Coefficient;
inline bool isPrime(const Coefficient n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (Coefficient p = 3; p <= n / p; p += 2)
		if (!(n % p)) return false;
	return true;
}

inline void multInvrVect(const Coefficient m, std::vector<Coefficient>& inverse) {
	// std::vector<Coefficient> inverse(m);
    inverse.resize(m);
	inverse[1] = 1;
	// m = a * (m / a) + m % a
	// Multipying with inverse(a) * inverse(m % a):
	// 0 = inverse(m % a) * (m / a) + inverse(a)  (mod m)
	for (Coefficient a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
	return;
}


template<class Entry>
class sparseColMat
{
public:
    using matColumn = std::vector<Entry>;

    const matColumn& column(size_t i) const
    {
        return columns.at(i);
    }

    matColumn& column(size_t i)
    {
        return columns.at(i);
    }

    void resize(size_t size)
    {
        columns.resize(size);
    }

    void setColumn(size_t i, matColumn column)
    {
        columns[i] = std::move(column);
    }

    // Column[j] += coef * Column[i] using positional, dense entry arithmetic.
    // Both columns must have the same number of entries; no modulus is applied.
    // void addColumn(Entry coef, size_t i,size_t j)
    // {
    //     const matColumn& source = columns.at(i);
    //     matColumn& destination = columns.at(j);
    //     if(source.size() != destination.size())
    //     {
    //         throw std::invalid_argument("Cannot add matrix columns with different row counts.");
    //     }

    //     for(size_t row = 0; row < source.size(); ++row)
    //     {
    //         destination[row] += coef * source[row];
    //     }
    // }

    // Add or combine one sparse entry in the given column over Z/modulus.
    void addEntry(size_t columnIndex, Entry entry, Coefficient modulus)
    {
        if(!isPrime(modulus))
        {
            throw std::invalid_argument("The coefficient modulus must be prime.");
        }

        entry.coef = static_cast<Coefficient>(entry.coef % modulus);
        if(entry.coef == 0)
        {
            return;
        }

        matColumn& destination = columns.at(columnIndex);
        auto existing = std::find(destination.begin(), destination.end(), entry);
        if(existing == destination.end())
        {
            destination.push_back(std::move(entry));
            return;
        }

        existing->coef = static_cast<Coefficient>(
            (static_cast<uint64_t>(existing->coef) + entry.coef) % modulus);
        if(existing->coef == 0)
        {
            destination.erase(existing);
        }
    }

    // Column[j] += factor * Column[i] over Z/modulus, matching sparse row labels.
    // Only stored entries participate; implicit diagonal entries are not included.
    void addCol(Coefficient factor, size_t i, size_t j, Coefficient modulus)
    {
        if(!isPrime(modulus))
        {
            throw std::invalid_argument("The coefficient modulus must be prime.");
        }

        factor = static_cast<Coefficient>(factor % modulus);
        if(factor == 0)
        {
            return;
        }

        // A copy is needed for a self-add because addEntry may invalidate iterators.
        if(i == j)
        {
            const matColumn source = columns.at(i);
            for(Entry entry : source)
            {
                entry.coef = static_cast<Coefficient>(
                    static_cast<uint64_t>(factor) * entry.coef % modulus);
                addEntry(j, std::move(entry), modulus);
            }
            return;
        }

        if(i >= columns.size() || j >= columns.size())
        {
            throw std::out_of_range("Matrix column index is out of range.");
        }
        const matColumn& source = columns[i];
        for(Entry entry : source)
        {
            entry.coef = static_cast<Coefficient>(
                static_cast<uint64_t>(factor) * entry.coef % modulus);
            addEntry(j, std::move(entry), modulus);
        }
    }

    // Column[i] += factor * colToAdd over Z/modulus, using only stored entries.
    void addCol(Coefficient factor, size_t i, matColumn& colToAdd, Coefficient modulus)
    {
        if(!isPrime(modulus))
        {
            throw std::invalid_argument("The coefficient modulus must be prime.");
        }

        factor = static_cast<Coefficient>(factor % modulus);
        if(factor == 0)
        {
            return;
        }

        // Reuse the snapshot-based self-add if the supplied column is the target.
        if(&colToAdd == &columns.at(i))
        {
            addCol(factor, i, i, modulus);
            return;
        }

        for(Entry entry : colToAdd)
        {
            entry.coef = static_cast<Coefficient>(
                static_cast<uint64_t>(factor) * entry.coef % modulus);
            addEntry(i, std::move(entry), modulus);
        }
    }

    void scaleCol(Coefficient factor, size_t columnIndex, Coefficient modulus)
    {
        if(modulus < 2)
        {
            throw std::invalid_argument("The coefficient modulus must be at least two.");
        }

        matColumn& target = columns.at(columnIndex);
        factor = static_cast<Coefficient>(factor % modulus);
        if(factor == 0)
        {
            target.clear();
            return;
        }

        for(Entry& entry : target)
        {
            entry.coef = static_cast<Coefficient>(
                static_cast<uint64_t>(factor) * entry.coef % modulus);
        }
    }

    Coefficient coefficient(const Entry& row, size_t columnIndex) const
    {
        const matColumn& source = columns.at(columnIndex);
        const auto entry = std::find(source.begin(), source.end(), row);
        return entry == source.end() ? 0 : entry->coef;
    }

    void eraseEntry(const Entry& row, size_t columnIndex)
    {
        matColumn& target = columns.at(columnIndex);
        const auto entry = std::find(target.begin(), target.end(), row);
        if(entry != target.end())
        {
            target.erase(entry);
        }
    }

    void sortCol(size_t columnIndex)
    {
        matColumn& target = columns.at(columnIndex);
        std::sort(target.begin(), target.end());
    }

    // Description:
    // Move row[i] to row[j], and rows between i,j 1 row forward/backward
    void moveRow(size_t orig, size_t term)
    {
        for(const matColumn& column : columns)
        {
            if(orig >= column.size() || term >= column.size())
            {
                throw std::out_of_range("Matrix row index is out of range.");
            }
        }

        if(orig == term)
        {
            return;
        }

        for(matColumn& column : columns)
        {
            if(orig < term)
            {
                std::rotate(column.begin() + orig,
                            column.begin() + orig + 1,
                            column.begin() + term + 1);
            }
            else
            {
                std::rotate(column.begin() + term,
                            column.begin() + orig,
                            column.begin() + orig + 1);
            }
        }
    }

    void moveCol(size_t orig, size_t term)
    {
        if(orig >= columns.size() || term >= columns.size())
        {
            throw std::out_of_range("Matrix column index is out of range.");
        }
        if(orig == term)
        {
            return;
        }

        if(orig < term)
        {
            std::rotate(columns.begin() + orig,
                        columns.begin() + orig + 1,
                        columns.begin() + term + 1);
        }
        else
        {
            std::rotate(columns.begin() + term,
                        columns.begin() + orig,
                        columns.begin() + orig + 1);
        }
    }
private:
    std::vector<matColumn> columns;
};

    class Coordinate 
    {
        // coord: [63:60]=m_dim, [59:56]=gird size, [55:53]=extend direction,
        //        grid size = # of minimum edge length (length of the finest cube), 
        //        the edge length is 2^grid size
        //        [47:32]=x_min coordinates, [31:16]=y_min coordinates, [15:0]=z_min coordinates
        //        each coordinates can be divided to [15:8]=coordinate corr. to first grid, 
        //        [7:0]=coordinate corr. to subdivision grid
        // NOTATION: location is the coordinate of the lower-left corner of the node, it has
        // dimension same as the input mesh, here always assumes to be 3.
        uint64_t m_value;

        constexpr static uint64_t DimMask=0xFULL << 60;
        constexpr static uint64_t GridSizeMask=0xFULL << 56;
        constexpr static uint64_t ExtendDirMask=0x7ull << 53;
        constexpr static uint64_t LocationMask=0xFFFFFFFFFFFFULL;
        constexpr static uint64_t NullValue=0xF000000000000000ULL;
    public:
        // constructor
        Coordinate(uint64_t value) : m_value(value) {};
        Coordinate(int dim=-1,uint16_t gridSize=0, int* extendDir=nullptr, uint16_t *loc=nullptr)
        : m_value(encode(dim,gridSize,extendDir,loc)) {};
        ~Coordinate() = default;
        Coordinate(const Coordinate& other) : m_value(other.m_value) {};

        inline bool isNull() const
        {
            return m_value==NullValue;
        };

        inline uint64_t data() const
        {
            return m_value;
        };

        inline bool operator==(const Coordinate& other) const
        {
            return m_value == other.m_value;
        }

        inline int dim() const
        {
            return static_cast<int>((m_value >> 60) & 0xFu);
        };
        inline bool setDim(int dim)
        {
            if(dim<0||dim>3)
            {
                return false;
            }
            m_value &= ~DimMask;
            m_value |= (static_cast<uint64_t>(dim) & 0xFu) << 60;
            return true;
        }

        // calculate grid size
        inline uint16_t gridSize() const
        {
            int num=static_cast<int>((m_value >> 56) & 0xFu);
            if(num==0)
            {
                return static_cast<uint16_t>(0);
            }
            return static_cast<uint16_t>(1<<(num-1));
        };
        // =-1 if dim()==0, log2(gridSize)-1 otherwise
        inline int subdivNum() const
        {
            return dim()==0?-1:9-static_cast<int>((m_value >> 56) & 0xFu);
        };
        inline uint16_t encodeGridSize(uint16_t x) const
        {
            if(x==0)
            {
                return static_cast<uint16_t>(0);
            }
            uint16_t s=1;
            for(uint16_t g=x; g>1; g>>=1)
            {
                ++s;
            }
            return s;
        };
        inline void setGridSize(uint16_t x)
        {
            m_value &= ~GridSizeMask;
            m_value |= (static_cast<uint64_t>(encodeGridSize(x)) & 0xFu) << 56;
        }
        
        // calculate extend direction
        inline void getExtendDir(int* extendDir) const
        {
            int extendDirBin=(m_value >> 53) & 0x7u;
            extendDir[0]=(extendDirBin>>2)&0b1u;
            extendDir[1]=(extendDirBin>>1)&0b1u;
            extendDir[2]=extendDirBin&0b1u;
        };
        inline void setExtendDir(int* extendDir)
        {
            if(extendDir==nullptr)// clear extendDir
            {
                m_value &= ~ExtendDirMask;
                return;
            }
            uint64_t extendDirBin=static_cast<uint64_t>(extendDir[0]*4+extendDir[1]*2+extendDir[2]);
            m_value &= ~ExtendDirMask;
            m_value |= (static_cast<uint64_t>(extendDirBin) & 0x7u) << 53;
        };

        // location of the lower-left corner of the node
        inline void location(uint16_t *loc) const
        {
            for(int i=0;i<3;i++)
            {
                loc[i] = static_cast<uint16_t>((m_value >> (16 * (3-1-i))) & 0xFFFFu);
            }
        };

        inline void setLocation(uint16_t *loc)
        {
            if(loc==nullptr)// clear location
            {
                m_value &= ~LocationMask;
                return;
            }
            m_value &= ~LocationMask;
            for(int i=0;i<3;i++)
            {
                m_value |= (static_cast<uint64_t>(loc[i]) & 0xFFFFu) << (16 * (3-1-i));
            }
        };
        inline void setLocationI(int i, uint16_t loc)
        {
            m_value &= ~(0xFFFFULL << (16 * (3-1-i)));
            m_value |= (static_cast<uint64_t>(loc) & 0xFFFFULL) << (16 * (3-1-i));
        };
        inline uint64_t encode(int dim=-1,uint16_t gridSize=0, int* extendDir=nullptr, uint16_t *loc=nullptr)
        {
            uint64_t value=0;
            const uint64_t encodedDim = dim == -1 ? 15ULL : static_cast<uint64_t>(dim);
            value |= (encodedDim & 0xFULL) << 60;
            // encode gridSize
            uint16_t s=encodeGridSize(gridSize);
            value |= (static_cast<uint64_t>(s) & 0xFu) << 56;
            if(extendDir!=nullptr)
            {
                uint64_t extendDirBin=static_cast<uint64_t>(extendDir[0]*4+extendDir[1]*2+extendDir[2]);
                value |= (static_cast<uint64_t>(extendDirBin) & 0x7u) << 53;
            }
            if(loc!=nullptr)
            {
                for(int i=0;i<3;i++)
                {
                value |= (static_cast<uint64_t>(loc[i]) & 0xFFFFu) << (16 * (3-1-i));
                }
            }
            return value;
        };

        // get the coordinate of children
        // children are sorted by the lexicographic order of the coordinates,
        // e.g. in dim=3, children[0]={xmin,ymin,zmin}, children[1]={xmin,ymin,zmax}, ...
        // we should work on geo information rather than code directly
        bool getChildrenCoord(std::vector<Coordinate>& children) const
        {
            const int d = dim();
            const uint16_t gs = gridSize();
            if(d == 0 || gs <= 1)
            {
                return false;
            }

            const int numChildren = 1 << d;
            children.resize(numChildren);
            int extendDir[3];
            getExtendDir(extendDir);
            uint16_t loc[3];
            location(loc);
            const uint16_t halfGrid = static_cast<uint16_t>(gs / 2);
            Coordinate coordHeader = m_value;
            coordHeader.setGridSize(halfGrid);
            coordHeader.setLocation(nullptr);
            // children[4] -- {1,0,0} -- {xmax,ymin,zmin}
            for(int idx = 0; idx < numChildren; ++idx)
            {
                children[idx] = coordHeader;
                int activeIdx = 0;
                for(int i = 0; i < 3; ++i)
                {
                    uint16_t childLoc = loc[i];
                    if(extendDir[i] == 1)
                    {
                        if((idx >> (d - 1 - activeIdx)) & 1)
                        {
                            childLoc += halfGrid;
                        }
                        ++activeIdx;
                    }
                    children[idx].setLocationI(i, childLoc);
                }
            }
            return true;
        };

        Coordinate getParentCoord() const
        {
            Coordinate coordParent = m_value;
            uint16_t gs = gridSize();
            uint16_t gridSizeBin = encodeGridSize(gs);
            coordParent.setGridSize(static_cast<uint16_t>(gs * 2));
            uint16_t loc[3];
            location(loc);
            int extendDir[3];
            getExtendDir(extendDir);
            for(int i = 0; i < 3; ++i)
            {
                // Preserve the fixed coordinates of edges and faces.
                if(extendDir[i] == 0) continue;
                uint16_t parentLoc =
                    static_cast<uint16_t>((static_cast<uint32_t>(loc[i]) >> gridSizeBin) << gridSizeBin);
                coordParent.setLocationI(i, parentLoc);
            }
            return coordParent;
        }

        // get idx-th facet coordinate of a simplex. 
        // A simplex with given dim have 2*dim facets, each eliminating a extendDir,
        // and set the corr. coordinate to min/max
        // idx is set as follows: idx set idx/2-th extendDir to 0,
        // and if idx is even, set the corr. coordinate to min; odd to max.
        bool getFacetICoord(int idx, Coordinate& coord)
        {
            // check if idx is valid
            if(idx<0||idx>=2*dim()||dim()==0)
            {
                return false;
            }

            // set the corr. extendDir to 0
            int extendDir[3];
            getExtendDir(extendDir);
            uint16_t loc[3];
            location(loc);
            uint16_t gs=gridSize();

            int count=0;
            for(int i=0;i<3;i++)
            {
                if(extendDir[i]==1)
                {
                    if(count==idx/2)
                    {
                        extendDir[i]=0;
                        if(idx%2==1)
                        {
                            loc[i]+=gs;
                        }
                        break;
                    }
                    else
                    {
                        count++;
                    }
                }    
            }

            // encode the new coordinate
            uint16_t newGridSize=dim()==1?0:gs;
            coord = Coordinate(dim()-1,newGridSize,extendDir,loc);
            return true;
        };

        int getAllFacetsCoord(std::vector<Coordinate>& coords) const
        {
            const int currentDim = dim();
            if(currentDim == 0)
            {
                coords.clear();
                return 0;
            }

            int extendDir[3];
            getExtendDir(extendDir);
            uint16_t loc[3];
            location(loc);
            const uint16_t gs = gridSize();
            const uint16_t facetGridSize = currentDim == 1 ? 0 : gs;

            coords.resize(static_cast<size_t>(2 * currentDim));
            int facetIdx = 0;
            for(int axis = 0; axis < 3; ++axis)
            {
                if(extendDir[axis] == 0)
                {
                    continue;
                }

                extendDir[axis] = 0;
                coords[facetIdx++] = Coordinate(currentDim - 1, facetGridSize, extendDir, loc);

                loc[axis] = static_cast<uint16_t>(loc[axis] + gs);
                coords[facetIdx++] = Coordinate(currentDim - 1, facetGridSize, extendDir, loc);

                loc[axis] = static_cast<uint16_t>(loc[axis] - gs);
                extendDir[axis] = 1;
            }

            return facetIdx;
        };
        // get i-th cofacet coordinate of a simplex.
        // A simplex of dimension dim has 2*(3-dim) cofacets, each adding an active axis,
        // and set the corr. coordinate to loc-gridSize or loc.
        // Return the incidence coefficient (+/-1), or 0 if no valid cofacet exists.
        int getCofacetICoord(int idx, Coordinate& coord)
        {
            const int currentDim = dim();
            if(idx < 0 || currentDim >= 3 || idx >= 2 * (3 - currentDim))
            {
                return 0;
            }

            int extendDir[3];
            getExtendDir(extendDir);
            uint16_t loc[3];
            location(loc);
            const uint16_t gs = currentDim == 0 ? 1 : gridSize();
            constexpr uint32_t maxResolution = 1u << 16;

            int missingAxisRank = 0;
            for(int axis = 0; axis < 3; ++axis)
            {
                if(extendDir[axis] != 0)
                {
                    continue;
                }

                if(missingAxisRank == idx / 2)
                {
                    const int k = axis - missingAxisRank;
                    const int coefficient = ((k + idx) % 2 == 0) ? 1 : -1;
                    extendDir[axis] = 1;
                    if(idx % 2 == 0)
                    {
                        if(loc[axis] < gs)
                        {
                            return 0;
                        }
                        loc[axis] = static_cast<uint16_t>(loc[axis] - gs);
                        coord = Coordinate(currentDim + 1, gs, extendDir, loc);
                        return coefficient;
                    }

                    if(static_cast<uint32_t>(loc[axis]) + gs > maxResolution)
                    {
                        return 0;
                    }

                    coord = Coordinate(currentDim + 1, gs, extendDir, loc);
                    return coefficient;
                }

                ++missingAxisRank;
            }

            return 0;
        };
        int getAllCofacetCoord(std::vector<Coordinate>& coords)
        {
            const int currentDim = dim();
            if(currentDim >= 3)
            {
                coords.clear();
                return 0;
            }

            int extendDir[3];
            getExtendDir(extendDir);
            uint16_t loc[3];
            location(loc);
            const uint16_t gs = currentDim == 0 ? 1 : gridSize();
            constexpr uint32_t maxResolution = 1u << 16;

            coords.clear();
            coords.reserve(static_cast<size_t>(2 * (3 - currentDim)));
            for(int axis = 0; axis < 3; ++axis)
            {
                if(extendDir[axis] != 0)
                {
                    continue;
                }

                extendDir[axis] = 1;

                // Negative-side cofacet: its lower corner is one grid size below this simplex.
                if(loc[axis] >= gs)
                {
                    loc[axis] = static_cast<uint16_t>(loc[axis] - gs);
                    coords.emplace_back(currentDim + 1, gs, extendDir, loc);
                    loc[axis] = static_cast<uint16_t>(loc[axis] + gs);
                }

                // Positive-side cofacet has the same lower corner as this simplex.
                if(static_cast<uint32_t>(loc[axis]) + gs <= maxResolution)
                {
                    coords.emplace_back(currentDim + 1, gs, extendDir, loc);
                }
                extendDir[axis] = 0;
            }

            return static_cast<int>(coords.size());
        };
    };

    class CoordinateHash
    {
    public:
        size_t operator()(const Coordinate& coord) const
        {
            return std::hash<uint64_t>{}(coord.data());
        }
    };



class CubicalCplxSubd
{
public:
    using real=float;

    struct CellData
    {
        Coordinate domiVert;
        Coordinate gradPair;
    };
    struct VertData
    {
        real filtration;
        Coordinate gradPair;
    };
    
    struct CellOrder
    {
        CellOrder()=default;
        CellOrder(real filt,Coordinate vert, Coordinate coor, Coefficient c=0)
        : filtration(filt), domiVert(vert), coord(coor), coef(c){};
        real filtration;
        Coordinate domiVert;
        Coordinate coord;
        Coefficient coef;
        inline bool operator==(const CellOrder& other) const
        {
            return coord == other.coord;// if coordinate are equal, all are equal
        };
        inline bool operator > (const CellOrder b) const
        {
            return filtration>b.filtration||
            (filtration==b.filtration&& domiVert.data()>b.domiVert.data())||
            (filtration==b.filtration && domiVert.data()==b.domiVert.data()&& coord.data()>b.coord.data());

        };
        inline bool operator < (const CellOrder b) const
        {
            return filtration<b.filtration||
            (filtration==b.filtration&& domiVert.data()<b.domiVert.data())||
            (filtration==b.filtration && domiVert.data()==b.domiVert.data()&& coord.data()<b.coord.data());
        };
    };
    class CellOrderHash
    {
    public:
        size_t operator()(const CellOrder& cell) const
        {
            return std::hash<uint64_t>{}(cell.coord.data());
        }

    };
public:
    CubicalCplxSubd() = default;

    CubicalCplxSubd(const nsSF::OctreeSF<real>* scalarField, uint16_t prime=2)
    {
        if(!isPrime(prime))
        {
            throw "input coefficient field is not a prime";
        }
        m_modulus=prime;
        multInvrVect(prime, m_multInvr);
        buildVFilt(scalarField);
    }

    // Description:
    // Build simplex filtration from octree leaves:
    // filt(splx)=max{filt(facet)|facet<splx, dim(facet)=sim(splx)-1}
    // Parameters:
    // scalarField: I, supplies a leaf iterator and point-value queries
    void buildVFilt(const nsSF::OctreeSF<real>* scalarField);

    using ExistFacets = std::array<std::unordered_set<Coordinate, CoordinateHash>, 3>;

    bool insertVertex(const Coordinate& coord, real filtration, ExistFacets& existFacets);
    bool checkAndInst(int dim, const Coordinate& coord, const Coordinate& dominatingVertex, ExistFacets& existFacets);
    Coordinate facetDomiVert(const Coordinate& cell);
    void addNode(const nsSF::NodeInfo<real>& node, const nsSF::OctreeSF<real>* scalarField, ExistFacets& existFacets);

    // Description:
    // just subdivide a given node, do not update its gradient pair and persistent homology
    // we suggest to use it only before ph calculation.  
    // Parameters:
    // coord: I, node to be subdivided, only top-dimension nodes allowed
    bool combSubdiv(Coordinate coord);

//private:
public:
    // part 1: combinatorics

    // Description:
    inline bool contain(Coordinate coord) const
    {
        if(coord.dim()==0)
        {
            return(m_vertFilt.find(coord)!=m_vertFilt.end());
        }
        return(m_splxFilt.find(coord)!=m_splxFilt.end());
    };

    // Description:
    // calculate the i-th facets of simplex. a facet might be the sum of several simplices
    // Parameters: 
    // idx: I, 
    // simplex: I,
    // facet: O,
    // Return:
    // int: number of facets found
    int getFacetI(int idx, Coordinate simplex, std::vector<Coordinate>& facet)
    {
        facet.clear();
        Coordinate directFacet;
        if(!simplex.getFacetICoord(idx, directFacet))
        {
            return -1;
        }

        if(contain(directFacet))
        {
            facet.push_back(directFacet);
            return 1;
        }
        std::stack<Coordinate> facetCand;
        int count=0;
        std::vector<Coordinate> children;
        if(!directFacet.getChildrenCoord(children))
        {
            return -1;
        }
        for(auto child : children)
        {
            facetCand.push(child);
        }
        while (!facetCand.empty())
        {
            Coordinate currentFacet=facetCand.top();
            facetCand.pop();
            if(contain(currentFacet))
            {
                facet.push_back(currentFacet);
                count++;
            }
            else
            {
                if(!currentFacet.getChildrenCoord(children))
                {
                    facet.clear();
                    return -1;
                }
                for(auto child : children)
                {
                    facetCand.push(child);
                }
            }

        }
        return count;
    };
    // Description:
    // calculate the i-th cofacet of simplex. There is only 1 possible cofacet
    // Parameters: 
    // idx: I, 
    // simplex: I,
    // cofacet: O,
    // Return:
    // int: incidence coefficient (+/-1), or 0 for invalid input or no stored cofacet
    int getCofacetI(int idx, Coordinate simplex, Coordinate& cofacet)
    {
        Coordinate cofacetCand;
        const int coefficient = simplex.getCofacetICoord(idx,cofacetCand);
        if(coefficient == 0)
        {
            return 0;
        }
        while(m_splxFilt.find(cofacetCand)==m_splxFilt.end())
        {
            Coordinate parentCofacet = cofacetCand.getParentCoord();
            if(parentCofacet.gridSize() >( 1<< g_secondSubdivNum))
            {
                return 0;
            }
            cofacetCand = parentCofacet;
        }
        cofacet = cofacetCand;
        return coefficient;
    };

    inline CellOrder getCellOrder(const Coordinate coord,const Coefficient coef=0)
    {
        real filtration;
        Coordinate domiVert;
        if(coord.dim()==0)
        {
            filtration=m_vertFilt[coord].filtration;
            domiVert=coord;
        }
        else
        {
            domiVert=m_splxFilt[coord].domiVert;
            filtration=m_vertFilt[domiVert].filtration;
        }
        return CellOrder(filtration,domiVert,coord,coef);
    }

    inline real getFiltration(const Coordinate& coord) const
    {
        if(coord.dim()==0)
        {
            return m_vertFilt.at(coord).filtration;
        }
        const Coordinate& dominatingVertex = m_splxFilt.at(coord).domiVert;
        return m_vertFilt.at(dominatingVertex).filtration;
    }

    inline CellOrder getCellOrder(const Coordinate& coord) const
    {
        const Coordinate owner = coord.dim() == 0 ? coord : m_splxFilt.at(coord).domiVert;
        return CellOrder(m_vertFilt.at(owner).filtration, owner, coord);
    }
    // part 2: persistent homology calculation

    bool oldestCofacet(Coordinate sigma, Coordinate& tao);

    bool youngestFacet(Coordinate tao, Coordinate& sigma);

    // Description:
    // find and set the apparent pair of a given vertex or simplex
    // an apparent pair is (sigma,tao), such that sigma is a facet with dim -1 of tao, and 
    // sigma is the youngest facet of tao; tao is the oldest cofacet of sigma.
    // Parameters:
    // sigma: I, lower dimensional cell in an apparent pair
    // Return:
    // bool: true if an apparent pair is found, false otherwise
    bool findApntPair(Coordinate sigma);

    void setApntPair(Coordinate sigma, Coordinate tao);

    void markCrit(Coordinate coord);

    void computeApntPair();

    
    using ReduceCol=std::priority_queue<CellOrder, std::vector<CellOrder>,std::greater<CellOrder>>;
    // Description:
    // initialize a column Rk=D Vk, similar to init_coboundary_and_get_pivot in Ripser
    CellOrder initCol(Coordinate splx, ReduceCol& Rk);
    
    // Combine duplicate entries and discard zeros; an empty result has coefficient 0.
    inline CellOrder popPivot(ReduceCol& column);
    // Return the combined pivot while retaining it in the column.
    inline CellOrder getPivot(ReduceCol& Rk);

    // Coefficient of alpha in beta's coboundary, modulo m_modulus.
    // Assumes beta is a valid facet (or refined facet piece) of alpha.
    inline Coefficient getCobdCoef(Coordinate beta, Coordinate alpha);
    // Description:
    // Add factor times the original cell and its stored corrections to Vk,
    // and their coboundaries to Rk.
    void addCobd(int index, Coefficient factor, int dim, ReduceCol& Vk, ReduceCol& Rk);

    // Add splx with weight factor to Vk and its weighted coboundary to Rk.
    void addSplxCobd(Coordinate splx, Coefficient factor, int dim, ReduceCol& Vk, ReduceCol& Rk);
    // Description:
    // compute persistent homology of cubical complex refinement
    // we compute the matrix decomposition R=DV, where D is the filtration coboundary matrix of 
    // relative cohomology C^d(K,K.)->C^{d+1}(K,K.)
    // This matrix is given as transpose of the filtration boundary matrix with columns and rows
    // in reverse filtration order. 
    // Parameters: 
    // dim: I, 
    // persPair: O,
    void computePersPair(int dim, std::vector<std::pair<real,real>>& persPair);

    void computePH(std::vector<std::pair<real, real>> &persPair);
    
   
    // part 3: persistent homology update

    void moveLeft(int dim, int j,int i);

    void moveRight(int dim, int i, int j);

    void MoveUp(int dim, int j,int i);
    void moveDown(int dim,int i,int j);

    void subdivide();

//private:
public:
    std::unordered_map<Coordinate, CellData, CoordinateHash> m_splxFilt;
    std::unordered_map<Coordinate, VertData, CoordinateHash> m_vertFilt;
    int m_size[3]={0,0,0};

    Coefficient m_modulus = 2;
    std::vector<Coefficient> m_multInvr;

    bool phFlag=false;
    std::array<std::vector<CellOrder>,4> m_columnToReduce;
    // m_pivotCol[dim][k] is the pivot of m_columnToReduce[dim][k].
    // A null CellOrder represents a zero column.
    std::array<std::vector<CellOrder>,4> m_pivotCol;

    std::array<sparseColMat<CellOrder>,4> m_Vk;
};
}
