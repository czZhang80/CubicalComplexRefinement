#include "CubicalComplexRefinement.h"
#include <algorithm>
#include <initializer_list>
#include <iterator>

bool nsCCR::CubicalComplexRefinement::insertVertex(
    const Coordinate& coord, real filtration, ExistFacets& existFacets)
{
    auto result = existFacets[0].emplace(coord);
    if(!result.second)
    {
        return true;
    }

    m_vertFilt.emplace(coord, VertData{filtration});
    return true;
}

bool nsCCR::CubicalComplexRefinement::checkAndInst(
    int dim, const Coordinate& coord, const Coordinate& dominatingVertex,
    ExistFacets& existFacets)
{
    if(dim <= 0 || dim >= static_cast<int>(existFacets.size()) + 1)
    {
        return false;
    }

    if(dim == static_cast<int>(existFacets.size()))
    {
        auto result = m_splxFilt.emplace(coord, CellData{dominatingVertex});
        return result.second;
    }

    auto result = existFacets[dim].emplace(coord);
    if(!result.second)
    {
        return true;
    }

    m_splxFilt.emplace(coord, CellData{dominatingVertex});

    Coordinate parent = coord.getParentCoord();
    for(int depth = coord.subdivNum(); depth > 0; --depth)
    {
        auto parentResult = existFacets[dim].emplace(parent);
        if(!parentResult.second)
        {
            break;
        }
        parent = parent.getParentCoord();
    }

    return true;
}

nsCCR::Coordinate nsCCR::CubicalComplexRefinement::facetDomiVert(
    const Coordinate& cell)
{
    std::vector<Coordinate> facets;
    Coordinate dominatingVertex;
    bool hasCandidate = false;
    const auto lessDomi = [this](const Coordinate& lhs, const Coordinate& rhs)
    {
        const real lhsFilt = m_vertFilt.at(lhs).filtration;
        const real rhsFilt = m_vertFilt.at(rhs).filtration;
        return lhsFilt<rhsFilt || (lhsFilt==rhsFilt && lhs.data() < rhs.data());
        // (lhsFilt-rhsFilt)<g_zero violate strict ordering, and break e.g. std::priority_queue assumptions(suggested by ai)
    };
    for(int idx = 0; idx < 2 * cell.dim(); ++idx)
    {
        const int facetCount = getFacetI(idx, cell, facets);
        if(facetCount <= 0)
        {
            throw std::logic_error("A cell was inserted before all of its facets.");
        }

        for(const Coordinate& facet : facets)
        {
            const Coordinate candidate = facet.dim() == 0? facet : m_splxFilt.at(facet).domiVert;
            if(!hasCandidate)
            {
                dominatingVertex = candidate;
                hasCandidate = true;
                continue;
            }
            dominatingVertex = std::max(dominatingVertex, candidate, lessDomi);
        }
    }
    return dominatingVertex;
}

void nsCCR::CubicalComplexRefinement::addNode(
    const nsSF::NodeInfo<real>& node,
    const nsSF::OctreeSF<real>* scalarField,
    ExistFacets& existFacets)
{
    const uint16_t gs = static_cast<uint16_t>(
        1u << (g_secondSubdivNum - node.m_numSubDiv));
    uint16_t vertLocLF[3] = {
        node.m_location[0], node.m_location[1], node.m_location[2]};
    int extendDir[3] = {1, 1, 1};

    for(int x = 0; x < 2; ++x)
    {
        for(int y = 0; y < 2; ++y)
        {
            for(int z = 0; z < 2; ++z)
            {
                int vertED[3] = {0, 0, 0};
                uint16_t vertLoc[3] = {
                    static_cast<uint16_t>(vertLocLF[0] + x * gs),
                    static_cast<uint16_t>(vertLocLF[1] + y * gs),
                    static_cast<uint16_t>(vertLocLF[2] + z * gs)};
                Coordinate vert(0, 0, vertED, vertLoc);
                insertVertex(vert, scalarField->getValue(vertLoc), existFacets);
            }
        }
    }

    for(int i = 0; i < 3; ++i)
    {
        int edgeED[3] = {0, 0, 0};
        edgeED[i] = 1;
        for(int j = 0; j < 2; ++j)
        {
            for(int k = 0; k < 2; ++k)
            {
                uint16_t vertLoc[3] = {
                    vertLocLF[0], vertLocLF[1], vertLocLF[2]};
                vertLoc[(i + 1) % 3] += j * gs;
                vertLoc[(i + 2) % 3] += k * gs;
                Coordinate edge(1, gs, edgeED, vertLoc);
                checkAndInst(1, edge, facetDomiVert(edge), existFacets);
            }
        }
    }

    for(int i = 0; i < 6; ++i)
    {
        int faceED[3] = {extendDir[0], extendDir[1], extendDir[2]};
        faceED[i / 2] = 0;
        uint16_t vertLoc[3] = {
            vertLocLF[0], vertLocLF[1], vertLocLF[2]};
        vertLoc[i / 2] += gs * (i % 2);
        Coordinate face(2, gs, faceED, vertLoc);
        checkAndInst(2, face, facetDomiVert(face), existFacets);
    }

    Coordinate cube(3, gs, extendDir, vertLocLF);
    checkAndInst(3, cube, facetDomiVert(cube), existFacets);
}

void nsCCR::CubicalComplexRefinement::buildVFilt(const nsSF::OctreeSF<real>* scalarField)
{
    if(scalarField == nullptr)
    {
        throw std::invalid_argument("scalarField must not be null.");
    }

    scalarField->rootNodeSize(m_size);

    // Insert simplices with gridsize from small to large, so either a simplex does not have children,
    // or all its children must exist. 
    // In each step we add a simplex, we need to check if its children already exists
    // step1: sort the simplex from small to large
    // step2: insert all simplices; construct a hash table that stores all possible parents of existing splx

    // step1: 
    std::array<std::vector<nsSF::NodeInfo<real>>, g_secondSubdivNum + 1> leafBuckets;
    size_t leafCount = 0;

    std::unique_ptr<nsSF::OctreeSF<real>::LeafIterator> leafIterator =
        scalarField->createLeafIterator();
    if(!leafIterator)
    {
        throw std::runtime_error("scalarField returned a null leaf iterator.");
    }

    nsSF::NodeInfo<real> nodeInfo;
    while(leafIterator->next(nodeInfo))
    {
        if(nodeInfo.m_numSubDiv < 0 || nodeInfo.m_numSubDiv > g_secondSubdivNum)
        {
            throw std::invalid_argument("An octree leaf exceeds the supported refinement depth.");
        }

        // A larger subdivision count means a smaller leaf.  Reversing the
        // bucket index preserves the existing small-to-large insertion order.
        const int level = g_secondSubdivNum - nodeInfo.m_numSubDiv;
        leafBuckets[level].push_back(nodeInfo);
        ++leafCount;
    }

    m_splxFilt.clear();
    m_vertFilt.clear();
    // In a cubical complex, all simplex corr. to a cube is approximately 2^3=8
    m_splxFilt.reserve(leafCount * 8);
    m_vertFilt.reserve(leafCount * 2);

    // step2
    ExistFacets existFacets;
    for(int i=0;i<g_secondSubdivNum+1;i++)
    {
        for(auto node : leafBuckets[i])
        {
            addNode(node, scalarField, existFacets);
        }
    }
}

bool nsCCR::CubicalComplexRefinement::oldestCofacet(Coordinate sigma, Coordinate &tao)
{
    bool found = false;
    CellOrder oldest;
    for(int i = 0; i < 2 * (3 - sigma.dim()); ++i)
    {
        Coordinate cofacet;
        if(getCofacetI(i, sigma, cofacet) == 0) continue;

        const CellOrder candidate = getCellOrder(cofacet);
        if(!found || candidate < oldest)
        {
            oldest = candidate;
            found = true;
        }
    }
    if(found) tao = oldest.coord;
    return found;
}

bool nsCCR::CubicalComplexRefinement::youngestFacet(Coordinate tao, Coordinate &sigma)
{
    bool found = false;
    CellOrder youngest;
    std::vector<Coordinate> facets;
    for(int i = 0; i < 2 * tao.dim(); ++i)
    {
        if(getFacetI(i, tao, facets) <= 0) return false;
        for(const Coordinate& facet : facets)
        {
            const CellOrder candidate = getCellOrder(facet);
            if(!found || youngest < candidate)
            {
                youngest = candidate;
                found = true;
            }
        }
    }
    if(found) sigma = youngest.coord;
    return found;
}

bool nsCCR::CubicalComplexRefinement::findApntPair(Coordinate sigma)
{
    if(sigma.dim()==0)
    {
        if(!m_vertFilt[sigma].gradPair.isNull())
        {
            return true;
        }
    }
    else
    {
        if(!m_splxFilt[sigma].gradPair.isNull())
        {
            return true;
        }
    }
    Coordinate tao, youngest;
    if(oldestCofacet(sigma, tao) && youngestFacet(tao, youngest) && youngest == sigma)
    {
        setApntPair(sigma, tao);
        return true;
    }
    return false;
}

void nsCCR::CubicalComplexRefinement::setApntPair(Coordinate sigma, Coordinate tao)
{
    const auto taoIt = m_splxFilt.find(tao);
    if(sigma.dim() == 0)
    {
        const auto sigmaIt = m_vertFilt.find(sigma);
        sigmaIt->second.gradPair = taoIt->first;
        taoIt->second.gradPair = sigmaIt->first;
    }
    else
    {
        const auto sigmaIt = m_splxFilt.find(sigma);
        sigmaIt->second.gradPair = taoIt->first;
        taoIt->second.gradPair = sigmaIt->first;
    }
}

void nsCCR::CubicalComplexRefinement::markCrit(Coordinate coord)
{
    m_columnToReduce[coord.dim()].push_back(getCellOrder(coord));
}

void nsCCR::CubicalComplexRefinement::computeApntPair()
{
    for(int i=0;i<m_columnToReduce.size();i++)
    {
        m_columnToReduce[i].clear();
    }
    // Determine every pair before collecting unmatched cells.
    for(const auto& vertEntry : m_vertFilt)
    {
        if(!findApntPair(vertEntry.first))
        {
            markCrit(vertEntry.first);
        }
    }
    for(const auto& splxEntry : m_splxFilt)
    {
        if(splxEntry.first.dim() < 3)
        {
            if(!findApntPair(splxEntry.first))
            {
                continue;
            }
        }
    }
    for(const auto& splxEntry : m_splxFilt)
    {
        if(splxEntry.second.gradPair.isNull())
        {
            markCrit(splxEntry.first);
        }
    }
    for(auto& columns : m_columnToReduce)
    {
        std::sort(columns.begin(), columns.end(), std::greater<CellOrder>());
    }
}


nsCCR::CubicalComplexRefinement::CellOrder nsCCR::CubicalComplexRefinement::initCol(Coordinate splx, ReduceCol& Rk)
{
    int d=splx.dim();
    for(int i=0;i<(3-d)*2;i++)
    {
        Coordinate cofacet;
        const int incidence=getCofacetI(i,splx,cofacet);
        if(incidence==0) continue;
        Coefficient coef=incidence>0?1:(m_modulus-1);
        Rk.push(getCellOrder(cofacet,coef));
    }
    return getPivot(Rk);
}

inline nsCCR::CubicalComplexRefinement::CellOrder nsCCR::CubicalComplexRefinement::popPivot(ReduceCol& column)
{
    while(!column.empty())
    {
        CellOrder pivot=column.top();
        column.pop();
        while(!column.empty() && column.top().coord==pivot.coord)
        {
            pivot.coef=(pivot.coef+column.top().coef)%m_modulus;
            column.pop();
        }
        if(pivot.coef!=0) return pivot;
    }
    return CellOrder(0.0,Coordinate(),Coordinate(),0.0);
}

inline nsCCR::CubicalComplexRefinement::CellOrder nsCCR::CubicalComplexRefinement::getPivot(ReduceCol& Rk)
{
    CellOrder pivot=popPivot(Rk);
    if(pivot.coef!=0) Rk.push(pivot);
    return pivot;
}

inline nsCCR::Coefficient nsCCR::CubicalComplexRefinement::getCobdCoef(Coordinate beta, Coordinate alpha)
{
    int betaAxes[3], alphaAxes[3];
    beta.getExtendDir(betaAxes);
    alpha.getExtendDir(alphaAxes);

    int axis=0, k=0;
    while(betaAxes[axis]==alphaAxes[axis])
    {
        k+=alphaAxes[axis];
        ++axis;
    }

    uint16_t betaLoc[3], alphaLoc[3];
    beta.location(betaLoc);
    alpha.location(alphaLoc);
    const bool lowerSide=betaLoc[axis]==alphaLoc[axis];
    return ((k+lowerSide)%2==0)?1:(m_modulus-1);
}

void nsCCR::CubicalComplexRefinement::addCobd(int index, Coefficient factor, int dim, ReduceCol &Vk, ReduceCol &Rk)
{
    // The original coefficient-1 cell is implicit in the stored reduction column.
    addSplxCobd(m_columnToReduce[dim][index].coord,factor,dim,Vk,Rk);
    for(const CellOrder& simplex : m_Vk[dim].column(index))
    {
        Coefficient coef=simplex.coef*factor%m_modulus;
        addSplxCobd(simplex.coord,coef,dim,Vk,Rk);
    }
}

void nsCCR::CubicalComplexRefinement::addSplxCobd(Coordinate splx, Coefficient factor, int dim, ReduceCol &Vk, ReduceCol &Rk)
{
    Vk.push(getCellOrder(splx,factor));
    for(int i=0;i<2*(3-dim);i++)
    {
        Coordinate cofacet;
        const int incidence=getCofacetI(i,splx,cofacet);
        if(incidence==0) continue;
        Coefficient coef=incidence>0?1:(m_modulus-1);
        coef=coef*factor%m_modulus;
        Rk.push(getCellOrder(cofacet,coef));
    }
}

void nsCCR::CubicalComplexRefinement::computePersPair(int dim, std::vector<std::pair<real, real>> &persPair)
{
    for(int i=0;i<m_columnToReduce[dim].size();i++)
    {
        if(!m_Vk[dim].column(i).empty())// a cleared column cannot have itself as Vk entry only
        {
            continue;
        }
        // m_Vk[dim].appendColumn({});
        ReduceCol Rk,Vk;
        CellOrder pivot=initCol(m_columnToReduce[dim][i].coord,Rk);
        while(true)
        {
            if(!(pivot==CellOrder(0.0,Coordinate(),Coordinate(),0.0)))
            {
                auto pair = m_pivotCol.find(pivot);
                if(pair!=m_pivotCol.end())
                {
                    CellOrder otherPivot=pair->first;
                    uint64_t indexColToAdd=pair->second;
                    Coefficient factor=
                    m_modulus - pivot.coef*m_multInvr[otherPivot.coef] % m_modulus;
                    addCobd(indexColToAdd,factor, dim, Vk,Rk);
                    pivot=getPivot(Rk);
                }
                else if(!m_splxFilt[pivot.coord].gradPair.isNull())
                {
                    Coordinate apparentFacet=m_splxFilt[pivot.coord].gradPair;
                    Coefficient coef=getCobdCoef(apparentFacet, pivot.coord);
                    // Incidence coefficients are +/-1, so coef is its own inverse.
                    Coefficient factor=(m_modulus-pivot.coef*coef%m_modulus)%m_modulus;
                    addSplxCobd(apparentFacet,factor,dim,Vk,Rk);
                    pivot=getPivot(Rk);
                }
                else// persistent pair
                {
                    persPair.push_back({m_columnToReduce[dim][i].filtration,pivot.filtration});
                    m_pivotCol.insert({pivot,i});
                    while(true)
                    {
                        CellOrder e=popPivot(Vk);
                        if(e.coef==0) break;
                        m_Vk[dim].column(i).push_back(e);
                    }

                    // clearing
                    if(dim==m_columnToReduce.size()-1)// max dimension
                    {
                        continue;
                    }
                    auto pivotIt = std::find(m_columnToReduce[dim+1].begin(),m_columnToReduce[dim+1].end(),pivot);
                    int pivotPos=std::distance(m_columnToReduce[dim+1].begin(),pivotIt);
                    popPivot(Rk);
                    while(true)
                    {
                        CellOrder e=popPivot(Rk);
                        if(e.coef==0) break;
                        m_Vk[dim+1].column(pivotPos).push_back(e);
                    }
                    break;
                }
            }
            else
            {
                //Rk==0, a infinite column
                persPair.push_back({m_columnToReduce[dim][i].filtration,g_infty});
                break;
            }
        }
    }
}

void nsCCR::CubicalComplexRefinement::computePH(std::vector<std::pair<real, real>> &persPair)
{
    computeApntPair();
    // init m_Vk
    for(int i=0;i<m_Vk.size();i++)
    {
        m_Vk[i].resize(m_columnToReduce[i].size());
    }
    for(int i=0;i<m_columnToReduce.size();i++)
    {
        computePersPair(i,persPair);
    }
}
