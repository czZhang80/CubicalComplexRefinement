#include "CubicalComplexRefinement.h"
#include <algorithm>
#include <initializer_list>
#include <iterator>

bool nsCCR::CubicalCplxSubd::insertVertex(
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

bool nsCCR::CubicalCplxSubd::checkAndInst(
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

nsCCR::Coordinate nsCCR::CubicalCplxSubd::facetDomiVert(
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

void nsCCR::CubicalCplxSubd::addNode(
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

void nsCCR::CubicalCplxSubd::buildVFilt(const nsSF::OctreeSF<real>* scalarField)
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

bool nsCCR::CubicalCplxSubd::oldestCofacet(Coordinate sigma, Coordinate &tao)
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

bool nsCCR::CubicalCplxSubd::youngestFacet(Coordinate tao, Coordinate &sigma)
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

bool nsCCR::CubicalCplxSubd::findApntPair(Coordinate sigma)
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

void nsCCR::CubicalCplxSubd::setApntPair(Coordinate sigma, Coordinate tao)
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

void nsCCR::CubicalCplxSubd::markCrit(Coordinate coord)
{
    m_columnToReduce[coord.dim()].push_back(getCellOrder(coord));
}

void nsCCR::CubicalCplxSubd::computeApntPair()
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


nsCCR::CubicalCplxSubd::CellOrder nsCCR::CubicalCplxSubd::initCol(Coordinate splx, ReduceCol& Rk)
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

inline nsCCR::CubicalCplxSubd::CellOrder nsCCR::CubicalCplxSubd::popPivot(ReduceCol& column)
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

inline nsCCR::CubicalCplxSubd::CellOrder nsCCR::CubicalCplxSubd::getPivot(ReduceCol& Rk)
{
    CellOrder pivot=popPivot(Rk);
    if(pivot.coef!=0) Rk.push(pivot);
    return pivot;
}

inline nsCCR::Coefficient nsCCR::CubicalCplxSubd::getCobdCoef(Coordinate beta, Coordinate alpha)
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

void nsCCR::CubicalCplxSubd::addCobd(int index, Coefficient factor, int dim, ReduceCol &Vk, ReduceCol &Rk)
{
    // The original coefficient-1 cell is implicit in the stored reduction column.
    addSplxCobd(m_columnToReduce[dim][index].coord,factor,dim,Vk,Rk);
    for(const CellOrder& simplex : m_Vk[dim].column(index))
    {
        Coefficient coef=simplex.coef*factor%m_modulus;
        addSplxCobd(simplex.coord,coef,dim,Vk,Rk);
    }
}

void nsCCR::CubicalCplxSubd::addSplxCobd(Coordinate splx, Coefficient factor, int dim, ReduceCol &Vk, ReduceCol &Rk)
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

void nsCCR::CubicalCplxSubd::computePersPair(int dim, std::vector<std::pair<real, real>> &persPair)
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
                const auto pivotColumn = std::find(
                    m_pivotCol[dim].begin(), m_pivotCol[dim].end(), pivot);
                if(pivotColumn != m_pivotCol[dim].end())
                {
                    const uint64_t indexColToAdd = static_cast<uint64_t>(
                        std::distance(m_pivotCol[dim].begin(), pivotColumn));
                    const CellOrder& otherPivot=*pivotColumn;
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
                    m_pivotCol[dim][i]=pivot;
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

void nsCCR::CubicalCplxSubd::computePH(std::vector<std::pair<real, real>> &persPair)
{
    computeApntPair();
    // Initialize V and the column-aligned pivot vectors.
    for(int i=0;i<m_Vk.size();i++)
    {
        m_Vk[i].resize(m_columnToReduce[i].size());
        m_pivotCol[i].assign(m_columnToReduce[i].size(), CellOrder());
    }
    for(int i=0;i<m_columnToReduce.size();i++)
    {
        computePersPair(i,persPair);
    }
}

void nsCCR::CubicalCplxSubd::moveLeft(int dim, int j, int i)
{
    if(dim < 0 || dim >= static_cast<int>(m_Vk.size()) || i < 0 || i >= j ||
       j >= static_cast<int>(m_columnToReduce[dim].size()))
    {
        return;
    }

    auto& cells = m_columnToReduce[dim];
    auto& pivots = m_pivotCol[dim];

    for(int k = j - 1; k >= i; --k)
    {
        // Since k < j, V[k,j] is an explicitly stored off-diagonal entry.
        const Coefficient coefficient = m_Vk[dim].coefficient(cells[k],j);
        if(coefficient == 0) continue;

        const Coefficient factor = static_cast<Coefficient>(m_modulus - coefficient);


        //Check if pivot(R[k])>pivot(R[j])
        if(!pivots[k].coord.isNull() && (pivots[j].coord.isNull() || pivots[k] < pivots[j]))
        {
            // store donor column
            sparseColMat<CellOrder>::matColumn donorCol = m_Vk[dim].column(j);
            CellOrder diagJ=cells[j];
            diagJ.coef=1;

            // V[j] += factor * V[k].
            m_Vk[dim].addCol(factor, k, j, m_modulus);
            // addRow sees only the stored part of V[k]; account for its implicit V[k,k]=1.
            m_Vk[dim].eraseEntry(cells[k],j);

            // V[k]=donor column
            m_Vk[dim].setColumn(k, std::move(donorCol));
            m_Vk[dim].addEntry(k, diagJ, m_modulus);
            Coefficient invrCoef=m_multInvr[coefficient];
            m_Vk[dim].scaleCol(invrCoef, k, m_modulus); // set diagonal entry to 1
            m_Vk[dim].eraseEntry(cells[k],k);

            // swap their pivots and update pivot coefficients
            std::swap(pivots[j], pivots[k]);
            pivots[j].coef = static_cast<Coefficient>(
                static_cast<uint64_t>(factor) * pivots[j].coef % m_modulus);

            if(!pivots[k].coord.isNull())
            {
                pivots[k].coef = static_cast<Coefficient>(
                    static_cast<uint64_t>(invrCoef) * pivots[k].coef % m_modulus);
            }
        }
        else
        {
            // V[j] += factor * V[k].
            m_Vk[dim].addCol(factor, k, j, m_modulus);
            // addRow sees only the stored part of V[k]; account for its implicit V[k,k]=1.
            m_Vk[dim].eraseEntry(cells[k],j);
        }
    }

    // Move column j to i in D, R and V.  V's row move is represented by the
    // cell coordinates used as sparse row labels.
    std::rotate(cells.begin() + i, cells.begin() + j, cells.begin() + j + 1);
    m_Vk[dim].moveCol(j, i);
    std::rotate(pivots.begin() + i, pivots.begin() + j, pivots.begin() + j + 1);
}

void nsCCR::CubicalCplxSubd::moveRight(int dim, int i, int j)
{
    // assemble K={column indexes of V, s.t. V[i,k]!=0}={k1,...ks} //here k1=i
    CellOrder cellI=m_columnToReduce[dim][i];
    std::vector<int> K;
    for(int k=i+1;k<=j;k++)
    {
        if(m_Vk[dim].coefficient(cellI,k)!=0)
        {
            K.push_back(k);
        }
    }
    auto donorV=m_Vk[dim].column(i);
    // A saved donor must include the diagonal that is implicit in matrix storage.
    CellOrder donorDiagonal=cellI;
    donorDiagonal.coef=1;
    donorV.push_back(donorDiagonal);
    Coefficient donorRowCoef=1; // coef of donor column at celli 
    CellOrder donorPiv=m_pivotCol[dim][i];
    for(auto ki:K)
    {
        const Coefficient rowCoef=m_Vk[dim].coefficient(cellI,ki);
        const Coefficient coef=static_cast<Coefficient>(m_modulus -
            static_cast<uint64_t>(rowCoef)*m_multInvr[donorRowCoef]%m_modulus);
        const auto& targetPiv = m_pivotCol[dim][ki];
        if(donorPiv.coord.isNull() ||
        (!targetPiv.coord.isNull() && targetPiv < donorPiv))        
        {
            // V[ki]+=coef*donor_V
            m_Vk[dim].addCol(coef, ki, donorV, m_modulus);
        }
        else
        {
            auto nextDonorV=m_Vk[dim].column(ki);
            CellOrder nextDonorDiagonal=m_columnToReduce[dim][ki];
            nextDonorDiagonal.coef=1;
            nextDonorV.push_back(nextDonorDiagonal);
            const CellOrder nextDonorPiv=m_pivotCol[dim][ki];

            m_Vk[dim].addCol(coef,ki,donorV,m_modulus);
            // R is implicit through D*V; this branch gives it the donor's pivot.
            m_pivotCol[dim][ki]=donorPiv;
            if(!donorPiv.coord.isNull())
            {
                m_pivotCol[dim][ki].coef=static_cast<Coefficient>(
                    static_cast<uint64_t>(coef)*donorPiv.coef%m_modulus);
            }

            donorV=std::move(nextDonorV);
            donorPiv=nextDonorPiv;
            donorRowCoef=rowCoef;
        }
    }
    // Install the last saved donor, normalized to the implicit diagonal one.
    const Coefficient donorRowInv=m_multInvr[donorRowCoef];
    m_Vk[dim].setColumn(i,std::move(donorV));
    m_Vk[dim].scaleCol(donorRowInv,i,m_modulus);
    m_Vk[dim].eraseEntry(cellI,i);
    // R updates implicitly through D*V; normalize its cached pivot as well.
    m_pivotCol[dim][i]=donorPiv;
    if(!donorPiv.coord.isNull())
    {
        m_pivotCol[dim][i].coef=static_cast<Coefficient>(
            static_cast<uint64_t>(donorRowInv)*donorPiv.coef%m_modulus);
    }
    // Move column i to j in D, R and V. V's row permutation follows the
    // cell coordinates used as sparse row labels.
    auto& cells=m_columnToReduce[dim];
    auto& pivots=m_pivotCol[dim];
    std::rotate(cells.begin()+i,cells.begin()+i+1,cells.begin()+j+1);
    m_Vk[dim].moveCol(i,j);
    std::rotate(pivots.begin()+i,pivots.begin()+i+1,pivots.begin()+j+1);
}

void nsCCR::CubicalCplxSubd::MoveUp(int dim, int j, int i)
{
    // Suppose we already changed the filtration of simplex corresponding to row j.
    // But its order in m_columnToReduce did not update yet.
    // We update m_Vk to make it still reduced, but do not update m_columnToRed yet.

    // Step 1: reconstruct DV[k] for the column whose pivot is cellJ.
    CellOrder cellJ=m_columnToReduce[dim][j];
    auto cellK0Itr=std::find(m_pivotCol[dim-1].begin(),m_pivotCol[dim-1].end(),cellJ);
    if(cellK0Itr==m_pivotCol[dim-1].end())
    {
        return;
    }
    int k0=static_cast<int>(std::distance(m_pivotCol[dim-1].begin(),cellK0Itr));
    ReduceCol Vk,Rk;
    addCobd(k0,1,dim-1,Vk,Rk);
    // popPivot(Rk);

    //step 2: find k1: piv(R[k1])=piv-j(R[k0])
    CellOrder pivK=getPivot(Rk);
    if(pivK==cellJ)// pivot after move is the same as before
    {
        return;
    }

    auto cellK1Ptr = std::find(m_pivotCol[dim-1].begin(),m_pivotCol[dim-1].end(),pivK);
    m_pivotCol[dim-1][k0]=pivK;
    if(cellK1Ptr==m_pivotCol[dim-1].end())
    {
        return;
    }
    int k1=static_cast<int>(std::distance(m_pivotCol[dim-1].begin(),cellK1Ptr));
    // Coefficients at the shared row pivK, before ordering the two columns.
    Coefficient pivotCoef0=pivK.coef;
    Coefficient pivotCoef1=cellK1Ptr->coef;
    while(true)
    {
        // k0=min(K), k1=max(K)
        if(k0>k1)
        {
            std::swap(k0,k1);
            std::swap(pivotCoef0,pivotCoef1);
        }

        // Cancel the shared row: pivotCoef1 + factor*pivotCoef0 = 0 (mod p).
        const Coefficient factor=static_cast<Coefficient>(m_modulus -
            static_cast<uint64_t>(pivotCoef1)*m_multInvr[pivotCoef0]%m_modulus);
        m_Vk[dim-1].addCol(factor,k0,k1,m_modulus);
        CellOrder sourceDiagonal=m_columnToReduce[dim-1][k0];
        sourceDiagonal.coef=factor;
        m_Vk[dim-1].addEntry(k1,sourceDiagonal,m_modulus);

        // R is implicit: reconstruct the updated target column as D*V[k1].
        Vk=ReduceCol();
        Rk=ReduceCol();
        addCobd(k1,1,dim-1,Vk,Rk);
        CellOrder pivotK1=getPivot(Rk);
        auto cellKItr=std::find(m_pivotCol[dim-1].begin(),m_pivotCol[dim-1].end(),pivotK1);
        m_pivotCol[dim-1][k1]=pivotK1;
        if(cellKItr==m_pivotCol[dim-1].end())
        {
            break;
        }
        k0=static_cast<int>(std::distance(m_pivotCol[dim-1].begin(),cellKItr));
        pivotCoef0=cellKItr->coef;
        pivotCoef1=pivotK1.coef;
    }
}

void nsCCR::CubicalCplxSubd::moveDown(int dim, int i, int j)
{
    // we update matrixV and pivot column before moving simplex at i down to j

    // step 1: K = {column indexes k of R, s.t. R[i,k]!=0, and piv[k]<j}
    std::vector<int> K;
    std::vector<Coefficient> coefK;//coef of R[i,k]

    const CellOrder cellJ=m_columnToReduce[dim][j],cellI=m_columnToReduce[dim][i];

    for(int k=0;k<static_cast<int>(m_columnToReduce[dim-1].size());++k)
    {
        const CellOrder& pivot=m_pivotCol[dim-1][k];
        if(pivot.coord.isNull() || pivot < cellJ) continue;
        ReduceCol Vk,Rk;
        addCobd(k,1,dim-1,Vk,Rk);
        while(true)
        {
            const CellOrder entry=popPivot(Rk);
            if(entry.coef==0) break;
            if(entry==cellI)
            {
                K.push_back(k);
                coefK.push_back(entry.coef);
                break;
            }
        }
    }

    if(K.empty()) return;

    auto donorCol=m_Vk[dim-1].column(K[0]);
    // Materialize the donor's implicit diagonal in the saved column.
    CellOrder donorDiagonal=m_columnToReduce[dim-1][K[0]];
    donorDiagonal.coef=1;
    donorCol.push_back(donorDiagonal);
    CellOrder donorPiv=m_pivotCol[dim-1][K[0]];
    Coefficient donorCoef=coefK[0];
    for (int index=1;index<K.size();index++)
    {
        int k=K[index];
        if(m_pivotCol[dim-1][k]<donorPiv)
        {
            // V[k]+=coef*donor_column_V to eliminate R[i,k]
            const Coefficient coef=static_cast<Coefficient>(m_modulus -
                static_cast<uint64_t>(coefK[index])*m_multInvr[donorCoef]%m_modulus);
            m_Vk[dim-1].addCol(coef,k,donorCol,m_modulus);
        }
        else
        {
            // Save the original target, including its implicit diagonal, as the next donor.
            auto nextDonorCol=m_Vk[dim-1].column(k);
            CellOrder nextDonorDiagonal=m_columnToReduce[dim-1][k];
            nextDonorDiagonal.coef=1;
            nextDonorCol.push_back(nextDonorDiagonal);
            const CellOrder nextDonorPiv=m_pivotCol[dim-1][k];
            const Coefficient nextDonorCoef=coefK[index];

            const Coefficient coef=static_cast<Coefficient>(m_modulus -
                static_cast<uint64_t>(nextDonorCoef)*m_multInvr[donorCoef]%m_modulus);
            m_Vk[dim-1].addCol(coef,k,donorCol,m_modulus);

            m_pivotCol[dim-1][k]=donorPiv;
            m_pivotCol[dim-1][k].coef = static_cast<Coefficient>(static_cast<uint64_t>(coef) * donorPiv.coef % m_modulus);

            donorCol=std::move(nextDonorCol);
            donorPiv=nextDonorPiv;
            donorCoef=nextDonorCoef;
        }
    }
    
    m_pivotCol[dim - 1][K.front()] = getCellOrder(cellI.coord, coefK.front());    
}
