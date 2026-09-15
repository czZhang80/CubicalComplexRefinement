// psedo-code for moves algorithms and insertion algorithm
// 1. MoveLeft
// psedo-code for move left: for decomposition R=DV, 
// Conduct column operation on R and V s.t. after moving column j of D to i,
// (this corresponds to move j-th column of R and V to i-th column, move j-th row of V to i-th row)
// restore R reduced and V invertible upper trianglar
// We write this psedo-code over field Z/2Z
// notes:
// columns between [i,j-1] cannot satisfy condition low(R[k])>low(R[j])
// MoveLeft(R,V,j,i) (moves column j to i, i<j)
// for(int k=j-1;k>=i;k--)
// {
//   if(V[k][j]!=0)
//   {
//       if(low(R[k])>low(R[j]))
//       {

//           r_donor_column=R[j]
//           R[j] = R[j]+R[k]
//           R[k] = r_donor_column

//           v_donor_column=V[j]
//           V[j]=V[j]+V[k]
//           V[k]=v_donor_column

//       }
//       else
//         {          
//             R[j] = R[j] + R[k]
//             V[j] = V[j]+V[k]
//         }      
//   }
// }
// move j-th column of R and V to i-th column, move j-th row of V to i-th row
// return (R,V,d_R,d_V)

// 2. MoveRight
// psedo-code for move right: for decomposition R=DV, 
// Conduct column operation on R and V s.t. after moving column i of D to j, (i<j)
// (this corresponds to move i-th column of R and V to j-th column, move i-th row of V to j-th row)
// restore R reduced and V invertible upper trianglar
// notes:
// This is exactly the 1st part of the original MoveRight algorithm
// For column l corr. to a 0pap in [i+1,j], V[l]=el, so l not in K
// MoveRight(R,V,i,j)
//{
//  K={column indexes of V, s.t. V[i,k]!=0}={k1,...ks} //here k1=i
//  donor_R=R[k1], donor_V=V[k1], donor_piv=piv(R[k1])
//  for(ki:K\{k1})
//  {
//      if(piv(R[ki])>donor_piv)
//      {
//          R[ki]+=donor_R
//          V[ki]+=donor_V
//      }
//      else
//      {
//          temp_R=R[ki], temp_V=V[ki], temp_piv=piv(R[ki])
//          R[ki]+=donor_R
//          V[ki]+=donor_V
//          donor_R=temp_R, donor_V=temp_V, donor_piv=temp_piv
//      }
//  }
//}

// 3. MoveUp
// for decomposition R=DV, move j-th row of D to i, (j>i).
// (This corresponds to moving j-th row of R to i)
// notes:
// this will not affect 0pap: 0pap have all its left entry of pivot 0, so it cannot appear at k1
// only 2 columns are not reduced in each stepof the process and k found in if(k: piv[k]=piv[k1] not found) is unique
// each iteration makes pivot up by at least one, and the iteration has at most j-1-i loops. 
//  MoveUp(R,V,j,i)
// {
//  K={k0, k1: piv(R[k0])=i, piv(R[k1])=piv-j(R[k0])}
//  while(true)
//  {
//      k0=min(K), k1=max(K)
//      R[k1]+=R[k0], V[k1]+=V[k0]
//      if(k: piv[k]=piv[k1] not found) break
//      K={k1,k: piv[k]=piv[k1]}
//  }
// }

// 4. MoveDown(i,j)
// notes: 
// for column l corr. to 0pap piv[l]=l-1>j, otherwise move will be illegal, so l not in K
// Let piv-i(k) be the pivot of column k without row i.
// K = {column indexes k of R, s.t. R[i,k]=1, and piv[k]<j}
// donor_column_V = V[k_min], donor_column_R = R[k_min], donor_piv-i=piv-i(kmin)
// for(k:K\{k_min} in ascending order)
//{
//  if(piv-i(k)>donor_piv-i)
//  {
//      R[k]+=donor_column_R
//      V[k]+=donor_column_V
//  }
//  else
//  {
//      temp_R=R[k], temp_V=V[k], donor_piv-i=piv-i(k)
//      R[k]+=donor_column_R
//      V[k]+=donor_column_V
//      donor_column_R=temp_R, donor_column_V=temp_V 
//  }
//}