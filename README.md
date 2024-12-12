# RNA_Triplets
BioInformatics project 

---

We aim at providing a C++ implementation of a dynamic programming scheme for the Strand soup Interaction model. This model is found in the following paper : 

RNA Triplet Repeats: Improved Algorithms for Structure Prediction and Interactions
Kimon Boehmer1, Sarah J. Berkemer1,2, Sebastian Will1, Yann Ponty1


___

Our first task is to provide an implementation for the nussinov algorithm. 
We will make a first implementation in python and then when the logic is well understood, change to a c++ implementation.

___

## Nussinov algorithm

This algorithm aims at predicting the secondary structure of the RNA. It does so by finding the optimal base-pairing that maximizes the number of base pairs in a sequence.

#### The problem

Given an RNA sequence S of length n. The goal is to find the structure tha maximizes the number of base pairing.
It is represented via a Matrix M[i,j] :
- i and j represent indices of bases in the RNA sequence
- M[i,j] is the maximum number of base pairing in the sequence S[i,j]

#### Solving

We use a dynamic programming approach. It is composed of two main functions : the filling of the energy matrix and then the reconstruction of the minimum energy structure based on this matrix.

### Function 1 : Filling the minimum energy matrix

#### **Input**: 
- `ω` – RNA of size n

#### **Output**:
- `m` – Minimum energy matrix m

``` 
    Fonction FillMatrix (ω):
        m ← EmptyMatrix(n × n)
        ## Initialize with 0 all the values of the diagonal up to θ.
        for i ← 1 to n do
            for j ← i to min(i + θ, n) do
                mᵢ,ⱼ ← 0
        for i ← n to 1 do
            for j ← i + θ + 1 to n do
                Case A: Position i left without partner
                    mᵢ,ⱼ ← mᵢ₊₁,ⱼ
                Case B: Positions i and j form a base pair
                    mᵢ,ⱼ ← min(mᵢ,ⱼ, mᵢ₊₁,ⱼ₋₁ + E^ωᵢⱼ)
                Case C: Position i paired to k < j
                    for k ← i + θ + 1 to j − 1 do
                        mᵢ,ⱼ ← min(mᵢ,ⱼ, mᵢ₊₁,ₖ₋₁ + mₖ₊₁,ⱼ + E^ωᵢ,ₖ)
        return m
```

### Function 2: Backtracking for the Minimum Energy Structure

#### **Input**:
- `[i, j]` – Region under consideration  
- `m` – Dynamic programming matrix, previously computed  
- `ω` – RNA sequence of length \( n \)  

#### **Output**:
- `S*` – Structure minimizing free energy  

```
Function Backtrack(i, j, m, w):
    if j - i <= theta then 
        return *....* #The empty structure has min energy
    else
        **Case A: Position i left without partner**
        if mᵢ,ⱼ = mᵢ₊₁,ⱼ then
            S*ᵢ ← Backtrack(i+1, j, m, w)
            return • S*ᵢ
        Case B : Positions i and j form a base pair
        if mi,j = mi+1,j−1 + Eω then
            Sᵢ,ⱼ⋆ ← Backtrack(i+1, j−1, m, w)
            return S*ᵢⱼ
        Case C: Position i paired to k < j
        for k ← i + θ + 1 to j − 1 do
            if mᵢⱼ = mᵢ₊₁,ₖ₋₁ + mₖ₊₁,ⱼ + Eωᵢₖ then
                S*₁ ← Backtrack(i+1, k−1, m, w)
                S*₂ ← Backtrack(k+1, j, m, w)
                return (S1⋆)S2⋆
```

---

Questions and or Ideas that needs to be discussed :

* the value of theta
* The value of the energy pairing
* Input format (what is the common input, for the moment we just have a string : 'ACGUC...')
* Output format (list of pairs of indices, Dot-Bracket Notation, adjacency matrix, graph representation for better visualization?)
* do we have standard RNA sequence that can be used for tests to verify correctness ?
* should add error handling such as input error etc for better readability
* We are only computing one minimal energy secondary structure S but there are multiple minimum structures. Do we care ? How do we compute ? 



---

Strand soup in the future...

We start with m strands of various sizes. We select a starting strand s and an ending strand r.

Then for each of the remaining strands we will recursively call the strands to extend the structure until every strand is used. 
