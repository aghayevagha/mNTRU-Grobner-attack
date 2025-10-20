# mNTRU-Grobner-attack
This repository contains source code for attacking multiple-key NTRU problem.
It builds a system of polynomial equations using Arora-Ge modeling and solves it using Grobner basis algorithms. 
In mNTRU problem, there is a fixed polynomial $g$ and $m$ number of of secret polynomials $f_i$ with a narrow support, and the public keys are computed as $h_i=f_i* g^{-1}$ for $i \in [m]$. Using public key polynomials $h_i$, we generate a set of polynomials equations using Arora-Ge modeling.
