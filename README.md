# mNTRU-Grobner-attack
This repository contains source code for attacking multiple-key NTRU problem.
It builds a system of polynomial equations using Arora-Ge modeling and solves it using Grobner basis algorithms. 
In mNTRU problem, there is a fixed polynomial $g$ and $m$ number of of secret polynomials $f_i$ with a narrow support, and the public keys are computed as $h_i=f_i* g^{-1}$ for $i \in [m]$. Using public key polynomials $h_i$, we generate a set of polynomials equations using Arora-Ge modeling.

## Code
We have code examples for SageMath and Magma systems. Sagemath seems to be inefficient due to the poor choice of underlying Grobner basis algorithm, and computation of variety is also slow. The variety function in Sage also includes Grobner basis computation as a subroutine. 

Magma algebra system is known for paralellized and efficient computation for Grobner bases and variety of ideals. binary.magma contains source code for attacking the generated system with or without binary constraints :
$x_i(x_i-1)=0$

ternary_general.magma contains the same attack style and also allows for generating samples with E=[-B,-B+1,..,B-1,B] support, and computes grobner basis with or without the following constraints:
$(x_i+B)(x_i+B-1)...x_i...(x_i-B+1)(x_i-B)=0.$
