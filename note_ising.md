# Ising Model

## Description

The main objective of this notebook is to present the famous Ising model in 2 dimensions. In particular, explore the numerical simulations and compare them with the theoretical predictions, well known in $1D$ and $2D$. The Ising model considers a uniform and equally spaced grid $[0, \cdots, N - 1] \times [0, \cdots, N - 1]$, where at each node $(i,j)$ is located an atom, with a spin quantum number $S_{i,j}$ assigned. One makes the assumption that this (normalized) spin can take only two values: $S = 1$ or $S = -1$, accordingly to the Quantum theory. This is simplified because in practice, an atom could have a much larger spin, possibly integer.


This toy model is quite useful in statistical physics and Quantum Field Theory, because, despite its simplistic approach to model magnetism, it contains an impressive amount of physical phenomenon, such as $2^{\text{nd}}$ order phase transition, renormalization groups, mean field... It is important to keep in mind that despite its popularity among physicists and age, a lot questions remain open. The partition function could only be calculated analytically for systems living in one or two spatial dimensions, with periodic boundary conditions. The latter case assumes no magnetic field is applied to the system.

For two dimensional systems, the very first calculation of the partition function is due to Onsager. Its demonstration is globally recognised as an impressive feat for its difficulty and length; it is accessible only for master-level students at least. 
Finally, the cases $\mathbb{R}^d, d \geq 4$ is approximated quite well by the Curie-Weiss model, showing its interest despite its unrealistic assumption. In a nutshell, only the case $d=3$ strangely resists to the physicist community, which is unfortunate. The world we live in is precisely $\mathbb{R}^3$!


In statistical physics, a problem is in general fully solved when you know how to compute analytically the partition function, usually noted $Z$, because most of the thermodynamical quantities of interest can be directly derived from it and its partial derivatives. The partition function is defined in the classical theory as $Z = \sum\limits_\alpha e^{-E_\alpha/kT}$, where the sum is made over all configurations $(\alpha)$ of the sysem in the thermal bath at temperature $T$. Indeed, it is known that the free energy $F$ of a system if directly related to the partition function by the reatlion $F = - kT\log(Z)$. In an Ising model, a configuration of the system is the knowledge of the spin's value at every node.

The hamiltonian/Energy of an Ising system is $E = -g\mu\_B B \sum\limits_{i, j} S_{i,j} - \sum\limits_{i, j} \sum\limits_{(u,v) \sim \text{NN}(i, j)}J_{i,j} S_{i,j}S_{u,v} = \sum_{(i,j)} E_{(i,j)}$, where $E_{(i,j)}$ is the hamiltonian at node $(i,j)$.
The Index $\text{NN}(i, j)$ denotes the set of nearest neighbours of the node $(i,j)$, this is the term that is truly dimension-dependant. In $1D$, an inside-node has $2$ nearest neighbours. No phase transition occurs at all. In $2D$ however, an inside-node has $4$ nearest neighbours, and a $2^{\text{nd}}$ order phase transition actually occurs. Hence, this slight change in the equations fundamentally change the behaviours of the system, and induces a discontinuity around a critical temperature $T\_c$.

Let's break down the different terms. The first term contains the interaction with the environment. Recall that in classical theory of magnetism, one can show that $E = -m\cdot B$, where $m$ is the magnetization vector of the system. There is a kind of similar term for quantum systems: $H = -g\mu\_B S \cdot B$.
1. $g$ is the g-factor. It is a scalar that relates the angular momentum with the magnetic moment. In first approximation, for quantum spin, $g \approx 2$. 
2. $B$ is the exterior magnetic field imposed by the operator. If $B \neq 0$, then no phase transition occurs.
3. $\mu\_B$ is the Bohr magneton.
4. $S_{i,j}$ is the spin at node $(i,j)$.

The second term describes the self-interaction of the system, via a coupling parameter $J_{i,j}$ called [exchange interaction](https://en.wikipedia.org/wiki/Exchange_interaction). Its justifications lies in quantum physics and the nature of the particles (bosons or fermions), it is rooted to the indistinguishable nature of the quantum components, a key feature too in classical statistical physics. It affects the enumeration of the possible configurations of a system, which mechanically impacts the computation of the partition function $Z$. When $J > 0$, the system is said ferromagnetic, when $J < 0$, it is said paramagnetic.

### Exact solution in 2D

Like I said above, the 2D model admits an exact solution, written below:
$Z = \bigg(2\cosh\large(\frac{2J}{k\_BT}\large)\exp{I}\bigg)^N,$
where $N$ is the number of spin in the lattice, $I = \frac{1}{2\pi}\displaystyle\int_{0}^\pi d\varphi \log(1/2[1 + (1 - x(T)^2\sin^2(\varphi))^{1/2}])$ is an elliptical integral and $x(T) = \frac{2\sinh(2J/k\_BT)}{\cosh^2(2J/k\_BT)}$. Define $\beta = \frac{1}{k\_B T}$
In 2D Ising Model, the critical temperature $T\_c$ verifies $\sinh(2J/kT\_c) = 1$ according to Onsager's exact model. It can be equivalently expressed as:
$T\_c = \frac{2J}{k\_B\log(1 + \sqrt{2})}.$

Onsager demonstrated the such system exhibits a second order phase transition, with order parameter $M = \sum\limits_{(i, j) \in L} S_{(i,j)}$. It means that whenever $T > T\_c$, $M = 0$, no macroscopic magnetization is observed, but if $T < T\_c$, the material acquires a spontaneous magnetization $M = (1 - \sinh(2\beta J)^{-1})^{1/8}.$ Note that the transition is continuous, i.e. $M = 0$ at $T = T\_c$. That is a specific behaviour to second order phase transition. However, if a magnetic field is applied to the system, no phase transition occurs at all, for any temperature. It is theoretically very simple to prove. If we note $\phi = F - HM$, the free energy of the system, then at equilibrium, one should have $\partial\_M\phi = \partial\_M F - H = 0$, implying that $\partial\_M F \neq 0$.

Let's now more discuss about the simulation itself.

## About the simulation

The first remark comes from the fact that computing the partition function by actually summing over all possible configurations $(\alpha)$ is infeasible in practice. If there is $N^2$ nodes in the grid, then there is a total of $2^{N^2}$ possible configurations. In practice, one sums only over the most probable configurations, that is, those that are not too 'far' from the current one. Any macroscopic quantity at equilibrium can be computed as $\langle A \rangle = \displaystyle\frac{\sum_\alpha A_\alpha e^{-E_\alpha/kT}}{Z}$, where $A_\alpha$ is the value taken by $A$ at the particular configuration $\alpha$. This quantity is an average because the system explores many different configurations on a macroscopic time scale, so measuring $A_\alpha$ would not make actually sense. 

Also, $J_{i,j}$ will be assumed constant at all node.
The algorithm implemented as described below.

1. Choice of the network of size $N \times N$ and temperature $T$
2. Sample an initial state
3. Choice of a site $(i, j)$ at random.
4. Evolution of the site.
5. Iterate steps 3. and 4.
6. Compute the average quantities over the multiple configurations visited.

To make evolve a site $(i,j)$ taken at random, first compute its energy: $H_{i,j} = -g\mu\_B B S_{i,j} - J S_{i,j}\sum\limits_{(i, j) \sim (m,n)}S_{m,n}$.
Then flip its sign $\bar E_{(i,j)} = -E_{(i,j)}$. If the flipped spin is energetically less favoured, do not change it. Otherwise, flip it with probability $\exp(-(\bar E_{(i,j)}  - E_{(i,j)})/kT) = \exp(2E_{(i,j)}/kT)$. 

In a physically exact model, $T\_c = 657 K$ for $J = 4e-21$ Joule. To simulate a ferromagnetic material, simply take the opposite sign.

Note: The `exp` method is computationally very expensive, so avoiding this calculation at every loop would save a significat amount of time. Noticing that only five energy states are possible if $B=0$, the exponential factor takes value in the set $\{1, e^{2\beta J}, e^{4\beta J}, e^{-2\beta J}, e^{-4\beta J} \}$, computing all of them once is enough. We can actually restrict to the set $\{1, e^{2\beta J}\}$ and get the three other values by taking the square and/or the inverse. Finally, since a flip might occurs if and only if $E_{i,j} < 0$, in practice, one considers only the set: $\{e^{-2\times 2\beta J}, e^{-2 \times 4\beta J}\}$.
I did not observe a significant increase in speed with this technique.
