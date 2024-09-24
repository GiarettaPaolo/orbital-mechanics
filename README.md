# Orbital Mechanics: Optimizing Satellite Transfers

_How to efficiently transfer a satellite between two orbits?_

## Objective
The goal is to transfer a satellite from an initial orbit to a final orbit, identified by the following orbital parameters:

\[
\begin{array}{|c|c|c|}
\hline
\textbf{Parameter} & \textbf{Initial Orbit (r_1, v_1)} & \textbf{Final Orbit (a_2, e_2, i_2, \Omega_2, \omega_2, \theta_2)} \\
\hline
r_1 & \begin{bmatrix} 2254.3254 \\ -8092.3126 \\ -4199.8027 \end{bmatrix} & \text{—} \\
v_1 & \begin{bmatrix} 5.6120 \\ 2.4220 \\ -1.7020 \end{bmatrix} & \text{—} \\
a_2 & \text{—} & 16410.0000 \\
e_2 & \text{—} & 0.2678 \\
i_2 & \text{—} & 0.5612 \\
\Omega_2 & \text{—} & 0.4075 \\
\omega_2 & \text{—} & 1.0700 \\
\theta_2 & \text{—} & 1.3420 \\
\hline
\end{array}
\]

The challenge is to optimize the transfer using impulsive burns by minimizing the required **delta-v** and transit time. The project explores not only the standard Hohmann transfer but also **six additional transfer methods** to improve efficiency.

---

## Initial Orbit Visualization
Below are some plots showing the initial orbital conditions:

<div align="center">
  <img src="Images/orbits.png" width="250"/>
  <img src="Images/initialPerifocal.png" width="250"/>
  <img src="Images/finalPerifocal.png" width="250"/>
</div>

---

## Example Transfers
Here are examples of different transfer methods, including the standard and additional proposed techniques:

<div align="center">
  <img src="Images/standard.png" width="250"/>
  <img src="Images/alternative2.png" width="250"/>
  <img src="Images/alternative6.png" width="250"/>
</div>

---

## Optimization Approach
The transfer optimization is based on:
- **Delta-v minimization**: Reducing the total velocity change required.
- **Time minimization**: Balancing fuel efficiency with transit time.

---

Feel free to explore the details of the transfer methods, code, and simulations in the repository!
