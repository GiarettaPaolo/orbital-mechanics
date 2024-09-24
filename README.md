# Orbital Mechanics: Optimizing Satellite Transfers

_How to efficiently transfer a satellite between two orbits?_

## Objective
The goal is to transfer a satellite from an initial orbit to a final orbit, identified by the following orbital parameters:
<div align="center">
  <img src="Images/tables.png" width="600"/>
</div>

The challenge is to optimize the transfer using impulsive burns by minimizing the required **delta-v** and transit time. The project explores the standard Hohmann transfer with periapsis and plane change maneuvers and **six additional transfer methods** to improve efficiency.

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
