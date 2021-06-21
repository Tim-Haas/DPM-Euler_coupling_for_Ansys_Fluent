# DPM-Euler coupling for Ansys Fluent

Lagrange Particle Tracking (LPT) is commonly implemented with the particle source in a cell method by Clift et al. This method was originally developed for sprays, thus, is based on the assumption that the particle is small compared to the Eulerian grid. This assumption is violated for larger particles like bubbles, especially on fine Eulerian grids such as those for LES simulations.

In those situations, momentum coupling with a single cell causes unphysical results like pseudo-turbulence (particle velocity fluctuations then the particle crosses a cell face) and an overestimate of both, the fluid and particle velocity (see Fig 1-2). Finally, it causes divergence. 

This UDF assembly for unstructured hexahedral grid in Ansys Fluent makes LPT fully independent of the Eulerian grid. Cells within the particle‘s influencing sphere are identified by an face loop cascade. The Eulerian fluid properties are mapped to the Lagrange particle by the Deen mapping polynom and the particle's momentum is distributed back by the same approach. Different drag and lift closures can be chosen.

 <img src="Coupling_schemes.png" alt="Coupling Schemes: (a) Forward (b) Backward coupling"> 

<b>Available drag coefficient models</b>
<br> -Tomiyama drag correlation for pure systems
<br>-Tomiyama drag correlation for slightly contamined systems
<br>-Tomiyama drag correlation for fully contamined systems
<br>-Bozzano drag correlation
<br>-Dijkhuizen drag correlation
<br><i>Choose the model by setting the DRAG_FLAG</i>

<b>Available lift coefficient models</b>
<br>-No lift force
<br>-Tomiyama lift correlation
<br>-Ziegenheim lift correlation
<br>-Constant lift coefficient 0.5
<br><i>Choose the model by setting the LIFT_FLAG</i>
