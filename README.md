# Nori
Volumetric path tracer based on the educational path tracer Nori, by Wenzel Jakob. Developed as the final project for the course "Modeling and Simulation of Appearance" within the MRGCV program at Universidad de Zaragoza, during the academic year 2022-2023.
The program allows for rendering all kinds participating media, both homogeneous and heterogeneous. It allows one global, infinite homogeneous medium and an infinite amount of bounded homogeneous/heterogeneous mediums.

It supports ```.vdb``` file loading in an indirect manner, by converting them to a Mitsuba-friendly format, ```.vol``` and loading those. Only supports the Henyey-Greenstein phase function (Rayleigh is coded but not tested).

The integrator also incorporates null-collision techniques such as Woodcock Tracking for the distance sampling, and Ratio Tracking for the transmittance estimation.

# Report
The report with some details on the implemented techniques can be found [here](https://github.com/PedroPerez14/MSoA-participating-media/master/report_compressed.pdf)!

# Compilation details
Since it's an academic project there are not many details about how the program works, or how to compile and run it. If you need them feel free to contact me.

# Renders
<p align="center">
  <img src="https://github.com/PedroPerez14/Nori2/assets/56193395/3c6aad43-c295-41b4-a313-2eb5cfda3ae9" width="40%" height="40%">
  <img src="https://github.com/PedroPerez14/Nori2/assets/56193395/97e0c894-c59e-4a30-a41c-c93ef9e1c3ab" width="70%" height="70%">
</p>
