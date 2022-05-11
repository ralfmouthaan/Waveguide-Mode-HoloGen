# Waveguide-Mode-HoloGen

Code for generating holograms capable of exciting discrete modes in waveguides. 

The target profile of the waveguide mode must be known, and can be calculated:
- Using analytical techniques for some waveguide geometries. For example, the modes of a step index fibre can be calculated using the mathematics described in Gloge's [paper](https://doi.org/10.1364/AO.10.002252), and as implemented in Waveguide-Analytical.
- Calculated using an FDFD model for any waveguide geometry, as described in my [paper](https://doi.org/10.1109/JLT.2021.3124469) and as implemented in Waveguide-FDFD.

Extremely high-purity modes can be coupled into with large amounts of power if direct search (DS) algorithms. are used to generate the hologram. The DS algorithms implemented here for this purpose have been extensively used, and should work reasonably well.

Other approaches for generating holograms capable of modal excitation an run faster, but yield lower-quality modes (mostly with less power coupled into the target mode). Some options include algorithms that exploit an experimental geometry with a spatial filter, such as those by [Arrizon](https://doi.org/10.1364/JOSAA.24.003500) and that by [Bolduc](https://doi.org/10.1364/OL.38.003546). An attempt has been made to implement these algorithms, but my efforts should be treated with caution - I think only Arrizon's third algorithm gives reasonable results.

Alternatively, iterative phase retrieval algorithms such as Gerchberg-Saxton can be used with a region of interest defined in the replay field plane. This will be implemented at some point.
