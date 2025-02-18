# Rotational-Symmetrical-DNA-origami-design-code
This git contains information about the rotational symmetry DNA origami design and code for calculation on the bridge staples distance.

# The first square tile
The DNA origami is a type of programmable DNA nanostructures with nucleotide level accuracy. The design of origami itself follows different rules. The basic rule to design classic 2D and 3D origamis can be found here: https://doi.org/10.6028/jres.126.001. Here I am focused on a special type of DNA origmai that is made to be rotational symmetrical, first published by Tikhomirov et al. in 2011 in the paper "Programmable disorder in random DNA tilings" (https://www.nature.com/articles/nnano.2016.256). Shortly followed the publication, Tikhomirov et al. drew the world's first microscale Mona Lisa painting with this DNA origami tile in the paper "Fractal assembly of micrometre-scale DNA origami arrays with arbitrary patterns" (https://www.nature.com/articles/nature24655). The micro Mona Lisa painting became the journal cover: https://www.nature.com/nature/volumes/552/issues/7683.
![image](https://github.com/user-attachments/assets/7037190b-37c6-492f-8a8b-86cfa7f5fa1c)

The sequence and design file is publicly available: https://nanobase.org/structure/132. The AFM image of the 4 fold rotational symmetrical DNA origami tile:
![08191454_0_00001_spm](https://github.com/user-attachments/assets/d8c8e566-4835-4a89-a315-e339d53d01d2)

# Polygon tiling expansion
The rotational symmetrical origamis continue to be developed into, Tikhomirov et al. published the triangle version of the rotational symmetrical DNA origami triangle tile in the paper "Triangular DNA Origami Tilings" (https://pubs.acs.org/doi/10.1021/jacs.8b10609). Recently (at the time of this repository), the Tang et al. published the "DNA Origami Tessellations" (https://pubs.acs.org/doi/10.1021/jacs.3c03044), where the idea of rotational symmetry was further developed to enable different design such as the hexagon, and the kite shape as well. In which they introduced the individual modular approach, I described it as pizza slicing, as polygons can always be sliced into a triangle, just like a pizza.
The AFM image of the 3 fold rotational symmetrical DNA origami tile:
![08191526_0_00001_spm](https://github.com/user-attachments/assets/8e09f481-6b20-4829-a340-7b726540fb31)

# Design
The design of this rotational symmetrical tile is not straight forward, first of all, it does not quite follow the standard DNA origami design rules that was outlined in the tutorial document above, it was quite overwhelming when I first opened the Cadnano design file for the 4 fold rotational symmetrical tile. The calculation of scaffold route length for each triangle sub-unit must be accounted for and connecting staples' lengths (the bridge staple) as well. This is to make sure that the structure wount be too tight but it is "stitched" into the shape of the tile, be that square or triangle or hexagon. Tikhomirov et al. and Tang et al. both published the calculation to perform the calculation in Euclidean space (Tang et al. provided a very detail description of the calculation procedure in the supporting information of the paper).

# Design a pentagon
Here, I will design a pentagon, partially inspired by the yorkshire white rose in the UK and the property of pentagon, for what makes pentagon so unique: "The Infinite Pattern That Never Repeats" by Veritasium
(https://www.youtube.com/watch?v=48sCx-wBs34&t=1s)





# Disclaimer
The images presented in this repository (github) are original and generated from my laboratory research.
