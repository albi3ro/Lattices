# Lattices
## Module to initialize common material lattices in Julia
### Written by Christina Lee
### Graduate Student, Okinawa Institute of Science and Technology

Useful for prototyping computational condensed matter physicists

Types of Lattices:
<ul>
<li>Square
<li>Chain
<li>To be added

Members of Type
<ul>
<li>name  , a string
<li>l , length in number of unit cells
<li>dim, dimension of lattice
<li>a, array containing the basis vectors by which positions are generated
<li>unitA. array of positions inside a single unit
<li> N , number of total sites
<li> X, array of positions
<li> nnei, number of nearest neighbors
<li> neigh, Array of nearest neighbors [i][j], where i is site and j is 1:nnei

To Initialize a new lattice:

  Call `lt=MakeLattice("Name",sizeoflattice)`
  
To plot the bonds

  call `PlotNeighbors(lt)`
  
To use:
  
  Put this file on your path
  
  Call `using Lattices`
