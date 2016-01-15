# Lattices
## Module to initialize common material lattices in Julia
### Written by Christina Lee
### Graduate Student, Okinawa Institute of Science and Technology

Useful for prototyping computational condensed matter physicists

### Types of Lattices:
<ul>
<li>Square
<li>Triangular
<li>Chain
<li>Honeycomb
<li>To be added

### Members of Type
<ul>
<li>name  , a string
<li>l , length in number of unit cells
<li>dim, dimension of lattice
<li>a, array containing the basis vectors by which positions are generated
<li>unit. array of positions inside a single unit
<li> N , number of total sites
<li> X, array of positions
<li> nnei, number of nearest neighbors
<li> neigh, Array of nearest neighbors [i][j], where i is site and j is 1:nnei


### Calling
| To:   | Do:  |
|------------------------------|---------------------------------------------|
| Use: | Put this file on your path, call `using Lattices`  |
| Initialize a new lattice: | Call `lt=MakeLattice("Name",sizeoflattice)` |
| Access a member: | ex. `lt.name` or  `lt.nnei` |
| Plot the bonds | call `PlotNeighbors(lt)` |



To put a file on your path:
<ul>
<li>each time execute the code `push!(LOAD_PATH,"/some/location")`
<li>or add to to .juliarc.jl in your home directory the line `@everywhere push!(LOAD_PATH,"/some/location")`
