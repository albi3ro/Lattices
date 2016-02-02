module Lattices
export Lattice, MakeLattice, PlotNeighbors, MakeX, MakeSquare, MakeCube, MakeChain, MakeHoneycomb,
    MakeTriangular, MakeCheckerboard, MakeKagome
    using PyPlot

    type Lattice
        name::ASCIIString
        l::Int
        dim::Int
        a::Array
        unit::Array
        N::Int
        X::Array
        nnei::Int
        neigh::Array
    end

    function armod(x,y)
        return mod(x-1+y,y)+1
    end

    function parity(X)
        return 1-2*mod(X,2)
    end


    function MakeLattice(name::ASCIIString,l::Int)
        if(name=="Square")
            return MakeSquare(l);
        elseif(name=="Chain")
            return MakeChain(l);
        elseif(name=="Honeycomb")
            return MakeHoneycomb(l);
        elseif(name=="Triangular")
            return MakeTriangular(l);
        elseif(name=="Checkerboard")
            return MakeCheckerboard(l);
            #elseif(name=="Pyrochlore")
        elseif(name=="Kagome")
            return MakeKagome(l);
        else
            println("I didnt give a good name, making a 2x2 square")
            return MakeSquare(2)
        end

    end

    function MakeSquare(l::Int)
        N::UInt16=l^2;
        d::Int8=2;
        X=Array{Int16}(N,2);
        nnei::UInt8=4;
        neigh=Array{Int16}(N,4);

        X[:,1]=armod(collect(1:N),l);
        X[:,2]=ceil(collect(1:(N))/l);
        a=[[1 0]
            [0 1]];
        unit=[0 0];

        neigh[:,1]=armod(X[:,1]+1,l)+l*(X[:,2]-1);
        neigh[:,2]=armod(X[:,1]+l-1,l)+l*(X[:,2]-1);
        neigh[:,3]=X[:,1]+mod(X[:,2],l)*l;
        neigh[:,4]=X[:,1]+mod(X[:,2]+l-2,l)*l;

        return Lattice("Square",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeCube(l::Int)
        N::UInt16=l^3;
        d::Int8=3;
        nnei::UInt8=6;
        neigh=Array{Int16}(N,6);

        a=[[1 0 0 ]
           [0 1 0]
           [0 0 1]];
        unit=[0 0 0];
        X=MakeX(a,unit,l,d);

        neigh[:,1]=armod(X[:,1]+2,l)+l*X[:,2]+X[:,3]*l^2;
        neigh[:,2]=armod(X[:,1],l)+l*X[:,2]+X[:,3]*l^2;
        neigh[:,3]=X[:,1]+1+l*mod(X[:,2]+2,l)+X[:,3]*l^2;
        neigh[:,4]=X[:,1]+1+l*mod(X[:,2],l)+X[:,3]*l^2
        neigh[:,5]=X[:,1]+1+l*X[:,2]+mod(X[:,3]+2,l)*l^2;
        neigh[:,6]=X[:,1]+1+l*X[:,2]+mod(X[:,3],l)*l^2;

        return Lattice("Cube",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeCheckerboard(l::Int)
      N::UInt16=l^2;
      d::Int8=2;
      X=Array{Int16}(N,2);
      nnei::UInt8=6;
      neigh=Array{Int16}(N,6);

      X[:,1]=armod(collect(1:N),l);
      X[:,2]=ceil(collect(1:(N))/l);
      a=[[1 0]
          [0 1]];
      unit=[0 0];

      neigh[:,1]=armod(X[:,1]+1,l)+l*(X[:,2]-1);
      neigh[:,2]=armod(X[:,1]-1,l)+l*(X[:,2]-1);
      neigh[:,3]=X[:,1]+mod(X[:,2],l)*l;
      neigh[:,4]=X[:,1]+mod(X[:,2]+l-2,l)*l;

      neigh[:,5]=armod(X[:,1]+1,l)+mod(X[:,2]+l-1+parity(X[:,1]).*parity(X[:,2]),l)*l;
      neigh[:,6]=armod(X[:,1]-1,l)+mod(X[:,2]+l-1-parity(X[:,1]).*parity(X[:,2]),l)*l;

      return Lattice("Checkerboard",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeChain(l::Int)
        N::UInt16=l
        d::UInt8=1
        nnei::UInt8=2
        neigh=Array{Int8}(N,2)
        a=[1];
        unit=[0];

        X=collect(1:l);

        neigh[:,1]=armod(collect(2:(l+1)),l);
        neigh[:,2]=armod(collect(l:(2*l-1)),l);

        return Lattice("Chain",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeTriangular(l::Int)
        N::UInt16=l^2
        d=2
        nnei::UInt8=6
        neigh=Array{Int8}(N,6)
        a=[[1 0]
           [cos(π/3) sin(π/3)]];
        unit=[0 0];

        X=MakeX(a,unit,l,d);

        for ii in 1:N
            column=armod(ii,l)
            row=floor((ii-1)/l)

            neigh[ii,1]=armod(column+1,l)+l*row;
            neigh[ii,2]=armod(column-1,l)+l*row;
            neigh[ii,3]=column+l*mod(row+1,l);
            neigh[ii,4]=armod(column-1,l)+l*mod(row+1,l);
            neigh[ii,5]=column+l*mod(row-1,l);
            neigh[ii,6]=armod(column+1,l)+l*mod(row-1,l);
        end

        return Lattice("Triangular",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeHoneycomb(l::Int)
        d::Int8=2;
        nnei::UInt8=3;
        N::UInt16=2*l^2;
        a=[[2*cos(pi/6) 0]
            [cos(pi/6) 1+sin(pi/6)]];
        unit=[[0 0]
                [cos(π/6) sin(π/6)]];

        X=MakeX(a,unit,l,d);

        neigh=Array{Int16}(N,3);
        for i in 1:l
            neigh[2*i-1,1]=2*i;
            neigh[2*i,1]=2*i-1;

            neigh[2*i+1,2]=2*i;
            neigh[2*i,2]=2*i+1;
        end

        for j in 1:l
            for i in 1:l
                neigh[2*i-1+2*l*(j-1),3]=armod(2*i+2*l*(j-2)+N,N);
                neigh[2*i+2*l*(j-1),3]=armod(2*i+2*l*(j)+N-1,N);
            end
        end

        neigh[1,2]=2*l;
        neigh[2*l,2]=1;

        for j in 2:l
            neigh[(2*l*(j-1)+1):(2*l*j),1:2]=neigh[1:(2*l),1:2]+2*l*(j-1);
        end

        return Lattice("Honeycomb",l,d,a,unit,N,X,nnei,neigh);
    end

    function MakeKagome(l::Int)
      d::Int8=2;
      nnei::UInt8=4;
      N::UInt16=3*l^2;
      a=[[2 0]
         [1 sqrt(3)]];
      unit=[[0 0]
            [1 0]
            [.5 .5*sqrt(3)]];

      X=MakeX(a,unit,l,d);

      ind=collect(1:N);
      ti=armod(ind,3);
      tri=Array{Int64}(N);
      for ii in 1:N
        tri[ii]=convert(Int64,floor((ind[ii]-1)/3))
      end

      neigh=Array{Int16}(N,4);
      neigh[:,1]=mod(ti+1,3)+3*tri+1;
      neigh[:,2]=mod(ti+2,3)+3*tri+1;
      neigh[:,3]=armod(armod(ti+2,3)+(ti-2)*3*l+3*(mod(ti,3)-1)+tri*3,N);
      neigh[:,4]=armod(armod(ti+1,3)+(1-mod(ti,3))*3*l+3-3*mod(ti+1,3)+tri*3,N);

      return Lattice("Kagome",l,d,a,unit,N,X,nnei,neigh);
    end

    function Bravis2D(l::Int,a1::Float64,a2::Float64,ϕ::Float64)
          if (a1!=a2)
              if (ϕ!=π/2)
                  if ( isapprox(2*a1*cos(ϕ),a2) )
                    name="rectangular"
                  else
                    name="oblique"
                  end
              else
                name="rectangular"
              end
          elseif (isapprox(ϕ,2π/3))
              name="hexagonal"
          elseif (isapprox(ϕ,π/2))
              name="square"
          end


          N::UInt16=l^2;
          d::Int8=2;
          nnei::UInt8=4;
          neigh=Array{Int16}(N,4);

          a=[[a1 0]
              [a2*cos(ϕ) a2*sin(ϕ)]];
          unit=[0 0];
          X=MakeX(a,unit,l,d);

          ind=collect(1:N);
          y=floor(ind/l);
          neigh[:,1]=armod(ind+1,l)+l*(X[:,2]-;
          neigh[:,2]=armod(ind-1,l)+l*(X[:,2]-1);
          neigh[:,3]=armod(ind,l)+mod(X[:,2],l)*l;
          neigh[:,4]=armod(ind,l)+mod(X[:,2]+l-2,l)*l;

          return Lattice("Square",l,d,a,unit,N,X,nnei,neigh);
    end

    function Bravais3D(l::Int,a1::Float64,a2::Float64,a3::Float64,α::Float64,β::Float64,γ::Float64)

        N::UInt16=l^3;
        d::Int8=3;
        nnei=6;
        neigh=Array{Int16}(N,6);

        a=[[a1 0 0]
           [a2*cos(α) 0 a2*sin(α)]
           []]




    end

    function MakeX(a::Array,unit::Array,l::Int,d::Int8)
        ncell=size(unit)[1];
        N=ncell*l^d;

        a1=repeat(a[1,:],outer=[ncell,1]);
        a2=repeat(a[2,:],outer=[ncell*l,1]);
        if d==3
            a3=repeat(a[3,:],outer=[ncell*l^2,1]);
        end


        X=Array{Float64}(N,d);
        # Here we are actually calculating the positions for every site
        for i in 1:l    #for the first row
            X[ncell*i-ncell+1:ncell*i,:]=unit+(i-1)*a1;
        end

        for j in 2:l    #copying the first row into the first layer
            X[(ncell*l*(j-1)+1):(ncell*l*j),:]=X[1:ncell*l,:]+(j-1)*a2;
        end

        if d==3
            for k in 2:l    #copying the first row into the first layer
                X[(ncell*l^2*(k-1)+1):(ncell*l^2*k),:]=X[1:ncell*l^2,:]+(k-1)*a3;
            end
        end


        return X
    end



    function PlotNeighbors(lt::Lattice)
        fig=gcf()
        if(lt.dim==2)
            for i in 1:lt.N
                for j in 1:lt.nnei
                    xx=[lt.X[i,1], lt.X[lt.neigh[i,j],1] ]
                    yx=[lt.X[i,2], lt.X[lt.neigh[i,j],2] ]
                    plot(xx,yx)
                end
            end
        elseif(lt.dim==3)
            for i in 1:lt.N
                for j in 1:lt.nnei
                    xx=[lt.X[i,1], lt.X[lt.neigh[i,j],1] ]
                    yx=[lt.X[i,2], lt.X[lt.neigh[i,j],2] ]
                    zx=[lt.X[i,3], lt.X[lt.neigh[i,j],3] ]
                    plot3D(xx,yx)
                end
            end
        else
            println("Dimension not 2 or 3")
            println("Not plotting")
        end
    end

end
