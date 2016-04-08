module MyColors
export color, ColorSet, Minerals

  type color
    name::ASCIIString
    n::Int
    hex::ASCIIString
    RGB::Array{Float64
  end

  type ColorSet
    name::ASCIIString
    parts::Array{color}
  end

  function HextoRGB(hexi::ASCIIString)
    R=parse(Int,hexi[2:3],16)/255.;
    G=parse(Int,hexi[4:5],16)/255.;
    B=parse(Int,hexi[6:7],16)/255.;

    RGB=[R,G,B];

    return RGB
  end

  function Minerals()
    parts=Array{color}(6);

    greyblue="#d85bb75";
    color1=color("greyblue",greyblue,HextoRGB(greyblue));
    parts[1]=color1;

    brown="#9c8d85";
    color2=color("brown",brown,HextoRGB(brown));
    parts[2]=color2;

    darkpurple="#41364a";
    color3=color("darkpurple",darkpurple,HextoRGB(darkpurple));
    parts[3]=color3;

    purple="#7a5b75";
    color4=color("purple",purple,HextoRGB(purple));
    parts[4]=color4;

    bisque="#eb8c84";
    color5=color("bisque",bisque,HextoRGB(bisque));
    parts[5]=color5;

    tan="#edd4b6";
    color6=color("tan",tan,HextoRGB(tan));
    parts[6]=color6;

    Minerals=ColorSet("Minerals",6,parts);
  end

end
