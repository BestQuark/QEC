import syndrome_v2 as syn

g1 = "X Y X 1 X Y X 1 Y Z Y 1 1 1 1 1 Y Z Y 1;"
g2 = "Y Z Y 1 Y Z Y 1 Z X Z 1 1 1 1 1 Z X Z 1;"
g3 = "Y Z Y 1 X Y X 1 X Y X 1 Y Z 1 1 1 1 1 1;"
g4 = "Z X Z 1 Y Z Y 1 Y Z Y 1 Z X Z 1 1 1 1 1;"
g5 = "Z Y Y Z 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;"
g6 = "X Z Z X 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;"
g7 = "1 1 1 1 Z Y Y Z 1 1 1 1 1 1 1 1 1 1 1 1;"
g8 = "1 1 1 1 X Z Z X 1 1 1 1 1 1 1 1 1 1 1 1;"
g9 = "1 1 1 1 1 1 1 1 Z Y Y Z 1 1 1 1 1 1 1 1;"
g10= "1 1 1 1 1 1 1 1 X Z Z X 1 1 1 1 1 1 1 1;"
g11= "1 1 1 1 1 1 1 1 1 1 1 1 Z Y Y Z 1 1 1 1;"
g12= "1 1 1 1 1 1 1 1 1 1 1 1 X Z Z X 1 1 1 1;"
g13= "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 Z Y Y Z;"
g14= "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 X Z Z X"
# g15= "1 1 Z X 1 1 Z X X X X X X 1 Y Y X 1 Y Y;"
# g16= "1 X Z Z X 1 X Z X 1 X Z 1 X Z Z 1 1 1 1;"
# g17= "X 1 Y Y X X X X 1 1 Z X 1 1 Z X X 1 Y Y;"
# g18= "Y Y X Z Y Z Y 1 Y Y X Z Z Z Z Z Y Z Y 1;"
# g19= "X X X X X X X X X X X X X X X X X X X X;"
# g20= "Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z"

stabs = g1+g2+g3+g4+g5+g6+g7+g8+g9+g10+g11+g12+g13+g14 #+g15+g16+g17+g18+g19+g20

errorsHolography = syn.classify_errors(stabs)
syn.save_syndromes("holographic", errors)