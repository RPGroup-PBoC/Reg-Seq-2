N=9
(
#for gc in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 28 29 30 31 32 33 34 35 36 37 38 39 40 41; do
for gc in 41; do 	
   ((i=i%N)); ((i++==0)) && wait
   julia hypothesis_growth_condition.jl "$gc" & 
done
)
