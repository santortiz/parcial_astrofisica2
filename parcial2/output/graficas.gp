

lim = 4000
retardo = 0.001
f = 0.01

set size square
set xrange [-5e-5/2.0:1e-3+5e-5/2.0]
set yrange [-5e-5/2.0:1e-3+5e-5/2.0]

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("state_%.4d",i)
    plot file every ::0::1599 u 2:3 w p ps 1 pt 7 lc rgb "blue" not,\
	 "" every ::1600::1680 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1681::1759 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1760::1840 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1841::1919 u 2:3 w lp ps 1 pt 7 lc rgb "black" not
    
  pause retardo
}

do for [i=0:lim]{
  titulo = sprintf("paso = %.4d - tiempo = %.5f",i,i*5e-5)  
  set title titulo
  file = sprintf("state_%.4d",i)
    plot file every ::0::1599 u 2:3:($4*f):($5*f) w vectors lc rgb "blue" not,\
	 "" every ::1600::1680 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1681::1759 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1760::1840 u 2:3 w lp ps 1 pt 7 lc rgb "black" not,\
	 "" every ::1841::1919 u 2:3 w lp ps 1 pt 7 lc rgb "black" not
    
  pause retardo
}