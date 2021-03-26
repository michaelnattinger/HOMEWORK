pause on
use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\FRED-QD", clear

gen y = 100*((gdpc1/L.gdpc1)^4-1)
gen ylag = L.y
lpoly y ylag, ci nosc title("NW Estimator")
graph export "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\pings\11p1.png", as(png) name("Graph")
lpoly y ylag, degree(1) ci nosc title("LL Estimator")
graph export "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\pings\11p2.png", as(png) name("Graph")

