pause on
use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\Invest1993", clear

gen I = inva
gen Q = vala
lpoly I Q if Q<=5, ci nosc title("Nadaraya-Watson 95p CI")
graph export "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\pings\9p1.png", as(png) name("Graph")
lpoly I Q if Q<=5, degree(1) ci nosc title("LL 95p CI")
graph export "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS2\pings\9p2.png", as(png) name("Graph")