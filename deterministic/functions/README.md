# Functions

This folder contains all the code needed to run the deterministic simulations. 

## The `main.py` file

Most of the simulations are launched from the file `main.py`. All the parameters can be adjusted at the top of this file and the value of the variable `what_to_do` determine the type simulation ran:

1) `what_to_do = "evolution"` simulates system (4) or (5) in article <https://doi.org/10.1007/s00285-023-01926-4> depending on the conversion timing chosen. It plots screenshots of the traveling wave at different times in 1D space. The function used to simulate the system are written in `evolution.py` and the functions used to plot the graph are written in the `graph.py`.

2) `what_to_do = "evolution 2D"` simulates system (4) or (5) in article <https://doi.org/10.1007/s00285-023-01926-4> depending on the conversion timing chosen. It plots screenshots of the traveling wave at different times in 2D space. The function used to simulate the system are written in `evolution.py` and the functions used to plot the graph are written in the `graph.py`.

3) `what_to_do = "tanaka fraction"` simulates equation (13) in article <https://doi.org/10.1007/s00285-023-01926-4> or equivalently equation (5) with reaction term (6) in article <https://doi.org/10.1073/pnas.1705868114>. It plots screenshots of the traveling wave at different times in 1D space. The functions used to simulate and plot the graphs are written in `tanaka.py`.

4) `what_to_do = "tanaka cubic"` simulates equation (14) in article <https://doi.org/10.1007/s00285-023-01926-4> or equivalently equation (5) with reaction term (7) in article <https://doi.org/10.1073/pnas.1705868114>. It plots screenshots of the traveling wave at different times in 1D space. The functions used to simulate and plot the graphs are written in `tanaka.py`.

5) `what_to_do = "KPP"` simulates the following equation: partial_t p - partial_xx p = p * (1-p) * (1-2*s) where p is the proportion of drive alleles and s the drive fitness disadvantage. It plots screenshots of the traveling wave at different times in 1D space. The functions used to simulate and plot the graphs are written in `tanaka.py`.

6) `what_to_do = "pulled pushed"` illustrates the concept of pulled and pushed waves. It colors two neutral populations in the wave, one at the very front and the other containing all the individuals behind. The offspring of each population is colored the same way to show the affiliation among the generations. Once again, we plot screenshots of the wave at different times in 1D space. The functions used to simulate and plot the graphs are written in `pulled_pushed_wave.py`.

7) `what_to_do = "speed function of time"` compares the speed of the traveling wave function of time for cases 1) 3) 4) and 5).  

8) `what_to_do = "speed function of s"` compares the speed of the traveling wave function of s (the drive fitness disadvantage) for cases 1) 3) 4) and 5).  

9) `what_to_do = "speed function of r"` compares the speed of the traveling wave function of r (the intrinsic growth rate) for cases 1) 3) 4) and 5).  

10) `what_to_do = "heatmap"` plot an heatmap representing the speed of the traveling wave function of two variables denoted `x` and `y` in the code. The possibilities are: i) `x = 's'` and `y == 'r'`, ii) `x = 'h'` and `y == 'r'` and iii) `x = 'h'` and `y == 's'`. The functions used to simulate and plot the heatmap are written in `heatmap.py`.


## Other files `.py` that can be run independently

`barton.py` illustrates the comparison between the minimum width of the spatial step stopping the wave, obtained by Barton (1979) and ours for different values of alpha. This figure is presented only in the introduction of my PhD thesis (Figure 1.32). 

`illustrations_bio.py` plots various figures from the soon submitted article "The speed of advance of a gene drive is
affected by density dependence" . If `what_to_do = "density_nd"`, the code plots the final drive densities, if `what_to_do = "ext_bist_line"` it plots the persistence/bistability/eradication lines.

`fig_rode_debarre_2019.py` plots Figures 13 and 14 from article <https://doi.org/10.1007/s00285-023-01926-4>.

