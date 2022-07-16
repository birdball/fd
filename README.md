# fd
A Rust implementation of the pricer explained in this paper https://www.researchgate.net/publication/315615243_Excel_implementation_of_finite_difference_methods_for_option_pricing


## fd - "Finite Differences"
Financial derivatives are often evaluated with numerical analysis. A standard way of assessing their value is via the finite difference method. The paper I rewrote this method from gives a solution for an 8x10 grid that evaluates option values with relatively large gaps between possible underlying prices. It's quite difficult to modify the model used in the paper to account for longer timespans or finer grained underlying prices, without making Excel sweat. This Rust implementation is far more flexible.

## How it works
Our target is assessing option prices today based on what they may be worth at payoff. To do that, we estimate boundaries on the top and rightmost parts of the grid and work backwards to get values today using matrix multiplication.

## Todo
Missing comparison to analytic solution
Missing american payoff
Make more command line friendly
