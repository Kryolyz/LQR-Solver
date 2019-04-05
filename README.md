# LQR-Solver
A program that creates feedback-controllers using a linear-quadratic cost function. The function is solved by turning it into a infinite-horizon continous riccati equation and using Newton's method. It's just the solving algorithm though, there is no way to directly change the input state-space data or weights matrices without changing and recompiling the code (there is no GUI), which means the code is nothing but an example of how to realize an algorithm for LQR.

# LQR-I
More accurately it's LQR-I because integral action to prevent steady state errors is included, but I consider that an absolute necessity for a functional controller and added it right after I was done.

# Efficiency and stability
Keep in mind that the code is basically the most simple way to implement such a solver. I used the Kronecker-Product to solve the Lyapunov-equation resulting from Newton's method, which is EXTREMELY inefficient and even numerically unstable, not to mention Newton's method isn't exactly the safest solver either. I highly recommend looking into the Bartels-Stewart algorithm in case you are interested in coding a similar solver. Its almost just as simple, much more efficient and safer. To replace Newton's method I recommend looking into optimization algorithms in general since the correct choice often depends on the problem at hand. For example a problem that one expects to have many local minima may best be solved using particle swarm optimization or other methods that make use of multiple pathfinders to weed out suboptimal solutions.

