# Numerical-Python-Solving-Laplace-equation
In this project we try to implement Finite Difference method in iterative way to solve Laplace equation. The first derivative can be approximated as

![CodeCogsEqn](https://user-images.githubusercontent.com/97800241/181074646-990baee6-4fd8-4937-b3dc-0598396549f3.gif)

similarly

![CodeCogsEqn (1)](https://user-images.githubusercontent.com/97800241/181074962-e25d84a4-8481-4ef2-ae6e-f1114b7c07f8.gif)

Now Laplace equation:

![CodeCogsEqn (2)](https://user-images.githubusercontent.com/97800241/181075458-a3cee75f-7627-43c3-8939-3e9f7f3a5072.gif)

writing using finite difference method :

![CodeCogsEqn (3)](https://user-images.githubusercontent.com/97800241/181076726-55988399-747e-499f-bee6-0c51e1a4616a.gif)

We have been given Boundary conditions as:

![CodeCogsEqn (4)](https://user-images.githubusercontent.com/97800241/181077830-19726ff6-c41f-47c1-8edd-a7916cb201fd.gif)

we have to use this equation to get values of $\phi$.

We initially take a matrix for $\phi$ of dimension 11x11 of values 0. which divides the unit square bounded by x=0,x=1,y=0,y=1 into a total of 20 squares of edge 0.1.

We then update the matrix as per Boundary condition. Then used the algorith iteratively to get the desired non zero matrix which keeps the values of function $\phi(x,y)=x^2-y^2$ in 0<x<1, 0<y<1.
