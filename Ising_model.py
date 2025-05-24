import scipy.stats as stat
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, Slider
from collections.abc import Callable
import argparse

class grid(object):
    def __init__(self, Nx, Ny):
        """
        Nx: Number of nodes along the x-axis. Points toward the right.
        Ny: Number of nodes along the y-axis. Points toward the bottom.
        """
        self.Nx = Nx
        self.Ny = Ny

        #This grid stores the spin' values at all nodes.
        self.grid = None

    def initialize_grid(self, mode):
        if mode == "random":
            self.grid = 2*(stat.bernoulli.rvs(1/2, size = (self.Nx, self.Ny)) - 1/2)
        elif mode == "+":
            self.grid = np.ones((self.Nx, self.Ny))
        elif mode == "-":
            self.grid = -np.ones((self.Nx, self.Ny))
        elif mode == "half":
            self.grid = np.empty((self.Nx, self.Ny))
            self.grid[:self.Nx//2, :] = 1
            self.grid[self.Nx//2:, :] = -1
        else:
            print("Initialization mode should be 'random', '+', '-' or 'half'.")
            return 1
        return 0

    def is_border(self, i, j):
        if i == 0 or j == 0:
            return True
        elif i == self.Nx - 1 or j == self.Ny - 1:
            return True
        return False

class Ising(object):
    def __init__(self, Nx: int, Ny: int, B: float | Callable, T: float, bound_cond: str, uniform_mag: bool = True):
        """
        Args:
            B: float | Callable. Exterior magnetic field in Tesla. If B is callable, then it should take as argument the node's position (i,j) at least. If B is a float, then
            it is uniform over the whole lattice.
            T: float. Temperature of the thermal bath in Kelvin.
            Nx: int. Number of nodes along x-axis
            Ny: int. Number of nodes along y-axis
            bound_cond: float. Boundary conditions. Should be "isolated" or "periodic".
            uniform_mag: bool. States if the magnetic field is uniform over the lattice.
        """
        self.Nx = Nx
        self.Ny = Ny
        if uniform_mag and not isinstance(B, float) and not isinstance(B, int):
            raise ValueError(f"if Magnetic field is Uniform, B should be a float or an int but is instead {type(B)}.")
        elif not uniform_mag and isinstance(B, Callable):
            raise ValueError(f"If magnetic field is not uniform, B should be callable but is instead {type(B)}.")
        self.B = B
        self.T = T
        self.k = 1.38e-23
        self.g = 2
        self.mu_B = 1e-23
        self.J = 4e-21
        self.grid = grid(self.Nx, self.Ny)
        self.bound_cond = bound_cond

        self.uniform_mag = uniform_mag

        self.beta = 1/(self.k*self.T)
        self.Set = {-4*self.beta*self.J: np.exp(-4*self.beta*self.J), -8*self.beta*self.J: np.exp(-8*self.beta*self.J)}

    
    def initialize_Ising(self, mode):
        """
        Initialize the grid.
        """
        self.grid.initialize_grid(mode)
    
    def B_at_node(self, i: int, j: int):
        if self.uniform_mag:
            return self.B
        else:
            return self.B(i,j)

    def ratio(self, Energy: float, *args):
        """
        Compute the configuration probability associated to this energy.
        Args: 
            Energy (in J). It should be the difference of energy between the opposite and actual spin direction at one node. 
            If the difference of energy is negative, you accept the new position with 100% probability. 
        Returns:
            Compute the ratio of the Metropolis-Hasting steps.
        """
        if Energy >= 0:
            return 1
        else:
            if not uniform_mag or (uniform_mag and self.B != 0):
                return np.exp(-Energy/(self.k*self.T))
            else:
                return self.Set[Energy*self.beta]

    def get_ppv(self, bound_cond: float, i: int, j: int):
        """
        Return the list of nearest neighbours of node (i,j). Results depend on the boundary conditions.
        """
        if bound_cond == "isolated":
            if self.grid.is_border(i, j):
                if i == 0:
                    ppv = [(i + 1, 0)]
                    if j == 0:
                        ppv += [(i, j + 1)]
                    elif j == self.Nx - 1:
                        ppv += [(i, j - 1)]
                    else:
                        ppv += [(i, 1), (i, j - 1)]
                elif i == self.Nx - 1:
                    ppv = [(i - 1, j)]
                    if j == 0:
                        ppv += [(i, j + 1)]
                    elif j == self.Ny - 1:
                        ppv += [(i, j - 1)]
                    else:
                        ppv += [(i, j + 1), (i, j - 1)]
                elif j == 0:
                    ppv = [(i - 1, j), (i + 1, j), (i, j + 1)]
                elif j == self.Ny - 1:
                    ppv = [(i - 1, j), (i + 1, j), (i, j - 1)]
            else:
                ppv = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
        elif bound_cond == "periodic":
            ppv = [((i - 1)%self.Nx, j), (i, (j - 1)%self.Ny), ((i + 1)%self.Nx, j), (i, (j + 1)%self.Ny)]
        return ppv

    def Hamiltonian(self, i: int, j: int):
        assert i >= 0 and i < self.Nx
        assert j >= 0 and j < self.Ny
        H = - self.g*self.mu_B*self.grid.grid[i,j]*self.B_at_node(i,j)
        ppv = self.get_ppv(self.bound_cond, i, j)
        ppv = np.array(ppv)
        H -= self.J*self.grid.grid[i,j]*np.sum(self.grid.grid[ppv[:,0], ppv[:,1]])
        return H

    def evolve(self):
        """
        Flip one spin at a time.
        """
        i, j = np.random.randint(0, self.Nx), np.random.randint(0, self.Ny)
        E = self.Hamiltonian(i, j)
        ratio = self.ratio(2*E)
        Flip = np.random.rand() < ratio
        return Flip, (i, j)

if __name__ == '__main__':

    #I should add a parser
    parser = argparse.ArgumentParser()

    parser.add_argument("--Lx", default=50, action='store', help="Length of the lattice along the x-axis")
    parser.add_argument("--Ly", default=50, action='store', help="Length of the lattice along the y-axis")
    parser.add_argument("-T", default=300, action='store', help="Temperature of the thermal bath")
    parser.add_argument("-B", default=0, action='store', help="Temperature of the thermal bath")
    parser.add_argument("-J", default=4e-21, action='store', help="Exchange interaction of the system")

    parser.add_argument("--bound_cond", default="periodic", action='store', help="Periodic condition", choices=["periodic", "isolated"])
    parser.add_argument("--Initialization", default="random", action='store', help="Initialize the lattice", choices=["uniform", "+", "-", "half"])
    parser.add_argument("--uniform_mag",default=True, action='store', help="Is the magnetic field uniform ?")
    parser.add_argument("--cmap",default="viridis", action='store', help="Is the magnetic field uniform ?")

    args = parser.parse_args()

    #Initial parameters
    B = args.B
    T = args.T #Room temperature
    uniform_mag = args.uniform_mag

    I = Ising(Nx = args.Lx, Ny = args.Ly, B = B, T = T, bound_cond=args.bound_cond, uniform_mag=uniform_mag)
    I.initialize_Ising(mode=args.Initialization)

    fig, axes = plt.subplots(1, 2, figsize = (12, 6))

    im = axes[0].imshow(I.grid.grid, animated=True, cmap = args.cmap, vmin = -1/2, vmax = 1./2)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_title("Lattice")

    # Make a vertically oriented slider to control the Temperature
    axamp = fig.add_axes(rect=[0.05, 0.25, 0.0225, 0.33])
    Temp_slider = Slider(
        ax=axamp,
        label="Temperature [°Kelvin]",
        valmin=1e1,
        valmax=2e3,
        valinit=T,
        orientation="vertical"
    )

    def update_Temp(val):
        I.T = Temp_slider.val
        I.beta = 1/(I.k*I.T)
        I.Set = {-4*I.beta*I.J: np.exp(-4*I.beta*I.J), -8*I.beta*I.J: np.exp(-8*I.beta*I.J)}

    Temp_slider.on_changed(update_Temp)

    if uniform_mag == True:
        # Make a horizontal slider to control the magnitude of magnetic field.
        axmag = fig.add_axes(rect=(0.18, 0.05, 0.25, 0.03))
        mag_slider = Slider(
            ax=axmag,
            label='B [Tesla]',
            valmin=-10,
            valmax=10,
            valinit=B,
        )

        ##Callback functions of sliders
        def update_mag(val):
            I.B = mag_slider.val

        # register the update function with each slider
        mag_slider.on_changed(update_mag)

    line, = axes[1].plot([], [], 'b-', lw=2, label = "At an instant")
    axes[1].set_xlim(0, 10)
    axes[1].set_ylim(-100, 100)
    axes[1].set_title("Magnetization")

    line2, = axes[1].plot([], [], 'r-', lw=2, label = "On average")
    axes[1].legend()

    i = 0
    j = 0

    history = []
    history_mean = []

    def updatefig(*args):
        #Update the lattice
        Flip, pos = I.evolve()
        if Flip:
            I.grid.grid[pos[0], pos[1]] *= -1
        im.set_array(I.grid.grid)

        #Update the magnetization plot at time t
        if len(history) >= 1000:
            history.pop(0) #Keep the last hundreds elements
        history.append(np.sum(I.grid.grid))
        t = np.linspace(0,len(history), len(history))
        line.set_data(t, history)

        #Update the average magnetization
        if len(history_mean) >= 1000:
            history_mean.pop(0)
        history_mean.append(np.mean(history))
        t = np.linspace(0,len(history_mean), len(history_mean))
        line2.set_data(t, history_mean)
        return im, line, line2

    #Declare an animation
    ani = FuncAnimation(fig=fig, func=updatefig,  blit=True, interval = 1, cache_frame_data=False)

    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')


    def reset(event):
        mag_slider.reset()
        Temp_slider.reset()
    button.on_clicked(reset)

    plt.show()