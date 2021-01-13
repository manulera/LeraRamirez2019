import numpy as np
import math
import random
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# This was done before I was really comfortable with object-oriented programming, so there is definitely room for
# improving it in that sense, especially to make it the code more clear

# Function that calculates the deterministic solution, and returns "solution", an array of size(len(t),len(init)),
# containing the values of the concentrations of the chemical species in init for every time point in t
# This function is called by the function ReAct, look at it to understand the arguments. In fact, it's better looking by
# that function
def DetSol(elements, init, t, Gamma, k, decay):
    # Function to calculate V, a np.array that contains for every reaction Ri, the reaction velocity Vi
    def calcV(init, Gamma, k):

        # Create an empty list where we will append each Vi
        V = list()

        for i in range(0, Gamma.shape[1]):  # Iterate along the reactions
            # Initialize Vi
            Vi = 1
            Ri = Gamma[:, i]

            for j in range(0, Gamma.shape[0]):  # Iterate along the chemical species in the stoichiometry matrix

                # If the current species,Xj, is a reactant, the velocity of the reaction will be proportional to the
                # concentration of Xj, init[j], to the power of abs(Ri[j]), where abs(Ri[j]) is the stoichiometry index
                # of Xj in Ri
                if Ri[j] < 0:
                    Vi *= (init[j] ** abs(Ri[j]))

                # If the current species,Xj, is a catalizer/enzyme (explained in the notebook), basically the same as
                # if it was a reactant
                if decay:
                    for ele in decay:
                        if (ele[0], ele[1]) == (j, i):
                            Vi *= (init[j] ** ele[2])

            # Multiply Vi times the corresponding reaction constant, and append it to the list V
            V.append(Vi * k[i])
        return V

    # ODE system (in format used for odeint), this function calculates V and multiplies it by the stoichiometry matrix
    # to calculate "rates", a np.array which contains the values of d(init/dt)
    # t is not used inside the function, but is required as an input for the calculations when we call myODE with odeint
    def myODE(init, t, Gamma, k):

        V = calcV(init, Gamma, k)
        rates = np.dot(Gamma, V)
        return tuple(rates)

    # Solve using odeint
    solution = odeint(myODE, init, t, args=(Gamma, k,))
    # "solution" is a np.array containing the concentration values of all species at the different time points specified
    # in "t"
    return solution


# Performs several rounds of gillespie calculations, the number of rounds can be indicated with the parameter "rounds".
# The output is several np.arrays.
# "tgill_all" is a np.array of length "rounds" containing arrays with the times where Gillespie reaction steps happened,
# for every round of Gillespie.
# "valsgill_all" is a np.array of length "rounds" that contains arrays of size(len(init),len(tgill)). They store
# the quantities of all chemical species in the system at the times in "tgill_all" for every round of Gillespie.
# "mus_all" is a np.array of length "rounds" that contains arrays of length len(tgill) containing the value of mu (as
# named by Gillespie, the index of the reaction that happened at every step, this can be useful)
# "taus_all" is a np.array of length "rounds" that contains arrays of length len(tgill) containing the value of tau (as
# named by Gillespie, the index of the reaction that happened at every step, this can be useful)


def Gillespy(elements, init, t, Gamma, k, decay, rounds=0):
    # Returns the number of possible collisions * index!(factorial), we will see why
    def Comby(quant, index):
        """
        Input is quant (number of molecules) and index (the reaction index for that particular molecule in the reaction).
        A very interesting question is why this function will always return an int, think about it! rounds
        """
        result = 1
        for i in range(index):
            result *= (quant - i)
            # result /= (i + 1)
            # The line above would be needed if the kinetic constant did not have this information
            # inside already (we will speak about this in the class, what are the units of k in a chemical reaction?)
        return result

    # Function to calculate P, a np.array that contains the cumulative sum of a vector Px,
    # Px contains for every reaction Ri, the reaction propensity Pi
    def calcP(init, Gamma, k):
        """
        Will calculate the cummulative sum of propensities for all reactions
        """
        Px = list()

        for i in range(0, Gamma.shape[1]):  # Iterate along the reactions

            Pi = 1
            Ri = Gamma[:, i]
            for j in range(0, Gamma.shape[0]):  # Iterate along the chemical species in the system
                # This is to solve some very particular situations in which one of the species acts both as a "normal"
                # reactant, and as a catalizer, for instance: 2A + B -> C + A
                counted = 0
                if decay:
                    for ele in decay:
                        if (ele[0], ele[1]) == (j, i):
                            # In case one of the elements is catalizer and also "normal reactant. Probably it's somehow
                            # redundant, one should check it.
                            Pi *= Comby(init[j], ele[2] - Ri[j] * (Ri[j] < 0))
                            counted = 1
                if not counted and Ri[j] < 0:
                    # Ri[j] is the stochiometry index of the element j
                    # init[j] is the number of molecules of element j
                    Pi *= Comby(init[j], -Ri[j])
                    counted = -Ri[j]

            Px.append(Pi * k[i])
        return np.cumsum(Px)

    # Perform a step of gillespie algorithm
    def gilleStep(init, tcount, Gamma, k, tend):

        # This method is not optimized, we do not need to recalculate all the propensities, only the ones that are
        # affected by the previous reaction
        P = calcP(init, Gamma, k)

        # If all propensities are equal to 0, the system is dead!
        # Sometimes, you can reach a situation where the propensity of all reactions is 0; stop
        # your algorithm then, otherwise it will return an error when it tries to divide by P
        if P[-1] == 0:
            return (init, tend, None, None)
            # The system is dead! Sometimes, you can reach a situation where the propensity of all reactions is 0,
            # stop your algorithm then, otherwise it will return an error when it tries to divide by P
        else:
            # P[-1] is what gillespie calls a_0 in the paper
            # Calculate the time step, as described in the paper
            tau_step = (math.log(1) - math.log(random.random())) / (P[-1])

            # Find the first non-false element in the logical array resulting from P > random.random() * P[-1]
            # This will select one of the reactions according to their propensities.
            mu = next(i for i, x in enumerate(P > random.random() * P[-1]) if x)

            # For gillespie, one round of the reaction results in the addition of Gamma[:, mu] (the column of
            # the stoichiometry matrix corresponding to the selected reaction mu) to the init values
            init += Gamma[:, mu]
            # We add to the counter of the time the tau_step
            tcount += tau_step
            return (init, tcount, mu, tau_step)

    # These lists will contain as many elements as indicated by the number "rounds", corresponding to several Gillespie
    # calculations.
    tgill_all = list()
    valsgill_all = list()
    old_init = init
    mus_all = list()
    taus_all = list()
    for i in range(rounds):
        init = old_init
        # tcount increases by tau_step every time gilleStep is called
        tcount = 0

        # When tcount reaches the time limit indicated by the user, stop the algorithm
        tend = t[-1]

        # tgill stores all the values of tcount, so we can plot the progression of the system in time afterwards
        tgill = [0]

        # valsgill is a np.array of size(len(init),numsteps), where numsteps is the number of gillespie steps performed.
        # It stores the quantities of all chemical species in the system at the times in tgill
        # It increases size every time, which is not very efficient, but if we want the calculation to stop when we reach a
        # certain value of tcount, we cannot know how many steps we will make, so we cannot prelocate the memory
        valsgill = np.array(init)
        mus = list()
        taus = list()
        while tcount < tend:
            (init, tcount, mu, tau) = gilleStep(init, tcount, Gamma, k, tend)
            # Append the values to valsgill
            valsgill = np.c_[valsgill, np.array(init)]

            tgill.append(tcount)
            mus.append(mu)
            taus.append(tau)

        taus_all.append(np.array(taus))
        mus_all.append(np.array(mus))
        valsgill_all.append(valsgill)
        tgill_all.append(tgill)

    return np.array(tgill_all), valsgill_all, mus_all, taus_all


# This can also be used to solve a linear system in which all reactions are A becomes B with a constant rate.
def LinearSystem(Gamma,k):

    system = np.zeros([Gamma.shape[0],Gamma.shape[0]])
    for s,species in enumerate(Gamma):
        for r,reaction in enumerate(species):
            if reaction==-1:
                system[s,s]-= k[r]
            elif reaction==1:
                # Look for the positive value
                s2 = np.where(Gamma[:,r]==-1)[0]
                system[s,s2] += k[r]
    return system




# "user_input" is a list of pairs Xi, amount of Xi. Xi being a string with the number of a chemical species. All the
# chemical species that appear in the reactions must be declared in this list, and given an initial concentration.
# "reactions" is a tuple with the reactions, see the examples and the notebook to see how to use them.
# "t" is an np.array with the timepoints for the deterministic calculation. The first and last values are used for
# gillespie algorithm
# "mode": 0: Deterministic and Stochastic solutions; 1: Deterministic only;2: Stochastic only

def ReAct(user_input, reactions, t, mode=0, rounds=1):
    # Create "elements" and "rows"
    # "elements" is a list containing the names of the chemical species in "user_input", in that same order. The rows of the
    # stoichiometry matrix, "Gamma", correspond to the chemical species in "elements". For instance, the second row in
    # "Gamma" will correspond to elements[2]
    # In order to access this relationship i <-> elements[i] fast, we create "row", a dictionary that returns i when we
    # ask rows[elements[i]]
    elements = list()
    row = dict()
    # "init" is a list containing the concentrations of the different chemical species at the beginning, the index
    # corresponds to the index of the species in the list "elements"
    init = list()
    for i in range(0, len(user_input), 2):
        elements.append(user_input[i])
        row.update({user_input[i]: i / 2})
        init.append(user_input[i + 1])

    # Initialize "Gamma" and "decay
    # "Gamma is the stoichiometry matrix. As explained in the notebook, if we use only this matrix to calculate the
    # velocities, species acting as enzymes/catalizers for which we do not include an E+S<->ES step cannot be included.
    # To overcome this problem, we introduce "decay", a tuple that contains the coordinates of the stoichiometry matrix
    # (species, reaction), where a chemical species behaves as an enzyme/catalizer as described above.

    Gamma = np.zeros((len(elements), len(reactions) / 3), int)
    decay = tuple()

    # k is a tuple that contains the kinetic constants for the reactions, the index corresponds to the index of the
    # reactions in the reactions tuple provided by the user
    k = tuple()

    # Iterate through the reactions tuple: j +0 :Reactants/ +1 Products/ +2 Kinetic constants
    # This loop constructs "Gamma" and "k"
    for j in range(0, len(reactions), 3):
        # In each reaction:
        # 1. Iterate through the tuple the reactants, get the value for Gamma, and maybe decay
        for ind in range(0, len(reactions[j]), 2):
            # the odd elements in the tuple are the stoichiometric indexes
            sto = reactions[j][ind]
            # the even elements in the tuple are the names of the species
            name = reactions[j][ind + 1]

            # "normal" indexes
            if sto > 0:
                Gamma[row[name], j / 3] = -sto
            # negative indexes, decay
            else:
                decay += ((row[name], j / 3, abs(sto)),)
        # 2. Iterate along the products
        for ind in range(0, len(reactions[j + 1]), 2):  # In the products, get the value for Gamma
            sto = reactions[j + 1][ind]
            name = reactions[j + 1][ind + 1]
            if sto > 0:
                Gamma[row[name], j / 3] = sto
        # 3. Get the value of k for this reaction
        k += (reactions[j + 2],)

    # Print the reaction in a "textbook" form
    for i in range(0, Gamma.shape[1]):  # Along the reactions
        reactants = list()
        products = list()
        catalizers = list()  # List of the catalizers in "decay"
        Ri = Gamma[:, i]
        for j in range(0, Gamma.shape[0]):  # Along the reactants
            if Ri[j] < 0:
                reactants.append(str(abs(Ri[j])) + elements[j])
            elif Ri[j] > 0:
                products.append(str(abs(Ri[j])) + elements[j])
            if decay:
                for ele in decay:
                    if (ele[0], ele[1]) == (j, i):
                        catalizers.append(str(ele[2]) + elements[j])
                        # reactants.append(str(ele[2]) + elements[j])
                        # products.append(str(ele[2]) + elements[j])
        #print reactants, '--k' + str(i) + str(catalizers) + '-->', products

    # Different return depending on "mode"
    if mode == 0:  # mode 0: Gillespie and Deterministic solution
        return (
        DetSol(elements, init, t, Gamma, k, decay), Gillespy(elements, init, t, Gamma, k, decay, rounds), row, mode)
    elif mode == 1:  # mode 1: Deterministic solution only
        return (DetSol(elements, init, t, Gamma, k, decay), (None, None, None, None), row, mode)
    elif mode == 2:  # mode 2: Gillespie only
        return (None, Gillespy(elements, init, t, Gamma, k, decay, rounds), row, mode)
    elif mode ==3: # mode 3: steady state

        linsys = LinearSystem(Gamma,k)

        # TODO This could be done but one has to know how to identify them
        # # WE combine the last two just to make it linearly independent
        # linsys[-2,:] += linsys[-1,:]
        #
        # linsys = linsys[:-1,:]

        # we add the condition that the sum gives the total

        sum_total = np.ones([1,linsys.shape[1]])

        linsys = np.concatenate((linsys, sum_total), axis=0)


        # All equations are equal to zero except for the sum of total

        solution = np.zeros(linsys.shape[0])

        solution[-1] = np.sum(init)

        return np.linalg.lstsq(linsys,solution)


def Gillesplot(solution, t, tgill, valsgill, rows, mode, which2plot=False):
    # which2plot is a list of strings containing the names of the chemical species that we want to plot
    # If there is not a specific list, plot all of them
    if not which2plot:
        which2plot = rows.keys()
    # fig, ax = plt.subplots()
    ax = plt.gca()
    cmap = plt.get_cmap('jet')
    line_colors = cmap(np.linspace(0, 1, len(which2plot)))
    colorind = 0
    for i in which2plot:
        # color0=np.random.rand(3)
        # color0 = color0*0.8 / np.linalg.norm(color0)
        # color0=tuple(color0) #This method is supposed to generate random colors, but colormaps are better

        if mode != 2:
            ax.plot(t, solution[:, rows[i]], color=line_colors[colorind], label=i)
        if mode != 1:
            for j in range(len(tgill)):
                ax.step(tgill[j], valsgill[j][rows[i], :], color=line_colors[colorind])
        colorind += 1
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=5)
    plt.draw()
