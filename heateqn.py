# libraries
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.animation as animation
import click

print("Uniform Initial Temperature, Single Dimensional Heat Equation Modeling Tool")

debugging = True

if not debugging:
    # initial conditions
    l = eval(input("Rod length: "))
    T1 = eval(input("u(0,t) = "))
    T2 = eval(input(f"u({l},t) = "))
    T0 = eval(input("u(x,0) = "))
    dt = eval(input("Time step: "))
    tsteps = int(eval(input("Total time steps: ")))
    xdivs = int(eval(input("# of x-divisions: ")))
else:
    # debugging initialization
    l = 1.0
    T1 = 100.0
    T2 = 40.0
    T0 = 20.0
    dt = 0.5
    tsteps = 10
    xdivs = 5
dx = l/(xdivs)
# get decimal places to round to based on decimal places in dx
xplaces = len(str(dx).split(".")[1])
tplaces = len(str(dt).split(".")[1])

# saved values
#u = {"x,t": u(x,t)}
U = {}
def u(x,t):
    if x==0:
        return T1
    if x==l:
        return T2
    if t==0:
        return T0
    return U[f"{round(x,xplaces)},{round(t,tplaces)}"]

# differentiation functions
    
def getB(x,t):
    b = 0
    if x-dx == 0:
        b += (T1/np.square(dx))/dt
    if x+dx == l:
        b += (T2/np.square(dx))/dt
    b += u(x,t-dt)/dt

    if type(b)==np.float64:
        b = round(b.item(),10)
    return round(b,tplaces)

def getSystem(t):
    A = []
    for i in range(xdivs-1):
        A.append([])
        for j in range(xdivs-1):
            if i==j:
                A[i].append(round(((2/np.square(dx))+1)*dt,tplaces))
            elif i-1==j or i+1==j:
                A[i].append(round((-1/np.square(dx))*dt,tplaces))
            else:
                A[i].append(0)

    b = []
    x = dx
    while round(x,2)<l:
        b.append([getB(round(x,xplaces),round(t,tplaces))])
        x+=dx
    
    return (np.array(A),np.array(b))



def graph(dict):
    xvals = []
    yvals = []
    for x in range(xdivs+1):
        xvals.append(round(x*dx,xplaces))
        yvals.append(u(round(x*dx,xplaces),0))

    fig, ax = plot.subplots()
    line, = ax.plot(xvals, yvals)
    def animate(i):
        yvals = []
        for x in range(xdivs+1):
            yvals.append(u(round(x*dx,xplaces),round(i*dt,tplaces)))
        line.set_ydata(yvals)  # update the data.
        return line,
    ani = animation.FuncAnimation(
        fig, animate, frames=tsteps, interval=200, blit=True)
    plot.ylabel("Temperature (C)")
    plot.xlabel("Position Along Rod (m)")
    plot.show()

# main
def main():
    t = 0
    for iteration in range(tsteps):
        t += dt
        sys = getSystem(t)
        soln = np.matmul(np.linalg.inv(sys[0]),sys[1])
        print(soln)
        for entry in range(len(soln)):
            U[f"{round(dx*(entry+1),xplaces)},{round(t,tplaces)}"] = soln[entry][0]
    graph(U)
  
# run
main()