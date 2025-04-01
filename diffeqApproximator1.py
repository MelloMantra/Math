#1st Order DiffEQ Numeric Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sympy
from sympy import *
import asyncio
import re

func = None
x,y = symbols("x y")
dx = None

def euler(func, dx, condition):
    return [condition[0]+dx,condition[1] + dx*(func.evalf(subs={x:condition[0], y:condition[1]}))]

def taylor2(func, dx, condition):
    return [condition[0]+dx, euler(func, dx, condition)[1]+0.5*dx*dx*(diff(func,x)+func*diff(func,y)).evalf(subs={x:condition[0], y:condition[1]})]

def taylor3(func, dx, condition):
    return [condition[0]+dx, taylor2(func, dx, condition)[1]+(1/6)*dx*dx*dx*(diff(func,x,x)+(2*func*diff(func,x,y))+(func*func*diff(func,y,y))+(diff(func,x)*diff(func,y))+(func*diff(func,y,y)*diff(func,y,y))).evalf(subs={x:condition[0], y:condition[1]})]

def rungekutta(func, dx, condition):
    return [condition[0]+dx, condition[1]+dx*func.evalf(subs={x:(condition[0]+(dx*0.5)), y:(condition[1]+0.5*dx*func.evalf(subs={x:condition[0],y:condition[1]}))})]


def graph(eu,t2,t3,rk):

    methods = [eu,t2,t3,rk]
    colors = ["red","blue","green","purple"]
    labels = ["Euler", "Taylor O(2)", "Taylor O(3)", "Runge-Kutta"]

    for i in range(len(methods)):
        xvals = []
        yvals = []
        for pt in methods[i]:
            xvals.append(pt[0])
            yvals.append(pt[1])

        plt.plot(xvals,yvals,colors[i], label=labels[i])
        print(labels[i]+":")
        print(methods[i])

    plt.title(f"dy/dx = {func}")
    plt.legend()
    plt.show()

async def main():

    global func
    global dx
    x,y = symbols("x y")
    func = sympify(input("dy/dx = "))
    ic = re.split(",",input("Initial condition (x,y): "))
    dx = eval(input("Step size: "))
    span = re.split(",",input("Interval [a,b]: "))

    for i in range(2): ic[i] = eval(ic[i])
    for i in range(2): span[i] = eval(span[i])

    eu = [ic]
    t2 = [ic]
    t3 = [ic]
    rk = [ic]

    for i in range(int((span[1]-span[0])/dx)):
        eu.append(euler(func, dx, eu[i]))
        t2.append(taylor2(func,dx,t2[i]))
        t3.append(taylor3(func,dx,t3[i]))
        rk.append(rungekutta(func,dx,rk[i]))
    
    graph(eu, t2, t3, rk)


asyncio.run(main())