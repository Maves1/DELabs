# Euler
def eulerCalcNext(x_prev, y_prev, step):
    y_curr = y_prev + step * (y_prev ** 2 + x_prev * y_prev - x_prev ** 2) / (x_prev ** 2)
    return y_curr

# Improved Euler
def wowEulerCalcNext(x_prev, y_prev, step):
    f = (y_prev ** 2 + x_prev * y_prev - x_prev ** 2) / (x_prev ** 2)
    x = x_prev + step / 2
    y = y_prev + step / 2 * f

    y_curr = y_prev + step * (y ** 2 + x * y - x ** 2) / (x ** 2)
    return y_curr

# Runge-Kutta
def rungeKuttaCalcNext(x_prev, y_prev, step):
    k1 = (y_prev ** 2 + x_prev * y_prev - x_prev ** 2) / (x_prev ** 2)

    x = x_prev + step / 2
    y = y_prev + step * k1 / 2
    k2 = (y ** 2 + x * y - x ** 2) / (x ** 2)

    y = y_prev + step * k2 / 2
    k3 = (y ** 2 + x * y - x ** 2) / (x ** 2)

    x = x_prev + step
    y = y_prev + step * k3
    k4 = (y ** 2 + x * y - x ** 2) / (x ** 2)

    y_curr = y_prev + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return y_curr

def calcExactSolution(x):
    y = (x * (1 + x ** 2 / 3) ) / (1 - x ** 2 / 3)
    return y

def calcGTE(exactYs, eulerYs):
    gtes = []
    for i in range(len(exactYs)):
        gtes.append(abs(eulerYs[i] - exactYs[i]))
    return gtes

def calcLTE(exactYs, x_init, step, epsilon, method):
    ltes = [0]
    x_curr = x_init
    if method == "euler":
        for i in range(1, len(exactYs)):
            curr_error = abs(eulerCalcNext(x_curr, exactYs[i - 1], step) - exactYs[i])
            ltes.append(curr_error)
            x_curr += step
    elif method == "wowEuler":
        for i in range(1, len(exactYs)):
            curr_error = abs(wowEulerCalcNext(x_curr, exactYs[i - 1], step) - exactYs[i])
            ltes.append(curr_error)
            x_curr += step
    else:
        for i in range(1, len(exactYs)):
            curr_error = abs(rungeKuttaCalcNext(x_curr, exactYs[i - 1], step) - exactYs[i])
            ltes.append(curr_error)
            x_curr += step
    return ltes


# Variables shared by all methods
epsilon = 0.1

h1 = 0.1
h2 = 0.05
h3 = 0.01

leftBoundary = 1.0
rightBoundary = 1.5

x_init = 1
y_init = 2

# Euler method
exactYs = [y_init]

eulerApproximationXs = [x_init]
eulerApproximationYs = [y_init]

x_prev = x_init
x_curr = x_init + h1
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    exactYs.append(calcExactSolution(x_curr))
    y_curr = rungeKuttaCalcNext(x_prev, y_prev, h1)

    eulerApproximationXs.append(x_curr)
    eulerApproximationYs.append(y_curr)

    x_prev = x_curr
    x_curr += h1
    y_prev = y_curr

eulerGtes = calcGTE(exactYs, eulerApproximationYs)
print(*eulerGtes)

eulerLtes = calcLTE(exactYs, x_init, h1, epsilon, "rk")
print(*eulerLtes)

# h2
exactYs2 = [y_init]
eulerApproximationXs2 = [x_init]
eulerApproximationYs2 = [y_init]

x_prev = x_init
x_curr = x_init + h2
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    exactYs2.append(calcExactSolution(x_curr))
    y_curr = rungeKuttaCalcNext(x_prev, y_prev, h2)

    eulerApproximationXs2.append(x_curr)
    eulerApproximationYs2.append(y_curr)

    x_prev = x_curr
    x_curr += h2
    y_prev = y_curr

eulerGtes2 = calcGTE(exactYs2, eulerApproximationYs2)
print(*eulerGtes2)

eulerLtes2 = calcLTE(exactYs2, x_init, h2, epsilon, "rk")
print(*eulerLtes2)

for i in range(1, len(eulerLtes)):
    index = eulerApproximationXs2.index(eulerApproximationXs[i])
    # print(len(eulerLtes), len(eulerLtes2))
    print(eulerApproximationYs2[index])
    print("x =", eulerApproximationXs[i], "\nLTE1 / LTE2:", eulerLtes[i] / eulerLtes2[index])
    print("x =", eulerApproximationXs[i], "\nGTE1 / GTE2:", eulerGtes[i] / eulerGtes2[index])
    print()

