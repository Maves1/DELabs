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

def calcEulerLTE(exactYs, x_init, step, epsilon):
    ltes = [0]
    x_curr = x_init
    for i in range(1, len(exactYs)):
        curr_error = abs(eulerCalcNext(x_curr, exactYs[i - 1], step) - exactYs[i])
        ltes.append(curr_error)
        x_curr += step
    return ltes

def calcWowEulerLTE(exactYs, x_init, step, epsilon):
    ltes = [0]
    x_curr = x_init
    for i in range(1, len(exactYs)):
        curr_error = abs(wowEulerCalcNext(x_curr, exactYs[i - 1], step) - exactYs[i])
        ltes.append(curr_error)
        x_curr += step
    return ltes

def calcRKLTE(exactYs, x_init, step, epsilon):
    ltes = [0]
    x_curr = x_init
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

exactYs = [y_init]

# Euler method
eulerApproximationXs = [x_init]
eulerApproximationYs = [y_init]

x_prev = x_init
x_curr = x_init + h1
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    exactYs.append(calcExactSolution(x_curr))

    eulerApproximationXs.append(x_curr)
    y_curr = eulerCalcNext(x_prev, y_prev, h1)
    eulerApproximationYs.append(y_curr)
    x_prev = x_curr
    x_curr += h1
    y_prev = y_curr

eulerGtes = calcGTE(exactYs, eulerApproximationYs)
print(*eulerGtes)

eulerLtes = calcEulerLTE(exactYs, x_init, h1, epsilon)
print(*eulerLtes)

# Improved Euler method
wowEulerApproximationXsh1 = [x_init]
wowEulerApproximationYsh1 = [y_init]

x_prev = x_init
x_curr = x_init + h1
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    wowEulerApproximationXsh1.append(x_curr)
    y_curr = wowEulerCalcNext(x_prev, y_prev, h1)
    wowEulerApproximationYsh1.append(y_curr)
    x_prev = x_curr
    x_curr += h1
    y_prev = y_curr

print()
print("Improved Euler ys:")
print(*wowEulerApproximationYsh1)

wowEulerGtes = calcGTE(exactYs, wowEulerApproximationYsh1)
print(*wowEulerGtes)

wowEulerLtes = calcWowEulerLTE(exactYs, x_init, h1, epsilon)
print(*wowEulerLtes)

wowEulerApproximationXsh2 = [x_init]
wowEulerApproximationYsh2 = [y_init]

x_prev = x_init
x_curr = x_init + h2
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    wowEulerApproximationXsh2.append(x_curr)
    y_curr = wowEulerCalcNext(x_prev, y_prev, h2)
    wowEulerApproximationYsh2.append(y_curr)
    x_prev = x_curr
    x_curr += h2
    y_prev = y_curr

print()
print("Improved Euler ys:")
print(*wowEulerApproximationYsh1)

wowEulerGtesh2 = calcGTE(exactYs, wowEulerApproximationYsh2)
print(*wowEulerGtes)

wowEulerLtesh2 = calcWowEulerLTE(exactYs, x_init, h2, epsilon)
print(*wowEulerLtes)

# Runge-Kutta
rgApproximationXs = [x_init]
rgApproximationYs = [y_init]

x_prev = x_init
x_curr = x_init + h1
y_prev = y_init
while x_curr <= rightBoundary + epsilon:
    rgApproximationXs.append(x_curr)
    y_curr = rungeKuttaCalcNext(x_prev, y_prev, h1)
    rgApproximationYs.append(y_curr)
    x_prev = x_curr
    x_curr += h1
    y_prev = y_curr

print()
print("Runge Kutta ys:")
print(*rgApproximationYs)
print("Exact ys")
print(*exactYs)

rgGtes = calcGTE(exactYs, rgApproximationYs)
print(*rgGtes)

rgLtes = calcRKLTE(exactYs, x_init, h1, epsilon)
print(*rgLtes)