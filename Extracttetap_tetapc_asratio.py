# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 12:45:46 2023

@author: User
"""

# Modify the values as needed
d = 21.5
bf = 11.8
tf = 2.3
hw = 16.9
tw = 1.28
Radius = 0.4
ASec = 76.0
AS2 = 27.52
AS3 = 45.23
J = 103.0
I22 = 628.0
I33 = 5510.0
S22 = 106.44
S33 = 512.56
Z22 = 166.0
Z33 = 611.0
R22 = 2.8746
R33 = 8.5147

#material properties
Fy = 67.93
Es = 31186.21

nFac = 10
inchToCurrUnit = 1
kipsToCurrUnit = 1
isA992Gr50 = 1
# Additional parameters
LbBeam = 356.7  # Example value, replace with the actual value


# Original code
def computeHingeRBSBeam(d, tw, bf, tf, ISec, ZSec, ry, Lmem, Lb, Es, Fy, nFac, inchToCurrUnit, kipsToCurrUnit, isA992Gr50, cRBS, LRatio):
    ksiToCurrUnit = kipsToCurrUnit / (inchToCurrUnit ** 2)
    # calculate thetaP, thetaPC and lambda
    hw = d - (2.0 * tf)
    LShear = Lmem / 2.0
    L = LShear
    if isA992Gr50 == 1:
        thetaP = 0.09 * (hw / tw) ** -0.3 * (bf / tf / 2.0) ** -0.1 * (L / d) ** 0.1 * (d / 21. / inchToCurrUnit) ** -0.8
        thetaP = thetaP / LRatio
        thetaPC = 6.5 * (hw / tw) ** -0.5 * (bf / tf / 2.0) ** -0.9
        thetaPC = thetaPC / LRatio
        Lambda = 49.0 * (hw / tw) ** -1.14 * (bf / tf / 2.0) ** -0.632 * (Lb / ry) ** -0.205 * (Es / Fy) ** 0.391
    else:
        thetaP = 0.19 * (hw / tw) ** -0.314 * (bf / tf / 2.0) ** -0.1 * (Lb / ry) ** -0.1185 * (L / d) ** 0.113 * (
                    d / 21. / inchToCurrUnit) ** -0.76 * (Fy / 50.0 / ksiToCurrUnit) ** -0.07
        thetaP = thetaP / LRatio
        thetaPC = 9.62 * (hw / tw) ** -0.513 * (bf / tf / 2.0) ** -0.863 * (Lb / ry) ** -0.108 * (
                    Fy / 50.0 / ksiToCurrUnit) ** -0.36
        thetaPC = thetaPC / LRatio
        Lambda = 592.0 * (hw / tw) ** -1.138 * (bf / tf / 2.0) ** -0.632 * (Lb / ry) ** -0.205 * (
                    Fy / 50.0 / ksiToCurrUnit) ** -0.391
    # calculate effective yield moment(My)
    beta = 1.1
    ZRBS = ZSec - (2.0 * cRBS * tf * (d - tf))
    Myp = ZRBS * Fy
    My = LRatio * beta * Myp

    # define peak moment (Mu) to effective yield moment (My) ratio

    MuMyFac = 1.1

    # define some parameters

    c = 1.0
    D = 1.0

    KResidual = 0.4
    thetaU = 0.2 / LRatio

    # calculate member stiffness and strain hardening coefficient
    ke = 6.0 * Es * (ISec / Lmem)
    asRatio = My * ((MuMyFac - 1) / (ke * thetaP))
    return thetaP, thetaPC, asRatio

# Calculate shRBS and Lmem based on provided formulas
aRBS = 0.625 * bf
bRBS = 0.75 * d
cRBS = 0.25 * bf
shRBS = aRBS + bRBS / 2
LRatio = 1.0
Lmem = LbBeam - 2 * shRBS



# Call computeHingeRBSBeam with the section properties and additional parameters
beamHingeParams = computeHingeRBSBeam(d, tw, bf, tf, I33, Z33, R22, Lmem, LbBeam, Es, Fy, nFac, inchToCurrUnit,
                                     kipsToCurrUnit, isA992Gr50, cRBS, LRatio)

# Print or use the computed parameters as needed
print(beamHingeParams)

#%%
d = 27.1
bf = 13.4
tf = 2.28
hw = 22.540000000000003
tw = 1.26
Radius = 0.5
ASec = 89.7
AS2 = 34.15
AS3 = 50.92
J = 117.0
I22 = 919.0
I33 = 10700.0
S22 = 137.16
S33 = 789.67
Z22 = 214.0
Z33 = 922.0
R22 = 3.2008
R33 = 10.9218

#material properties
Fy = 67.93
Es = 31186.21

nFac = 10
inchToCurrUnit = 1
kipsToCurrUnit = 1
isA992Gr50 = 1
HStory = 164.03
deckSpacing = 20
deckThick = 12
LbCol = HStory-deckThick
LbCol = HStory-0.2
Lmem=HStory

#Pg
Pg = 245.48


def computeHingeWColumn(d, tw, bf, tf, I33, Z33, R22, Lmem, Lb, Es, Fy, nFac, inchToCurrUnit, kipsToCurrUnit, isA992Gr50, ASec, Pg):

    ksiToCurrUnit = kipsToCurrUnit / inchToCurrUnit ** 2

    # calculate thetaP, thetaPC and lambda
    hw = d - 2 * tf
    LShear = Lmem / 2
    L = LShear
    Pg = abs(Pg)
    Py = Fy * ASec

    thetaP = min((294 * (hw / tw) ** -1.7 * (Lb / R22) ** -0.7 * (1 - Pg / Py) ** 1.6), 0.2)
    thetaPC = min((90 * (hw / tw) ** -0.8 * (Lb / R22) ** -0.8 * (1 - Pg / Py) ** 2.5), 0.3)

    if isA992Gr50:
        Lambda = 85 * (hw / tw) ** -1.26 * (bf / tf / 2) ** -0.525 * (Lb / R22) ** -0.13 * (Es / Fy) ** 0.291
    else:
        Lambda = 500 * (hw / tw) ** -1.34 * (bf / tf / 2) ** -0.595 * (Fy / 50 / ksiToCurrUnit) ** -0.36

    if Pg / Py > 0.6:
        print('Warning: Pg/Pye > 0.6: current hinge model for columns may not be suitable')
        print('columns need to be treated as force-controlled elements')
        print('see computeHingeWColumn.tcl file and ASCE/SEI 41-17')

    # calculate effective yield moment(My)
    beta = 1.2
    Myp = Z33 * Fy
    My = beta * Myp

    # calculate effective yield moment reduced by the applied load (My)
    # it is showed by My* in the references
    if Pg / Py <= 0.2:
        My = 1.15 * Z33 * Fy * (1 - Pg / Py)
    else:
        My = 1.15 * Z33 * Fy * (9 / 8) * (1 - Pg / Py)

    # calculate peak moment (Mu) to effective yield moment (My) ratio
    alpha = 12.5 * (hw / tw) ** -0.2 * (Lb / R22) ** -0.4 * (1 - Pg / Py) ** 0.4
    if alpha < 1:
        alpha = 1
    elif alpha > 1.3:
        alpha = 1.3

    MuMyFac = alpha
    # define some parameters
    c = 1
    D = 1
    KResidual = 0.5 - 0.4 * Pg / Py
    thetaU = 0.15

    # calculate member stiffness and strain hardening coefficient
    ke = 6 * Es * I33 / Lmem
    asRatio = My * (MuMyFac - 1) / (ke * thetaP)

    
    return thetaP, thetaPC, asRatio



# Call computeHingeWColumn with the section properties and additional parameters
beamHingeParams = computeHingeWColumn(d, tw, bf, tf, I33, Z33, R22, Lmem, LbCol, Es, Fy, nFac,
                                      inchToCurrUnit, kipsToCurrUnit, isA992Gr50, ASec, Pg)

# Print or use the computed parameters as needed
print(beamHingeParams)

#%%
from sklearn.tree import plot_tree
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier

iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(iris.data, iris.target, random_state=42)
clf = DecisionTreeClassifier()
clf.fit(X_train, y_train)

plt.figure(figsize=(12, 8))
plot_tree(clf, filled=True, feature_names=iris.feature_names, class_names=iris.target_names, rounded=True)
plt.show()
