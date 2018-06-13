import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, datasets


def normal_data(set, config):
    data_vector = []

    f_tm = (set["F_tm"]-int(config["Tm"]["Min"]))/(int(config["Tm"]["Max"])-int(config["Tm"]["Min"]))
    data_vector.append(f_tm)
    r_tm = (set["R_tm"]-int(config["Tm"]["Min"]))/(int(config["Tm"]["Max"])-int(config["Tm"]["Min"]))
    data_vector.append(r_tm)

    f_gc = (set["F_gc"] - int(config["GC Percent"]["Min"])) / (int(config["GC Percent"]["Max"]) - int(config["GC Percent"]["Min"]))
    data_vector.append(f_gc)
    r_gc = (set["R_gc"] - int(config["GC Percent"]["Min"])) / (int(config["GC Percent"]["Max"]) - int(config["GC Percent"]["Min"]))
    data_vector.append(r_gc)

    f_len = (set["F_len"] - int(config["Length"]["Min"])) / (int(config["Length"]["Max"]) - int(config["Length"]["Min"]))
    data_vector.append(f_len)
    r_len = (set["R_len"] - int(config["Length"]["Min"])) / (int(config["Length"]["Max"]) - int(config["Length"]["Min"]))
    data_vector.append(r_len)

    tm_dif = (set["Tm_dif"] - int(config["Temperature difference"]["Min"])) / (int(config["Temperature difference"]["Max"])
                                                                              - int(config["Temperature difference"]["Min"]))
    data_vector.append(tm_dif)

    amp_len = (set["Amp_len"] - int(config["Amplicon Length"]["Min"])) / (int(config["Amplicon Length"]["Max"]) - int(config["Amplicon Length"]["Min"]))
    data_vector.append(r_len)

    return data_vector


#def init():
config = json.load(open("config.json"))

# import some data to play with
#iris = datasets.load_iris()
#X = iris.data[:, :2]  # we only take the first two features.
#print(X)
#Y = iris.target

data = json.load(open("Primers_data.json"))
X_sets_vector = []
Y_scores = []
for set in data["Sets"]:
    X_sets_vector.append(normal_data(set,config))
    Y_scores.append(set["score"])

h = .02  # step size in the mesh

logreg = linear_model.LogisticRegression(C=1e5)

# we create an instance of Neighbours Classifier and fit the data.
logreg.fit(X_sets_vector,Y_scores)
print("HERE")

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, x_max]x[y_min, y_max].
#x_min, x_max = X_sets_vector[:, 0].min() - .5, X_sets_vector[:, 0].max() + .5
#y_min, y_max = X_sets_vector[:, 1].min() - .5, X_sets_vector[:, 1].max() + .5
#xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = logreg.predict([X_sets_vector[1]])
print(Z)
# Put the result into a color plot
#Z = Z.reshape(xx.shape)'''




